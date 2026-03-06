function generate_dashboard_direct()
% GENERATE_DASHBOARD_DIRECT  Unified pipeline to generate all dashboard images.
%
%   GENERATE_DASHBOARD_DIRECT() reads the experiment log, filters to T4T5
%   Summer/Autumn 2025 cells (excluding test genotype), and generates all
%   dashboard PNG images by:
%     1. Converting Protocol 1 GridPlots (.fig → PNG)
%     2. Loading saved bar analysis results → regenerating polar plots as PNGs
%     3. Running flash RF analysis from raw data (4px) and loading saved
%        results (6px) → generating timeseries, heatmap, and Gaussian contour PNGs
%     4. Loading saved bar flash results (Autumn only) → generating bar flash PNGs
%
%   Also extracts analysis metrics and writes them into cell_index.json.
%
%   Cells that fail are logged to failed_cells.txt in the output directory.
%
%   PREREQUISITES:
%     >> addpath(genpath('/Users/burnettl/Documents/GitHub/nested_RF_stimulus/src'))
%     >> generate_dashboard_direct
%
%   See also PROCESS_BARS_P2, PROCESS_FLASH_P2, PROCESS_BAR_FLASHES_P2,
%            BUILD_CELL_INDEX, EXPORT_DASHBOARD_IMAGES

%% Configuration
PROJECT_ROOT = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2';
EXCEL_PATH   = '/Users/burnettl/Documents/GitHub/nested_RF_stimulus/reiser-lab-ephys-experiment-log.xlsx';
P1_BASE      = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol1/bkg4/protocol1_10kHz_4reps_12px_6px_LHS_2sbkg_200msfl_150msint_07-10-25_15-45-60/';
P2_SUMMER    = fullfile(PROJECT_ROOT, 'data', 'T4T5_Summer2025');
P2_AUTUMN    = fullfile(PROJECT_ROOT, 'data', 'T4T5_autumn2025');
OUTPUT_DIR   = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/dashboard_images/';

FLASH_RESULTS_DIR = fullfile(PROJECT_ROOT, 'results', 'flash_results');
BAR_RESULTS_DIR   = fullfile(PROJECT_ROOT, 'results', 'bar_results');
BAR_FLASH_DIR     = fullfile(PROJECT_ROOT, 'results', 'bar_flash_results');
BAR_FLASH_1DRF    = fullfile(PROJECT_ROOT, 'results', '1DRF', 'bar_flash_data');
MARCM_BAR_DIR     = fullfile(PROJECT_ROOT, 'results', 'MARCM_ttl_2025', 'protocol2_bar');

EXPORT_OPTS  = {'Resolution', 300, 'BackgroundColor', 'white'};
THUMB_WIDTH  = 400;
FAILURE_LOG  = fullfile(OUTPUT_DIR, 'failed_cells.txt');

%% Create output subdirectories
img_dirs = {'gridplots', 'bar_polar', 'polar_arrow', ...
            'flash_rf_4px', 'flash_heatmap_4px', ...
            'flash_rf_6px', 'flash_heatmap_6px', ...
            'gaussian_rf', 'bar_flash_slow', 'bar_flash_fast'};
for k = 1:numel(img_dirs)
    d = fullfile(OUTPUT_DIR, img_dirs{k});
    if ~isfolder(d), mkdir(d); end
    t = fullfile(OUTPUT_DIR, 'thumbs', img_dirs{k});
    if ~isfolder(t), mkdir(t); end
end

%% Read and filter experiment log
T = readtable(EXCEL_PATH, 'Sheet', 'Form responses 1', 'VariableNamingRule', 'preserve');
T.Properties.VariableNames = strrep(T.Properties.VariableNames, ' ', '_');
T.Properties.VariableNames = strrep(T.Properties.VariableNames, '?', '_');
T.Properties.VariableNames = strrep(T.Properties.VariableNames, '/', '_');

mask = contains(T.project_type, 'Summer2025') | contains(T.project_type, 'Autumn2025');
mask = mask & ~contains(T.strain, 'test', 'IgnoreCase', true);
cells = T(mask, :);
n_cells = height(cells);
fprintf('Processing %d cells (Summer + Autumn 2025, excluding test)\n', n_cells);

%% Initialize metrics storage
cell_metrics = cell(n_cells, 1);
cell_errors_all = cell(n_cells, 1);

%% Open failure log
fid_fail = fopen(FAILURE_LOG, 'w');
fprintf(fid_fail, 'cell_id\texp_folder\terror_message\n');

%% Main processing loop
orig_dir = pwd;
n_success = 0;
n_fail = 0;

for i = 1:n_cells
    date_str = safe_str(cells.date(i));
    time_str = safe_str(cells.time(i));
    cell_id  = sprintf('%s_%s', date_str, time_str);
    strain   = safe_str(cells.strain(i));
    pc       = safe_str(cells.preferred_contrast(i));
    cell_errors = {};  % Accumulate per-section errors for this cell

    if strcmpi(pc, 'ON')
        on_off = 'on';
    else
        on_off = 'off';
    end

    % Determine genotype
    if ismember('genotype', cells.Properties.VariableNames) && ...
            is_valid_str(cells.genotype(i))
        genotype = safe_str(cells.genotype(i));
    else
        genotype = strain;
    end

    % Determine if this is a Summer cell (no bar flash stimulus)
    is_summer = contains(string(safe_str(cells.project_type(i))), 'Summer2025');

    fprintf('\n[%d/%d] %s (%s, %s, %s)\n', i, n_cells, cell_id, strain, pc, ...
        ternary(is_summer, 'Summer', 'Autumn'));

    % Find experiment folder (with fuzzy ±2 min matching)
    exp_folder = find_experiment_folder({P2_SUMMER, P2_AUTUMN}, cell_id, date_str, time_str);

    if isempty(exp_folder)
        fprintf('  WARNING: No experiment folder found for %s\n', cell_id);
        cell_errors{end+1} = 'No experiment folder found';
    elseif ~strcmp(exp_folder(end-length(cell_id)+1:end), cell_id)
        fprintf('  NOTE: Matched to folder %s (fuzzy time match)\n', exp_folder);
    end

    % Create results subfolder for new saves
    results_folder = fullfile(PROJECT_ROOT, 'results', genotype, on_off);
    if ~isfolder(results_folder)
        mkdir(results_folder);
    end

    % Initialize metrics for this cell
    m = struct();
    m.strain = strain;
    m.preferred_contrast = pc;
    m.genotype = genotype;

    % Metadata from currentExp.mat
    metadata = struct();  % default in case currentExp.mat is missing
    if ~isempty(exp_folder)
        meta_path = fullfile(exp_folder, 'currentExp.mat');
        if isfile(meta_path)
            S = load(meta_path, 'metadata');
            metadata = S.metadata;
            if isfield(metadata, 'Frame'), m.frame = metadata.Frame; end
            if isfield(metadata, 'Side'), m.side = metadata.Side; end
            if isfield(metadata, 'Age'), m.age = metadata.Age; end
        end
    end

    try
        cleanup = onCleanup(@() cd(orig_dir));

        %% ====== Section 1: GridPlots (.fig → 4 PNGs) ======
        gridplot_1_path = fullfile(OUTPUT_DIR, 'gridplots', [cell_id '_gridplot_1.png']);
        if ~isfile(gridplot_1_path)
            fig_path = find_matching_p1_gridplot(P1_BASE, date_str, time_str);
            if ~isempty(fig_path) && isfile(fig_path)
                try
                    fprintf('  Exporting GridPlots...\n');
                    figs = openfig(fig_path, 'invisible');
                    for fi = 1:numel(figs)
                        figs(fi).Color = 'w';
                        gp_out = fullfile(OUTPUT_DIR, 'gridplots', ...
                            sprintf('%s_gridplot_%d.png', cell_id, fi));
                        exportgraphics(figs(fi), gp_out, EXPORT_OPTS{:});
                        make_thumbnail(gp_out, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'gridplots', ...
                            sprintf('%s_gridplot_%d.png', cell_id, fi)), THUMB_WIDTH);
                    end
                    close(figs);
                catch ME_gp
                    warning('  GridPlot failed: %s', ME_gp.message);
                    cell_errors{end+1} = sprintf('GridPlot: %s', ME_gp.message);
                end
                close all force
            else
                fprintf('  No matching P1 GridPlot found\n');
                cell_errors{end+1} = 'No matching P1 GridPlot';
            end
        else
            fprintf('  GridPlots already exist, skipping\n');
        end

        %% ====== Section 2: Bar Analysis (bar_polar + polar_arrow) ======
        bar_polar_path  = fullfile(OUTPUT_DIR, 'bar_polar', [cell_id '_bar_polar.png']);
        polar_arrow_path = fullfile(OUTPUT_DIR, 'polar_arrow', [cell_id '_polar_arrow.png']);

        resultant_angle = 0;  % default

        try  % --- Section 2 try/catch (isolated from flash/bar_flash) ---

        if ~isfile(bar_polar_path) || ~isfile(polar_arrow_path)
            % Search for saved peak_vals
            pv_path = find_peak_vals(BAR_RESULTS_DIR, strain, on_off, date_str, time_str, MARCM_BAR_DIR);

            if ~isempty(pv_path)
                fprintf('  Loading saved bar results: %s\n', pv_path);
                PV = load(pv_path, 'bar_results', 'data');

                params_bar = struct('date', date_str, 'time', time_str, ...
                    'strain', strain, 'on_off', on_off);

                % Check speed count: plot functions hardcode 3 speeds
                n_data_rows = size(PV.data, 1);
                n_speeds = n_data_rows / 16;

                if n_speeds == 3
                    % Bar polar timeseries (3-speed only)
                    if ~isfile(bar_polar_path)
                        [~, ~] = plot_timeseries_polar_bars(PV.data, ...
                            PV.bar_results.median_voltage, params_bar, false, '');
                        fig_bp = gcf;
                        export_and_thumb(fig_bp, bar_polar_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'bar_polar', [cell_id '_bar_polar.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(fig_bp);
                    end

                    % Polar arrow (3-speed only)
                    if ~isfile(polar_arrow_path)
                        max_v = PV.bar_results.max_v{1};
                        resultant_angle = plot_polar_with_arrow(max_v, ...
                            PV.bar_results.median_voltage, params_bar, false, '');
                        fig_pa = gcf;
                        export_and_thumb(fig_pa, polar_arrow_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'polar_arrow', [cell_id '_polar_arrow.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(fig_pa);
                    else
                        resultant_angle = PV.bar_results.resultant_angle;
                    end
                elseif n_speeds == 2
                    % Bar polar timeseries (2-speed, Summer data)
                    if ~isfile(bar_polar_path)
                        [~, ~] = plot_timeseries_polar_bars_flex(PV.data, ...
                            PV.bar_results.median_voltage, params_bar, false, '');
                        fig_bp = gcf;
                        export_and_thumb(fig_bp, bar_polar_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'bar_polar', [cell_id '_bar_polar.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(fig_bp);
                    end

                    % Polar arrow (2-speed, Summer data)
                    if ~isfile(polar_arrow_path)
                        max_v = PV.bar_results.max_v{1};
                        resultant_angle = plot_polar_with_arrow_flex(max_v, ...
                            PV.bar_results.median_voltage, params_bar, false, '');
                        fig_pa = gcf;
                        export_and_thumb(fig_pa, polar_arrow_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'polar_arrow', [cell_id '_polar_arrow.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(fig_pa);
                    else
                        resultant_angle = PV.bar_results.resultant_angle;
                    end
                else
                    fprintf('  Skipping bar plots: data has %d speeds (expected 2 or 3)\n', n_speeds);
                    resultant_angle = PV.bar_results.resultant_angle;
                end

                % Extract bar metrics (works regardless of speed count)
                m.resultant_angle_deg = rad2deg(PV.bar_results.resultant_angle);
                if isfield(PV.bar_results, 'slow')
                    m.bar_slow = extract_bar_speed_metrics(PV.bar_results, 'slow');
                end
                if isfield(PV.bar_results, 'fast')
                    m.bar_fast = extract_bar_speed_metrics(PV.bar_results, 'fast');
                end
                if isfield(PV.bar_results, 'vfast')
                    m.bar_vfast = extract_bar_speed_metrics(PV.bar_results, 'vfast');
                end

                close all force
            else
                fprintf('  No saved bar results found — running analysis from raw data\n');
                if isempty(exp_folder)
                    cell_errors{end+1} = 'Bar: no peak_vals found and no exp_folder';
                end
                if ~isempty(exp_folder)
                    resultant_angle = process_bars_p2(exp_folder, metadata, PROJECT_ROOT);
                    close all force

                    % Re-find saved peak_vals to generate PNGs
                    pv_path = find_peak_vals(BAR_RESULTS_DIR, strain, on_off, date_str, time_str, MARCM_BAR_DIR);
                    if ~isempty(pv_path)
                        PV = load(pv_path, 'bar_results', 'data');
                        params_bar = struct('date', date_str, 'time', time_str, ...
                            'strain', strain, 'on_off', on_off);

                        n_data_rows = size(PV.data, 1);
                        n_speeds = n_data_rows / 16;
                        if n_speeds == 3
                            [~, ~] = plot_timeseries_polar_bars(PV.data, ...
                                PV.bar_results.median_voltage, params_bar, false, '');
                            fig_bp = gcf;
                            export_and_thumb(fig_bp, bar_polar_path, ...
                                fullfile(OUTPUT_DIR, 'thumbs', 'bar_polar', [cell_id '_bar_polar.png']), ...
                                EXPORT_OPTS, THUMB_WIDTH);
                            close(fig_bp);

                            max_v = PV.bar_results.max_v{1};
                            plot_polar_with_arrow(max_v, PV.bar_results.median_voltage, params_bar, false, '');
                            fig_pa = gcf;
                            export_and_thumb(fig_pa, polar_arrow_path, ...
                                fullfile(OUTPUT_DIR, 'thumbs', 'polar_arrow', [cell_id '_polar_arrow.png']), ...
                                EXPORT_OPTS, THUMB_WIDTH);
                            close(fig_pa);
                        elseif n_speeds == 2
                            [~, ~] = plot_timeseries_polar_bars_flex(PV.data, ...
                                PV.bar_results.median_voltage, params_bar, false, '');
                            fig_bp = gcf;
                            export_and_thumb(fig_bp, bar_polar_path, ...
                                fullfile(OUTPUT_DIR, 'thumbs', 'bar_polar', [cell_id '_bar_polar.png']), ...
                                EXPORT_OPTS, THUMB_WIDTH);
                            close(fig_bp);

                            max_v = PV.bar_results.max_v{1};
                            plot_polar_with_arrow_flex(max_v, PV.bar_results.median_voltage, params_bar, false, '');
                            fig_pa = gcf;
                            export_and_thumb(fig_pa, polar_arrow_path, ...
                                fullfile(OUTPUT_DIR, 'thumbs', 'polar_arrow', [cell_id '_polar_arrow.png']), ...
                                EXPORT_OPTS, THUMB_WIDTH);
                            close(fig_pa);
                        else
                            fprintf('  Skipping bar plots: data has %d speeds (expected 2 or 3)\n', n_speeds);
                        end

                        m.resultant_angle_deg = rad2deg(PV.bar_results.resultant_angle);
                        if isfield(PV.bar_results, 'slow')
                            m.bar_slow = extract_bar_speed_metrics(PV.bar_results, 'slow');
                        end
                        if isfield(PV.bar_results, 'fast')
                            m.bar_fast = extract_bar_speed_metrics(PV.bar_results, 'fast');
                        end
                        if isfield(PV.bar_results, 'vfast')
                            m.bar_vfast = extract_bar_speed_metrics(PV.bar_results, 'vfast');
                        end
                    end
                    close all force
                end
            end
        else
            fprintf('  Bar plots already exist, skipping\n');
            % Still load resultant_angle for flash RF
            pv_path = find_peak_vals(BAR_RESULTS_DIR, strain, on_off, date_str, time_str, MARCM_BAR_DIR);
            if ~isempty(pv_path)
                PV = load(pv_path, 'bar_results');
                resultant_angle = PV.bar_results.resultant_angle;
                m.resultant_angle_deg = rad2deg(resultant_angle);
                if isfield(PV.bar_results, 'slow')
                    m.bar_slow = extract_bar_speed_metrics(PV.bar_results, 'slow');
                end
                if isfield(PV.bar_results, 'fast')
                    m.bar_fast = extract_bar_speed_metrics(PV.bar_results, 'fast');
                end
                if isfield(PV.bar_results, 'vfast')
                    m.bar_vfast = extract_bar_speed_metrics(PV.bar_results, 'vfast');
                end
            end
        end

        catch ME_bar
            warning('  Bar analysis failed: %s', ME_bar.message);
            cell_errors{end+1} = sprintf('Bar analysis: %s', ME_bar.message);
            close all force
        end  % --- End Section 2 try/catch ---

        %% ====== Section 3: Flash RF (4px + 6px + Gaussian) ======
        try  % --- Section 3 try/catch (isolated from bar/bar_flash) ---

        flash_rf_4px_path     = fullfile(OUTPUT_DIR, 'flash_rf_4px', [cell_id '_flash_rf_4px.png']);
        flash_heatmap_4px_path = fullfile(OUTPUT_DIR, 'flash_heatmap_4px', [cell_id '_flash_heatmap_4px.png']);
        flash_rf_6px_path     = fullfile(OUTPUT_DIR, 'flash_rf_6px', [cell_id '_flash_rf_6px.png']);
        flash_heatmap_6px_path = fullfile(OUTPUT_DIR, 'flash_heatmap_6px', [cell_id '_flash_heatmap_6px.png']);
        gaussian_rf_path      = fullfile(OUTPUT_DIR, 'gaussian_rf', [cell_id '_gaussian_rf.png']);

        % Check if we have saved rf_results
        rf_path = find_rf_results(FLASH_RESULTS_DIR, cell_id, date_str, on_off);
        has_saved_rf = ~isempty(rf_path);

        % Determine idx filtering mode: Summer OFF cells need unfiltered idx
        % (process-summer-2025-data branch removes the [1,2,5,6,9,10] filter)
        use_autumn_idx_filter = ~is_summer;

        all_flash_exist = isfile(flash_rf_4px_path) && isfile(flash_heatmap_4px_path) && ...
                          isfile(flash_rf_6px_path) && isfile(flash_heatmap_6px_path) && ...
                          isfile(gaussian_rf_path);

        if ~all_flash_exist && ~isempty(exp_folder)
            fprintf('  Loading raw data for flash RF analysis...\n');

            % Load raw TDMS data
            [d_str, t_str, Log, params, ~] = load_protocol2_data(exp_folder);
            f_data = Log.ADC.Volts(1, :);
            v_data = Log.ADC.Volts(2, :) * 10;
            median_v = median(v_data);
            v2_data = v_data - median_v;
            params.date = d_str;
            params.time = t_str;
            params.strain = strain;
            params.resultant_angle = resultant_angle;

            % Compute flash start indices (Summer-aware)
            diff_f_data = diff(f_data);
            n_flashes_4px = 196;
            n_flashes_6px = 100;
            if on_off == "off"
                idx = find(diff_f_data == 1 & f_data(2:end) == 1);
                if use_autumn_idx_filter
                    idx = idx([1,2,5,6,9,10]);  % Autumn: filter to 6 rep-boundary transitions
                end
            elseif on_off == "on"
                % ON-cell logic is identical for Summer and Autumn
                idx_4 = find(diff_f_data == 1 + n_flashes_4px & f_data(2:end) == 1 + n_flashes_4px);
                idx_4(:, [2,4,6]) = [];
                idx_6 = find(diff_f_data == 1 + n_flashes_6px & f_data(2:end) == 1 + n_flashes_6px);
                idx = sort(horzcat(idx_4, idx_6));
            end

            % Voltage variance
            filtered_voltage_data = movmean(v_data, 20000);
            var_filtered_v = var(filtered_voltage_data);

            % Initialize rf_results struct for new saves
            rf_results = struct();
            rf_results.Date = d_str;
            rf_results.Time = t_str;
            rf_results.Strain = strain;
            rf_results.Type = on_off;

            % --- 4px analysis (always from raw data) ---
            if ~isfile(flash_rf_4px_path) || ~isfile(flash_heatmap_4px_path)
                fprintf('  Running 4px flash analysis from raw data...\n');
                try
                    [data_comb_4, cmap_id_4, var_across_4, var_within_4, diff_mean_4, max_data_4, min_data_4] = ...
                        parse_flash_data(f_data, v_data, on_off, "slow", 4, PROJECT_ROOT, use_autumn_idx_filter);
                    data_comb2_4 = rescale(data_comb_4, 0, 1);

                    % Timeseries
                    if ~isfile(flash_rf_4px_path)
                        f_ts = plot_rf_estimate_timeseries_line(data_comb2_4, cmap_id_4, ...
                            f_data, v2_data, "slow", 4, idx, params);
                        export_and_thumb(f_ts, flash_rf_4px_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'flash_rf_4px', [cell_id '_flash_rf_4px.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_ts);
                    end

                    % Heatmap
                    if ~isfile(flash_heatmap_4px_path)
                        f_hm = plot_heatmap_flash_responses(data_comb2_4);
                        export_and_thumb(f_hm, flash_heatmap_4px_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'flash_heatmap_4px', [cell_id '_flash_heatmap_4px.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_hm);
                    end

                    % Store 4px metrics
                    med_var_X_4 = median(reshape(var_across_4, [1, n_flashes_4px]));
                    med_var_W_4 = var(reshape(var_within_4, [1, n_flashes_4px]));
                    rf_results.px4.data_comb = {data_comb_4};
                    rf_results.px4.cmap_id = {cmap_id_4};
                    rf_results.px4.max_data = {max_data_4};
                    rf_results.px4.min_data = {min_data_4};

                    % Gaussian fit on 4px data
                    try
                        exc_data_4 = data_comb2_4;
                        inh_data_4 = data_comb_4;
                        inh_data_4(cmap_id_4 ~= 2) = 0;
                        [optEx4, R2_4, optInh4, R2i_4, f1_4, ~] = gaussian_RF_estimate(exc_data_4, inh_data_4);
                        close(f1_4);
                        close all force  % clean up Gaussian figures

                        rf_results.px4.R_squared = R2_4;
                        rf_results.px4.R_squaredi = R2i_4;
                        rf_results.px4.sigma_x_exc = optEx4(4);
                        rf_results.px4.sigma_y_exc = optEx4(5);
                        rf_results.px4.sigma_x_inh = optInh4(4);
                        rf_results.px4.sigma_y_inh = optInh4(5);
                        rf_results.px4.optExc = {optEx4};
                        rf_results.px4.optInh = {optInh4};

                        m.rf_4px = struct('R_squared', R2_4, 'R_squaredi', R2i_4, ...
                            'sigma_x_exc', optEx4(4), 'sigma_y_exc', optEx4(5), ...
                            'sigma_x_inh', optInh4(4), 'sigma_y_inh', optInh4(5));
                    catch ME_g4
                        warning('  4px Gaussian fit failed: %s', ME_g4.message);
                        cell_errors{end+1} = sprintf('4px Gaussian: %s', ME_g4.message);
                        close all force
                    end
                catch ME_4px
                    warning('  4px flash analysis failed: %s', ME_4px.message);
                    cell_errors{end+1} = sprintf('4px flash: %s', ME_4px.message);
                    close all force
                end
            end

            % --- 6px analysis ---
            if has_saved_rf
                % Load 6px data from saved rf_results
                fprintf('  Loading saved 6px rf_results: %s\n', rf_path);
                RF = load(rf_path, 'rf_results');
                rf_saved = RF.rf_results;

                % Generate 6px PNGs from saved data
                if ~isfile(flash_rf_6px_path) || ~isfile(flash_heatmap_6px_path)
                    data_comb_6 = rf_saved.slow.data_comb{1};
                    cmap_id_6 = rf_saved.slow.cmap_id{1};
                    data_comb2_6 = rescale(data_comb_6, 0, 1);

                    % Auto-detect px_size from data dimensions
                    grid_sz = size(data_comb_6, 1);
                    if grid_sz == 14
                        px_sz_saved = 4;
                    elseif grid_sz == 10
                        px_sz_saved = 6;
                    else
                        px_sz_saved = 6;  % default
                    end

                    if ~isfile(flash_rf_6px_path)
                        f_ts6 = plot_rf_estimate_timeseries_line(data_comb2_6, cmap_id_6, ...
                            f_data, v2_data, "slow", px_sz_saved, idx, params);
                        export_and_thumb(f_ts6, flash_rf_6px_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'flash_rf_6px', [cell_id '_flash_rf_6px.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_ts6);
                    end

                    if ~isfile(flash_heatmap_6px_path)
                        f_hm6 = plot_heatmap_flash_responses(data_comb2_6);
                        export_and_thumb(f_hm6, flash_heatmap_6px_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'flash_heatmap_6px', [cell_id '_flash_heatmap_6px.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_hm6);
                    end
                end

                % Gaussian contour from saved parameters
                if ~isfile(gaussian_rf_path)
                    try
                        data_comb_6 = rf_saved.slow.data_comb{1};
                        data_comb2_6 = rescale(data_comb_6, 0, 1);
                        cmap_id_6 = rf_saved.slow.cmap_id{1};
                        exc_data_6 = data_comb2_6;
                        inh_data_6 = data_comb_6;
                        inh_data_6(cmap_id_6 ~= 2) = 0;

                        if isfield(rf_saved.slow, 'optExc') && isfield(rf_saved.slow, 'optInh')
                            % Re-plot using saved parameters
                            optEx6 = rf_saved.slow.optExc{1};
                            optInh6 = rf_saved.slow.optInh{1};
                            f_gauss = plot_gaussian_contour(data_comb2_6, optEx6, optInh6);
                        else
                            % Run Gaussian fit
                            [optEx6, ~, optInh6, ~, f1_6, f_gauss] = gaussian_RF_estimate(exc_data_6, inh_data_6);
                            close(f1_6);
                        end
                        export_and_thumb(f_gauss, gaussian_rf_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'gaussian_rf', [cell_id '_gaussian_rf.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_gauss);
                    catch ME_g6
                        warning('  Gaussian contour failed: %s', ME_g6.message);
                        cell_errors{end+1} = sprintf('Gaussian contour: %s', ME_g6.message);
                    end
                end

                % Extract 6px metrics from saved rf_results
                m.rf_6px = struct( ...
                    'R_squared', rf_saved.slow.R_squared, ...
                    'R_squaredi', rf_saved.slow.R_squaredi, ...
                    'sigma_x_exc', rf_saved.slow.sigma_x_exc, ...
                    'sigma_y_exc', rf_saved.slow.sigma_y_exc, ...
                    'sigma_x_inh', rf_saved.slow.sigma_x_inh, ...
                    'sigma_y_inh', rf_saved.slow.sigma_y_inh);

                % Copy saved rf_results fields into our struct for saving
                rf_results.slow = rf_saved.slow;

                close all force
            else
                % No saved rf_results — run full 6px analysis from raw data
                fprintf('  Running 6px flash analysis from raw data...\n');
                try
                    [data_comb_6, cmap_id_6, var_across_6, var_within_6, diff_mean_6, max_data_6, min_data_6] = ...
                        parse_flash_data(f_data, v_data, on_off, "slow", 6, PROJECT_ROOT, use_autumn_idx_filter);
                    data_comb2_6 = rescale(data_comb_6, 0, 1);

                    if ~isfile(flash_rf_6px_path)
                        f_ts6 = plot_rf_estimate_timeseries_line(data_comb2_6, cmap_id_6, ...
                            f_data, v2_data, "slow", 6, idx, params);
                        export_and_thumb(f_ts6, flash_rf_6px_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'flash_rf_6px', [cell_id '_flash_rf_6px.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_ts6);
                    end

                    if ~isfile(flash_heatmap_6px_path)
                        f_hm6 = plot_heatmap_flash_responses(data_comb2_6);
                        export_and_thumb(f_hm6, flash_heatmap_6px_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'flash_heatmap_6px', [cell_id '_flash_heatmap_6px.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_hm6);
                    end

                    % Gaussian fit on 6px data
                    if ~isfile(gaussian_rf_path)
                        exc_data_6 = data_comb2_6;
                        inh_data_6 = data_comb_6;
                        inh_data_6(cmap_id_6 ~= 2) = 0;
                        [optEx6, R2_6, optInh6, R2i_6, f1_6, f_gauss] = gaussian_RF_estimate(exc_data_6, inh_data_6);
                        close(f1_6);
                        export_and_thumb(f_gauss, gaussian_rf_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'gaussian_rf', [cell_id '_gaussian_rf.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(f_gauss);
                    else
                        % Still run Gaussian for metrics
                        exc_data_6 = data_comb2_6;
                        inh_data_6 = data_comb_6;
                        inh_data_6(cmap_id_6 ~= 2) = 0;
                        [optEx6, R2_6, optInh6, R2i_6, ~, ~] = gaussian_RF_estimate(exc_data_6, inh_data_6);
                        close all force
                    end

                    % Store 6px results
                    med_var_X_6 = median(reshape(var_across_6, [1, n_flashes_6px]));
                    med_var_W_6 = var(reshape(var_within_6, [1, n_flashes_6px]));
                    rf_results.slow.data_comb = {data_comb_6};
                    rf_results.slow.cmap_id = {cmap_id_6};
                    rf_results.slow.max_data = {max_data_6};
                    rf_results.slow.min_data = {min_data_6};
                    rf_results.slow.diff_mean = {diff_mean_6};
                    rf_results.slow.var_within_reps = {var_within_6};
                    rf_results.slow.var_across_reps = {var_across_6};
                    rf_results.slow.var_filtered_v = var_filtered_v;
                    rf_results.slow.med_var_X_reps = med_var_X_6;
                    rf_results.slow.med_var_W_reps = med_var_W_6;
                    rf_results.slow.R_squared = R2_6;
                    rf_results.slow.sigma_x_exc = optEx6(4);
                    rf_results.slow.sigma_y_exc = optEx6(5);
                    rf_results.slow.optExc = {optEx6};
                    rf_results.slow.R_squaredi = R2i_6;
                    rf_results.slow.optInh = {optInh6};
                    rf_results.slow.sigma_x_inh = optInh6(4);
                    rf_results.slow.sigma_y_inh = optInh6(5);

                    m.rf_6px = struct('R_squared', R2_6, 'R_squaredi', R2i_6, ...
                        'sigma_x_exc', optEx6(4), 'sigma_y_exc', optEx6(5), ...
                        'sigma_x_inh', optInh6(4), 'sigma_y_inh', optInh6(5));
                catch ME_6px
                    warning('  6px flash analysis failed: %s', ME_6px.message);
                    cell_errors{end+1} = sprintf('6px flash: %s', ME_6px.message);
                end
                close all force
            end

            % Save rf_results to genotype/on_off folder (only if we ran new analysis)
            if ~has_saved_rf && isfield(rf_results, 'slow')
                rf_save_path = fullfile(results_folder, ...
                    sprintf('rf_results_%s_%s_%s_%s.mat', d_str, t_str, strain, on_off));
                save(rf_save_path, 'rf_results');
                fprintf('  Saved rf_results: %s\n', rf_save_path);
            end

        elseif all_flash_exist
            fprintf('  Flash RF plots already exist, skipping\n');
            % Load metrics from saved rf_results if available
            rf_path = find_rf_results(FLASH_RESULTS_DIR, cell_id, date_str, on_off);
            if ~isempty(rf_path)
                RF = load(rf_path, 'rf_results');
                rf_saved = RF.rf_results;
                m.rf_6px = struct( ...
                    'R_squared', rf_saved.slow.R_squared, ...
                    'R_squaredi', rf_saved.slow.R_squaredi, ...
                    'sigma_x_exc', rf_saved.slow.sigma_x_exc, ...
                    'sigma_y_exc', rf_saved.slow.sigma_y_exc, ...
                    'sigma_x_inh', rf_saved.slow.sigma_x_inh, ...
                    'sigma_y_inh', rf_saved.slow.sigma_y_inh);
            end
        else
            fprintf('  No experiment folder for flash RF\n');
            cell_errors{end+1} = 'Flash RF: no exp_folder and no saved results';
        end

        catch ME_flash
            warning('  Flash RF analysis failed: %s', ME_flash.message);
            cell_errors{end+1} = sprintf('Flash RF: %s', ME_flash.message);
            close all force
        end  % --- End Section 3 try/catch ---

        %% ====== Section 4: Bar Flash (Autumn only) ======
        try  % --- Section 4 try/catch (isolated from bar/flash) ---

        bar_flash_slow_path = fullfile(OUTPUT_DIR, 'bar_flash_slow', [cell_id '_bar_flash_slow.png']);
        bar_flash_fast_path = fullfile(OUTPUT_DIR, 'bar_flash_fast', [cell_id '_bar_flash_fast.png']);

        if is_summer
            fprintf('  Skipping bar flash (Summer — stimulus not used)\n');
        elseif isfile(bar_flash_slow_path) && isfile(bar_flash_fast_path)
            fprintf('  Bar flash plots already exist, skipping\n');
        else
            % Search for saved bar_flash_results
            bf_path = find_bar_flash(BAR_FLASH_DIR, BAR_FLASH_1DRF, cell_id, strain, on_off, date_str, time_str);

            if ~isempty(bf_path)
                fprintf('  Loading saved bar flash results: %s\n', bf_path);
                BF = load(bf_path, 'data_slow', 'data_fast', 'mean_slow', 'mean_fast');

                if ~isfile(bar_flash_slow_path)
                    n_int_slow = 7500;  % 0.75 * 10000
                    fig_slow = plot_bar_flash_data(BF.data_slow, BF.mean_slow, n_int_slow);
                    fig_slow.Color = 'w';
                    export_and_thumb(fig_slow, bar_flash_slow_path, ...
                        fullfile(OUTPUT_DIR, 'thumbs', 'bar_flash_slow', [cell_id '_bar_flash_slow.png']), ...
                        EXPORT_OPTS, THUMB_WIDTH);
                    close(fig_slow);
                end

                if ~isfile(bar_flash_fast_path)
                    n_int_fast = 3750;  % 0.75 * 5000
                    fig_fast = plot_bar_flash_data(BF.data_fast, BF.mean_fast, n_int_fast);
                    fig_fast.Color = 'w';
                    export_and_thumb(fig_fast, bar_flash_fast_path, ...
                        fullfile(OUTPUT_DIR, 'thumbs', 'bar_flash_fast', [cell_id '_bar_flash_fast.png']), ...
                        EXPORT_OPTS, THUMB_WIDTH);
                    close(fig_fast);
                end
                close all force
            else
                fprintf('  No saved bar flash results — running from raw data\n');
                if ~isempty(exp_folder)
                    process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT);
                    close all force

                    % Re-find and generate PNGs
                    bf_path = find_bar_flash(BAR_FLASH_DIR, BAR_FLASH_1DRF, cell_id, strain, on_off, date_str, time_str);
                    if ~isempty(bf_path)
                        BF = load(bf_path, 'data_slow', 'data_fast', 'mean_slow', 'mean_fast');
                        n_int_slow = 7500;
                        fig_slow = plot_bar_flash_data(BF.data_slow, BF.mean_slow, n_int_slow);
                        fig_slow.Color = 'w';
                        export_and_thumb(fig_slow, bar_flash_slow_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'bar_flash_slow', [cell_id '_bar_flash_slow.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(fig_slow);

                        n_int_fast = 3750;
                        fig_fast = plot_bar_flash_data(BF.data_fast, BF.mean_fast, n_int_fast);
                        fig_fast.Color = 'w';
                        export_and_thumb(fig_fast, bar_flash_fast_path, ...
                            fullfile(OUTPUT_DIR, 'thumbs', 'bar_flash_fast', [cell_id '_bar_flash_fast.png']), ...
                            EXPORT_OPTS, THUMB_WIDTH);
                        close(fig_fast);
                    end
                    close all force
                end
            end
        end

        catch ME_bf
            warning('  Bar flash analysis failed: %s', ME_bf.message);
            cell_errors{end+1} = sprintf('Bar flash: %s', ME_bf.message);
            close all force
        end  % --- End Section 4 try/catch ---

        % Store metrics and errors
        cell_metrics{i} = m;
        cell_errors_all{i} = cell_errors;
        n_success = n_success + 1;
        fprintf('  Done.\n');

    catch ME
        fprintf('  FAILED: %s\n', ME.message);
        fprintf(fid_fail, '%s\t%s\t%s\n', cell_id, exp_folder, strrep(ME.message, newline, ' '));
        cell_errors{end+1} = sprintf('FATAL: %s', ME.message);
        n_fail = n_fail + 1;
        cell_metrics{i} = m;
        cell_errors_all{i} = cell_errors;
        cd(orig_dir);
    end

    close all force
end

fclose(fid_fail);

%% Write diagnostic report
DIAG_PATH = fullfile(OUTPUT_DIR, 'diagnostic_report.txt');
fprintf('\nWriting diagnostic report...\n');

% Column definitions: {display_name, subdirectory, filename_suffix}
plot_cols = {
    'bar_polar',        'bar_polar',        'bar_polar'
    'polar_arrow',      'polar_arrow',      'polar_arrow'
    'flash_rf_4px',     'flash_rf_4px',     'flash_rf_4px'
    'flash_heatmap_4px','flash_heatmap_4px','flash_heatmap_4px'
    'flash_rf_6px',     'flash_rf_6px',     'flash_rf_6px'
    'flash_heatmap_6px','flash_heatmap_6px','flash_heatmap_6px'
    'gaussian_rf',      'gaussian_rf',      'gaussian_rf'
    'bar_flash_slow',   'bar_flash_slow',   'bar_flash_slow'
    'bar_flash_fast',   'bar_flash_fast',   'bar_flash_fast'
    'gridplot_1',       'gridplots',        'gridplot_1'
    'gridplot_2',       'gridplots',        'gridplot_2'
    'gridplot_3',       'gridplots',        'gridplot_3'
    'gridplot_4',       'gridplots',        'gridplot_4'
};
n_plot_types = size(plot_cols, 1);

fid_diag = fopen(DIAG_PATH, 'w');

% Header row
fprintf(fid_diag, 'cell_id\tstrain\tcontrast\tseason\texp_folder');
for k = 1:n_plot_types
    fprintf(fid_diag, '\t%s', plot_cols{k, 1});
end
fprintf(fid_diag, '\ttotal_plots\terrors\n');

% Track per-plot-type totals for console summary
plot_totals = zeros(1, n_plot_types);

% One row per cell
for i = 1:n_cells
    d_str = safe_str(cells.date(i));
    t_str = safe_str(cells.time(i));
    cid = sprintf('%s_%s', d_str, t_str);
    strain_i = safe_str(cells.strain(i));
    pc_i = safe_str(cells.preferred_contrast(i));
    pt_i = safe_str(cells.project_type(i));
    if contains(string(pt_i), 'Summer2025')
        season = 'Summer';
    else
        season = 'Autumn';
    end

    % Re-find experiment folder
    ef = find_experiment_folder({P2_SUMMER, P2_AUTUMN}, cid, d_str, t_str);

    fprintf(fid_diag, '%s\t%s\t%s\t%s\t%s', cid, strain_i, pc_i, season, ef);

    total = 0;
    for k = 1:n_plot_types
        fname = sprintf('%s_%s.png', cid, plot_cols{k, 3});
        exists = isfile(fullfile(OUTPUT_DIR, plot_cols{k, 2}, fname));
        fprintf(fid_diag, '\t%d', exists);
        total = total + exists;
        plot_totals(k) = plot_totals(k) + exists;
    end

    % Error column
    if i <= numel(cell_errors_all) && ~isempty(cell_errors_all{i})
        err_str = strjoin(cell_errors_all{i}, '; ');
        err_str = strrep(err_str, sprintf('\t'), ' ');
        err_str = strrep(err_str, newline, ' ');
    else
        err_str = '';
    end

    fprintf(fid_diag, '\t%d\t%s\n', total, err_str);
end

fclose(fid_diag);
fprintf('Diagnostic report: %s\n', DIAG_PATH);

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('  Processed: %d\n', n_success);
fprintf('  Failed:    %d\n', n_fail);
fprintf('  Total:     %d\n', n_cells);
if n_fail > 0
    fprintf('  See: %s\n', FAILURE_LOG);
end

fprintf('\n=== Plot Completeness ===\n');
for k = 1:n_plot_types
    fprintf('  %-20s: %d/%d\n', plot_cols{k, 1}, plot_totals(k), n_cells);
end

%% Build cell_index.json with metrics
fprintf('\nBuilding cell index with metrics...\n');
build_cell_index_with_metrics(cells, OUTPUT_DIR, cell_metrics);

fprintf('\nDone! Dashboard images saved to:\n  %s\n', OUTPUT_DIR);

end


%% ========================================================================
%  Helper Functions
%  ========================================================================

function export_and_thumb(fig, png_path, thumb_path, export_opts, thumb_width)
% EXPORT_AND_THUMB  Export a figure as PNG and create a thumbnail.
    fig.Color = 'w';
    exportgraphics(fig, png_path, export_opts{:});
    make_thumbnail(png_path, thumb_path, thumb_width);
end


function make_thumbnail(src_path, thumb_path, max_width)
% MAKE_THUMBNAIL  Create a resized thumbnail of a PNG image.
    try
        img = imread(src_path);
        [h, w, ~] = size(img);
        if w > max_width
            scale = max_width / w;
            new_h = round(h * scale);
            thumb = imresize(img, [new_h, max_width]);
            imwrite(thumb, thumb_path);
        else
            copyfile(src_path, thumb_path);
        end
    catch
        copyfile(src_path, thumb_path);
    end
end


function found_path = find_peak_vals(bar_results_dir, strain, on_off, date_str, time_str, marcm_bar_dir)
% FIND_PEAK_VALS  Search for peak_vals .mat file recursively (fuzzy ±2 min).
%   Searches bar_results_dir first, then marcm_bar_dir as fallback.
    if nargin < 6
        marcm_bar_dir = '';
    end

    found_path = '';
    time_variants = fuzzy_time_variants(time_str);

    % Build list of directories to search (primary first, then MARCM fallback)
    search_dirs = {bar_results_dir};
    if ~isempty(marcm_bar_dir) && isfolder(marcm_bar_dir)
        search_dirs{end+1} = marcm_bar_dir;
    end

    for sd = 1:numel(search_dirs)
        search_root = search_dirs{sd};

        for tv = 1:numel(time_variants)
            ts = time_variants{tv};

            % Exact strain pattern
            pattern = fullfile(search_root, '**', ...
                sprintf('peak_vals_%s_%s_%s_%s.mat', strain, on_off, date_str, ts));
            files = dir(pattern);
            if ~isempty(files)
                found_path = fullfile(files(1).folder, files(1).name);
                return;
            end

            % Broader pattern (strain may differ in case or format)
            pattern2 = fullfile(search_root, '**', ...
                sprintf('peak_vals_*_%s_%s_%s.mat', on_off, date_str, ts));
            files = dir(pattern2);
            if ~isempty(files)
                found_path = fullfile(files(1).folder, files(1).name);
                return;
            end
        end
    end
end


function found_path = find_rf_results(flash_results_dir, cell_id, date_str, on_off)
% FIND_RF_RESULTS  Search for rf_results .mat file (fuzzy ±2 min).
    found_path = '';

    % Parse time from cell_id (format: YYYY_MM_DD_HH_MM)
    parts = strsplit(cell_id, '_');
    if numel(parts) >= 5
        time_str = sprintf('%s_%s', parts{4}, parts{5});
        time_variants = fuzzy_time_variants(time_str);
    else
        time_variants = {};
    end

    % Try cell_id pattern (exact)
    pattern = fullfile(flash_results_dir, sprintf('rf_results_%s_*_%s.mat', cell_id, on_off));
    files = dir(pattern);
    if ~isempty(files)
        found_path = fullfile(flash_results_dir, files(1).name);
        return;
    end

    % Broader date-only pattern
    pattern2 = fullfile(flash_results_dir, sprintf('rf_results_%s*%s*.mat', date_str, on_off));
    files = dir(pattern2);
    if ~isempty(files)
        found_path = fullfile(flash_results_dir, files(1).name);
        return;
    end

    % Fuzzy: try alternative time offsets
    for tv = 2:numel(time_variants)  % skip tv=1, already tried as exact
        fuzzy_cell_id = sprintf('%s_%s', date_str, time_variants{tv});
        pattern3 = fullfile(flash_results_dir, sprintf('rf_results_%s_*_%s.mat', fuzzy_cell_id, on_off));
        files = dir(pattern3);
        if ~isempty(files)
            found_path = fullfile(flash_results_dir, files(1).name);
            return;
        end
    end
end


function found_path = find_bar_flash(primary_dir, alt_dir, cell_id, strain, on_off, date_str, time_str)
% FIND_BAR_FLASH  Search for bar_flash_results .mat file (fuzzy ±2 min).
    found_path = '';
    time_variants = fuzzy_time_variants(time_str);

    % Try each time variant with exact-strain patterns first
    for tv = 1:numel(time_variants)
        ts = time_variants{tv};
        fuzzy_cell_id = sprintf('%s_%s', date_str, ts);

        % Search primary directory (exact pattern)
        pattern1 = fullfile(primary_dir, '**', ...
            sprintf('bar_flash_results_%s_%s_%s_%s.mat', date_str, ts, strain, on_off));
        files = dir(pattern1);
        if ~isempty(files)
            found_path = fullfile(files(1).folder, files(1).name);
            return;
        end

        % Search alternative directory (exact pattern)
        if isfolder(alt_dir)
            pattern2 = fullfile(alt_dir, ...
                sprintf('bar_flash_results_%s_%s_%s.mat', fuzzy_cell_id, strain, on_off));
            files = dir(pattern2);
            if ~isempty(files)
                found_path = fullfile(alt_dir, files(1).name);
                return;
            end
        end
    end

    % Broad fallback: date + on_off only (no specific time or strain)
    pattern_broad = fullfile(primary_dir, '**', ...
        sprintf('bar_flash_results_%s*%s*.mat', date_str, on_off));
    files = dir(pattern_broad);
    if ~isempty(files)
        found_path = fullfile(files(1).folder, files(1).name);
        return;
    end

    if isfolder(alt_dir)
        pattern_broad2 = fullfile(alt_dir, ...
            sprintf('bar_flash_results_%s*%s*.mat', date_str, on_off));
        files = dir(pattern_broad2);
        if ~isempty(files)
            found_path = fullfile(alt_dir, files(1).name);
        end
    end
end


function sm = extract_bar_speed_metrics(bar_results, speed_name)
% EXTRACT_BAR_SPEED_METRICS  Pull metrics for one speed from bar_results struct.
    sm = struct();
    if isfield(bar_results, speed_name)
        s = bar_results.(speed_name);
        if isfield(s, 'magnitude'), sm.magnitude = s.magnitude; end
        if isfield(s, 'DSI_vector'), sm.DSI_vector = s.DSI_vector; end
        if isfield(s, 'DSI_pdnd'), sm.DSI_pdnd = s.DSI_pdnd; end
        if isfield(s, 'fwhm'), sm.fwhm = s.fwhm; end
        if isfield(s, 'cv'), sm.cv = s.cv; end
        if isfield(s, 'angle_rad'), sm.angle_rad = s.angle_rad; end
    end
end


function f = plot_gaussian_contour(response, optEx, optInh)
% PLOT_GAUSSIAN_CONTOUR  Re-create Gaussian RF contour plot from saved parameters.
%   Used when optExc/optInh are available from saved rf_results.
    theta_vals = linspace(0, 2*pi, 100);

    R_exc = [cos(optEx(6)), -sin(optEx(6)); sin(optEx(6)), cos(optEx(6))];
    R_inh = [cos(optInh(6)), -sin(optInh(6)); sin(optInh(6)), cos(optInh(6))];

    exc_ellipse = R_exc * [1.5 * optEx(4) * cos(theta_vals); 1.5 * optEx(5) * sin(theta_vals)];
    exc_x = optEx(2) + exc_ellipse(1, :);
    exc_y = optEx(3) + exc_ellipse(2, :);

    inh_ellipse = R_inh * [1.5 * optInh(4) * cos(theta_vals); 1.5 * optInh(5) * sin(theta_vals)];
    inh_x = optInh(2) + inh_ellipse(1, :);
    inh_y = optInh(3) + inh_ellipse(2, :);

    figure;
    imagesc(response); hold on;
    colormap redblue; colorbar;
    title('Receptive Field with 1.5\sigma Contours');

    med_val = median(response(:));
    max_val = prctile(response(:), 98);
    clim([med_val - max_val, med_val + max_val]);

    plot(exc_x, exc_y, 'r', 'LineWidth', 2);
    plot(inh_x, inh_y, 'k', 'LineWidth', 2);
    axis square;

    f = gcf;
    f.Position = [712, 576, 560, 420];
end


function result = ternary(condition, true_val, false_val)
% TERNARY  Inline conditional.
    if condition
        result = true_val;
    else
        result = false_val;
    end
end


function exp_folder = find_experiment_folder(base_dirs, cell_id, date_str, time_str)
% FIND_EXPERIMENT_FOLDER  Find experiment folder with fuzzy ±2 min matching.
%   Tries exact cell_id match first across all base directories, then
%   searches with ±1 and ±2 minute offsets on the time component.
    exp_folder = '';
    time_variants = fuzzy_time_variants(time_str);

    for tv = 1:numel(time_variants)
        fuzzy_id = sprintf('%s_%s', date_str, time_variants{tv});
        for d = 1:numel(base_dirs)
            candidate = fullfile(base_dirs{d}, fuzzy_id);
            if isfolder(candidate)
                exp_folder = candidate;
                return;
            end
        end
    end
end


function variants = fuzzy_time_variants(time_str)
% FUZZY_TIME_VARIANTS  Return cell array of time strings: exact, then ±1, ±2 min.
%   Given time_str in 'HH_MM' format, returns up to 5 variants with hour rollover.
    variants = {time_str};
    parts = strsplit(time_str, '_');
    if numel(parts) < 2, return; end
    hh = str2double(parts{1});
    mm = str2double(parts{2});
    if isnan(hh) || isnan(mm), return; end

    for offset = [-1, 1, -2, 2]
        new_mm = mm + offset;
        new_hh = hh;
        if new_mm < 0
            new_mm = new_mm + 60;
            new_hh = new_hh - 1;
        end
        if new_mm >= 60
            new_mm = new_mm - 60;
            new_hh = new_hh + 1;
        end
        variants{end+1} = sprintf('%02d_%02d', new_hh, new_mm); %#ok<AGROW>
    end
end


function tf = is_valid_str(val)
% IS_VALID_STR  Return true if val is a non-empty, non-missing string/char.
%   Handles cell arrays, string arrays, char, categorical, and NaN gracefully.
    if iscell(val)
        val = val{1};
    end
    if isnumeric(val)
        tf = false;
        return;
    end
    if isstring(val)
        tf = strlength(val) > 0 && ~ismissing(val);
    elseif iscategorical(val)
        tf = ~isundefined(val);
    elseif ischar(val)
        tf = ~isempty(val);
    else
        tf = false;
    end
end


function val = safe_str(tbl_cell)
% SAFE_STR  Extract a char/string from a table cell, returning '' if missing.
    if iscell(tbl_cell)
        val = tbl_cell{1};
    else
        val = tbl_cell;
    end
    if isnumeric(val) && isnan(val)
        val = '';
        return;
    end
    if isstring(val)
        if ismissing(val)
            val = '';
        else
            val = char(val);
        end
    elseif iscategorical(val)
        if isundefined(val)
            val = '';
        else
            val = char(val);
        end
    elseif ~ischar(val)
        val = '';
    end
end


function build_cell_index_with_metrics(cells_table, output_dir, cell_metrics)
% BUILD_CELL_INDEX_WITH_METRICS  Write cell_index.json with image paths and metrics.

n = height(cells_table);
cells_json = cell(n, 1);

% Image type definitions: {key, subdirectory}
image_types = {
    'bar_polar',        'bar_polar'
    'polar_arrow',      'polar_arrow'
    'flash_rf_4px',     'flash_rf_4px'
    'flash_heatmap_4px','flash_heatmap_4px'
    'flash_rf_6px',     'flash_rf_6px'
    'flash_heatmap_6px','flash_heatmap_6px'
    'gaussian_rf',      'gaussian_rf'
    'bar_flash_slow',   'bar_flash_slow'
    'bar_flash_fast',   'bar_flash_fast'
};

for i = 1:n
    row = cells_table(i, :);
    d_str = safe_str(row.date(1));
    t_str = safe_str(row.time(1));
    cell_id = sprintf('%s_%s', d_str, t_str);

    c = struct();
    c.cell_id = cell_id;

    % Parse date for display
    date_parts = strsplit(d_str, '_');
    c.date = sprintf('%s-%s-%s', date_parts{1}, date_parts{2}, date_parts{3});
    time_parts = strsplit(t_str, '_');
    c.time = sprintf('%s:%s', time_parts{1}, time_parts{2});

    c.strain = safe_str(row.strain(1));

    if isfield(row, 'genotype') && is_valid_str(row.genotype(1))
        c.genotype = safe_str(row.genotype(1));
    else
        c.genotype = '';
    end

    c.preferred_contrast = safe_str(row.preferred_contrast(1));

    if isfield(row, 'project_type') && is_valid_str(row.project_type(1))
        c.project_type = safe_str(row.project_type(1));
    else
        c.project_type = '';
    end

    c.age = row.age(1);

    % Sweep / bar flash inclusion flags
    if isfield(row, 'Incl_SweepAnalysis_') && is_valid_str(row.Incl_SweepAnalysis_(1))
        c.incl_sweep = safe_str(row.Incl_SweepAnalysis_(1));
    else
        c.incl_sweep = '';
    end
    if isfield(row, 'Incl_BarFlashAnalysis_') && is_valid_str(row.Incl_BarFlashAnalysis_(1))
        c.incl_bar_flash = safe_str(row.Incl_BarFlashAnalysis_(1));
    else
        c.incl_bar_flash = '';
    end

    % Notes
    if isfield(row, 'NotesStart') && is_valid_str(row.NotesStart(1))
        c.notes_start = safe_str(row.NotesStart(1));
    else
        c.notes_start = '';
    end
    if isfield(row, 'NotesEnd') && is_valid_str(row.NotesEnd(1))
        c.notes_end = safe_str(row.NotesEnd(1));
    else
        c.notes_end = '';
    end

    % Gridplot images (4 per cell)
    for gi = 1:4
        gp_fname = sprintf('%s_gridplot_%d.png', cell_id, gi);
        c.(sprintf('has_gridplot_%d', gi)) = isfile(fullfile(output_dir, 'gridplots', gp_fname));
        c.(sprintf('gridplot_%d_path', gi)) = ['gridplots/' gp_fname];
        c.(sprintf('gridplot_%d_thumb', gi)) = ['thumbs/gridplots/' gp_fname];
    end

    % Analysis image types
    for k = 1:size(image_types, 1)
        key = image_types{k, 1};
        subdir = image_types{k, 2};
        fname = sprintf('%s_%s.png', cell_id, key);
        c.(sprintf('has_%s', key)) = isfile(fullfile(output_dir, subdir, fname));
        c.(sprintf('%s_path', key)) = [subdir '/' fname];
        c.(sprintf('%s_thumb', key)) = ['thumbs/' subdir '/' fname];
    end

    % Display label
    c.display_label = sprintf('%s-%s-%s-%s', c.strain, c.preferred_contrast, d_str, t_str);

    % Metrics
    if i <= numel(cell_metrics) && ~isempty(cell_metrics{i})
        c.metrics = cell_metrics{i};
    else
        c.metrics = struct();
    end

    cells_json{i} = c;
end

% Build JSON structure
index = struct();
index.generated_at = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
index.images_dir = output_dir;
index.cells = cells_json;

json_str = jsonencode(index, 'PrettyPrint', true);
json_path = fullfile(output_dir, 'cell_index.json');
fid = fopen(json_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

fprintf('Cell index written to %s (%d cells)\n', json_path, n);

end
