function results = analyze_voltage_histograms(data_root, opts)
% ANALYZE_VOLTAGE_HISTOGRAMS  Voltage distribution analysis for 28 dps bar sweeps.
%
%   RESULTS = ANALYZE_VOLTAGE_HISTOGRAMS(DATA_ROOT) processes all 1DRF
%   experiments, extracts the contiguous voltage block during 28 dps
%   (slowest speed) bar sweeps, and generates voltage histograms comparing
%   control vs TTL conditions grouped by ON/OFF cell type.
%
%   RESULTS = ANALYZE_VOLTAGE_HISTOGRAMS(DATA_ROOT, OPTS) uses options to
%   override defaults.
%
%   INPUTS:
%     data_root - Path to 1DRF data directory containing experiment folders
%     opts      - (Optional) structure with fields:
%                   .on_threshold - Frame threshold for ON/OFF (default: 129)
%                   .bin_width    - Histogram bin width in mV (default: 0.1)
%                   .save_figs    - Save figures as PDF (default: true)
%                   .save_dir     - Output directory for population PDFs
%                                   (default: <data_root>/population_results)
%
%   OUTPUT:
%     results - Structure array with per-cell fields:
%       .folder        - Experiment folder name
%       .group         - 'on_control', 'on_ttl', 'off_control', or 'off_ttl'
%       .strain        - Strain name from metadata
%       .frame         - Frame number from metadata
%       .is_on         - Logical, true if ON cell
%       .is_ttl        - Logical, true if TTL condition
%       .v_28dps       - Concatenated voltage during 28 dps epochs (3 reps)
%       .boundaries    - 3x2 matrix of [start, end] sample indices per rep
%
%   FIGURES GENERATED:
%     Per cell: voltage histogram + verification plot (frame data with
%     extraction boundaries). Population: ON and OFF histograms with
%     control (black/grey) vs TTL (red/light-red) overlay.
%
%   See also BATCH_ANALYZE_1DRF, PARSE_BAR_DATA, LOAD_PROTOCOL2_DATA

    if nargin < 2, opts = struct(); end
    opts = set_defaults(opts, data_root);

    %% Discover experiment folders
    d = dir(data_root);
    d = d([d.isdir]);
    d = d(~startsWith({d.name}, '.'));

    valid = false(numel(d), 1);
    for i = 1:numel(d)
        valid(i) = isfile(fullfile(data_root, d(i).name, 'currentExp.mat'));
    end
    d = d(valid);
    n_exp = numel(d);
    fprintf('Found %d experiment folders in %s\n', n_exp, data_root);

    %% Process each experiment and generate verification plots
    results = struct([]);

    for exp_idx = 1:n_exp
        folder = d(exp_idx).name;
        exp_folder = fullfile(data_root, folder);
        fprintf('\n[%d/%d] Processing %s...\n', exp_idx, n_exp, folder);

        try
            [r, f_data_full] = process_single_cell(exp_folder, opts);
            r.folder = folder;

            % Generate verification plot immediately (no need to store f_data)
            fig_verif = plot_verification(f_data_full, r.boundaries, ...
                sprintf('Verification: %s', folder));
            if opts.save_figs
                exportgraphics(fig_verif, ...
                    fullfile(exp_folder, 'verification_28dps.pdf'), ...
                    'ContentType', 'vector', 'BackgroundColor', 'none');
                close(fig_verif);
            end

            if isempty(results)
                results = r;
            else
                results(end + 1) = r; %#ok<AGROW>
            end

            fprintf('  -> %s | n_samples: %d\n', r.group, numel(r.v_28dps));

        catch ME
            fprintf('  ERROR: %s\n', ME.message);
            continue;
        end
    end

    fprintf('\n=== Processing Complete ===\n');
    fprintf('Total cells: %d\n', numel(results));

    groups = {results.group};
    for g = ["on_control", "on_ttl", "off_control", "off_ttl"]
        fprintf('  %s: %d\n', g, sum(strcmp(groups, g)));
    end

    %% Determine shared bin edges from global voltage range
    all_v = cellfun(@(x) [min(x), max(x)], {results.v_28dps}, 'UniformOutput', false);
    all_v = vertcat(all_v{:});
    global_min = floor(min(all_v(:, 1)) * 10) / 10;  % Round down to 0.1 mV
    global_max = ceil(max(all_v(:, 2)) * 10) / 10;    % Round up to 0.1 mV
    bin_edges = global_min:opts.bin_width:global_max;
    fprintf('\nHistogram range: [%.1f, %.1f] mV, %d bins\n', ...
        global_min, global_max, numel(bin_edges) - 1);

    %% Generate per-cell histogram figures
    fprintf('\nGenerating per-cell histograms...\n');
    for i = 1:numel(results)
        r = results(i);
        exp_folder = fullfile(data_root, r.folder);

        fig_hist = plot_individual_histogram(r.v_28dps, r.group, ...
            sprintf('%s — %s', r.folder, r.group), bin_edges);

        if opts.save_figs
            exportgraphics(fig_hist, ...
                fullfile(exp_folder, 'voltage_histogram_28dps.pdf'), ...
                'ContentType', 'vector', 'BackgroundColor', 'none');
            close(fig_hist);
        end
    end

    %% Generate population figures
    fprintf('\nGenerating population figures...\n');

    if ~isfolder(opts.save_dir)
        mkdir(opts.save_dir);
    end

    for on_off_label = ["ON", "OFF"]
        if on_off_label == "ON"
            mask = [results.is_on];
        else
            mask = ~[results.is_on];
        end

        ctrl_mask = mask & ~[results.is_ttl];
        ttl_mask  = mask & [results.is_ttl];

        v_ctrl = {results(ctrl_mask).v_28dps};
        v_ttl  = {results(ttl_mask).v_28dps};

        % Population histogram: median +/- MAD of per-cell PDFs
        pop_title = sprintf('%s Cells — 28 dps Voltage Distribution', on_off_label);
        fig_pop = plot_population_histogram(v_ctrl, v_ttl, bin_edges, pop_title);

        if opts.save_figs
            fname = sprintf('%s_voltage_histogram_ctrl_vs_ttl.pdf', lower(char(on_off_label)));
            exportgraphics(fig_pop, fullfile(opts.save_dir, fname), ...
                'ContentType', 'vector', 'BackgroundColor', 'none');
            fprintf('  Saved: %s\n', fname);
        end
    end

    fprintf('\nDone.\n');

end


%% ========================= Core Processing ============================

function [r, f_data] = process_single_cell(exp_folder, opts)
% PROCESS_SINGLE_CELL  Load data, classify cell, extract 28 dps voltage.
%   Returns f_data separately for verification plots (not stored in struct).

    % Save/restore working directory (load_protocol2_data uses cd)
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);

    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    ce = load(fullfile(exp_folder, 'currentExp.mat'), 'metadata');
    metadata = ce.metadata;

    % Classify cell
    r.strain = metadata.Strain;
    r.frame  = metadata.Frame;
    r.is_on  = metadata.Frame > opts.on_threshold;
    r.is_ttl = contains(metadata.Strain, 'ttl');

    if r.is_on && ~r.is_ttl
        r.group = 'on_control';
    elseif r.is_on && r.is_ttl
        r.group = 'on_ttl';
    elseif ~r.is_on && ~r.is_ttl
        r.group = 'off_control';
    else
        r.group = 'off_ttl';
    end

    % Extract contiguous 28 dps voltage
    [r.v_28dps, r.boundaries] = extract_28dps_voltage(f_data, v_data);

end


function [v_28, boundaries] = extract_28dps_voltage(f_data, v_data)
% EXTRACT_28DPS_VOLTAGE  Extract contiguous voltage during 28 dps bar sweeps.
%
%   Finds the 28 dps stimulus blocks across all 3 repetitions by detecting
%   grey screen intervals (>=30,000 samples of frame==0) in the frame data.
%
%   Extraction window per rep:
%     Start: first sample after the pre-28dps grey gap ends (first sweep onset)
%     End:   last sample of the post-28dps grey gap (end of interval before
%            56 dps sweeps begin)
%
%   Returns concatenated voltage and boundary indices for verification.

    % Find all grey screen gaps >= 3s (30,000 samples at 10 kHz)
    zero_mask = f_data == 0;
    dd = diff([0, zero_mask, 0]);
    start_idx = find(dd == 1);
    end_idx   = find(dd == -1) - 1;
    len = end_idx - start_idx + 1;

    mask_3s = len >= 30000;
    idx_3_start = start_idx(mask_3s);  % Start of each >= 3s gap
    idx_3_end   = end_idx(mask_3s);    % End of each >= 3s gap

    % 28 dps block boundaries per rep:
    %   Gap indices:  idx_3(1) = 10s pre-rep gap
    %                 idx_3(2) = 3s gap before 28 dps sweeps
    %                 idx_3(3) = 3s gap between 28 dps and 56 dps
    %   Start: end of gap(2) + 1  (first sweep sample)
    %   End:   end of gap(3)      (last sample of post-28dps interval)
    % For reps 2,3: offset by 6 gaps per rep (k = 2, 8, 14)
    rep_offsets = [2, 8, 14];
    boundaries = zeros(3, 2);

    v_28 = [];
    for rep = 1:3
        k = rep_offsets(rep);
        b_start = idx_3_end(k) + 1;      % First sample after pre-sweep grey gap
        b_end   = idx_3_end(k + 1);      % Last sample of post-sweep grey gap
        boundaries(rep, :) = [b_start, b_end];
        v_28 = [v_28, v_data(b_start:b_end)]; %#ok<AGROW>
    end

end


%% ========================= Plotting Functions ==========================

function fig = plot_individual_histogram(v_data, group, title_str, bin_edges)
% PLOT_INDIVIDUAL_HISTOGRAM  Single-cell voltage histogram.

    if contains(group, 'control')
        line_col = [0 0 0];
        fill_col = [0.80 0.80 0.80];
    else
        line_col = [1 0 0];
        fill_col = [1 0.70 0.70];
    end

    [counts, ~] = histcounts(v_data, bin_edges, 'Normalization', 'pdf');
    bin_centers = bin_edges(1:end-1) + diff(bin_edges) / 2;

    fig = figure('Name', title_str);
    hold on;

    % Shaded area fill
    fill([bin_centers, fliplr(bin_centers)], [counts, zeros(size(counts))], ...
        fill_col, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

    % Line on top
    plot(bin_centers, counts, '-', 'Color', line_col, 'LineWidth', 1.2);

    xlabel('Voltage (mV)');
    ylabel('Probability Density');
    title(title_str, 'Interpreter', 'none');
    box off;
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.LineWidth = 1.2;
    ax.FontSize = 12;

    hold off;

end


function fig = plot_population_histogram(v_ctrl, v_ttl, bin_edges, title_str)
% PLOT_POPULATION_HISTOGRAM  Population voltage histogram with median +/- MAD.
%
%   Computes per-cell normalized histograms (PDF), then takes the median
%   across cells with MAD spread bands.

    bin_centers = bin_edges(1:end-1) + diff(bin_edges) / 2;
    n_bins = numel(bin_centers);

    % Colors
    col_ctrl_line = [0 0 0];
    col_ctrl_fill = [0.80 0.80 0.80];
    col_ttl_line  = [1 0 0];
    col_ttl_fill  = [1 0.70 0.70];
    alpha_val     = 0.4;

    fig = figure('Name', title_str);
    hold on;

    % Compute and plot control
    n_ctrl = numel(v_ctrl);
    if n_ctrl > 0
        hist_matrix_ctrl = zeros(n_ctrl, n_bins);
        for k = 1:n_ctrl
            hist_matrix_ctrl(k, :) = histcounts(v_ctrl{k}, bin_edges, 'Normalization', 'pdf');
        end
        med_ctrl = median(hist_matrix_ctrl, 1);
        mad_ctrl = mad(hist_matrix_ctrl, 1, 1);

        % Shaded spread band
        upper = med_ctrl + mad_ctrl;
        lower = max(med_ctrl - mad_ctrl, 0);
        fill([bin_centers, fliplr(bin_centers)], [upper, fliplr(lower)], ...
            col_ctrl_fill, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
        plot(bin_centers, med_ctrl, '-', 'Color', col_ctrl_line, 'LineWidth', 1.5);
    end

    % Compute and plot TTL
    n_ttl = numel(v_ttl);
    if n_ttl > 0
        hist_matrix_ttl = zeros(n_ttl, n_bins);
        for k = 1:n_ttl
            hist_matrix_ttl(k, :) = histcounts(v_ttl{k}, bin_edges, 'Normalization', 'pdf');
        end
        med_ttl = median(hist_matrix_ttl, 1);
        mad_ttl = mad(hist_matrix_ttl, 1, 1);

        upper = med_ttl + mad_ttl;
        lower = max(med_ttl - mad_ttl, 0);
        fill([bin_centers, fliplr(bin_centers)], [upper, fliplr(lower)], ...
            col_ttl_fill, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
        plot(bin_centers, med_ttl, '-', 'Color', col_ttl_line, 'LineWidth', 1.5);
    end

    xlabel('Voltage (mV)');
    ylabel('Probability Density');
    title(title_str);

    % Legend
    h = [];
    labs = {};
    if n_ctrl > 0
        h(end+1) = plot(NaN, NaN, '-', 'Color', col_ctrl_line, 'LineWidth', 1.5);
        labs{end+1} = sprintf('control (n=%d)', n_ctrl);
    end
    if n_ttl > 0
        h(end+1) = plot(NaN, NaN, '-', 'Color', col_ttl_line, 'LineWidth', 1.5);
        labs{end+1} = sprintf('ttl (n=%d)', n_ttl);
    end
    if ~isempty(h)
        legend(h, labs, 'Location', 'northwest', 'FontSize', 10);
    end

    box off;
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.LineWidth = 1.2;
    ax.FontSize = 12;

    hold off;

end


function fig = plot_verification(f_data, boundaries, title_str)
% PLOT_VERIFICATION  Frame data with red lines at 28 dps extraction boundaries.
%
%   Plots Log.ADC.Volts(1,:) and overlays red vertical lines showing the
%   start and end of each 28 dps block across the 3 repetitions.

    fig = figure('Name', title_str);
    fig.Position = [18 714 1749 244];

    plot(f_data, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5);
    hold on;

    y_lim = [min(f_data) - 5, max(f_data) + 5];

    for rep = 1:3
        b_start = boundaries(rep, 1);
        b_end   = boundaries(rep, 2);
        plot([b_start, b_start], y_lim, 'r-', 'LineWidth', 1.5);
        plot([b_end, b_end], y_lim, 'r-', 'LineWidth', 1.5);
    end

    ylim(y_lim);
    xlabel('Sample');
    ylabel('Frame');
    title(title_str, 'Interpreter', 'none');

    box off;
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.LineWidth = 1.2;
    ax.FontSize = 12;

    hold off;

end


%% ========================= Default Options ============================

function opts = set_defaults(opts, data_root)
% SET_DEFAULTS  Fill in default values for histogram analysis options.

    if ~isfield(opts, 'on_threshold')
        opts.on_threshold = 129;
    end
    if ~isfield(opts, 'bin_width')
        opts.bin_width = 0.1;   % mV
    end
    if ~isfield(opts, 'save_figs')
        opts.save_figs = true;
    end
    if ~isfield(opts, 'save_dir')
        opts.save_dir = fullfile(data_root, 'population_results');
    end

end
