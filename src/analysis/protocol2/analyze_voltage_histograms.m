function results = analyze_voltage_histograms(data_root, opts)
% ANALYZE_VOLTAGE_HISTOGRAMS  Voltage distribution analysis for 28 dps bar sweeps.
%
%   RESULTS = ANALYZE_VOLTAGE_HISTOGRAMS(DATA_ROOT) processes all 1DRF
%   experiments, extracts the contiguous voltage block during 28 dps
%   (slowest speed) bar sweeps, and generates voltage histograms comparing
%   control vs TTL conditions grouped by ON/OFF cell type.
%
%   Includes six analyses:
%     0. Baseline voltage comparison (box + whisker with Mann-Whitney)
%     1. Raw voltage histograms (PDF, median +/- MAD)
%     2. Baseline-subtracted voltage histograms
%     3. Cumulative distribution functions with K-S testing
%     4. Quantile comparison (5th, 10th percentiles)
%     5. Per-direction histogram breakdown (4x4 grid)
%     6. Time-resolved voltage analysis (1s windows)
%
%   RESULTS = ANALYZE_VOLTAGE_HISTOGRAMS(DATA_ROOT, OPTS) uses options to
%   override defaults.
%
%   INPUTS:
%     data_root - Path to 1DRF data directory containing experiment folders
%     opts      - (Optional) structure with fields:
%                   .on_threshold  - Frame threshold for ON/OFF (default: 129)
%                   .bin_width     - Histogram bin width in mV (default: 0.1)
%                   .save_figs     - Save figures as PDF (default: true)
%                   .save_dir      - Output directory for population PDFs
%                                    (default: <data_root>/population_results)
%                   .plot_order    - 1x16 data row to angle mapping
%                   .window_size   - Samples per time-resolved window (default: 10000)
%
%   OUTPUT:
%     results - Structure array with per-cell fields:
%       .folder          - Experiment folder name
%       .group           - 'on_control', 'on_ttl', 'off_control', or 'off_ttl'
%       .strain, .frame, .is_on, .is_ttl - Classification info
%       .v_28dps         - Concatenated voltage during 28 dps epochs (3 reps)
%       .v_28dps_bl_sub  - Baseline-subtracted version
%       .boundaries      - 3x2 matrix of [start, end] sample indices per rep
%       .baseline_median - Median voltage during 10s pre-stimulus grey screen
%       .pct_5, .pct_10  - 5th and 10th percentiles of v_28dps
%       .bar_data_mean   - 16x1 cell of mean traces per direction
%       .dir_angles      - 16x1 direction angles (degrees, ordered by plot_order)
%
%   See also BATCH_ANALYZE_1DRF, PARSE_BAR_DATA, LOAD_PROTOCOL2_DATA

    if nargin < 2, opts = struct(); end
    opts = set_defaults(opts, data_root);

    %% Load LUT (shared across experiments)
    S_lut = load(opts.lut_path, 'Tbl');
    Tbl = S_lut.Tbl;

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
            [r, f_data_full] = process_single_cell(exp_folder, Tbl, opts);
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

            fprintf('  -> %s | n_samples: %d | baseline: %.1f mV\n', ...
                r.group, numel(r.v_28dps), r.baseline_median);

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
    global_min = floor(min(all_v(:, 1)) * 10) / 10;
    global_max = ceil(max(all_v(:, 2)) * 10) / 10;
    bin_edges = global_min:opts.bin_width:global_max;
    fprintf('\nHistogram range: [%.1f, %.1f] mV, %d bins\n', ...
        global_min, global_max, numel(bin_edges) - 1);

    % Baseline-subtracted bin edges
    all_v_bl = cellfun(@(x) [min(x), max(x)], {results.v_28dps_bl_sub}, 'UniformOutput', false);
    all_v_bl = vertcat(all_v_bl{:});
    bl_min = floor(min(all_v_bl(:, 1)) * 10) / 10;
    bl_max = ceil(max(all_v_bl(:, 2)) * 10) / 10;
    bin_edges_bl = bl_min:opts.bin_width:bl_max;

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

    %% Ensure output directory exists
    if ~isfolder(opts.save_dir)
        mkdir(opts.save_dir);
    end

    %% Generate all population figures (loop over ON/OFF)
    fprintf('\nGenerating population figures...\n');

    for on_off_label = ["ON", "OFF"]
        if on_off_label == "ON"
            mask = [results.is_on];
        else
            mask = ~[results.is_on];
        end

        ctrl_mask = mask & ~[results.is_ttl];
        ttl_mask  = mask & [results.is_ttl];
        prefix = lower(char(on_off_label));

        % --- Raw voltage data ---
        v_ctrl = {results(ctrl_mask).v_28dps};
        v_ttl  = {results(ttl_mask).v_28dps};

        % --- Baseline-subtracted voltage data ---
        v_ctrl_bl = {results(ctrl_mask).v_28dps_bl_sub};
        v_ttl_bl  = {results(ttl_mask).v_28dps_bl_sub};

        % ==================== Analysis 0: Baseline comparison ====================
        bl_ctrl = [results(ctrl_mask).baseline_median];
        bl_ttl  = [results(ttl_mask).baseline_median];

        fig_bl = plot_baseline_comparison(bl_ctrl, bl_ttl, ...
            sprintf('%s Cells — Baseline Voltage', on_off_label));
        save_fig(fig_bl, opts, sprintf('%s_baseline_voltage_ctrl_vs_ttl.pdf', prefix));

        % ==================== Analysis 1: Raw voltage histograms =================
        pop_title = sprintf('%s Cells — 28 dps Voltage Distribution', on_off_label);
        fig_pop = plot_population_histogram(v_ctrl, v_ttl, bin_edges, pop_title);
        save_fig(fig_pop, opts, sprintf('%s_voltage_histogram_ctrl_vs_ttl.pdf', prefix));

        % ==================== Analysis 2: Baseline-subtracted histograms =========
        bl_title = sprintf('%s Cells — 28 dps Baseline-Subtracted Distribution', on_off_label);
        fig_bl_hist = plot_population_histogram(v_ctrl_bl, v_ttl_bl, bin_edges_bl, bl_title);
        save_fig(fig_bl_hist, opts, sprintf('%s_voltage_histogram_bl_sub_ctrl_vs_ttl.pdf', prefix));

        % ==================== Analysis 3: CDFs with K-S testing ==================
        cdf_title = sprintf('%s Cells — 28 dps Voltage CDF', on_off_label);
        fig_cdf = plot_population_cdf(v_ctrl, v_ttl, bin_edges, cdf_title);
        save_fig(fig_cdf, opts, sprintf('%s_cdf_ctrl_vs_ttl.pdf', prefix));

        % ==================== Analysis 4: Quantile comparison ====================
        pct5_ctrl  = [results(ctrl_mask).pct_5];
        pct5_ttl   = [results(ttl_mask).pct_5];
        pct10_ctrl = [results(ctrl_mask).pct_10];
        pct10_ttl  = [results(ttl_mask).pct_10];

        fig_p5 = plot_quantile_comparison(pct5_ctrl, pct5_ttl, '5th', ...
            sprintf('%s Cells — 5th Percentile', on_off_label));
        save_fig(fig_p5, opts, sprintf('%s_5th_percentile_ctrl_vs_ttl.pdf', prefix));

        fig_p10 = plot_quantile_comparison(pct10_ctrl, pct10_ttl, '10th', ...
            sprintf('%s Cells — 10th Percentile', on_off_label));
        save_fig(fig_p10, opts, sprintf('%s_10th_percentile_ctrl_vs_ttl.pdf', prefix));

        % ==================== Analysis 5: Per-direction breakdown ================
        traces_ctrl = {results(ctrl_mask).bar_data_mean};
        traces_ttl  = {results(ttl_mask).bar_data_mean};

        % Get direction angles (same for all cells)
        if any(ctrl_mask)
            angles = results(find(ctrl_mask, 1)).dir_angles;
        elseif any(ttl_mask)
            angles = results(find(ttl_mask, 1)).dir_angles;
        else
            angles = linspace(0, 337.5, 16)';
        end

        dir_title = sprintf('%s Cells — Per-Direction Voltage', on_off_label);
        fig_dir = plot_per_direction_histograms(traces_ctrl, traces_ttl, ...
            angles, bin_edges, dir_title);
        save_fig(fig_dir, opts, sprintf('%s_per_direction_histogram.pdf', prefix));

        % ==================== Analysis 6: Time-resolved analysis =================
        tr_title = sprintf('%s Cells — Time-Resolved Voltage', on_off_label);
        fig_tr = plot_time_resolved_voltage(v_ctrl, v_ttl, opts.window_size, tr_title);
        save_fig(fig_tr, opts, sprintf('%s_time_resolved_voltage.pdf', prefix));
    end

    fprintf('\nDone.\n');

end


%% ========================= Core Processing ============================

function [r, f_data] = process_single_cell(exp_folder, Tbl, opts)
% PROCESS_SINGLE_CELL  Load data, classify cell, extract all per-cell data.
%   Returns f_data separately for verification plots (not stored in struct).

    % Save/restore working directory (load_protocol2_data uses cd)
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);

    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    ce = load(fullfile(exp_folder, 'currentExp.mat'), ...
        'metadata', 'pattern_order', 'func_order');
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

    % Extract contiguous 28 dps voltage and baseline
    [r.v_28dps, r.boundaries, r.baseline_median] = ...
        extract_28dps_voltage(f_data, v_data);

    % Baseline-subtracted voltage
    r.v_28dps_bl_sub = r.v_28dps - r.baseline_median;

    % Percentiles
    r.pct_5  = prctile(r.v_28dps, 5);
    r.pct_10 = prctile(r.v_28dps, 10);

    % Per-direction traces
    bar_data = parse_bar_data(f_data, v_data);

    [lut_directions, ~, ~, ~] = ...
        verify_lut_directions(Tbl, ce.pattern_order, ce.func_order, opts.plot_order);

    % Store mean trace per direction (ordered by plot_order = sequential angles)
    r.bar_data_mean = cell(16, 1);
    r.dir_angles = zeros(16, 1);
    for subplot_idx = 1:16
        data_row = opts.plot_order(subplot_idx);
        r.bar_data_mean{subplot_idx} = bar_data{data_row, 4};
        r.dir_angles(subplot_idx) = lut_directions(data_row);
    end

end


function [v_28, boundaries, baseline_median] = extract_28dps_voltage(f_data, v_data)
% EXTRACT_28DPS_VOLTAGE  Extract contiguous voltage during 28 dps bar sweeps.
%
%   Extraction window per rep:
%     Start: first sample after the pre-28dps grey gap ends (first sweep onset)
%     End:   last sample of the post-28dps grey gap (end of interval before
%            56 dps sweeps begin)
%
%   Also computes baseline_median from the 10s pre-stimulus grey screen.

    % Find all grey screen gaps >= 3s (30,000 samples at 10 kHz)
    zero_mask = f_data == 0;
    dd = diff([0, zero_mask, 0]);
    start_idx = find(dd == 1);
    end_idx   = find(dd == -1) - 1;
    len = end_idx - start_idx + 1;

    mask_3s = len >= 30000;
    idx_3_start = start_idx(mask_3s);
    idx_3_end   = end_idx(mask_3s);

    % Baseline: median voltage during 10s pre-stimulus grey screen (gap 1)
    baseline_median = median(v_data(idx_3_start(1):idx_3_end(1)));

    % 28 dps block boundaries per rep
    rep_offsets = [2, 8, 14];
    boundaries = zeros(3, 2);

    v_28 = [];
    for rep = 1:3
        k = rep_offsets(rep);
        b_start = idx_3_end(k) + 1;
        b_end   = idx_3_end(k + 1);
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

    fill([bin_centers, fliplr(bin_centers)], [counts, zeros(size(counts))], ...
        fill_col, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    plot(bin_centers, counts, '-', 'Color', line_col, 'LineWidth', 1.2);

    xlabel('Voltage (mV)');
    ylabel('Probability Density');
    title(title_str, 'Interpreter', 'none');
    format_axes();
    hold off;

end


function fig = plot_population_histogram(v_ctrl, v_ttl, bin_edges, title_str)
% PLOT_POPULATION_HISTOGRAM  Population voltage histogram with median +/- MAD.

    bin_centers = bin_edges(1:end-1) + diff(bin_edges) / 2;
    n_bins = numel(bin_centers);
    [col_ctrl_line, col_ctrl_fill, col_ttl_line, col_ttl_fill, alpha_val] = get_colors();

    fig = figure('Name', title_str);
    hold on;

    n_ctrl = numel(v_ctrl);
    if n_ctrl > 0
        hist_matrix = zeros(n_ctrl, n_bins);
        for k = 1:n_ctrl
            hist_matrix(k, :) = histcounts(v_ctrl{k}, bin_edges, 'Normalization', 'pdf');
        end
        plot_median_mad_band(bin_centers, hist_matrix, col_ctrl_line, col_ctrl_fill, alpha_val);
    end

    n_ttl = numel(v_ttl);
    if n_ttl > 0
        hist_matrix = zeros(n_ttl, n_bins);
        for k = 1:n_ttl
            hist_matrix(k, :) = histcounts(v_ttl{k}, bin_edges, 'Normalization', 'pdf');
        end
        plot_median_mad_band(bin_centers, hist_matrix, col_ttl_line, col_ttl_fill, alpha_val);
    end

    xlabel('Voltage (mV)');
    ylabel('Probability Density');
    title(title_str);
    add_legend(n_ctrl, n_ttl, col_ctrl_line, col_ttl_line);
    format_axes();
    hold off;

end


function fig = plot_population_cdf(v_ctrl, v_ttl, bin_edges, title_str)
% PLOT_POPULATION_CDF  Population CDF with median +/- MAD and K-S test.

    bin_centers = bin_edges(1:end-1) + diff(bin_edges) / 2;
    n_bins = numel(bin_centers);
    [col_ctrl_line, col_ctrl_fill, col_ttl_line, col_ttl_fill, alpha_val] = get_colors();

    fig = figure('Name', title_str);
    hold on;

    n_ctrl = numel(v_ctrl);
    if n_ctrl > 0
        cdf_matrix_ctrl = zeros(n_ctrl, n_bins);
        for k = 1:n_ctrl
            cdf_matrix_ctrl(k, :) = histcounts(v_ctrl{k}, bin_edges, 'Normalization', 'cdf');
        end
        plot_median_mad_band(bin_centers, cdf_matrix_ctrl, col_ctrl_line, col_ctrl_fill, alpha_val);
    end

    n_ttl = numel(v_ttl);
    if n_ttl > 0
        cdf_matrix_ttl = zeros(n_ttl, n_bins);
        for k = 1:n_ttl
            cdf_matrix_ttl(k, :) = histcounts(v_ttl{k}, bin_edges, 'Normalization', 'cdf');
        end
        plot_median_mad_band(bin_centers, cdf_matrix_ttl, col_ttl_line, col_ttl_fill, alpha_val);
    end

    % K-S test on pooled voltage
    if n_ctrl > 0 && n_ttl > 0
        pooled_ctrl = horzcat(v_ctrl{:});
        pooled_ttl  = horzcat(v_ttl{:});
        [~, ks_p] = kstest2(pooled_ctrl, pooled_ttl);
        text(0.02, 0.95, sprintf('K-S p = %.4g', ks_p), ...
            'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
    end

    xlabel('Voltage (mV)');
    ylabel('Cumulative Probability');
    title(title_str);
    add_legend(n_ctrl, n_ttl, col_ctrl_line, col_ttl_line);
    format_axes();
    hold off;

end


function fig = plot_baseline_comparison(bl_ctrl, bl_ttl, title_str)
% PLOT_BASELINE_COMPARISON  Box + whisker of baseline voltage, ctrl vs TTL.

    fig = figure('Name', title_str);
    hold on;

    col_ctrl = [0 0 0];
    col_ttl  = [1 0 0];
    jitter_w = 0.12;

    n_ctrl = numel(bl_ctrl);
    n_ttl  = numel(bl_ttl);

    % Boxcharts
    if n_ctrl > 0
        boxchart(ones(n_ctrl, 1), bl_ctrl(:), ...
            'BoxFaceColor', 'k', 'BoxFaceAlpha', 0.15, 'MarkerStyle', 'none');
    end
    if n_ttl > 0
        boxchart(2 * ones(n_ttl, 1), bl_ttl(:), ...
            'BoxFaceColor', col_ttl, 'BoxFaceAlpha', 0.15, 'MarkerStyle', 'none');
    end

    % Jittered scatter overlay
    rng('default');
    if n_ctrl > 0
        scatter(1 + jitter_w * randn(n_ctrl, 1), bl_ctrl, 45, 'o', ...
            'MarkerEdgeColor', col_ctrl, 'LineWidth', 1);
    end
    if n_ttl > 0
        scatter(2 + jitter_w * randn(n_ttl, 1), bl_ttl, 45, 'o', ...
            'MarkerEdgeColor', col_ttl, 'LineWidth', 1);
    end

    % Mann-Whitney U test
    if n_ctrl >= 2 && n_ttl >= 2
        [p_val, ~, ~] = ranksum(bl_ctrl(:), bl_ttl(:));
        yl = ylim;
        plot([1 2], [yl(2) yl(2)], 'k-', 'LineWidth', 1);
        text(1.5, yl(2) + 0.05 * range(yl), sprintf('p = %.4g', p_val), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
        ylim([yl(1), yl(2) + 0.12 * range(yl)]);
    end

    set(gca, 'XTick', [1 2], 'XTickLabel', ...
        {sprintf('control (n=%d)', n_ctrl), sprintf('ttl (n=%d)', n_ttl)});
    ylabel('Baseline Voltage (mV)');
    title(title_str);
    format_axes();
    hold off;

end


function fig = plot_quantile_comparison(pct_ctrl, pct_ttl, pct_label, title_str)
% PLOT_QUANTILE_COMPARISON  Box + whisker of a voltage percentile, ctrl vs TTL.

    fig = figure('Name', title_str);
    hold on;

    col_ctrl = [0 0 0];
    col_ttl  = [1 0 0];
    jitter_w = 0.12;

    n_ctrl = numel(pct_ctrl);
    n_ttl  = numel(pct_ttl);

    if n_ctrl > 0
        boxchart(ones(n_ctrl, 1), pct_ctrl(:), ...
            'BoxFaceColor', 'k', 'BoxFaceAlpha', 0.15, 'MarkerStyle', 'none');
    end
    if n_ttl > 0
        boxchart(2 * ones(n_ttl, 1), pct_ttl(:), ...
            'BoxFaceColor', col_ttl, 'BoxFaceAlpha', 0.15, 'MarkerStyle', 'none');
    end

    rng('default');
    if n_ctrl > 0
        scatter(1 + jitter_w * randn(n_ctrl, 1), pct_ctrl, 45, 'o', ...
            'MarkerEdgeColor', col_ctrl, 'LineWidth', 1);
    end
    if n_ttl > 0
        scatter(2 + jitter_w * randn(n_ttl, 1), pct_ttl, 45, 'o', ...
            'MarkerEdgeColor', col_ttl, 'LineWidth', 1);
    end

    if n_ctrl >= 2 && n_ttl >= 2
        [p_val, ~, ~] = ranksum(pct_ctrl(:), pct_ttl(:));
        yl = ylim;
        plot([1 2], [yl(2) yl(2)], 'k-', 'LineWidth', 1);
        text(1.5, yl(2) + 0.05 * range(yl), sprintf('p = %.4g', p_val), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
        ylim([yl(1), yl(2) + 0.12 * range(yl)]);
    end

    set(gca, 'XTick', [1 2], 'XTickLabel', ...
        {sprintf('control (n=%d)', n_ctrl), sprintf('ttl (n=%d)', n_ttl)});
    ylabel(sprintf('%s Percentile Voltage (mV)', pct_label));
    title(title_str);
    format_axes();
    hold off;

end


function fig = plot_per_direction_histograms(traces_ctrl, traces_ttl, angles, bin_edges, title_str)
% PLOT_PER_DIRECTION_HISTOGRAMS  4x4 grid of per-direction voltage histograms.
%
%   traces_ctrl/ttl: cell arrays, one per cell. Each element is a 16x1 cell
%   of mean voltage traces (one per direction, ordered by plot_order).

    bin_centers = bin_edges(1:end-1) + diff(bin_edges) / 2;
    n_bins = numel(bin_centers);
    [col_ctrl_line, col_ctrl_fill, col_ttl_line, col_ttl_fill, alpha_val] = get_colors();

    n_ctrl = numel(traces_ctrl);
    n_ttl  = numel(traces_ttl);

    fig = figure('Name', title_str);
    fig.Position = [50 50 1200 1000];
    tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    for dir_idx = 1:16
        nexttile;
        hold on;

        % Collect per-cell histograms for this direction
        if n_ctrl > 0
            hist_mat = zeros(n_ctrl, n_bins);
            for k = 1:n_ctrl
                v_dir = traces_ctrl{k}{dir_idx};
                if ~isempty(v_dir)
                    hist_mat(k, :) = histcounts(v_dir, bin_edges, 'Normalization', 'pdf');
                end
            end
            plot_median_mad_band(bin_centers, hist_mat, col_ctrl_line, col_ctrl_fill, alpha_val);
        end

        if n_ttl > 0
            hist_mat = zeros(n_ttl, n_bins);
            for k = 1:n_ttl
                v_dir = traces_ttl{k}{dir_idx};
                if ~isempty(v_dir)
                    hist_mat(k, :) = histcounts(v_dir, bin_edges, 'Normalization', 'pdf');
                end
            end
            plot_median_mad_band(bin_centers, hist_mat, col_ttl_line, col_ttl_fill, alpha_val);
        end

        title(sprintf('%.1f%s', angles(dir_idx), char(176)));

        if dir_idx == 13  % bottom-left
            xlabel('mV');
            ylabel('PDF');
        else
            set(gca, 'XTickLabel', [], 'YTickLabel', []);
        end

        box off;
        ax = gca;
        ax.TickDir = 'out';
        ax.FontSize = 8;
        hold off;
    end

    sgtitle(title_str, 'FontSize', 14);

end


function fig = plot_time_resolved_voltage(v_ctrl, v_ttl, window_size, title_str)
% PLOT_TIME_RESOLVED_VOLTAGE  Time-resolved median voltage in fixed windows.

    [col_ctrl_line, col_ctrl_fill, col_ttl_line, col_ttl_fill, alpha_val] = get_colors();
    fs = 10000;  % 10 kHz sampling rate

    fig = figure('Name', title_str);
    hold on;

    n_ctrl = numel(v_ctrl);
    if n_ctrl > 0
        % Find minimum number of complete windows across cells
        min_windows = Inf;
        for k = 1:n_ctrl
            n_win = floor(numel(v_ctrl{k}) / window_size);
            min_windows = min(min_windows, n_win);
        end

        win_matrix = zeros(n_ctrl, min_windows);
        for k = 1:n_ctrl
            v = v_ctrl{k};
            for w = 1:min_windows
                idx_start = (w - 1) * window_size + 1;
                idx_end   = w * window_size;
                win_matrix(k, w) = median(v(idx_start:idx_end));
            end
        end

        t_sec = ((1:min_windows) - 0.5) * window_size / fs;
        med_vals = median(win_matrix, 1);
        mad_vals = mad(win_matrix, 1, 1);
        plot_shaded_trace(t_sec, med_vals, mad_vals, ...
            col_ctrl_line, col_ctrl_fill, alpha_val, 1.5);
    end

    n_ttl = numel(v_ttl);
    if n_ttl > 0
        min_windows = Inf;
        for k = 1:n_ttl
            n_win = floor(numel(v_ttl{k}) / window_size);
            min_windows = min(min_windows, n_win);
        end

        win_matrix = zeros(n_ttl, min_windows);
        for k = 1:n_ttl
            v = v_ttl{k};
            for w = 1:min_windows
                idx_start = (w - 1) * window_size + 1;
                idx_end   = w * window_size;
                win_matrix(k, w) = median(v(idx_start:idx_end));
            end
        end

        t_sec = ((1:min_windows) - 0.5) * window_size / fs;
        med_vals = median(win_matrix, 1);
        mad_vals = mad(win_matrix, 1, 1);
        plot_shaded_trace(t_sec, med_vals, mad_vals, ...
            col_ttl_line, col_ttl_fill, alpha_val, 1.5);
    end

    xlabel('Time (s)');
    ylabel('Median Voltage (mV)');
    title(title_str);
    add_legend(n_ctrl, n_ttl, col_ctrl_line, col_ttl_line);
    format_axes();
    hold off;

end


function fig = plot_verification(f_data, boundaries, title_str)
% PLOT_VERIFICATION  Frame data with red lines at 28 dps extraction boundaries.

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
    format_axes();
    hold off;

end


%% ========================= Shared Helpers ==============================

function plot_median_mad_band(x, data_matrix, line_col, fill_col, alpha_val)
% PLOT_MEDIAN_MAD_BAND  Plot median line with MAD shaded band on current axes.

    med_vals = median(data_matrix, 1);
    mad_vals = mad(data_matrix, 1, 1);

    upper = med_vals + mad_vals;
    lower = max(med_vals - mad_vals, 0);

    fill([x(:)', fliplr(x(:)')], [upper, fliplr(lower)], ...
        fill_col, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
    plot(x, med_vals, '-', 'Color', line_col, 'LineWidth', 1.5);

end


function plot_shaded_trace(x, center_vals, spread_vals, line_col, fill_col, alpha_val, lw)
% PLOT_SHADED_TRACE  Plot center line with shaded spread band on current axes.

    x = x(:)';
    m = center_vals(:)';
    s = spread_vals(:)';

    upper = m + s;
    lower = m - s;

    fill([x, fliplr(x)], [upper, fliplr(lower)], fill_col, ...
        'FaceAlpha', alpha_val, 'EdgeColor', 'none');
    plot(x, m, '-', 'Color', line_col, 'LineWidth', lw);

end


function [col_ctrl_line, col_ctrl_fill, col_ttl_line, col_ttl_fill, alpha_val] = get_colors()
% GET_COLORS  Return the standard control/TTL color scheme.

    col_ctrl_line = [0 0 0];
    col_ctrl_fill = [0.80 0.80 0.80];
    col_ttl_line  = [1 0 0];
    col_ttl_fill  = [1 0.70 0.70];
    alpha_val     = 0.4;

end


function add_legend(n_ctrl, n_ttl, col_ctrl_line, col_ttl_line)
% ADD_LEGEND  Add control/TTL legend with sample sizes.

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

end


function format_axes()
% FORMAT_AXES  Apply standard axis formatting.

    box off;
    ax = gca;
    ax.TickDir = 'out';
    ax.TickLength = [0.015 0.015];
    ax.LineWidth = 1.2;
    ax.FontSize = 12;

end


function save_fig(fig, opts, fname)
% SAVE_FIG  Export figure as vector PDF and close.

    if opts.save_figs
        exportgraphics(fig, fullfile(opts.save_dir, fname), ...
            'ContentType', 'vector', 'BackgroundColor', 'none');
        close(fig);
        fprintf('  Saved: %s\n', fname);
    end

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
    if ~isfield(opts, 'plot_order')
        opts.plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    end
    if ~isfield(opts, 'lut_path')
        script_dir = fileparts(mfilename('fullpath'));
        opts.lut_path = fullfile(script_dir, 'bar_lut.mat');
    end
    if ~isfield(opts, 'window_size')
        opts.window_size = 10000;  % 1 second at 10 kHz
    end

end
