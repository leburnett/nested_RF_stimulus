%% analyze_bar_flash.m
% Analyze voltage recordings from bar flash stimuli experiments.
%
% Pipeline:
%   1. Load .mat files and group by condition (control/ttl x on/off)
%   2. Identify preferred orientation column per cell (two methods)
%   3. Reorder rows so maximum response is at position 6
%   4. Align time series by peak (centered at sample 10,000 in 20,000 array)
%   5. Average across cells per condition
%   6. Visualize: 1x11 subplots with individual traces + mean overlay
%   7. Generate summary report

clear; close all; clc;

%% ========================================================================
%  PARAMETERS
%  ========================================================================

data_dir = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/bar_flash_results/1DRF';
output_dir = data_dir;  % Save outputs alongside data

num_positions = 11;    % Rows in mean_slow
num_orientations = 8;  % Columns in mean_slow
center_position = 6;   % Target position for max response
aligned_length = 20000; % Length of peak-aligned arrays
peak_center = 10000;   % Sample index where peak is centered

%% ========================================================================
%  PHASE 1: FILE LOADING AND GROUPING
%  ========================================================================

fprintf('=== Phase 1: Loading and grouping files ===\n');

mat_files = dir(fullfile(data_dir, '*.mat'));
num_files = numel(mat_files);
fprintf('Found %d .mat files\n', num_files);

% Group files by condition
conditions = {'control_on', 'control_off', 'ttl_on', 'ttl_off'};
data_by_condition = struct();
for c = 1:numel(conditions)
    data_by_condition.(conditions{c}) = {};
end
file_conditions = cell(num_files, 1);

for f = 1:num_files
    fname = mat_files(f).name;
    fname_lower = lower(fname);

    % Determine condition
    if contains(fname_lower, 'control') && contains(fname_lower, '_on')
        cond = 'control_on';
    elseif contains(fname_lower, 'control') && contains(fname_lower, '_off')
        cond = 'control_off';
    elseif contains(fname_lower, 'ttl') && contains(fname_lower, '_on')
        cond = 'ttl_on';
    elseif contains(fname_lower, 'ttl') && contains(fname_lower, '_off')
        cond = 'ttl_off';
    else
        warning('Could not classify file: %s', fname);
        file_conditions{f} = 'unknown';
        continue;
    end
    file_conditions{f} = cond;

    % Load mean_slow
    fpath = fullfile(data_dir, fname);
    loaded = load(fpath, 'mean_slow');

    if ~isfield(loaded, 'mean_slow')
        warning('No mean_slow variable in %s, skipping.', fname);
        continue;
    end

    mean_slow = loaded.mean_slow;

    % Validate dimensions
    [nr, nc] = size(mean_slow);
    if nr ~= num_positions || nc ~= num_orientations
        warning('Unexpected dimensions [%d x %d] in %s, expected [%d x %d]. Skipping.', ...
            nr, nc, fname, num_positions, num_orientations);
        continue;
    end

    % Flatten any (1,N) arrays to (N,) for consistency
    for i = 1:nr
        for j = 1:nc
            ts = mean_slow{i,j};
            if size(ts, 1) == 1
                mean_slow{i,j} = ts(:);
            end
        end
    end

    data_by_condition.(cond){end+1} = struct('mean_slow', {mean_slow}, 'filename', fname);
end

% Print summary
fprintf('\nFile grouping:\n');
for c = 1:numel(conditions)
    fprintf('  %-12s: %d cells\n', conditions{c}, numel(data_by_condition.(conditions{c})));
end

%% ========================================================================
%  PHASE 2: COLUMN IDENTIFICATION (PER CELL)
%  ========================================================================

fprintf('\n=== Phase 2: Identifying preferred orientation column ===\n');

% Store results for all cells across all conditions
all_results = struct();

for c = 1:numel(conditions)
    cond = conditions{c};
    cells = data_by_condition.(cond);
    n_cells = numel(cells);

    col_method1 = zeros(n_cells, 1);
    col_method2 = zeros(n_cells, 1);
    val_method1 = zeros(n_cells, 1);
    val_method2 = zeros(n_cells, 1);

    for ci = 1:n_cells
        ms = cells{ci}.mean_slow;

        % Method 1: Maximum voltage response
        [col_method1(ci), val_method1(ci)] = find_preferred_column_method1(ms);

        % Method 2: Largest response difference
        [col_method2(ci), val_method2(ci)] = find_preferred_column_method2(ms);
    end

    all_results.(cond).col_method1 = col_method1;
    all_results.(cond).col_method2 = col_method2;
    all_results.(cond).val_method1 = val_method1;
    all_results.(cond).val_method2 = val_method2;

    % Print comparison
    fprintf('\n--- %s (%d cells) ---\n', cond, n_cells);
    fprintf('  Cell | Method1 (col) | Method2 (col) | Agree?\n');
    fprintf('  -----|---------------|---------------|-------\n');
    for ci = 1:n_cells
        agree = col_method1(ci) == col_method2(ci);
        if agree
            agree_str = 'Yes';
        else
            agree_str = '** NO **';
        end
        fprintf('  %4d | %13d | %13d | %s\n', ci, col_method1(ci), col_method2(ci), agree_str);
    end
    if n_cells > 0
        pct_agree = sum(col_method1 == col_method2) / n_cells * 100;
        fprintf('  Agreement: %.1f%%\n', pct_agree);
    end
end

%% ========================================================================
%  PHASE 3-4: DATA REORDERING AND PEAK ALIGNMENT
%  ========================================================================

fprintf('\n=== Phase 3-4: Reordering and peak alignment ===\n');

aligned_data = struct();

for c = 1:numel(conditions)
    cond = conditions{c};
    cells = data_by_condition.(cond);
    n_cells = numel(cells);
    col_idxs = all_results.(cond).col_method1;  % Use Method 1

    % Store aligned data: n_cells x num_positions x aligned_length
    cond_aligned = NaN(n_cells, num_positions, aligned_length);

    for ci = 1:n_cells
        ms = cells{ci}.mean_slow;
        best_col = col_idxs(ci);

        % Phase 3: Reorder rows to center max at position 6
        reordered = reorder_to_center_max(ms, best_col, center_position);

        % Phase 4: Align peaks
        for pos = 1:num_positions
            ts = reordered{pos};
            cond_aligned(ci, pos, :) = align_peak(ts, aligned_length, peak_center);
        end
    end

    aligned_data.(cond) = cond_aligned;
    fprintf('  %s: aligned %d cells\n', cond, n_cells);
end

%% ========================================================================
%  PHASE 5: CROSS-CELL AVERAGING
%  ========================================================================

fprintf('\n=== Phase 5: Cross-cell averaging ===\n');

mean_data = struct();
se_data = struct();

for c = 1:numel(conditions)
    cond = conditions{c};
    ad = aligned_data.(cond);  % n_cells x num_positions x aligned_length
    n_cells = size(ad, 1);

    if n_cells == 0
        mean_data.(cond) = [];
        se_data.(cond) = [];
        continue;
    end

    cond_mean = NaN(num_positions, aligned_length);
    cond_se = NaN(num_positions, aligned_length);

    for pos = 1:num_positions
        pos_data = squeeze(ad(:, pos, :));  % n_cells x aligned_length
        if n_cells == 1
            pos_data = reshape(pos_data, 1, []);
        end
        cond_mean(pos, :) = nanmean(pos_data, 1);
        cond_se(pos, :) = nanstd(pos_data, 0, 1) / sqrt(n_cells);
    end

    mean_data.(cond) = cond_mean;
    se_data.(cond) = cond_se;
    fprintf('  %s: computed mean and SE across %d cells\n', cond, n_cells);
end

%% ========================================================================
%  PHASE 6: VISUALIZATION
%  ========================================================================

fprintf('\n=== Phase 6: Generating figures ===\n');

for c = 1:numel(conditions)
    cond = conditions{c};
    ad = aligned_data.(cond);
    n_cells = size(ad, 1);

    if n_cells == 0
        fprintf('  %s: no cells, skipping figure.\n', cond);
        continue;
    end

    cm = mean_data.(cond);

    fig = figure('Position', [50 200 2000 400], 'Visible', 'off');
    t = tiledlayout(1, num_positions, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Determine shared y-axis limits
    all_vals = ad(:);
    all_vals = all_vals(~isnan(all_vals));
    if ~isempty(all_vals)
        y_lo = prctile(all_vals, 0.5);
        y_hi = prctile(all_vals, 99.5);
    else
        y_lo = -70; y_hi = -40;
    end

    for pos = 1:num_positions
        nexttile;
        hold on;

        % Plot individual cells as thin lines
        for ci = 1:n_cells
            ts = squeeze(ad(ci, pos, :));
            plot(ts, 'Color', [0.6 0.6 0.6 0.3], 'LineWidth', 0.5);
        end

        % Overlay mean as thick line
        plot(cm(pos, :), 'Color', [0.2 0.2 0.8], 'LineWidth', 2.5);

        ylim([y_lo y_hi]);
        xlim([1 aligned_length]);
        title(sprintf('Pos %d', pos));

        if pos == 1
            ylabel('Voltage');
        else
            set(gca, 'YTickLabel', []);
        end
        set(gca, 'XTick', []);

        % Highlight position 6
        if pos == center_position
            ax = gca;
            ax.Box = 'on';
            ax.LineWidth = 2;
            ax.XColor = [0.8 0 0];
            ax.YColor = [0.8 0 0];
        end

        grid on;

        % Legend in first subplot only
        if pos == 1
            legend({'Individual cells', 'Mean'}, 'Location', 'northwest', 'FontSize', 6);
        end
    end

    cond_title = strrep(cond, '_', ' ');
    title(t, sprintf('%s - Bar Flash Response by Position (n=%d)', cond_title, n_cells));

    % Save figure
    saveas(fig, fullfile(output_dir, sprintf('bar_flash_analysis_%s.fig', cond)));
    exportgraphics(fig, fullfile(output_dir, sprintf('bar_flash_analysis_%s.png', cond)), 'Resolution', 300);
    fprintf('  Saved figure for %s\n', cond);
    close(fig);
end

%% ========================================================================
%  PHASE 7: SUMMARY REPORT
%  ========================================================================

fprintf('\n=== Phase 7: Generating summary report ===\n');

report_path = fullfile(output_dir, 'analysis_summary.txt');
fid = fopen(report_path, 'w');

fprintf(fid, 'Bar Flash Analysis Summary Report\n');
fprintf(fid, '=================================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

% File counts
fprintf(fid, '1. FILE GROUPING\n');
fprintf(fid, '   Total files: %d\n', num_files);
for c = 1:numel(conditions)
    cond = conditions{c};
    fprintf(fid, '   %-12s: %d cells\n', cond, numel(data_by_condition.(cond)));
end

% Column selection comparison
fprintf(fid, '\n2. COLUMN SELECTION COMPARISON (Method 1 vs Method 2)\n');
for c = 1:numel(conditions)
    cond = conditions{c};
    n_cells = numel(data_by_condition.(cond));
    if n_cells == 0, continue; end

    m1 = all_results.(cond).col_method1;
    m2 = all_results.(cond).col_method2;

    fprintf(fid, '\n   --- %s ---\n', cond);
    fprintf(fid, '   Cell | Method1 | Method2 | Agree\n');
    for ci = 1:n_cells
        if m1(ci) == m2(ci)
            agree_str = 'Yes';
        else
            agree_str = 'NO';
        end
        fprintf(fid, '   %4d | %7d | %7d | %s\n', ci, m1(ci), m2(ci), agree_str);
    end
    fprintf(fid, '   Agreement: %.1f%%\n', sum(m1 == m2) / n_cells * 100);
end

% Per-condition statistics
fprintf(fid, '\n3. PER-CONDITION STATISTICS\n');
for c = 1:numel(conditions)
    cond = conditions{c};
    ad = aligned_data.(cond);
    n_cells = size(ad, 1);
    if n_cells == 0, continue; end

    % Peak voltage at center position
    center_traces = squeeze(ad(:, center_position, :));
    if n_cells == 1
        center_traces = reshape(center_traces, 1, []);
    end
    peak_voltages = nanmax(center_traces, [], 2);

    fprintf(fid, '\n   --- %s (n=%d) ---\n', cond, n_cells);
    fprintf(fid, '   Mean peak voltage (position %d): %.3f\n', center_position, mean(peak_voltages));
    fprintf(fid, '   Std peak voltage:  %.3f\n', std(peak_voltages));
    fprintf(fid, '   Preferred columns (Method 1): %s\n', mat2str(all_results.(cond).col_method1'));
end

% Time series length report
fprintf(fid, '\n4. TIME SERIES LENGTH CHECK\n');
for c = 1:numel(conditions)
    cond = conditions{c};
    cells = data_by_condition.(cond);
    for ci = 1:numel(cells)
        ms = cells{ci}.mean_slow;
        lengths = cellfun(@numel, ms);
        unique_lens = unique(lengths(:));
        if numel(unique_lens) > 1
            fprintf(fid, '   %s cell %d (%s): variable lengths %s\n', ...
                cond, ci, cells{ci}.filename, mat2str(unique_lens'));
        end
    end
end

fclose(fid);
fprintf('  Saved report to %s\n', report_path);

%% Save processed data
save_path = fullfile(output_dir, 'processed_data.mat');
save(save_path, 'aligned_data', 'mean_data', 'se_data', 'all_results', ...
    'data_by_condition', 'conditions', '-v7.3');
fprintf('  Saved processed data to %s\n', save_path);

fprintf('\n=== Analysis complete ===\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [col_idx, max_val] = find_preferred_column_method1(mean_slow)
% Find column with maximum voltage response across all positions.
    [num_pos, num_orient] = size(mean_slow);
    col_maxes = zeros(num_orient, 1);
    for j = 1:num_orient
        for i = 1:num_pos
            ts = mean_slow{i,j};
            col_maxes(j) = max(col_maxes(j), max(ts));
        end
    end
    [max_val, col_idx] = max(col_maxes);
end

function [col_idx, diff_val] = find_preferred_column_method2(mean_slow)
% Find column with largest response difference (max - min) across positions.
    [num_pos, num_orient] = size(mean_slow);
    col_diffs = zeros(num_orient, 1);
    for j = 1:num_orient
        pos_maxes = zeros(num_pos, 1);
        for i = 1:num_pos
            ts = mean_slow{i,j};
            pos_maxes(i) = max(ts);
        end
        col_diffs(j) = max(pos_maxes) - min(pos_maxes);
    end
    [diff_val, col_idx] = max(col_diffs);
end

function reordered = reorder_to_center_max(mean_slow, col_idx, center_pos)
% Circularly shift rows of the selected column so max response is at center_pos.
    num_pos = size(mean_slow, 1);

    % Find which row has the maximum voltage in the selected column
    row_maxes = zeros(num_pos, 1);
    for i = 1:num_pos
        ts = mean_slow{i, col_idx};
        row_maxes(i) = max(ts);
    end
    [~, max_row] = max(row_maxes);

    % Calculate circular shift
    shift_amount = center_pos - max_row;

    % Extract column data and circularly shift
    col_data = mean_slow(:, col_idx);
    reordered = circshift(col_data, shift_amount);
end

function aligned = align_peak(ts, output_length, peak_center_idx)
% Align a time series so its peak is at peak_center_idx within an array of output_length.
    aligned = NaN(output_length, 1);

    if isempty(ts) || all(isnan(ts))
        return;
    end

    [~, peak_idx] = max(ts);
    n = numel(ts);

    % Calculate where to insert the time series
    start_idx = peak_center_idx - peak_idx + 1;
    end_idx = start_idx + n - 1;

    % Handle edge cases: clamp to array bounds
    ts_start = 1;
    ts_end = n;

    if start_idx < 1
        ts_start = 2 - start_idx;
        start_idx = 1;
    end
    if end_idx > output_length
        ts_end = n - (end_idx - output_length);
        end_idx = output_length;
    end

    aligned(start_idx:end_idx) = ts(ts_start:ts_end);
end
