function analyze_single_experiment(exp_folder)
% ANALYZE_SINGLE_EXPERIMENT  Bar sweep and bar flash analysis for one experiment.
%
%   ANALYZE_SINGLE_EXPERIMENT(EXP_FOLDER) loads data from a single Protocol 2
%   experiment, creates a slow-speed (28 dps) bar sweep polar timeseries plot
%   with LUT-based direction verification, then creates bar flash plots.
%
%   INPUT:
%     exp_folder - Full path to experiment directory
%                  e.g. '/Users/.../protocol2/data/1DRF/2025_11_10_10_17'
%
%   OUTPUTS:
%     Figure 1: LUT verification — bar orientation/direction visual check
%     Figure 2: Slow bar sweep polar timeseries with vector sum arrow
%     Figure 3: Full 8x11 bar flash heatmap (all orientations)
%     Figure 4: 1x11 bar flash subplots along PD-ND axis
%
%   See also LOAD_PROTOCOL2_DATA, PARSE_BAR_DATA, PARSE_BAR_FLASH_DATA,
%            VECTOR_SUM_POLAR, ADD_ARROW_TO_POLARPLOT

    %% Step 1: Data Loading

    % Load bar lookup table
    lut_path = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/1DRF/bar_lut.mat';
    S_lut = load(lut_path, 'Tbl');
    Tbl = S_lut.Tbl;

    % Save current directory and restore on exit (load_protocol2_data uses cd)
    orig_dir = pwd;
    cleanup = onCleanup(@() cd(orig_dir));

    % Load experiment data
    [date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);

    f_data = Log.ADC.Volts(1, :);     % frame position data
    v_data = Log.ADC.Volts(2, :)*10;  % voltage data (scaled)
    median_v = median(v_data);

    % Load currentExp for pattern/function ordering
    ce = load(fullfile(exp_folder, 'currentExp.mat'), 'pattern_order', 'func_order', 'metadata');
    pattern_order = ce.pattern_order;
    func_order = ce.func_order;
    metadata = ce.metadata;

    %% Step 2: Parse bar sweep data

    bar_data = parse_bar_data(f_data, v_data);

    %% Step 3: LUT-Based Verification

    % Find which entries in pattern_order/func_order correspond to slow bars
    slow_func_mask = (func_order == 3 | func_order == 4);
    slow_pattern_nums = pattern_order(slow_func_mask);
    slow_func_nums = func_order(slow_func_mask);

    % For each of the 16 slow bar data rows, look up direction and orientation
    n_dir = 16;
    lut_directions = zeros(n_dir, 1);
    lut_orientations = zeros(n_dir, 1);
    lut_patterns = zeros(n_dir, 1);
    lut_functions = zeros(n_dir, 1);

    for i = 1:n_dir
        p = slow_pattern_nums(i);
        f = slow_func_nums(i);
        row_mask = (Tbl.pattern == p) & (Tbl.function == f);
        lut_directions(i) = Tbl.direction(row_mask);
        lut_orientations(i) = Tbl.orientation(row_mask);
        lut_patterns(i) = p;
        lut_functions(i) = f;
    end

    % What does plot_order assume?
    plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    assumed_angles = linspace(0, 360 - 22.5, 16)';

    % Build comparison table
    fprintf('\n=== LUT-Based Direction Verification ===\n');
    fprintf('%-8s %-8s %-8s %-12s %-12s %-16s %-6s\n', ...
        'DataRow', 'Pattern', 'Func', 'LUT_Dir', 'LUT_Orient', 'Assumed_Dir', 'Match');
    fprintf('%s\n', repmat('-', 1, 72));

    all_match = true;
    for subplot_idx = 1:n_dir
        data_row = plot_order(subplot_idx);
        assumed_dir = assumed_angles(subplot_idx);
        actual_dir = lut_directions(data_row);
        actual_orient = lut_orientations(data_row);
        match = abs(actual_dir - assumed_dir) < 0.1;
        if ~match
            all_match = false;
        end
        fprintf('%-8d %-8d %-8d %-12.1f %-12.1f %-16.1f %-6s\n', ...
            data_row, lut_patterns(data_row), lut_functions(data_row), ...
            actual_dir, actual_orient, assumed_dir, string(match));
    end

    if all_match
        fprintf('\nAll directions match plot_order assumptions.\n');
    else
        fprintf('\nWARNING: Mismatches found! Using LUT directions for plotting.\n');
    end

    % Check for duplicate directions in the LUT
    unique_dirs = unique(lut_directions);
    if numel(unique_dirs) < n_dir
        fprintf('\nWARNING: LUT contains duplicate directions!\n');
        for d_idx = 1:numel(unique_dirs)
            rows_at_dir = find(lut_directions == unique_dirs(d_idx));
            if numel(rows_at_dir) > 1
                fprintf('  Direction %.1f° appears in data rows: %s\n', ...
                    unique_dirs(d_idx), mat2str(rows_at_dir'));
            end
        end
        missing = setdiff(assumed_angles, lut_directions);
        if ~isempty(missing)
            fprintf('  Missing directions: %s\n', mat2str(missing'));
        end
    end

    %% Step 3b: Visual verification figure
    fig_verify = figure('Name', 'LUT Verification: Bar Orientation & Direction');
    set(fig_verify, 'Position', [50 50 1600 800]);

    numPlots = 16;
    theta_plot = linspace(0, 2*pi, numPlots+1);
    theta_plot = theta_plot(1:end-1);

    verify_radius = 0.32;
    sw = 0.08; sh = 0.08;

    for panel = 1:2
        if panel == 1
            cX = 0.25;
            panel_title = 'Current plot\_order (assumed angles)';
        else
            cX = 0.75;
            panel_title = 'LUT directions (actual angles)';
        end
        cY = 0.5;

        for subplot_idx = 1:numPlots
            data_row = plot_order(subplot_idx);

            if panel == 1
                angle_rad = theta_plot(subplot_idx);
            else
                angle_rad = deg2rad(lut_directions(data_row));
            end

            x_pos = cX + verify_radius * cos(angle_rad);
            y_pos = cY + verify_radius * sin(angle_rad);

            ax = axes('Position', [x_pos - sw/2, y_pos - sh/2, sw, sh]);
            hold on;

            orient_deg = lut_orientations(data_row);
            dir_deg = lut_directions(data_row);

            orient_rad = deg2rad(orient_deg);
            bar_len = 0.35;
            bx = bar_len * cos(orient_rad);
            by = bar_len * sin(orient_rad);
            plot([-bx, bx], [-by, by], 'k-', 'LineWidth', 3);

            dir_rad = deg2rad(dir_deg);
            arr_len = 0.3;
            ax_arr = arr_len * cos(dir_rad);
            ay_arr = arr_len * sin(dir_rad);
            quiver(0, 0, ax_arr, ay_arr, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.8);

            xlim([-0.5 0.5]); ylim([-0.5 0.5]);
            axis equal off;

            text(0, -0.45, sprintf('R%d D:%.0f°', data_row, dir_deg), ...
                'HorizontalAlignment', 'center', 'FontSize', 5);
        end

        annotation('textbox', [cX-0.15, 0.88, 0.3, 0.05], ...
            'String', panel_title, ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
            'FontSize', 10, 'FontWeight', 'bold');

        annotation('textbox', [cX-0.1, 0.43, 0.2, 0.08], ...
            'String', {'Black = bar orientation', 'Red arrow = motion direction'}, ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 7);
    end

    sgtitle(sprintf('LUT Verification — %s %s — Black=orientation, Red=direction', ...
        strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' ')));

    %% Step 4: Slow-speed bar sweep polar timeseries plot

    fig_polar = figure('Name', 'Slow Bar Sweep Polar Timeseries (28 dps)');

    centerX = 0.5;
    centerY = 0.5;
    radius = 0.35;

    max_v_slow = zeros(numPlots, 1);
    col = [0.2 0.4 0.7];

    % LUT directions for each subplot position (in order of plot_order)
    lut_dirs_ordered = lut_directions(plot_order);  % 16x1, degrees
    lut_dirs_ordered_rad = deg2rad(lut_dirs_ordered);

    for subplot_idx = 1:numPlots
        data_row = plot_order(subplot_idx);

        % Position subplot at the LUT direction for this data row
        angle_rad = lut_dirs_ordered_rad(subplot_idx);
        x_pos = centerX + radius * cos(angle_rad);
        y_pos = centerY + radius * sin(angle_rad);

        subW = 0.15; subH = 0.15;
        ax = axes('Position', [x_pos - subW/2, y_pos - subH/2, subW, subH]);
        hold on;

        n_reps = size(bar_data, 2) - 1;

        for r = 1:n_reps+1
            d2plot = bar_data{data_row, r};
            x_vals = 1:numel(d2plot);

            if r == 1
                plot([1 x_vals(end)], [median_v, median_v], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            end

            if r < n_reps+1
                plot(x_vals, d2plot, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
            else
                plot(x_vals, d2plot, 'Color', col, 'LineWidth', 1.2);
            end
        end

        ylim([-80 -10]);
        axis(ax, 'off');

        % Compute max_v: depolarization relative to pre-stimulus baseline
        d_mean = bar_data{data_row, n_reps+1};
        d_before = d_mean(1000:9000);
        mean_before = mean(d_before);
        d_stim = d_mean(9000:end-7000);
        max_v_slow(subplot_idx) = abs(diff([prctile(d_stim, 98), mean_before]));
    end

    % Central polar plot — use LUT directions as the theta values
    centralSize = (2 * radius) * 0.65;
    centralPosition = [centerX - centralSize/2, centerY - centralSize/2, centralSize, centralSize];

    % Sort by direction for smooth curve
    dir_and_val = sortrows([lut_dirs_ordered_rad, max_v_slow], 1);
    polar_theta = [dir_and_val(:,1); dir_and_val(1,1)];  % close loop
    polar_rho = [dir_and_val(:,2); dir_and_val(1,2)];

    axCentral = polaraxes('Position', centralPosition);
    hold on;
    polarplot(polar_theta, polar_rho, 'Color', col, 'LineWidth', 2, ...
        'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w');

    % Vector sum arrow using existing functions
    rho_for_vecsum = max_v_slow';  % 1×16
    theta_for_vecsum = lut_dirs_ordered_rad';  % 1×16
    [~, resultant_angle] = vector_sum_polar(rho_for_vecsum, theta_for_vecsum);

    arrow_magnitude = max(polar_rho) * 0.9;
    add_arrow_to_polarplot(arrow_magnitude, resultant_angle, [0.3 0.3 0.3]);

    axCentral.LineWidth = 1.2;
    axCentral.ThetaTick = 0:22.5:337.5;
    axCentral.ThetaTickLabel = {};

    set(gcf, 'Position', [303 78 961 969]);
    sgtitle(sprintf('28 dps — %s — %s — %s', ...
        strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' '), params.on_off));

    %% Step 5: Find PD via vector sum — closest LUT direction

    fprintf('\n=== Preferred Direction (PD via vector sum) ===\n');
    fprintf('Vector sum angle: %.1f°\n', rad2deg(resultant_angle));

    % Find the closest LUT direction to the vector sum angle
    angle_diffs = abs(lut_dirs_ordered_rad - resultant_angle);
    % Handle wraparound
    angle_diffs = min(angle_diffs, 2*pi - angle_diffs);
    [~, pd_subplot_idx] = min(angle_diffs);
    pd_data_row = plot_order(pd_subplot_idx);
    pd_direction = lut_directions(pd_data_row);
    pd_orientation = lut_orientations(pd_data_row);
    pd_pattern = lut_patterns(pd_data_row);
    pd_function = lut_functions(pd_data_row);

    fprintf('Closest LUT direction: %.1f°\n', pd_direction);
    fprintf('Bar orientation:       %.1f°\n', pd_orientation);
    fprintf('Pattern:               %d\n', pd_pattern);
    fprintf('Function:              %d\n', pd_function);
    fprintf('Data row:              %d\n', pd_data_row);

    %% Step 6: Map PD orientation to bar flash column

    % Experiment patterns 3-10 → full-field patterns 1-8
    pd_ff_pattern = pd_pattern - 2;
    bar_flash_col = pd_ff_pattern;

    fprintf('Full-field pattern index: %d\n', pd_ff_pattern);
    fprintf('Bar flash column:         %d\n', bar_flash_col);

    %% Step 7: Parse bar flash data

    [data_slow_bf, ~, mean_slow_bf, ~] = parse_bar_flash_data(f_data, v_data);

    pd_flash_data = data_slow_bf(:, bar_flash_col, :);
    pd_flash_mean = mean_slow_bf(:, bar_flash_col);

    % Orthogonal orientation: 4 columns away (wraps around 8 orientations)
    ortho_flash_col = mod(bar_flash_col - 1 + 4, 8) + 1;
    ortho_flash_data = data_slow_bf(:, ortho_flash_col, :);
    ortho_flash_mean = mean_slow_bf(:, ortho_flash_col);

    % Look up orthogonal orientation angle
    ortho_exp_pat = ortho_flash_col + 2;
    ortho_mask = Tbl.pattern == ortho_exp_pat & Tbl.function == 3;
    ortho_orientation = Tbl.orientation(ortho_mask);

    fprintf('Orthogonal bar flash column:  %d (orientation: %.1f°)\n', ...
        ortho_flash_col, ortho_orientation);

    %% Step 8: Determine spatial ordering (ND→PD, left to right)

    is_forward = mod(pd_function, 2) == 1;
    if is_forward
        pos_order = 1:11;
    else
        pos_order = 11:-1:1;
    end

    fprintf('PD function %d is %s → position order: %s\n', ...
        pd_function, ternary(is_forward, 'forward', 'reverse'), mat2str(pos_order));

    %% Step 9: Full 8×11 bar flash plot (all orientations)

    orient_labels = cell(8, 1);
    for k = 1:8
        exp_pat = k + 2;
        mask = Tbl.pattern == exp_pat & Tbl.function == 3;
        orient_labels{k} = sprintf('Orient: %.0f°', Tbl.orientation(mask));
    end

    % Find max across all mean data for consistent y-limits
    max_overall = -Inf;
    for i = 1:88
        d = mean_slow_bf{i};
        if ~isempty(d)
            max_d = prctile(d(ceil(numel(d)*0.5):ceil(numel(d)*0.75)), 98);
            if max_d > max_overall
                max_overall = max_d;
            end
        end
    end

    % Compute normalised heatmap background values
    max_vals = zeros(11, 8);
    for orient_idx = 1:8
        for pos_idx = 1:11
            d = mean_slow_bf{pos_idx, orient_idx};
            if ~isempty(d)
                n_points = numel(d);
                max_vals(pos_idx, orient_idx) = prctile(d(ceil(n_points*0.5):ceil(n_points*0.75)), 98);
            end
        end
    end
    max_vals_med = max_vals - median_v;
    max_vals_med(max_vals_med < 0) = 0;
    if max(max_vals_med(:)) > 0
        normalizedArray = 1 - max_vals_med / max(max_vals_med(:));
    else
        normalizedArray = ones(size(max_vals_med));
    end

    fig_all_flash = figure('Name', 'Bar Flashes — All Orientations (80ms)');
    tiledlayout(8, 11, 'TileSpacing', 'compact', 'Padding', 'compact');

    for orient_idx = 1:8
        for pos_idx_raw = 1:11
            % For the PD row, use pos_order; for all others, use 1:11
            if orient_idx == bar_flash_col
                pos_idx = pos_order(pos_idx_raw);
            else
                pos_idx = pos_idx_raw;
            end

            nexttile;
            hold on;

            d = mean_slow_bf{pos_idx, orient_idx};
            if ~isempty(d)
                rectangle('Position', [0, -70, numel(d), 40], ...
                    'FaceColor', [1, normalizedArray(pos_idx, orient_idx), normalizedArray(pos_idx, orient_idx)]*0.9);
            end

            for r = 1:3
                ts = data_slow_bf{pos_idx, orient_idx, r};
                if ~isempty(ts)
                    plot(1:numel(ts), ts, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7);
                end
            end

            if ~isempty(d)
                if orient_idx == bar_flash_col
                    plot(1:numel(d), d, 'Color', [0.8 0 0], 'LineWidth', 1.2);
                else
                    plot(1:numel(d), d, 'k', 'LineWidth', 1);
                end
            end

            ylim([-70 max_overall*0.9]);
            if ~isempty(d)
                xlim([0 numel(d)]);
            end
            set(gca, 'XTick', [], 'YTick', []);

            if pos_idx_raw == 1
                if orient_idx == bar_flash_col
                    ylabel(sprintf('%s *PD*', orient_labels{orient_idx}), ...
                        'FontSize', 7, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                else
                    ylabel(orient_labels{orient_idx}, 'FontSize', 7);
                end
            end

            if orient_idx == 1
                title(sprintf('Pos %d', pos_idx_raw), 'FontSize', 6);
            end

            if orient_idx == bar_flash_col
                set(gca, 'XColor', [0.8 0 0], 'YColor', [0.8 0 0], 'LineWidth', 1.5);
                box on;
            end
        end
    end

    sgtitle(sprintf('Bar Flashes 80ms — PD: Dir %.0f° Orient %.0f° (row %d) — %s — %s', ...
        pd_direction, pd_orientation, bar_flash_col, ...
        strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' ')));
    set(gcf, 'Position', [45 128 1675 902]);

    %% Step 10: 1×11 bar flash subplot along PD-ND axis (baseline-subtracted)

    % Baseline period: first 5000 samples (~500ms before flash onset)
    baseline_samples = 1:5000;

    fig_flash = figure('Name', 'Bar Flash PD-ND Axis');
    tiledlayout(1, 11, 'TileSpacing', 'compact', 'Padding', 'compact');

    for pos_idx = 1:11
        flash_pos = pos_order(pos_idx);
        nexttile;
        hold on;

        for r = 1:3
            ts = pd_flash_data{flash_pos, 1, r};
            if ~isempty(ts)
                bl = mean(ts(baseline_samples));
                plot(1:numel(ts), ts - bl, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            end
        end

        ts_mean = pd_flash_mean{flash_pos};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(baseline_samples));
            plot(1:numel(ts_mean), ts_mean - bl_mean, 'k', 'LineWidth', 1.5);
        end

        ylim([-15 25]);
        if pos_idx == 1
            ylabel('\DeltamV');
        else
            set(gca, 'YTickLabel', []);
        end

        set(gca, 'XTick', []);

        if pos_idx == 1
            title('ND');
        elseif pos_idx == 6
            title('Center');
        elseif pos_idx == 11
            title('PD');
        else
            title(sprintf('%d', pos_idx));
        end

        box on;
    end

    sgtitle(sprintf('Bar Flash PD-ND — Dir:%.0f° Orient:%.0f° — %s — %s', ...
        pd_direction, pd_orientation, strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' ')));

    set(gcf, 'Position', [50 400 1800 300]);

    %% Step 11: 1×11 bar flash subplot — orthogonal orientation (baseline-subtracted)

    fig_flash_ortho = figure('Name', 'Bar Flash Orthogonal Axis');
    tiledlayout(1, 11, 'TileSpacing', 'compact', 'Padding', 'compact');

    for pos_idx = 1:11
        flash_pos = pos_order(pos_idx);
        nexttile;
        hold on;

        for r = 1:3
            ts = ortho_flash_data{flash_pos, 1, r};
            if ~isempty(ts)
                bl = mean(ts(baseline_samples));
                plot(1:numel(ts), ts - bl, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            end
        end

        ts_mean = ortho_flash_mean{flash_pos};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(baseline_samples));
            plot(1:numel(ts_mean), ts_mean - bl_mean, 'k', 'LineWidth', 1.5);
        end

        ylim([-15 25]);
        if pos_idx == 1
            ylabel('\DeltamV');
        else
            set(gca, 'YTickLabel', []);
        end

        set(gca, 'XTick', []);

        if pos_idx == 1
            title('ND');
        elseif pos_idx == 6
            title('Center');
        elseif pos_idx == 11
            title('PD');
        else
            title(sprintf('%d', pos_idx));
        end

        box on;
    end

    sgtitle(sprintf('Bar Flash Orthogonal — Orient:%.0f° — %s — %s', ...
        ortho_orientation, strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' ')));

    set(gcf, 'Position', [50 100 1800 300]);

    %% Save figures as PDFs
    save_dir = fullfile(exp_folder, 'analysis_output');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end

    exportgraphics(fig_verify, fullfile(save_dir, 'lut_verification.pdf'), ...
        'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_polar, fullfile(save_dir, 'slow_bar_sweep_polar.pdf'), ...
        'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_all_flash, fullfile(save_dir, 'bar_flash_all_orientations.pdf'), ...
        'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_flash, fullfile(save_dir, 'bar_flash_PD_ND_axis.pdf'), ...
        'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_flash_ortho, fullfile(save_dir, 'bar_flash_orthogonal_axis.pdf'), ...
        'ContentType', 'image', 'Resolution', 300);
    fprintf('\nFigures saved to: %s\n', save_dir);

end

function out = ternary(cond, val_true, val_false)
    if cond
        out = val_true;
    else
        out = val_false;
    end
end
