function analyze_single_experiment(exp_folder)
% ANALYZE_SINGLE_EXPERIMENT  Bar sweep and bar flash analysis for one experiment.
%
%   ANALYZE_SINGLE_EXPERIMENT(EXP_FOLDER) loads data from a single Protocol 2
%   experiment, creates a slow-speed (28 dps) bar sweep polar timeseries plot
%   with LUT-based direction verification, then creates a 1x11 bar flash
%   subplot along the PD-ND axis.
%
%   INPUT:
%     exp_folder - Full path to experiment directory
%                  e.g. '/Users/.../protocol2/data/1DRF/2025_11_10_10_17'
%
%   OUTPUTS:
%     Figure 1: LUT verification — bar orientation/direction visual check
%     Figure 2: Slow bar sweep polar timeseries (using LUT-derived directions)
%     Figure 3: 1x11 bar flash subplots along PD-ND axis
%
%   See also LOAD_PROTOCOL2_DATA, PARSE_BAR_DATA, PARSE_BAR_FLASH_DATA

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
    % bar_data: Nx4 cell. Rows 1-16 = slow (28 dps). Cols 1-3 = reps, col 4 = mean.

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
    assumed_angles = linspace(0, 360 - 22.5, 16)'; % 0, 22.5, 45, ..., 337.5

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

    %% Step 3b: Visual verification figure — orientation + direction arrows
    % Two side-by-side circular layouts:
    %   Left: subplots positioned by current plot_order (assumed angles)
    %   Right: subplots positioned by LUT directions (actual angles)

    fig_verify = figure('Name', 'LUT Verification: Bar Orientation & Direction');
    set(fig_verify, 'Position', [50 50 1600 800]);

    numPlots = 16;
    theta_plot = linspace(0, 2*pi, numPlots+1);
    theta_plot = theta_plot(1:end-1);

    radius = 0.32;
    sw = 0.08; sh = 0.08;

    for panel = 1:2
        if panel == 1
            centerX = 0.25;  % Left panel: plot_order assumed positions
            panel_title = 'Current plot\_order (assumed angles)';
        else
            centerX = 0.75;  % Right panel: LUT-based positions
            panel_title = 'LUT directions (actual angles)';
        end
        centerY = 0.5;

        for subplot_idx = 1:numPlots
            data_row = plot_order(subplot_idx);

            if panel == 1
                % Position at assumed angle (plot_order convention)
                angle_rad = theta_plot(subplot_idx);
            else
                % Position at actual LUT direction
                angle_rad = deg2rad(lut_directions(data_row));
            end

            x_pos = centerX + radius * cos(angle_rad);
            y_pos = centerY + radius * sin(angle_rad);

            ax = axes('Position', [x_pos - sw/2, y_pos - sh/2, sw, sh]);
            hold on;

            orient_deg = lut_orientations(data_row);
            dir_deg = lut_directions(data_row);

            % Draw bar orientation as a thick line segment
            orient_rad = deg2rad(orient_deg);
            bar_len = 0.35;
            bx = bar_len * cos(orient_rad);
            by = bar_len * sin(orient_rad);
            plot([-bx, bx], [-by, by], 'k-', 'LineWidth', 3);

            % Draw direction arrow
            dir_rad = deg2rad(dir_deg);
            arr_len = 0.3;
            ax_arr = arr_len * cos(dir_rad);
            ay_arr = arr_len * sin(dir_rad);
            quiver(0, 0, ax_arr, ay_arr, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.8);

            xlim([-0.5 0.5]); ylim([-0.5 0.5]);
            axis equal off;

            % Label with data row and actual direction
            text(0, -0.45, sprintf('R%d D:%.0f°', data_row, dir_deg), ...
                'HorizontalAlignment', 'center', 'FontSize', 5);
        end

        % Panel title
        annotation('textbox', [centerX-0.15, 0.88, 0.3, 0.05], ...
            'String', panel_title, ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
            'FontSize', 10, 'FontWeight', 'bold');

        % Center legend
        annotation('textbox', [centerX-0.1, 0.43, 0.2, 0.08], ...
            'String', {'Black = bar orientation', 'Red arrow = motion direction'}, ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 7);
    end

    sgtitle(sprintf('LUT Verification — %s %s — Black=orientation, Red=direction', ...
        strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' ')));

    %% Step 4: Slow-speed bar sweep polar timeseries plot

    % Use LUT directions (not assumed angles) for subplot positioning
    fig_polar = figure('Name', 'Slow Bar Sweep Polar Timeseries (28 dps)');

    % Reset layout parameters for this figure
    centerX = 0.5;
    centerY = 0.5;
    radius = 0.35;

    % Compute max_v using LUT-verified direction positions
    max_v_slow = zeros(numPlots, 1);

    col = [0.2 0.4 0.7]; % Dark blue for slow bars

    for subplot_idx = 1:numPlots
        data_row = plot_order(subplot_idx);

        % Position subplot at LUT direction angle
        actual_dir_rad = deg2rad(lut_directions(data_row));
        x_pos = centerX + radius * cos(actual_dir_rad);
        y_pos = centerY + radius * sin(actual_dir_rad);

        sw = 0.15; sh = 0.15;
        ax = axes('Position', [x_pos - sw/2, y_pos - sh/2, sw, sh]);
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

    % Add central polar plot
    centralSize = (2 * radius) * 0.65;
    centralPosition = [centerX - centralSize/2, centerY - centralSize/2, centralSize, centralSize];
    angls_polar = deg2rad(lut_directions(plot_order));
    angls_polar = [angls_polar; angls_polar(1)]; % close the loop
    max_v_polar = [max_v_slow; max_v_slow(1)];

    axCentral = polaraxes('Position', centralPosition);
    hold on;
    polarplot(angls_polar, max_v_polar, 'Color', col, 'LineWidth', 2);

    set(gcf, 'Position', [303 78 961 969]);
    sgtitle(sprintf('28 dps — %s — %s — %s', ...
        strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' '), params.on_off));

    %% Step 5: Find Preferred Direction (PD)

    [~, pd_subplot_idx] = max(max_v_slow);
    pd_data_row = plot_order(pd_subplot_idx);
    pd_direction = lut_directions(pd_data_row);
    pd_orientation = lut_orientations(pd_data_row);
    pd_pattern = lut_patterns(pd_data_row);
    pd_function = lut_functions(pd_data_row);

    fprintf('\n=== Preferred Direction (PD) ===\n');
    fprintf('PD direction:   %.1f°\n', pd_direction);
    fprintf('Bar orientation: %.1f°\n', pd_orientation);
    fprintf('Pattern:        %d\n', pd_pattern);
    fprintf('Function:       %d\n', pd_function);
    fprintf('Data row:       %d\n', pd_data_row);

    %% Step 6: Map PD orientation to bar flash column

    % Experiment patterns 3-10 → full-field patterns 1-8
    % Bar flash column in parse_bar_flash_data output = full-field pattern index
    pd_ff_pattern = pd_pattern - 2;  % full-field pattern index (1-8)
    bar_flash_col = pd_ff_pattern;   % column in data_slow (11×8×3)

    fprintf('Full-field pattern index: %d\n', pd_ff_pattern);
    fprintf('Bar flash column:         %d\n', bar_flash_col);

    %% Step 7: Parse bar flash data and extract PD-ND axis

    [data_slow_bf, ~, mean_slow_bf, ~] = parse_bar_flash_data(f_data, v_data);

    % data_slow_bf: 11×8×3 cell (positions × orientations × reps)
    % Extract the column for PD orientation: 11 positions × 3 reps
    pd_flash_data = data_slow_bf(:, bar_flash_col, :);  % 11×1×3
    pd_flash_mean = mean_slow_bf(:, bar_flash_col);      % 11×1

    %% Step 8: Determine spatial ordering (ND→PD, left to right)

    % For forward functions (3, 5, 7): lower frame = earlier in sweep
    %   Lower position index in bar flash = earlier frame = bar approaching from one side
    %   For forward sweep, position 1 = ND side, position 11 = PD side → order 1:11
    % For reverse/flip functions (4, 6, 8): the sweep direction is flipped
    %   Position 1 in the bar flash pattern still references the same physical frame from the
    %   original (non-flipped) pattern. But the flip function reverses the playback direction.
    %   So position 1 = what was the end of the forward sweep = PD side in the flipped case
    %   → order 11:-1:1 to get ND→PD

    is_forward = mod(pd_function, 2) == 1;  % odd functions = forward
    if is_forward
        pos_order = 1:11;  % ND→PD for forward sweep
    else
        pos_order = 11:-1:1;  % reverse to get ND→PD for flipped sweep
    end

    fprintf('PD function %d is %s → position order: %s\n', ...
        pd_function, ternary(is_forward, 'forward', 'reverse'), mat2str(pos_order));

    %% Step 9: Create 1×11 bar flash subplot

    fig_flash = figure('Name', 'Bar Flash PD-ND Axis');
    tiledlayout(1, 11, 'TileSpacing', 'compact', 'Padding', 'compact');

    for pos_idx = 1:11
        flash_pos = pos_order(pos_idx);
        nexttile;
        hold on;

        % Plot 3 repetitions in light grey
        for r = 1:3
            ts = pd_flash_data{flash_pos, 1, r};
            if ~isempty(ts)
                plot(1:numel(ts), ts, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            end
        end

        % Plot mean in thick black
        ts_mean = pd_flash_mean{flash_pos};
        if ~isempty(ts_mean)
            plot(1:numel(ts_mean), ts_mean, 'k', 'LineWidth', 1.5);
        end

        % Formatting
        ylim([-80 -10]);
        if pos_idx == 1
            ylabel('mV');
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

    sgtitle(sprintf('Bar Flash PD-ND Axis — Dir:%.0f° Orient:%.0f° — %s — %s', ...
        pd_direction, pd_orientation, strrep(date_str,'_','-'), strrep(metadata.Strain,'_',' ')));

    set(gcf, 'Position', [50 400 1800 300]);

    %% Save figures as PDFs
    save_dir = fullfile(exp_folder, 'analysis_output');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end

    exportgraphics(fig_verify, fullfile(save_dir, 'lut_verification.pdf'), ...
        'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_polar, fullfile(save_dir, 'slow_bar_sweep_polar.pdf'), ...
        'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_flash, fullfile(save_dir, 'bar_flash_PD_ND_axis.pdf'), ...
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
