function result = process_1DRF_single(exp_folder, figures_root)
    % Process a single 1DRF experiment: run bar sweep analysis, determine
    % preferred direction, create stimulus GIF, and plot bar flash spatial
    % profile along the PD axis.
    %
    % Inputs
    % ------
    %   exp_folder : str
    %       Full path to the experiment data folder.
    %
    %   figures_root : str
    %       Root directory for saving figures. A subfolder named by
    %       date_time will be created inside.
    %
    % Outputs
    % -------
    %   result : struct
    %       Summary of analysis results for this experiment.

    original_dir = cd;
    cleanup = onCleanup(@() cd(original_dir));

    %% A: Load data
    loaded_exp = load(fullfile(exp_folder, 'currentExp.mat'));
    if isfield(loaded_exp, 'metadata')
        metadata = loaded_exp.metadata;
    else
        warning('No metadata found for %s. Using defaults.', exp_folder);
        metadata.Strain = 'unknown';
        metadata.Frame = '';
        metadata.Age = '';
        metadata.Side = '';
    end

    [date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);
    on_off = params.on_off;
    params.date = date_str;
    params.time = time_str;
    params.strain = metadata.Strain;

    f_data = Log.ADC.Volts(1, :); % frame data
    v_data = Log.ADC.Volts(2, :) * 10; % voltage data
    median_v = median(v_data);

    fprintf('Processing: %s_%s (%s, %s)\n', date_str, time_str, metadata.Strain, on_off);

    %% B: Parse bar sweep data
    data = parse_bar_data(f_data, v_data);

    %% C: Create per-cell figure folder
    cell_folder_name = strcat(date_str, '_', time_str);
    cell_fig_folder = fullfile(figures_root, cell_folder_name);
    if ~isfolder(cell_fig_folder)
        mkdir(cell_fig_folder);
    end

    %% D: Plot bar sweep polar figure
    save_fig = 0;
    [max_v, min_v] = plot_timeseries_polar_bars(data, median_v, params, save_fig, cell_fig_folder);
    f = gcf;
    f.Position = [303, 380, 688, 667];
    close(f);

    %% E: Determine preferred direction
    % max_v is 16 x 3 (directions x speeds). PD = direction/speed with
    % highest 98th percentile response.
    [~, linear_idx] = max(max_v(:));
    [pd_dir_idx, pd_speed_idx] = ind2sub(size(max_v), linear_idx);

    angls = linspace(0, 2*pi, 17);
    angls = angls(1:16);
    pd_angle_rad = angls(pd_dir_idx);

    speed_labels = {'28dps', '56dps', '168dps'};
    fprintf('  PD: %.1f deg (direction %d, speed %s, max_v = %.2f)\n', ...
        rad2deg(pd_angle_rad), pd_dir_idx, speed_labels{pd_speed_idx}, max_v(pd_dir_idx, pd_speed_idx));

    % Map to data row for pattern identification
    plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    pd_data_row = plot_order(pd_dir_idx);

    %% F: Identify pattern file for GIF
    % Pattern files 0003-0010 correspond to the 8 bar orientations.
    % Data rows 1-16 map in pairs (odd=forward, even=flip) to patterns 3-10.
    pattern_file_idx = ceil(pd_data_row / 2) + 2;
    is_flip = mod(pd_data_row, 2) == 0;

    pattern_folder = fullfile(exp_folder, 'Patterns');
    pattern_files = dir(fullfile(pattern_folder, '*.mat'));
    [~, sort_idx] = sort({pattern_files.name});
    pattern_files = pattern_files(sort_idx);

    if pattern_file_idx <= numel(pattern_files)
        pd_pattern_file = fullfile(pattern_folder, pattern_files(pattern_file_idx).name);
        fprintf('  Pattern file: %s (flip=%d)\n', pattern_files(pattern_file_idx).name, is_flip);

        %% G: Create GIF (wrapped in try/catch so failure doesn't block remaining analysis)
        try
            create_bar_sweep_gif(pd_pattern_file, is_flip, cell_fig_folder, date_str, time_str, pd_angle_rad);
        catch ME
            warning('GIF creation failed for %s_%s: %s. Continuing with remaining analysis.', ...
                date_str, time_str, ME.message);
        end
    else
        warning('Pattern file index %d exceeds available patterns (%d). Skipping GIF.', ...
            pattern_file_idx, numel(pattern_files));
    end

    %% H: Parse bar flash data
    [data_slow, ~, mean_slow, ~] = parse_bar_flash_data(f_data, v_data);

    %% I: Map PD to bar flash column
    [bf_col, positions_order] = map_PD_to_bar_flash_column(pd_angle_rad);

    %% J: Plot bar flash PD positions (raw + baseline-subtracted, closes figures internally)
    plot_bar_flash_PD_positions(data_slow, mean_slow, bf_col, positions_order, ...
        median_v, params, cell_fig_folder, pd_angle_rad, exp_folder);

    %% K: Compute additional metrics for summary
    max_v_polar_slow = vertcat(max_v(:, 1), max_v(1, 1));
    [d_slow, ~, magnitude_slow, angle_rad_slow, fwhm_slow, cv_slow, ~, ~] = find_PD_and_order_idx(max_v_polar_slow, median_v);
    [sym_ratio_slow, DSI_slow, DSI_pdnd_slow, ~] = compute_bar_response_metrics(d_slow);

    %% L: Package results
    result = struct();
    result.date = date_str;
    result.time = time_str;
    result.strain = metadata.Strain;
    result.on_off = on_off;
    result.pd_angle_rad = pd_angle_rad;
    result.pd_angle_deg = rad2deg(pd_angle_rad);
    result.pd_speed_idx = pd_speed_idx;
    result.pd_speed_label = speed_labels{pd_speed_idx};
    result.pd_max_v = max_v(pd_dir_idx, pd_speed_idx);
    result.max_v = max_v;
    result.min_v = min_v;
    result.median_v = median_v;
    result.bf_col = bf_col;
    result.positions_order = positions_order;
    result.vector_sum_magnitude_slow = magnitude_slow;
    result.vector_sum_angle_slow = angle_rad_slow;
    result.DSI_slow = DSI_slow;
    result.DSI_pdnd_slow = DSI_pdnd_slow;
    result.sym_ratio_slow = sym_ratio_slow;
    result.fwhm_slow = fwhm_slow;
    result.cv_slow = cv_slow;

    fprintf('  Done: %s_%s\n\n', date_str, time_str);

end
