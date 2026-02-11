function create_bar_sweep_gif(pattern_file, is_flip, save_folder, date_str, time_str, pd_angle_rad)
    % Create an animated GIF of the bar sweep stimulus for the preferred
    % direction. Loads the pattern .mat file, expands to 192x48, applies a
    % black-and-green colormap, and saves as a looping GIF.
    %
    % Inputs
    % ------
    %   pattern_file : str
    %       Full path to the pattern .mat file for the PD bar orientation.
    %
    %   is_flip : logical
    %       If true, play frames in reverse order (flip/null direction).
    %
    %   save_folder : str
    %       Directory in which to save the GIF.
    %
    %   date_str : str
    %       Date string for filename (YYYY_MM_DD).
    %
    %   time_str : str
    %       Time string for filename (HH_MM).
    %
    %   pd_angle_rad : float
    %       Preferred direction angle in radians (used in filename).

    % Load the pattern
    loaded = load(pattern_file, 'pattern');
    Pats = loaded.pattern.Pats;

    % Frame range for bar sweeps (px_crop=30 protocol)
    frame_start = 11;
    frame_end = 62;

    if is_flip
        frame_range = frame_end:-1:frame_start;
    else
        frame_range = frame_start:frame_end;
    end

    % Handle different pattern storage formats
    [dim1, dim2, n_frames] = size(Pats);

    if dim1 == 48 && dim2 == 192
        % Standard format: 48 x 192 x N - use as-is
        expand_needed = false;
    elseif dim1 == 192 && dim2 == 3
        % Compressed G4 format: 192 x 3 x N - expand height from 3 to 48
        expand_needed = true;
        expand_dim = 'height';
    elseif dim1 == 3 && dim2 == 192
        % Compressed G4 format: 3 x 192 x N - expand height from 3 to 48
        expand_needed = true;
        expand_dim = 'height_transposed';
    else
        warning('Unexpected pattern dimensions: %d x %d x %d. Attempting to use as-is.', dim1, dim2, n_frames);
        expand_needed = false;
    end

    % Black-to-green colormap (gs_val=4 -> 16 levels, pixel values 0-15)
    n_levels = 16;
    green_cmap = zeros(n_levels, 3);
    for k = 1:n_levels
        green_cmap(k, :) = [0, (k-1)/(n_levels-1), 0];
    end

    % Build GIF filename
    gif_filename = fullfile(save_folder, sprintf('%s_%s_bar_sweep_PD_%ddeg.gif', ...
        date_str, time_str, round(rad2deg(pd_angle_rad))));

    delay_time = 0.05; % 50ms between frames

    for idx = 1:numel(frame_range)
        frame = Pats(:, :, frame_range(idx));

        % Expand if needed to reach 48 x 192
        if expand_needed
            if strcmp(expand_dim, 'height')
                % 192 x 3 -> 192 x 48: replicate each column 16 times
                frame = repelem(frame, 1, 16); % 192 x 48
                frame = frame'; % transpose to 48 x 192
            elseif strcmp(expand_dim, 'height_transposed')
                % 3 x 192 -> 48 x 192: replicate each row 16 times
                frame = repelem(frame, 16, 1); % 48 x 192
            end
        end

        % Convert pixel values (0-15) to 1-indexed for colormap
        frame_idx = uint8(frame) + 1;

        % Convert to RGB using green colormap
        rgb_frame = ind2rgb(frame_idx, green_cmap);
        rgb_uint8 = uint8(rgb_frame * 255);

        % Scale up for better visibility (each pixel -> 4x4 block)
        rgb_uint8 = repelem(rgb_uint8, 4, 4, 1);

        if idx == 1
            imwrite(rgb_uint8, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', delay_time);
        else
            imwrite(rgb_uint8, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
    end

    fprintf('  GIF saved: %s\n', gif_filename);

end
