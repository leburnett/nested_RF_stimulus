function generate_stimulus_gif(pattern_file, function_file, output_filename, options)
% GENERATE_STIMULUS_GIF Creates an animated GIF from pattern and function files
%
% Inputs:
%   pattern_file    - Full path to .mat file containing 'pattern' struct with 'Pats' field
%   function_file   - Full path to .mat file containing 'pfnparam' struct with 'func' field
%   output_filename - Output GIF filename (e.g., 'stimulus.gif') OR full path
%                     If just a filename, saves in same directory as pattern_file
%   options         - (optional) Struct with fields:
%       .fps              - Frame rate for GIF (default: 10)
%       .fps_initial      - Initial playback speed in fps (default: same as fps)
%       .fps_fast         - Accelerated playback speed in fps (default: same as fps)
%       .switch_time_step - Time step to switch from initial to fast speed (default: inf)
%       .max_time_steps   - Maximum time steps to include (default: all)
%       .time_step_skip   - Only write every Nth time step (default: 1)
%       .colormap         - 'green' or 'gray' (default: 'green')
%       .crop_outline     - [x y width height] to draw rectangle overlay (default: [])
%       .show_title       - Show frame info title (default: false)
%       .crop_region      - [x1 x2 y1 y2] to crop display to this region (default: [])
%       .show_full_arena  - Show full arena with crop_region highlighted (default: false)
%       .speed_label      - Show speed label on frames (default: true)
%
% Example:
%   % Basic usage (same as before)
%   generate_stimulus_gif(pattern_file, function_file, 'my_stimulus.gif');
%
%   % With options
%   opts.fps_initial = 5;
%   opts.fps_fast = 50;
%   opts.switch_time_step = 100;
%   opts.max_time_steps = 500;
%   opts.colormap = 'green';
%   opts.crop_outline = [40 10 30 30];
%   generate_stimulus_gif(pattern_file, function_file, 'stimulus.gif', opts);

    % Default options
    if nargin < 4
        options = struct();
    end

    % Set defaults
    if ~isfield(options, 'fps'), options.fps = 10; end
    if ~isfield(options, 'fps_initial'), options.fps_initial = options.fps; end
    if ~isfield(options, 'fps_fast'), options.fps_fast = options.fps; end
    if ~isfield(options, 'switch_time_step'), options.switch_time_step = inf; end
    if ~isfield(options, 'max_time_steps'), options.max_time_steps = inf; end
    if ~isfield(options, 'time_step_skip'), options.time_step_skip = 1; end
    if ~isfield(options, 'colormap'), options.colormap = 'green'; end
    if ~isfield(options, 'crop_outline'), options.crop_outline = []; end
    if ~isfield(options, 'show_title'), options.show_title = false; end
    if ~isfield(options, 'crop_region'), options.crop_region = []; end
    if ~isfield(options, 'show_full_arena'), options.show_full_arena = false; end
    if ~isfield(options, 'speed_label'), options.speed_label = true; end

    % Load the pattern file
    fprintf('Loading pattern file: %s\n', pattern_file);
    pattern_data = load(pattern_file);
    if ~isfield(pattern_data, 'pattern')
        error('Pattern file does not contain a "pattern" struct');
    end
    if ~isfield(pattern_data.pattern, 'Pats')
        error('Pattern struct does not contain a "Pats" field');
    end
    Pats = pattern_data.pattern.Pats;
    fprintf('Pattern size: [%d x %d x %d frames]\n', size(Pats, 1), size(Pats, 2), size(Pats, 3));

    % Load the function file
    fprintf('Loading function file: %s\n', function_file);
    function_data = load(function_file);
    if ~isfield(function_data, 'pfnparam')
        error('Function file does not contain a "pfnparam" struct');
    end
    if ~isfield(function_data.pfnparam, 'func')
        error('pfnparam struct does not contain a "func" field');
    end
    func = function_data.pfnparam.func;
    fprintf('Function length: %d time steps\n', length(func));

    % Determine output path
    [output_dir, ~, ~] = fileparts(output_filename);
    if isempty(output_dir)
        [pattern_dir, ~, ~] = fileparts(pattern_file);
        output_path = fullfile(pattern_dir, output_filename);
    else
        output_path = output_filename;
    end

    % Get intensity range for consistent scaling
    min_intensity = min(Pats(:));
    max_intensity = max(Pats(:));
    fprintf('Intensity range: [%d, %d]\n', min_intensity, max_intensity);

    % Create colormap
    % G4 arena uses 4-bit intensity (0-15), with typical values:
    %   0 = OFF (black)
    %   4 = background (dark green)
    %   15 = ON (bright green)
    if strcmp(options.colormap, 'green')
        % Create a 16-level green colormap for G4 arena (4-bit = 0-15)
        % Map intensity levels to appropriate green brightness
        n_levels = 16;
        green_map = zeros(n_levels, 3);
        % Green channel intensity scales with level, but not linearly
        % Level 0 = black, Level 4 = dim green (background), Level 15 = bright green
        for i = 1:n_levels
            level = i - 1;  % 0 to 15
            % Use a mapping that makes background (level 4) visibly dark green
            % and ON flashes (level 15) bright green
            green_map(i, 2) = (level / 15)^0.7;  % Slight gamma for better visibility
        end
        cmap = green_map;
    else
        cmap = gray(16);
    end

    % Limit time steps if requested
    total_time_steps = min(length(func), options.max_time_steps);

    % Determine which time steps to write (subsampling)
    time_steps_to_write = 1:options.time_step_skip:total_time_steps;

    fprintf('Generating GIF...\n');
    fprintf('Total time steps: %d, Writing: %d frames\n', total_time_steps, length(time_steps_to_write));

    % Create figure for rendering
    fig = figure('Visible', 'off', 'Color', 'k');

    % Determine display dimensions
    if ~isempty(options.crop_region)
        x1 = options.crop_region(1);
        x2 = options.crop_region(2);
        y1 = options.crop_region(3);
        y2 = options.crop_region(4);
        display_height = y2 - y1 + 1;
        display_width = x2 - x1 + 1;
    else
        display_height = size(Pats, 1);
        display_width = size(Pats, 2);
    end

    % Set figure size based on display dimensions (scale up for visibility)
    scale_factor = max(4, floor(400 / display_height));
    fig_width = display_width * scale_factor;
    fig_height = display_height * scale_factor;

    if options.show_title || options.speed_label
        fig_height = fig_height + 40;  % Extra space for text
    end

    set(fig, 'Position', [100 100 fig_width fig_height]);

    ax = axes('Parent', fig);
    colormap(ax, cmap);

    % Progress tracking
    progress_step = max(1, floor(length(time_steps_to_write) / 20));
    frame_count = 0;

    for idx = 1:length(time_steps_to_write)
        t = time_steps_to_write(idx);
        frame_count = frame_count + 1;

        % Get the frame index from the function
        frame_idx = func(t);

        % Validate frame index
        if frame_idx < 1 || frame_idx > size(Pats, 3)
            continue;
        end

        % Extract the current frame
        current_frame = double(Pats(:, :, frame_idx));

        % Apply crop if requested
        if ~isempty(options.crop_region) && ~options.show_full_arena
            x1 = options.crop_region(1);
            x2 = options.crop_region(2);
            y1 = options.crop_region(3);
            y2 = options.crop_region(4);
            current_frame = current_frame(y1:y2, x1:x2);
        end

        % Display the frame
        imagesc(ax, current_frame);
        caxis(ax, [0, 15]);  % G4 arena uses 4-bit intensity (0-15)
        axis(ax, 'equal', 'tight', 'off');
        set(ax, 'Position', [0 0 1 1]);

        % Draw crop outline if requested (for full arena view)
        if ~isempty(options.crop_outline) && (isempty(options.crop_region) || options.show_full_arena)
            hold(ax, 'on');
            rect_x = options.crop_outline(1);
            rect_y = options.crop_outline(2);
            rect_w = options.crop_outline(3);
            rect_h = options.crop_outline(4);
            rectangle(ax, 'Position', [rect_x rect_y rect_w rect_h], ...
                      'EdgeColor', [0.2 0.6 1], 'LineWidth', 2, 'LineStyle', '-');
            hold(ax, 'off');
        end

        % Add title if requested
        if options.show_title
            title(ax, sprintf('Step: %d/%d, Frame: %d', t, total_time_steps, frame_idx), ...
                  'Color', 'w', 'FontSize', 10);
        end

        % Add speed label if requested
        if options.speed_label
            if t < options.switch_time_step
                speed_text = '1x speed';
            else
                speed_ratio = round(options.fps_fast / options.fps_initial);
                speed_text = sprintf('%dx speed', speed_ratio);
            end
            text(ax, 3, 3, speed_text, 'Color', 'w', 'FontSize', 12, ...
                 'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
                 'BackgroundColor', [0 0 0 0.5]);
        end

        % Determine delay time based on current position
        if t < options.switch_time_step
            delay_time = 1 / options.fps_initial;
        else
            delay_time = 1 / options.fps_fast;
        end

        % Capture the frame
        drawnow;
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        % Write to GIF
        if frame_count == 1
            imwrite(imind, cm, output_path, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(imind, cm, output_path, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end

        % Progress indicator
        if mod(idx, progress_step) == 0
            fprintf('Progress: %d%% (%d / %d frames)\n', ...
                    round(100 * idx / length(time_steps_to_write)), idx, length(time_steps_to_write));
        end
    end

    close(fig);

    fprintf('\nGIF saved successfully to:\n%s\n', output_path);
    fprintf('Total frames written: %d\n', frame_count);
end
