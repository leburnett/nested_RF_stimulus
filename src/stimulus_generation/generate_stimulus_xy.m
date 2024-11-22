function generate_stimulus_xy(x, y)
% Generate the patterns and functions for the 4 x 4 pixel flashes - on a
% 50% overlapping grid of size 30 x 30 pixels, centred on the 'peak' of 
% responses to Protocol 1 [x,y]. 
    
    [disp_y1, disp_y2, disp_x1, disp_x2] = centeredSquare(x, y, 30);
    
    % disp_x1 = 17;
    % disp_x2 = 46;
    % disp_y1 = 1;
    % disp_y2 = 30;
    params.x = x;
    params.y = y;
    params.px_intensity = [6, 0, 15]; % Intensity parameters of flash stimulus - [bkg_color, off_color, on_color]
    params.px_rng = [disp_y1, disp_y2, disp_x1, disp_x2]; % Pixel range of screen to present on - [row_start, row_end, col_start, col_end]
    params.flash_sz_px = 4; % Size of flash in pixels
    params.overlap = 0.5; % Overlap of flashes - between 0 and 1. 
    
    params.bkg_frame = 1; % The background frame = frame 1 in the pattern:
    params.interval_dur = 0.05; % duration of interval background screen in seconds.
    params.flash_dur = 0.2; % duration of flash in seconds.
    
    generate_stimulus(params)
    
end 
    
function [row_start, row_end, col_start, col_end] = centeredSquare(x, y, square_size)
    % Screen dimensions
    screen_width_start = 17;
    screen_width_end = 196;
    screen_height_start = 1;
    screen_height_end = 48;

    % Square size and half-size for centering
    % square_size = 30;
    half_size = square_size / 2;

    % Calculate row and column start and end positions
    row_start = max(screen_height_start, min(y - half_size + 1, screen_height_end - square_size + 1));
    row_end = row_start + square_size - 1;

    col_start = max(screen_width_start, min(x - half_size + 1, screen_width_end - square_size + 1));
    col_end = col_start + square_size - 1;

    % Output the square boundaries
    fprintf('Square bounds: [%d, %d, %d, %d]\n', row_start, row_end, col_start, col_end);
end