function generate_flash_stimulus_xy(x, y, px_intensity, px_crop, exp_folder)
% Generate the patterns and functions for the 4 x 4 pixel flashes - on a
% 50% overlapping grid of size 30 x 30 pixels, centred on the 'peak' of 
% responses to Protocol 1 [x,y]. 
    
    [disp_y1, disp_y2, disp_x1, disp_x2] = centeredSquare(x, y, px_crop);
    
    % disp_x1 = 17;
    % disp_x2 = 46;
    % disp_y1 = 1;
    % disp_y2 = 30;
    params.x = x;
    params.y = y;
    params.px_intensity = px_intensity; % Intensity parameters of flash stimulus - [bkg_color, off_color, on_color]
    params.px_rng = [disp_y1, disp_y2, disp_x1, disp_x2]; % Pixel range of screen to present on - [row_start, row_end, col_start, col_end]
    params.flash_sz_px = 4; % Size of flash in pixels
    params.overlap = 0.5; % Overlap of flashes - between 0 and 1. 
    params.protocol = 'protocol2';
    
    params.bkg_frame = 1; % The background frame = frame 1 in the pattern:
    params.interval_dur = 0.05; % duration of interval background screen in seconds.
    params.flash_dur = 0.2; % duration of flash in seconds.

    params.root_dir = exp_folder;
    
    generate_stimulus(params)

end 
    