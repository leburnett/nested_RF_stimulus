function generate_protocol2(peak_frame)
% Takes an array containing the frames of the patterns displayed in
% protocol 1 to which the cell had the strongest voltage response. It then
% uses these frames to determine the pixel on the screen on which to centre
% the stimuli that are presented in this following protocol. 

% This script calls the functions that generate the 4 pixel, 50%
% overlapping flash stimulus and also the moving bar stimulus.

% It creates the position functions for these stimuli and runs the
% protocol. 

% Inputs:
% ______
% 'peak_frames' - [n_condition, n_rep, [peak_frame, peak_voltage]]

% _________________________________________________________________________
    px_intensity = [6, 0, 15];
    px_crop_flash = 30;
    px_crop_bar = 45;

    % Pixel limits of the screen:
    screen_width_start = 17;
    screen_width_end = 192;
    screen_height_start = 1;
    screen_height_end = 48;

    [x, y, on_off] = patt_frame_to_coord(peak_frame, px_intensity(1));

    % Warning message if [x,y] is close to the edge of the screen.
    if x < screen_width_start+(px_crop_flash/2) || x > screen_width_end-(px_crop_flash/2)
        warning('x coordinate is close to the edge of the screen. The stimulus will not be centred on the x coordinate.')
    end 

    if y < screen_height_start+(px_crop_flash/2)
        warning('y coordinate is close to the bottom edge of the screen. The stimulus will not be centred on the y coordinate, consider moving the screen up.')
    end 

    if y > screen_height_end-(px_crop_flash/2)
        warning('y coordinate is close to the top edge of the screen. The stimulus will not be centred on the y coordinate, consider moving the screen down.')
    end
    
    % 1 - create experiment folder for protocol 2 
            exp_path = 'C:\matlabroot\G4_Protocols\nested_RF_protocol2';
            exp_name = string(datetime('now','TimeZone','local','Format','yyyy_MM_dd_HH_mm'));
            exp_folder = create_exp_dir_G4(exp_name, exp_path);
    
    % 2 - create 4 pixel flash stimuli with 50% overlapping grid, centred on [X,Y]
    % - this makes both the patterns and the functions for the flash.         
            generate_flash_stimulus_xy(x, y, px_intensity, px_crop_flash, on_off, exp_folder)
    
    % 3 - create patterns with bar stimulus centred on [x,y] 
            generate_bar_stimulus_xy(x, y, px_intensity, px_crop_bar, on_off, exp_folder)
    
    % 4 - generate cropped bar position functions.
            bar_pos_fn_dir = fullfile(exp_folder, 'Functions');
            generate_bar_pos_fns(bar_pos_fn_dir)
    
    % 5 - generate 'CurrentExp' from these components. 
           [pattern_order, func_order, trial_dur] = create_protocol2(exp_folder);

    % 6 - run the protocol:
            run_protocol2(exp_folder, pattern_order, func_order, trial_dur)

    % For the experiment design: 
    %       First, 4 pixel flashes, then bars.
    %       Runs through each pattern 'orientation' ON - forward and flip
    %       direction - then OFF pattern forward and flip - at one speed, then
    %       repeats through the two other speeds.
end 

