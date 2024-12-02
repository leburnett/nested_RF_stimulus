function generate_protocol2_stimuli(peak_frames)
    % check that [x,y] from lisa will be in the same register. 1,1 etc and 17,1.
    % x = 1;
    % y = 1;
    % [1,1] is the bottom left pixel on the screen from the fly's view.
    % x and y are the coordinates of the pixel of the screen with the
    % highest response. 
% _________________________________________________________________________
    px_intensity = [6, 0 , 15];
    px_crop_flash = 30;
    px_crop_bar = 45;

    % Pixel limits of the screen:
    screen_width_start = 17;
    screen_width_end = 192;
    screen_height_start = 1;
    screen_height_end = 48;

    x,y = frame_to_coord(peak_frames);

    % Warning message if [x,y] is close to the edge of the screen.
    if x < screen_width_start+px_crop_flash || x > screen_width_end-px_crop_flash
        warning('x coordinate is close to the edge of the screen. The stimulus will not be centred on the x coordinate.')
    end 

    if y < screen_height_start+px_crop_flash
        warning('y coordinate is close to the bottom edge of the screen. The stimulus will not be centred on the y coordinate, consider moving the screen up.')
    end 

    if y > screen_height_end-px_crop_flash
        warning('y coordinate is close to the top edge of the screen. The stimulus will not be centred on the y coordinate, consider moving the screen down.')
    end
    
    % 1 - create experiment folder for protocol 2 
            exp_path = 'C:\matlabroot\nested_RF_protocols\protocol2';
            exp_name = string(datetime('now','TimeZone','local','Format','yyyy_MM_dd_HH_mm'));
            exp_folder = create_exp_dir_G4(exp_name, exp_path);
    
    % 2 - create 4 pixel flash stimuli with 50% overlapping grid, centred on [X,Y]
    % - this makes both the patterns and the functions for the flash.         
            generate_flash_stimulus_xy(x, y, px_intensity, px_crop_flash, exp_folder)
    
    % 3 - create patterns with bar stimulus centred on [x,y] 
            generate_bar_stimulus_xy(x, y, px_intensity, px_crop_bar, exp_folder)
    
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

