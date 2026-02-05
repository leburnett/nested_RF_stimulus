function generate_protocol2_stimuli(peak_frame, on_off)
% GENERATE_PROTOCOL2_STIMULI  Legacy function for Protocol 2 generation.
%
%   GENERATE_PROTOCOL2_STIMULI(PEAK_FRAME, ON_OFF) generates and runs
%   Protocol 2 stimuli centered on the identified RF location.
%
%   NOTE: This is an earlier version of the protocol generation workflow.
%   For current usage, use GENERATE_PROTOCOL2 instead, which includes
%   improved parameter handling and metadata collection via GUI.
%
%   INPUTS:
%     peak_frame - Frame number from Protocol 1 with strongest response
%     on_off     - Contrast preference:
%                  1 or 'on'  = bright (ON) flashes
%                  0 or 'off' = dark (OFF) flashes
%
%   WORKFLOW:
%     1. Converts peak frame to screen coordinates [x,y]
%     2. Creates experiment folder with timestamp
%     3. Generates 4px overlapping flash stimuli (30px crop area)
%     4. Generates cropped bar stimuli (45px crop area)
%     5. Creates position functions for bars
%     6. Assembles and runs the protocol
%
%   DIFFERENCES FROM GENERATE_PROTOCOL2:
%     - No GUI for metadata input
%     - Different crop sizes (45px vs 30px for bars)
%     - Different intensity values (6,0,15 vs 4,0,15)
%     - No 6px flash grid
%
%   See also GENERATE_PROTOCOL2, PATT_FRAME_TO_COORD
% _________________________________________________________________________
    px_intensity = [6, 0, 15];
    px_crop_flash = 30;
    px_crop_bar = 45;

    % Pixel limits of the screen:
    screen_width_start = 17;
    screen_width_end = 192;
    screen_height_start = 1;
    screen_height_end = 48;

    if on_off == 1
        on_off = 'on';
    elseif on_off == 0
        on_off = 'off';
    end 

    [x, y] = patt_frame_to_coord(peak_frame, px_intensity(1));

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

