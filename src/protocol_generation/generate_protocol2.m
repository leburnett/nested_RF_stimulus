function generate_protocol2()
% GENERATE_PROTOCOL2  Main entry point for creating and running Protocol 2.
%
%   GENERATE_PROTOCOL2() creates a complete Protocol 2 experiment for
%   high-resolution receptive field (RF) and direction selectivity (DS)
%   mapping. This is the second phase of the two-protocol RF mapping system
%   for Drosophila T4/T5 motion-sensitive neurons using the G4 LED arena.
%
%   WORKFLOW:
%     1. Prompts user for input parameters (peak frame from Protocol 1,
%        arena side, fly age, strain)
%     2. Converts peak frame to screen coordinates [x,y] and determines
%        ON/OFF contrast preference
%     3. Creates experiment folder with timestamp
%     4. Generates 4px flash stimuli (196 flashes, 14x14 grid, 50% overlap)
%     5. Generates 6px flash stimuli (100 flashes, 10x10 grid, 50% overlap)
%     6. Generates cropped bar patterns centered on [x,y]
%     7. Creates position functions for flashes and moving bars
%     8. Assembles currentExp.mat with pattern/function ordering
%     9. Runs the protocol on the G4 arena
%    10. Displays inverse peak frame for opposite-contrast experiment
%
%   PREREQUISITES:
%     - Protocol 1 must have been run to identify the peak_frame
%     - G4 Display Tools must be installed and configured
%     - Arena must be connected and calibrated
%
%   OUTPUT:
%     Creates experiment folder at:
%       C:\matlabroot\G4_Protocols\nested_RF_protocol2\YYYY_MM_DD_HH_MM\
%     containing:
%       - Patterns/ folder with flash and bar patterns
%       - Functions/ folder with position functions
%       - currentExp.mat with experiment configuration
%
%   STIMULUS PARAMETERS:
%     Flash duration: 160ms
%     Inter-flash interval: 440ms (600ms total per flash)
%     4px flash area: 30px crop around center
%     6px flash area: 33px crop around center
%     Bar area: 30px crop around center
%     Pixel intensity: [bkg=4, off=0, on=15]
%
%   See also GET_INPUT_PARAMETERS, PATT_FRAME_TO_COORD,
%            GENERATE_FLASH_STIMULUS_XY, GENERATE_BAR_STIMULUS_XY,
%            CREATE_PROTOCOL2, RUN_PROTOCOL2

% _________________________________________________________________________
    
    metadata = get_input_parameters();
    peak_frame = metadata.Frame;
    screen_hemi = metadata.Side;

    px_intensity = [4, 0, 15]; % used to be [6 0 15]
    px_crop_flash = 30;
    px_crop_bar = 30;

    % Pixel limits of the entire screen - this will not change for screen_hemi:
    screen_width_start = 17;
    screen_width_end = 192;
    screen_height_start = 1;
    screen_height_end = 48;

    [x, y, on_off] = patt_frame_to_coord(peak_frame, screen_hemi);

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
            
        % 28 dps - 4 pixel flashes 
            flash_dur_slow = 0.16;
            int_dur_slow = 0.44; % total 600ms    % before, 0.34 - total = 500ms
            px_flash = 4;
            px_crop_flash = 30;
            generate_flash_stimulus_xy(x, y, px_intensity, px_crop_flash, px_flash, on_off, flash_dur_slow, int_dur_slow, exp_folder)
    
        % 28 dps - 6 pixel flashes
            flash_dur_slow = 0.16;
            int_dur_slow = 0.44; % total 600ms    % before, 0.34 - total = 500ms
            px_flash = 6;
            px_crop_flash = 33;
            generate_flash_stimulus_xy(x, y, px_intensity, px_crop_flash, px_flash, on_off, flash_dur_slow, int_dur_slow, exp_folder)
    
        % % 56 dps - faster
        %     flash_dur_fast = 0.08;
        %     int_dur_fast = 0.22; % total 300ms    % before, 0.17 - total = 250ms
        %     generate_flash_stimulus_xy(x, y, px_intensity, px_crop_flash, on_off, flash_dur_fast, int_dur_fast, exp_folder)
    
    % 3 - create patterns with bar stimulus centred on [x,y] 
            generate_bar_stimulus_xy(x, y, px_intensity, px_crop_bar, on_off, exp_folder)
    
    % 4 - generate cropped bar position functions.
            bar_pos_fn_dir = fullfile(exp_folder, 'Functions');
            generate_bar_pos_fns(bar_pos_fn_dir, px_crop_bar)
    
    % 5 - generate 'CurrentExp' from these components. 
           [pattern_order, func_order, trial_dur] = create_protocol2(exp_folder, metadata);

    % 6 - run the protocol:
            run_protocol2(exp_folder, pattern_order, func_order, trial_dur)

    % For the experiment design: 
    %       First, 4 pixel flashes, then bars.
    %       Runs through each pattern 'orientation' ON - forward and flip
    %       direction - then OFF pattern forward and flip - at one speed, then
    %       repeats through the two other speeds.

    % Find the peak_frame number to run the same protocol with the opposite
    % contrast stimuli. 
    find_inv_peak_frame(peak_frame);

end 

