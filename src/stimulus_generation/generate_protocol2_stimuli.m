function generate_protocol2_stimuli(x, y)
    % INPUTS = [X,Y]
    
    % check that [x,y] from lisa will be in the same register. 1,1 etc and 17,1.
    % x = 1;
    % y = 1;
    px_intensity = [6, 0 , 15];
    px_crop_flash = 30;
    px_crop_bar = 45;
    
    % For making bar protocol:
    exp_path = 'C:\matlabroot\nested_RF_protocols\protocol2';
    exp_name = string(datetime('now','TimeZone','local','Format','yyyy_MM_dd_HH_mm'));
    
    % 1 - create experiment folder for protocol 2 
            exp_folder = create_exp_dir_G4(exp_name, exp_path);
    
    % 2 - create 4 pixel flash stimuli with 50% overlapping grid, centred on [X,Y]
    % - this makes both the patterns and the functions for the flash.         
            generate_flash_stimulus_xy(x, y, px_intensity, px_crop_flash, exp_folder)
    
    % 3 - create patterns with bar stimulus centred on [x,y] 
            generate_bar_stimulus_xy(x, y, px_intensity, px_crop_bar, exp_folder)
    
    % 4 - functions are already made - these should always stay the same.
    % functions for the 3 different speeds and back and forth.
    % bar_function_folder = 'C:\matlabroot\GitHub\nested_RF_stimulus\results\Functions\protocol2\bar_45px';
            bar_pos_fn_dir = fullfile(exp_folder, 'Functions');
            
    % 5 - copy bar position functions into experiment folder. 
            generate_bar_pos_fns(bar_pos_fn_dir)
    
    % 5 - generate protocol from these components. 
            create_protocol2(exp_folder)
    % start with the 4 px - 50% overlap - centred on the [X,Y]
    
    %       16 patterns generated within the structure: 
    %       results/ patterns/ protocol2/bar_"45" / "X_Y"
    
    % For the experiment design: 
    %       Runs through each pattern 'orientation' ON - forward and flip
    %       direction - then OFF pattern forward and flip - at one speed, then
    %       repeats through the two other speeds.
    
    %       14 dps
    %       0001 pattern - 0001 function
    %       0001 pattern - 0002 function
    %       0002 pattern - 0001 function
    %       0002 pattern - 0002 function
    %       ...
    %       28 dps
    %       0001 pattern - 0003 function
    %       0001 pattern - 0004 function
    %       0002 pattern - 0003 function
    %       0002 pattern - 0004 function
    %       ...
    %       56 dps
    %       0001 pattern - 0005 function
    %       0001 pattern - 0006 function
    %       0002 pattern - 0005 function
    %       0002 pattern - 0006 function
end 

