function generate_bar_stimulus_xy(x, y, px_intensity, px_crop, on_off, exp_folder)
% GENERATE_BAR_STIMULUS_XY  Create cropped moving bar patterns for Protocol 2.
%
%   GENERATE_BAR_STIMULUS_XY(X, Y, PX_INTENSITY, PX_CROP, ON_OFF, EXP_FOLDER)
%   generates bar stimulus patterns cropped to a region centered on [x,y].
%   Uses pre-made full-field bar patterns and crops them to the RF location.
%
%   INPUTS:
%     x            - Horizontal pixel coordinate (column) of RF center
%     y            - Vertical pixel coordinate (row) of RF center
%     px_intensity - [bkg_color, off_color, on_color] intensity values (0-15)
%     px_crop      - Size of square region to display bars within (pixels)
%     on_off       - Contrast preference: 'on', 'off', or 'both'
%     exp_folder   - Experiment directory path for saving patterns
%
%   BAR ORIENTATIONS:
%     Generates 8 or 16 patterns depending on on_off setting:
%     - Patterns 1-8:  ON (bright) bars moving in 8 directions
%     - Patterns 9-16: OFF (dark) bars moving in 8 directions
%     Directions: 0, 45, 90, 135, 180, 225, 270, 315 degrees
%
%   SOURCE PATTERNS:
%     Loads full-field bar patterns from:
%     C:\matlabroot\G4_Protocols\nested_RF_stimulus\results\patterns\
%       protocol2\full_field_bars4\
%     These patterns contain 288 frames of bar movement across the full arena.
%
%   OUTPUT:
%     Creates pattern files in exp_folder/Patterns/:
%       - One .pat and .mat file per bar direction/contrast combination
%       - Filenames include crop size and center coordinates
%
%   See also GENERATE_BAR_PATTERNS_XY, GENERATE_PROTOCOL2, CENTEREDSQUARE
    
    save_dir = strcat(exp_folder, '\patterns');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    
    % Path to the folder with the full field bar stimuli.
    full_bar_patterns_dir = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\results\patterns\protocol2\full_field_bars4'; % July 2025 - changed to 4bkg
    cd(full_bar_patterns_dir)
    
    bar_patts = dir('*.mat');

    if on_off == "on"
        patt_ids = 1:1:8;
    elseif on_off == "off"
        patt_ids = 9:1:16;
    elseif on_off == "both"
        patt_ids = 1:1:16;
    end 
    
    % Generate 'cropped' bar stimuli centred on the [x,y] coordinates.
    for p = patt_ids
        fname = bar_patts(p).name;
        patt_file = fullfile(bar_patts(p).folder, fname);
        patName = strcat(fname(6:end-14), 'CROP', string(px_crop), '_X', string(x), '_Y', string(y));
        generate_bar_patterns_xy(x, y, px_crop, px_intensity, patt_file, patName, save_dir)
    end 
end 