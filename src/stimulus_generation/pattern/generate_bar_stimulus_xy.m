function generate_bar_stimulus_xy(x, y, px_intensity, px_crop, on_off, exp_folder)
% Generate inset bar patterns centred on [X, Y]
    
    % Make patterns of bars moving within 45px areas.
    % x = 1;
    % y = 1;
    % px_intensity = [6, 0, 15];
    % px_crop = 45;
    % coord_str = strcat(string(x), '_', string(y));
    
    save_dir = strcat(exp_folder, '\patterns');
    if ~isfolder(save_dir)
        mkdir(save_dir);
    end
    
    % Path to the folder with the full field bar stimuli.
    full_bar_patterns_dir = 'C:\matlabroot\GitHub\nested_RF_stimulus\results\Patterns\protocol2\full_field_bars';
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