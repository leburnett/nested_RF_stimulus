%% Generate nested flash protocol for determining the receptive field of T4 / T5 cells
function generate_stimulus(params)
    
    ROOT_DIR =  params.root_dir; %'C:\matlabroot\GitHub\nested_RF_stimulus';
    
    %% Generate pattern:

    px_intensity = params.px_intensity;

    % Range of pixels of which to display the stimulus
    px_rng = params.px_rng;
    % px_rng_formatted = arrayfun(@(x) sprintf('%02d', x), px_rng, 'UniformOutput', false);

    protocol = params.protocol;
    if protocol == "protocol1" 
        patt_save_dir = fullfile(ROOT_DIR, 'results', 'Patterns', protocol);
    elseif protocol == "protocol2" 
        patt_save_dir = fullfile(ROOT_DIR, 'Patterns');
    end

    flash_sz_px = params.flash_sz_px;
    overlap = params.overlap;

    % String used for saving pattern: 
    patName = strcat(string(flash_sz_px),'px_square_RF_ON_OFF_', string(px_rng(3)), '_', string(px_rng(4)), '_overlap', string(overlap*100));
    
    % stim_name = strcat('flash_', string(flash_sz_px));

    % Directory to save pattern:
    if ~isfolder(patt_save_dir)
        mkdir(patt_save_dir);
    end 
    
    % Generate and save flash pattern:
    [n_frames, fl_rows, fl_cols] = generate_flash_pattern(px_intensity, px_rng, flash_sz_px, overlap, patName, patt_save_dir);
    
    % _________________________________________________________________________
    
    %% Generate position function:
  
    bkg_frame = params.bkg_frame; 
    interval_dur = params.interval_dur; 
    flash_dur = params.flash_dur;

    if protocol == "protocol1" 
        func_save_dir = fullfile(ROOT_DIR, 'results', 'Functions', protocol);
    elseif protocol == "protocol2" 
        func_save_dir = fullfile(ROOT_DIR, 'Functions');
    end
    if ~isfolder(func_save_dir)
        mkdir(func_save_dir);
    end 
    
    generate_flash_function(flash_sz_px, fl_rows, fl_cols, n_frames, bkg_frame, interval_dur, flash_dur, func_save_dir);
    
    % _________________________________________________________________________
    
    %% Save parameters
    params.n_frames = n_frames; 
    params.patName = patName;
    params.fl_rows = fl_rows; 
    params.fl_cols = fl_cols;
   
    % Directory to save parameters:
    if protocol == "protocol1" 
        params_save_dir = fullfile(ROOT_DIR, 'results', 'params', protocol);
    elseif protocol == "protocol2" 
        params_save_dir = fullfile(ROOT_DIR, 'params');
    end
    if ~isfolder(params_save_dir)
        mkdir(params_save_dir);
    end 
    
    save(fullfile(params_save_dir, patName), 'params');

end 
%%

% time to run:
% If every element is 2ms
% run_t_s = (n_frames_func*2)/1000;
