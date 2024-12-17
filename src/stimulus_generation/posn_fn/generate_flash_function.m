function generate_flash_function(flash_sz_px, fl_rows, fl_cols, n_frames, bkg_frame, interval_dur, flash_dur, on_off, func_save_dir)

% Generate position function for nested RF flash stimuli.

% 'fl_rows' = the number of rows in the grid of flashes.
% 'fl_cols' = the number of columns in the grid of flashes.
% 'bkg_frame' = the index in the pattern of the background frame. Should be
% the first frame. 
% 'interval_dur' = duration of interval background screen in seconds.
% 'flash_dur' = duration of flash in seconds.
% 'on_off' = whether to present bright "on", dark "off" or both "both"
% flashes.
% ______________________________________________________________________

    % Determine the order in which to display the flashes.
    flash_seq = generate_flash_order(fl_rows, fl_cols);
    
    func = generate_func_for_flash(bkg_frame, interval_dur, flash_seq, flash_dur, on_off);
    
    % Output checking metrics:
    n_frames_func = numel(func);
    disp(strcat("Number of elements in 'func': ", string(n_frames_func)))
    
    n_unique_frame_pos = numel(unique(func));
    disp(strcat("Number of unique frames in 'func': ", string(n_unique_frame_pos)))

    % Generate 'param' struct for making position function:
    param.func = func;
    function_type = 'pfn';
    ID = get_function_ID(function_type, func_save_dir);
    param.ID = ID;

    param.frames = n_frames; % number of frames in the pattern.
    param.gs_val = 4;
    param.type = function_type;
    
    % String to use for function name.
    n_flashes = fl_rows*fl_cols;

    % add duration of the stimulus
    param.dur = (numel(func)*2)/1000; % func consists of 1 element every 2ms. n_elements*2 = time in ms. Divide by 100 to get time in seconds. %(n_flashes*flash_dur + (n_flashes)*interval_dur); 

    filename = strcat(string(flash_sz_px), 'px_flashes_', string(n_flashes), 'flashes_', string(flash_dur*1000), 'ms_', string(interval_dur*1000), 'ms_', on_off);
    
    save_function_G4(func, param, func_save_dir, filename);

end