function generate_flash_function(fl_rows, fl_cols, bkg_frame, interval_dur, flash_dur, func_save_dir)

% Generate position function for nested RF flash stimuli.

% 'fl_rows' = the number of rows in the grid of flashes.
% 'fl_cols' = the number of columns in the grid of flashes.
% 'bkg_frame' = the index in the pattern of the background frame. Should be
% the first frame. 
% 'interval_dur' = duration of interval background screen in seconds.
% 'flash_dur' = duration of flash in seconds.
% ______________________________________________________________________

    % Determine the order in which to display the flashes.
    flash_seq = generate_flash_order(fl_rows, fl_cols);
    
    func = generate_func_for_flash(bkg_frame, interval_dur, flash_seq, flash_dur);
    
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
    filename = strcat(string(flash_sz_px), 'px_flashes_', string(px_rng(3)), '_', string(px_rng(4)), '_', string(n_frames_func), 'frames_', string(n_flashes), 'flashes');
    
    save_function_G4(func, param, func_save_dir, filename);

end