function generate_flash_function(flash_sz_px, fl_rows, fl_cols, n_frames, bkg_frame, interval_dur, flash_dur, on_off, func_save_dir)
% GENERATE_FLASH_FUNCTION  Create complete position function for flash stimuli.
%
%   GENERATE_FLASH_FUNCTION(FLASH_SZ_PX, FL_ROWS, FL_COLS, N_FRAMES, ...
%       BKG_FRAME, INTERVAL_DUR, FLASH_DUR, ON_OFF, FUNC_SAVE_DIR)
%   generates and saves the position function file that controls the
%   temporal sequence of flash presentations for the G4 display system.
%
%   INPUTS:
%     flash_sz_px   - Size of each flash square in pixels (for filename)
%     fl_rows       - Number of rows in the flash grid
%     fl_cols       - Number of columns in the flash grid
%     n_frames      - Total number of frames in the associated pattern
%     bkg_frame     - Frame index for background (typically 1)
%     interval_dur  - Duration between flashes in seconds
%     flash_dur     - Duration of each flash in seconds
%     on_off        - Flash polarity: 'on', 'off', or 'both'
%     func_save_dir - Directory to save position function files
%
%   WORKFLOW:
%     1. Calls GENERATE_FLASH_ORDER to determine spatial presentation order
%     2. Calls GENERATE_FUNC_FOR_FLASH to build the timing array
%     3. Constructs parameter structure with metadata
%     4. Saves function files using SAVE_FUNCTION_G4
%
%   OUTPUT FILES:
%     - <flash_sz>px_flashes_<n>flashes_<dur>ms_<int>ms_<on_off>.pfn
%     - <flash_sz>px_flashes_<n>flashes_<dur>ms_<int>ms_<on_off>_G4.mat
%
%   CONSOLE OUTPUT:
%     Prints number of elements in func and number of unique frames
%     for verification.
%
%   TOTAL DURATION:
%     Calculated as: (n_elements * 2ms) / 1000 = seconds
%     Includes initial background + all flashes + inter-flash intervals
%
%   See also GENERATE_FLASH_ORDER, GENERATE_FUNC_FOR_FLASH,
%            GENERATE_FLASH_PATTERN, SAVE_FUNCTION_G4
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