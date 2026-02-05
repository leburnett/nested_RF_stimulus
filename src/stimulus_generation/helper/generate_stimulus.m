function generate_stimulus(params)
% GENERATE_STIMULUS  Create flash pattern and position function for Protocol 1.
%
%   GENERATE_STIMULUS(PARAMS) generates both the visual pattern (.pat file)
%   and position function (.pfn file) for Protocol 1 flash stimuli. This is
%   the core stimulus generation function called by GENERATE_PROTOCOL1_STIMULI.
%
%   INPUT:
%     params - Structure containing stimulus parameters:
%       .root_dir      - Base directory for saving files
%       .protocol      - Protocol identifier string ('protocol1')
%       .px_intensity  - [bkg, off, on] pixel intensity values (0-15)
%       .px_rng        - [row_start, row_end, col_start, col_end] display region
%       .flash_sz_px   - Flash square size in pixels (e.g., 6, 12)
%       .overlap       - Flash overlap fraction (0-1, typically 0 for Protocol 1)
%       .bkg_frame     - Frame number for background (typically 1)
%       .interval_dur  - Duration of inter-flash interval in seconds
%       .flash_dur     - Duration of flash presentation in seconds
%       .on_off        - 'on', 'off', or 'both' flash polarities to generate
%
%   WORKFLOW:
%     1. Creates Patterns/ directory if needed
%     2. Calls GENERATE_FLASH_PATTERN to create spatial patterns
%     3. Creates Functions/ directory if needed
%     4. Calls GENERATE_FLASH_FUNCTION to create temporal sequence
%     5. Saves params structure to params/ directory for reference
%
%   OUTPUT FILES:
%     - Patterns/<patName>.pat  - Binary pattern file for G4 display
%     - Patterns/<patName>.mat  - MATLAB pattern structure
%     - Functions/<funcName>.pfn - Binary position function file
%     - Functions/<funcName>.mat - MATLAB function structure
%     - params/<patName>.mat    - Parameter structure for reference
%
%   See also GENERATE_PROTOCOL1_STIMULI, GENERATE_FLASH_PATTERN,
%            GENERATE_FLASH_FUNCTION
    
    ROOT_DIR =  params.root_dir; %'C:\matlabroot\GitHub\nested_RF_stimulus';
    
    %% Generate pattern:

    px_intensity = params.px_intensity;

    % Range of pixels of which to display the stimulus
    px_rng = params.px_rng;
    % px_rng_formatted = arrayfun(@(x) sprintf('%02d', x), px_rng, 'UniformOutput', false);

    flash_sz_px = params.flash_sz_px;
    overlap = params.overlap;

    % String used for saving pattern: 
    patName = strcat(string(flash_sz_px),'px_square_RF_ON_OFF_', string(px_rng(3)), '_', string(px_rng(4)), '_overlap', string(overlap*100));
    
    % stim_name = strcat('flash_', string(flash_sz_px));

    % Directory to save pattern:
    patt_save_dir = fullfile(ROOT_DIR, 'Patterns');
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

    func_save_dir = fullfile(ROOT_DIR, 'Functions');
    if ~isfolder(func_save_dir)
        mkdir(func_save_dir);
    end 

    on_off = params.on_off;
    
    generate_flash_function(flash_sz_px, fl_rows, fl_cols, n_frames, bkg_frame, interval_dur, flash_dur, on_off, func_save_dir);
    
    % _________________________________________________________________________
    
    %% Save parameters
    params.n_frames = n_frames; 
    params.patName = patName;
    params.fl_rows = fl_rows; 
    params.fl_cols = fl_cols;
   
    % Directory to save parameters:
    params_save_dir = fullfile(ROOT_DIR, 'params');
    if ~isfolder(params_save_dir)
        mkdir(params_save_dir);
    end 
    
    save(fullfile(params_save_dir, patName), 'params');

end 
%%

% time to run:
% If every element is 2ms
% run_t_s = (n_frames_func*2)/1000;
