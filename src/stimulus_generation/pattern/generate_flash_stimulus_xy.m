function generate_flash_stimulus_xy(x, y, px_intensity, px_crop, px_flash, on_off, flash_dur, int_dur, exp_folder)
% GENERATE_FLASH_STIMULUS_XY  Create flash patterns and functions for Protocol 2.
%
%   GENERATE_FLASH_STIMULUS_XY(X, Y, PX_INTENSITY, PX_CROP, PX_FLASH, ...
%       ON_OFF, FLASH_DUR, INT_DUR, EXP_FOLDER)
%   generates both pattern and position function files for flash stimuli
%   centered on the receptive field location identified in Protocol 1.
%
%   INPUTS:
%     x            - Horizontal pixel coordinate (column) of RF center
%     y            - Vertical pixel coordinate (row) of RF center
%     px_intensity - [bkg_color, off_color, on_color] intensity values (0-15)
%                    Example: [4, 0, 15] for medium gray background
%     px_crop      - Size of square region for stimulus presentation (pixels)
%                    Example: 30 for 30x30 pixel display area
%     px_flash     - Size of individual flash squares (pixels)
%                    Example: 4 for 4x4 pixel flashes
%     on_off       - Contrast to present: 'on' (bright), 'off' (dark), or 'both'
%     flash_dur    - Duration of each flash in seconds (e.g., 0.16 for 160ms)
%     int_dur      - Duration between flashes in seconds (e.g., 0.44 for 440ms)
%     exp_folder   - Experiment directory for saving files
%
%   PROTOCOL 2 FLASH GRID:
%     Uses 50% overlap between adjacent flashes to increase spatial resolution.
%     For 4px flashes in 30px area: creates 14x14 grid = 196 flashes
%     For 6px flashes in 33px area: creates 10x10 grid = 100 flashes
%
%   OUTPUT:
%     Creates files in exp_folder/Patterns/ and exp_folder/Functions/:
%       - Pattern files (.pat, .mat) with all flash positions
%       - Position function files (.pfn, .mat) with presentation timing
%
%   WORKFLOW:
%     1. Calculates display bounds using CENTEREDSQUARE
%     2. Configures parameters structure
%     3. Calls GENERATE_STIMULUS to create pattern and function files
%
%   See also GENERATE_STIMULUS, GENERATE_FLASH_PATTERN, CENTEREDSQUARE
    
    [disp_y1, disp_y2, disp_x1, disp_x2] = centeredSquare(x, y, px_crop);
    fprintf('Square bounds of area to present stimulus within [y1, y2, x1, x2]: [%d, %d, %d, %d]\n', disp_y1, disp_y2, disp_x1, disp_x2);
    
    % disp_x1 = 17;
    % disp_x2 = 46;
    % disp_y1 = 1;
    % disp_y2 = 30;
    params.x = x;
    params.y = y;
    params.px_intensity = px_intensity; % Intensity parameters of flash stimulus - [bkg_color, off_color, on_color]
    params.px_rng = [disp_y1, disp_y2, disp_x1, disp_x2]; % Pixel range of screen to present on - [row_start, row_end, col_start, col_end]
    params.flash_sz_px = px_flash; % Size of flash in pixels
    params.overlap = 0.5; % Overlap of flashes - between 0 and 1. 
    params.protocol = 'protocol2';
    
    params.bkg_frame = 1; % The background frame = frame 1 in the pattern:
    params.interval_dur = int_dur; % duration of interval background screen in seconds.
    params.flash_dur = flash_dur; % duration of flash in seconds.
    params.on_off = on_off;

    params.root_dir = exp_folder;
    
    generate_stimulus(params)

end 
    