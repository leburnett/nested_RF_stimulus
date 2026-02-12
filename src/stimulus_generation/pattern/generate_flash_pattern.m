function [n_frames, fl_rows, fl_cols] = generate_flash_pattern(px_intensity, px_rng, flash_sz_px, overlap, patName, save_dir)
% GENERATE_FLASH_PATTERN  Create spatial pattern array for flash stimuli.
%
%   [N_FRAMES, FL_ROWS, FL_COLS] = GENERATE_FLASH_PATTERN(PX_INTENSITY, ...
%       PX_RNG, FLASH_SZ_PX, OVERLAP, PATNAME, SAVE_DIR)
%   generates a G4 pattern file containing a grid of flash stimuli, with
%   both ON (bright) and OFF (dark) versions.
%
%   INPUTS:
%     px_intensity - [bkg_color, off_color, on_color] intensity values (0-15)
%                    Example: [4, 0, 15] for medium gray background
%     px_rng       - [row_start, row_end, col_start, col_end] display region
%                    Example: [1, 48, 17, 112] for left half of arena
%     flash_sz_px  - Size of each flash square in pixels (e.g., 4, 6, 12)
%     overlap      - Fraction of overlap between adjacent flashes (0-1)
%                    0 = no overlap, 0.5 = 50% overlap
%     patName      - String name for saving pattern files
%     save_dir     - Directory path to save pattern files
%
%   OUTPUTS:
%     n_frames - Total number of frames in the pattern
%                (1 background + n_flashes OFF + n_flashes ON)
%     fl_rows  - Number of flash rows in the grid
%     fl_cols  - Number of flash columns in the grid
%
%   PATTERN STRUCTURE:
%     Frame 1:           Uniform background (bkg_color)
%     Frames 2 to N+1:   OFF flashes (off_color on bkg_color)
%     Frames N+2 to 2N+1: ON flashes (on_color on bkg_color)
%
%   FRAME ORDERING:
%     Flashes are ordered column-major (down columns first):
%     Frame 2 = position [1,1], Frame 3 = position [2,1], etc.
%     This ordering is compatible with sub2ind() for position function
%     generation.
%
%   OUTPUT FILES:
%     - <patName>.pat - Binary pattern file for G4 controller
%     - <patName>.mat - MATLAB structure with pattern and parameters
%
%   EXAMPLE:
%     % Create 6px flash grid with no overlap for left arena half
%     [n, rows, cols] = generate_flash_pattern([4,0,15], [1,48,17,112], ...
%                                              6, 0, 'test_pattern', './Patterns');
%
%   See also GENERATE_STIMULUS, GENERATE_FLASH_FUNCTION, SAVE_PATTERN_G4
% ______________________________________________________________________

    % Arena parameters - [n_rows, n_cols, n_pix_per_panel]
    arena_size = [3, 12, 16];    

    % Arena specs:
    n_rows = arena_size(1);
    n_cols = arena_size(2);  
    px_per_panel = arena_size(3);
    
    h_display = n_rows*px_per_panel;
    w_display = n_cols*px_per_panel; 
    
    % Light intensity values for the pattern:
    bkg_color = px_intensity(1);
    off_color = px_intensity(2);
    on_color = px_intensity(3);
    
    % Subregion of screen over which to present the flashes.
    start_pixel_w = px_rng(3);
    end_pixel_w = px_rng(4);
    
    start_pixel_h = px_rng(1);
    end_pixel_h = px_rng(2);
    
    % Stimulus size in pixels.
    flash_px = flash_sz_px;

    grid_step = flash_px * (1-overlap);
    
    % Find the specs of the grid of flshes.
    edge_st_w = start_pixel_w:grid_step:(end_pixel_w - flash_px +1);
    edge_end_w = edge_st_w + flash_px - 1;
    fl_cols = numel(edge_st_w);
    
    edge_st_h = start_pixel_h:grid_step:(end_pixel_h - flash_px +1);
    edge_end_h = edge_st_h + flash_px - 1;
    
    fl_rows = numel(edge_st_h);
    n_flashes = fl_rows*fl_cols;
    disp(strcat("Flashes of size ", string(flash_px), "px with a ", string(overlap*100), "% overlap form a ", string(fl_rows), "x", string(fl_cols), " grid of ", string(n_flashes), " flashes."))
    
    %% Generate the pattern: 
    
    % First frame = background. 
    Bkg = ones(h_display, w_display, 1)*bkg_color;
    
    % Start by making every frame a uniform gray background. 
    OFF_Pats = ones(h_display, w_display, n_flashes)*bkg_color;
    ON_Pats = ones(h_display, w_display, n_flashes)*bkg_color;
    
    idx = 1;
    % Start with the OFF flashes
    for c = 1:fl_cols
        y_st = edge_st_w(c);
        y_end = edge_end_w(c);
    
        for r = 1:fl_rows
            x_st = edge_st_h(r);
            x_end = edge_end_h(r);
        
            % Add flash
            OFF_Pats(x_st:x_end, y_st:y_end, idx) = off_color;
            ON_Pats(x_st:x_end, y_st:y_end, idx) = on_color;
        
            % Update frame number 
            idx = idx+1;
    
        end 
    end 
    
    Pats = cat(3, Bkg, OFF_Pats, ON_Pats);
    
    % Check if there are other parameters that are required.
    param.stretch = zeros(size(Pats, 3), 1);
    param.gs_val = 4;
    param.arena_pitch = 0;
    param.px_rng = px_rng;
    
    param.ID = get_pattern_ID(save_dir);
    
    save_pattern_G4(Pats, param, save_dir, patName);
    
    % Total number of frames in the pattern - both ON and OFF.
    n_frames = size(Pats, 3);
    
    % TEST % % % % % % % 
    % Visualise the pattern in order.
    % 
    % figure
    % for i = 1:n_frames
    % idx = flash_seq(i)+1;
    % aa = Pats(:, :, idx);
    % imagesc(aa)
    % pause(0.2)
    % end 
    
    % % 

end 