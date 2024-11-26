%% For protocol 2 - 4 pix - 50% overlap 

ROOT_DIR = 'C:\matlabroot\GitHub\nested_RF_stimulus';

% Arena parameters - [n_rows, n_cols, n_pix_per_panel]
arena_size = [3, 12, 16];

% _________________________________________________________________________

% Generate pattern:

% Intensity parameters of flash stimulus - [bkg_color, off_color, on_color]
px_intensity = [6, 0, 15];

% For 30 x 30 pixel area 
disp_x1 = 17;
disp_x2 = 46;
disp_y1 = 1;
disp_y2 = 30;

% Pixel range of screen to present on - [row_start, row_end, col_start, col_end]
px_rng_to_use = [disp_y1, disp_y2, disp_x1, disp_x2];

% Size of flash in pixels:
flash_sz_px = 4;

% Overlap of grid
overlap = 0.5;

% String used for saving pattern: 
patName = strcat(string(flash_sz_px),'px_square_RF_ON_OFF_', string(px_rng_to_use(3)), '_', string(px_rng_to_use(4)), '_overlap', string(overlap*100));

% Directory to save pattern:
patt_save_dir = fullfile(ROOT_DIR, 'results\patterns');
if ~isfolder(patt_save_dir)
    mkdir(patt_save_dir);
end 

% Generate and save flash pattern:
[n_frames, fl_rows, fl_cols] = generate_flash_pattern(arena_size, px_intensity, px_rng_to_use, flash_sz_px, overlap, patName, patt_save_dir);

% _________________________________________________________________________

% Generate position function:

% Determine the order in which to display the flashes.
% screen_size = [px_rng_to_use(2)-px_rng_to_use(1)+1, px_rng_to_use(4)-px_rng_to_use(3)+1];
flash_seq = generate_flash_order(fl_rows, fl_cols);

% Generate 'func' - linear array of what frame to present every 2ms.
bkg_frame = 1; % The background frame = frame 1 in the pattern:
bkg_dur = 0.05; % duration of interval background screen in seconds.
flash_dur = 0.2; % duration of flash in seconds.

func = generate_func_for_flash(bkg_frame, bkg_dur, flash_seq, flash_dur);

% Output checking metrics:
n_frames_func = numel(func);
disp(strcat("Number of elements in 'func': ", string(n_frames_func)))

n_unique_frame_pos = numel(unique(func));
disp(strcat("Number of unique frames in 'func': ", string(n_unique_frame_pos)))

func_save_dir = fullfile(ROOT_DIR, 'results\functions');
if ~isfolder(func_save_dir)
    mkdir(func_save_dir);
end 

% Generate 'param' struct for making position function:
param.func = func;
function_type = 'pfn';
ID = get_function_ID(function_type, func_save_dir);
param.ID = ID;
param.frames = 65; % n_frames; % number of frames in the pattern.
param.gs_val = 4;
param.type = function_type;

n_flashes = fl_rows*fl_cols;
% String to use for function name.
filename = strcat(string(flash_sz_px), 'px_flashes_', string(px_rng_to_use(3)), '_', string(px_rng_to_use(4)), '_', string(n_frames_func), 'frames_', string(n_flashes), 'flashes');

save_function_G4(func, param, func_save_dir, filename);










