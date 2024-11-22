% Generate both the pattern and position function for RF flash stimuli - protocol 1. 
% Non-overlapping grid of flashes that cover ~1/2 of the arena.
% Initial protocol to find the centre of the receptive field.

%% Region of arena on which to present the stimulus. 

%  1 - For the 1/2 screen area. 6 panels. RHS of screen when look from the
%  front. LHS from fly's perspective. For recording from left hemisphere. 

disp_x1 = 17;
disp_x2 = 112;
disp_y1 = 1;
disp_y2 = 48;

%  2 - For the 1/2 screen area. 6 panels. LHS of screen when look from the
%  front. RHS from fly's perspective. For recording from right hemisphere. 

% disp_x1 = 101;
% disp_x2 = 196;
% disp_y1 = 1;
% disp_y2 = 48;

%% Other parameters which don't change with the area of the screen that the
% pattern is being presented. 

params.px_intensity = [6, 0, 15]; % Intensity parameters of flash stimulus - [bkg_color, off_color, on_color]
params.px_rng = [disp_y1, disp_y2, disp_x1, disp_x2]; % Pixel range of screen to present on - [row_start, row_end, col_start, col_end]
params.flash_sz_px = 6; % Size of flash in pixels
params.overlap = 0; % Overlap of flashes - between 0 and 1. 

params.bkg_frame = 1; % The background frame = frame 1 in the pattern:
params.interval_dur = 0.05; % duration of interval background screen in seconds.
params.flash_dur = 0.2; % duration of flash in seconds.

generate_stimulus(params)

