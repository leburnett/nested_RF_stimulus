function generate_bar_patterns_xy(x, y, px_crop, px_intensity, patt_file, patName, save_dir)
% GENERATE_BAR_PATTERNS_XY  Crop full-field bar pattern to RF location.
%
%   GENERATE_BAR_PATTERNS_XY(X, Y, PX_CROP, PX_INTENSITY, PATT_FILE, ...
%       PATNAME, SAVE_DIR)
%   takes a pre-made full-field moving bar pattern and creates a cropped
%   version centered on the receptive field location [x,y].
%
%   INPUTS:
%     x            - Horizontal pixel coordinate (column) of RF center
%     y            - Vertical pixel coordinate (row) of RF center
%     px_crop      - Size of square display region in pixels (e.g., 30, 45)
%     px_intensity - [bkg_color, off_color, on_color] intensity values (0-15)
%     patt_file    - Full path to source full-field bar pattern .mat file
%     patName      - Name for the output pattern files
%     save_dir     - Directory to save cropped pattern files
%
%   ALGORITHM:
%     1. Loads full-field bar pattern (288 frames, full arena size)
%     2. Calculates display region centered on [x,y] using CENTEREDSQUARE
%     3. Creates new pattern array with gray background
%     4. Extracts central portion of each frame from source pattern
%     5. Places cropped bar at the [x,y] location in new pattern
%     6. Saves cropped pattern for G4 display
%
%   CROPPING LOGIC:
%     The function extracts a px_crop x px_crop region from the CENTER
%     of the full-field pattern (around arena center [24, 96]) and
%     repositions it at the RF location [x,y]. This ensures the bar
%     movement passes through the neuron's receptive field.
%
%   OUTPUT:
%     Creates in save_dir:
%       - <patName>.pat - Binary pattern for G4 controller
%       - <patName>.mat - MATLAB structure with pattern data
%
%   See also GENERATE_BAR_STIMULUS_XY, CENTEREDSQUARE, SAVE_PATTERN_G4
% _______________________________________________________________________
load(patt_file, 'pattern');

% Find region of the screen upon which to display the stimulus. 
[disp_y1, disp_y2, disp_x1, disp_x2] = centeredSquare(x, y, px_crop);

% Make empty pattern 
n_frames_pattern = 288;
bkg_color = px_intensity(1);
arena_px_h = 48;
arena_px_w = 192;

% Make empty pattern array purely with grey background.
Pats = ones(arena_px_h, arena_px_w, n_frames_pattern)*bkg_color;

% crop around the centre of the arena screen [24, 96];
arena_px_centre_h = arena_px_h/2;
arena_px_centre_w = arena_px_w/2;
% 45 pixel area total - therefore +/- 22 either side. 
% [2:46, 74:118]
px_side = (px_crop-1)/2;
crop_h_st = int8(arena_px_centre_h-px_side); 
crop_h_end = int8(arena_px_centre_h+px_side);
crop_w_st = int8(arena_px_centre_w-px_side);
crop_w_end = int8(arena_px_centre_w+px_side);
px_rng = [crop_h_st, crop_h_end, crop_w_st, crop_w_end];

% Fill in the location with the central part of pattern. 
for i = 1:n_frames_pattern
    Pats(int16(disp_y1):int16(disp_y2), int16(disp_x1):int16(disp_x2), i) = pattern.Pats(crop_h_st:crop_h_end, crop_w_st:crop_w_end, i);
end

%
param = pattern.param;
param.stretch = zeros(size(Pats, 3), 1);
param.gs_val = 4;
param.arena_pitch = 0;
param.px_rng = px_rng;

param.ID = get_pattern_ID(save_dir);
    
save_pattern_G4(Pats, param, save_dir, patName);
% Generate protocol from these and test! % % % % % %

% TESTS % % % % % % % 

% Visualise the pattern in order.
% 
% figure
% for i = 1:288 %n_frames
% aa = Pats(:, :, i);
% imagesc(aa)
% xlim([0 193])
% ylim([0 49])
% clim([0 15])
% pause(0.1)
% disp(i)
% end 

% Check bar is not visible in first frame after bkg
% figure; imagesc(Pats(:, :, 2))

end 




















