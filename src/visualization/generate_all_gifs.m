% GENERATE_ALL_GIFS - Generate all GIFs for the nested RF protocol visualization
%
% This script generates the animated GIFs needed for the visualization webpage.
% Run this script from MATLAB after ensuring the paths are correct.
%
% Output files will be saved to: src/visualization/assets/

%% Setup paths
base_path = fileparts(mfilename('fullpath'));
assets_path = fullfile(base_path, 'assets');

% Ensure assets directory exists
if ~exist(assets_path, 'dir')
    mkdir(assets_path);
end

% Path to protocol files
repo_root = fileparts(fileparts(base_path));
protocol1_path = fullfile(repo_root, 'protocols', 'LHS', ...
    'protocol1_4reps_12px_6px_LHS_2sbkg_200msfl_50msint_12-03-24_15-11-40');

% Path to actual Protocol 2 data (from a real experiment run)
protocol2_data_path = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/RNAi_ttl_2025/control/ON/2024_12_18_15_07';

%% Parameters for demonstration
% Using center of the tested region as the example peak location
% For LHS protocol: x range is 17-112 (96 pixels wide), y range is 1-48
% Center would be approximately x=64, y=24
demo_peak_x = 64;  % Center of x range (17 + 95/2)
demo_peak_y = 24;  % Center of y range (1 + 47/2)

% Protocol 2 crop region (30x30 pixels centered on peak)
crop_size = 30;
crop_x1 = demo_peak_x - floor(crop_size/2);
crop_x2 = demo_peak_x + floor(crop_size/2) - 1;
crop_y1 = demo_peak_y - floor(crop_size/2);
crop_y2 = demo_peak_y + floor(crop_size/2) - 1;

fprintf('Demo peak location: [%d, %d]\n', demo_peak_x, demo_peak_y);
fprintf('Protocol 2 crop region: x=[%d,%d], y=[%d,%d]\n', crop_x1, crop_x2, crop_y1, crop_y2);

%% Generate Protocol 1 - 12px flashes GIF
fprintf('\n=== Generating Protocol 1 - 12px flashes ===\n');

pattern_file_12px = fullfile(protocol1_path, 'Patterns', '0001_12px_square_RF_ON_OFF_17_112_overlap0_G4.mat');
function_file_12px = fullfile(protocol1_path, 'Functions', '0001_12px_flashes_32flashes_200ms_50ms_both_G4.mat');
output_file_12px = fullfile(assets_path, 'protocol1_12px.gif');

opts_12px = struct();
opts_12px.colormap = 'green';
opts_12px.fps_initial = 5;           % Start at real-ish speed
opts_12px.fps_fast = 50;             % Then speed up
opts_12px.switch_time_step = 500;    % Switch after ~1 second of real time
opts_12px.max_time_steps = 4000;     % Limit total frames (show subset of protocol)
opts_12px.time_step_skip = 10;       % Subsample to reduce file size
opts_12px.speed_label = true;
opts_12px.show_title = false;
opts_12px.crop_outline = [crop_x1, crop_y1, crop_size, crop_size];  % Show where P2 will be

generate_stimulus_gif(pattern_file_12px, function_file_12px, output_file_12px, opts_12px);

%% Generate Protocol 1 - 6px flashes GIF
fprintf('\n=== Generating Protocol 1 - 6px flashes ===\n');

pattern_file_6px = fullfile(protocol1_path, 'Patterns', '0002_6px_square_RF_ON_OFF_17_112_overlap0_G4.mat');
function_file_6px = fullfile(protocol1_path, 'Functions', '0002_6px_flashes_128flashes_200ms_50ms_both_G4.mat');
output_file_6px = fullfile(assets_path, 'protocol1_6px.gif');

opts_6px = struct();
opts_6px.colormap = 'green';
opts_6px.fps_initial = 5;
opts_6px.fps_fast = 50;
opts_6px.switch_time_step = 500;
opts_6px.max_time_steps = 8000;      % More frames since there are more flashes
opts_6px.time_step_skip = 10;
opts_6px.speed_label = true;
opts_6px.show_title = false;
opts_6px.crop_outline = [crop_x1, crop_y1, crop_size, crop_size];

generate_stimulus_gif(pattern_file_6px, function_file_6px, output_file_6px, opts_6px);

%% Generate Protocol 2 region highlight GIF
% This shows the full arena with the Protocol 2 region outlined
fprintf('\n=== Generating Protocol 2 region highlight ===\n');

output_file_region = fullfile(assets_path, 'protocol2_region.gif');

opts_region = struct();
opts_region.colormap = 'green';
opts_region.fps_initial = 2;
opts_region.fps_fast = 2;
opts_region.switch_time_step = inf;
opts_region.max_time_steps = 100;    % Just a few frames to show the region
opts_region.time_step_skip = 25;
opts_region.speed_label = false;
opts_region.show_title = false;
opts_region.crop_outline = [crop_x1, crop_y1, crop_size, crop_size];

generate_stimulus_gif(pattern_file_12px, function_file_12px, output_file_region, opts_region);

%% Generate Protocol 2 flash GIFs using actual Protocol 2 patterns
fprintf('\n=== Generating Protocol 2 flash GIFs ===\n');

% Protocol 2 - 4px flashes (from actual Protocol 2 data)
p2_pattern_4px = fullfile(protocol2_data_path, 'Patterns', '0001_4px_square_RF_ON_OFF_66_95_overlap50_G4.mat');
p2_function_4px = fullfile(protocol2_data_path, 'Functions', '0001_4px_flashes_196flashes_160ms_340ms_on_G4.mat');

if exist(p2_pattern_4px, 'file') && exist(p2_function_4px, 'file')
    output_file_p2_4px = fullfile(assets_path, 'protocol2_4px.gif');

    opts_p2_4px = struct();
    opts_p2_4px.colormap = 'green';
    opts_p2_4px.fps_initial = 5;
    opts_p2_4px.fps_fast = 30;
    opts_p2_4px.switch_time_step = 300;
    opts_p2_4px.max_time_steps = 5000;
    opts_p2_4px.time_step_skip = 10;
    opts_p2_4px.speed_label = true;
    opts_p2_4px.show_title = false;
    % No crop needed - P2 patterns are already cropped to the relevant region

    generate_stimulus_gif(p2_pattern_4px, p2_function_4px, output_file_p2_4px, opts_p2_4px);
else
    fprintf('Protocol 2 4px flash files not found. Skipping.\n');
    fprintf('Expected pattern: %s\n', p2_pattern_4px);
    fprintf('Expected function: %s\n', p2_function_4px);
end

% For 6px flashes, check if there's a 6px pattern in Protocol 2
% If not, we can use the 4px pattern with different label or skip
p2_pattern_6px = fullfile(protocol2_data_path, 'Patterns', '0002_4px_square_RF_ON_OFF_66_95_overlap50_G4.mat');
p2_function_6px = fullfile(protocol2_data_path, 'Functions', '0002_4px_flashes_196flashes_80ms_170ms_on_G4.mat');

if exist(p2_pattern_6px, 'file') && exist(p2_function_6px, 'file')
    output_file_p2_6px = fullfile(assets_path, 'protocol2_6px.gif');

    opts_p2_6px = struct();
    opts_p2_6px.colormap = 'green';
    opts_p2_6px.fps_initial = 5;
    opts_p2_6px.fps_fast = 30;
    opts_p2_6px.switch_time_step = 300;
    opts_p2_6px.max_time_steps = 5000;
    opts_p2_6px.time_step_skip = 10;
    opts_p2_6px.speed_label = true;
    opts_p2_6px.show_title = false;

    generate_stimulus_gif(p2_pattern_6px, p2_function_6px, output_file_p2_6px, opts_p2_6px);
else
    fprintf('Protocol 2 6px flash files not found. Using 4px pattern as substitute.\n');
end

%% Generate Protocol 2 bars GIF using actual Protocol 2 bar patterns
fprintf('\n=== Generating Protocol 2 bars GIF ===\n');

% Use the actual cropped bar pattern from Protocol 2 data
% Using pattern 0003 which is a bright bar at pi orientation
p2_bar_pattern = fullfile(protocol2_data_path, 'Patterns', '0003_4pix_bar_15_6_pi_39shift_1-25step_288frames_CROP30_X80_Y40_G4.mat');
p2_bar_function = fullfile(protocol2_data_path, 'Functions', '0003_28dps_2-328s_11-62_288fr_1-25step_G4.mat');

if exist(p2_bar_pattern, 'file') && exist(p2_bar_function, 'file')
    output_file_bars = fullfile(assets_path, 'protocol2_bars.gif');

    opts_bars = struct();
    opts_bars.colormap = 'green';
    opts_bars.fps_initial = 10;
    opts_bars.fps_fast = 10;
    opts_bars.switch_time_step = inf;
    opts_bars.max_time_steps = inf;     % Show full bar sweep
    opts_bars.time_step_skip = 3;       % Subsample for file size
    opts_bars.speed_label = false;
    opts_bars.show_title = false;

    generate_stimulus_gif(p2_bar_pattern, p2_bar_function, output_file_bars, opts_bars);
else
    fprintf('Protocol 2 bar files not found. Skipping.\n');
    fprintf('Expected pattern: %s\n', p2_bar_pattern);
    fprintf('Expected function: %s\n', p2_bar_function);
end

%% Generate a composite bar GIF showing multiple orientations (optional)
% This cycles through several bar orientations to demonstrate the variety
fprintf('\n=== Generating multi-orientation bar GIF ===\n');

% List of bar patterns for different orientations (bright bars, patterns 3-10)
bar_orientations = {
    '0003_4pix_bar_15_6_pi_39shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',      % 0/8 pi
    '0004_4pix_bar_15_6_1-8pi_44shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',   % 1/8 pi
    '0005_4pix_bar_15_6_2-8pi_46shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',   % 2/8 pi
    '0006_4pix_bar_15_6_3-8pi_46shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',   % 3/8 pi
    '0007_4pix_bar_15_6_4-8pi_38shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',   % 4/8 pi
    '0008_4pix_bar_15_6_5-8pi_42shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',   % 5/8 pi
    '0009_4pix_bar_15_6_6-8pi_44shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',   % 6/8 pi
    '0010_4pix_bar_15_6_7-8pi_42shift_1-25step_288frames_CROP30_X80_Y40_G4.mat',   % 7/8 pi
};

% Generate individual bar GIFs for a few orientations
bar_function_28dps = fullfile(protocol2_data_path, 'Functions', '0003_28dps_2-328s_11-62_288fr_1-25step_G4.mat');

for i = 1:min(4, length(bar_orientations))  % Generate first 4 orientations
    bar_pat_file = fullfile(protocol2_data_path, 'Patterns', bar_orientations{i});
    if exist(bar_pat_file, 'file') && exist(bar_function_28dps, 'file')
        output_file = fullfile(assets_path, sprintf('protocol2_bar_orient%d.gif', i));

        opts_bar_i = struct();
        opts_bar_i.colormap = 'green';
        opts_bar_i.fps_initial = 15;
        opts_bar_i.fps_fast = 15;
        opts_bar_i.switch_time_step = inf;
        opts_bar_i.max_time_steps = inf;
        opts_bar_i.time_step_skip = 5;
        opts_bar_i.speed_label = false;
        opts_bar_i.show_title = false;

        fprintf('Generating bar orientation %d/%d...\n', i, min(4, length(bar_orientations)));
        generate_stimulus_gif(bar_pat_file, bar_function_28dps, output_file, opts_bar_i);
    end
end

%% Summary
fprintf('\n=== GIF Generation Complete ===\n');
fprintf('Output directory: %s\n', assets_path);
fprintf('\nGenerated files:\n');

generated_files = dir(fullfile(assets_path, '*.gif'));
for i = 1:length(generated_files)
    file_info = dir(fullfile(assets_path, generated_files(i).name));
    fprintf('  %s (%.1f KB)\n', generated_files(i).name, file_info.bytes/1024);
end

fprintf('\nNext steps:\n');
fprintf('1. Copy the provided images to the assets folder:\n');
fprintf('   - arena_setup.png\n');
fprintf('   - analysis_dark.png\n');
fprintf('   - analysis_light.png\n');
fprintf('2. Open index.html in Chrome to view the visualization\n');
