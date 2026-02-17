function process_protocol2()
% PROCESS_PROTOCOL2  Main analysis entry point for Protocol 2 data.
%
%   PROCESS_PROTOCOL2() processes electrophysiology data recorded during
%   Protocol 2, performing both bar stimulus (direction selectivity) and
%   flash stimulus (receptive field) analyses.
%
%   IMPORTANT: Run this function from WITHIN the experiment directory.
%   Example: cd('C:\matlabroot\G4_Protocols\nested_RF_protocol2\data\2025_05_09_15_46')
%            process_protocol2()
%
%   PREREQUISITES:
%     - currentExp.mat must exist in the experiment folder (contains metadata)
%     - Log Files/ folder with G4_TDMS_Log*.mat files (voltage and frame data)
%     - Patterns/ and Functions/ folders from protocol generation
%
%   WORKFLOW:
%     1. Loads metadata (Frame, Age, Strain, Side) from currentExp.mat
%     2. Processes bar responses via PROCESS_BARS_P2:
%        - Extracts responses to 16 bar directions at 3 speeds (28, 56, 168 dps)
%        - Computes direction selectivity index (DSI)
%        - Calculates preferred direction via vector sum
%        - Generates polar plots and heatmaps
%     3. Processes flash responses via PROCESS_FLASH_P2:
%        - Extracts responses to 4px and 6px flash grids
%        - Fits 2D Gaussian model to receptive field
%        - Separates excitatory and inhibitory components
%        - Generates heatmaps and RF contour plots
%     4. Processes bar flash responses via PROCESS_BAR_FLASHES_P2:
%        - Extracts responses to bar flashes at 2 speeds (80ms and 14ms)
%        - 8 orientations x 11 positions per orientation
%
%   OUTPUT:
%     Creates in PROJECT_ROOT/results/:
%       - bar_results/bar_results_*.mat with DS metrics
%       - flash_results/rf_results_*.mat with RF parameters
%     Creates in PROJECT_ROOT/figures/:
%       - bar_stimuli/*.pdf with polar plots and timeseries
%       - flash_stimuli/*.pdf with RF heatmaps and Gaussian fits
%
%   See also PROCESS_BARS_P2, PROCESS_FLASH_P2, LOAD_PROTOCOL2_DATA

exp_folder = cd;

PROJECT_ROOT = "/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2";

% Load metadata - contains 'Frame', 'Age', 'Strain' and 'Side'. Side is the
% side of the arena upon which protocol 1 was run. 
load('currentExp.mat', 'metadata');

resultant_angle = process_bars_p2(exp_folder, metadata, PROJECT_ROOT);

process_flash_p2(exp_folder, metadata, PROJECT_ROOT, resultant_angle)

process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT)

cd(exp_folder)

end 

