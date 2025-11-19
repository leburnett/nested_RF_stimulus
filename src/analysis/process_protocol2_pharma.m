function process_protocol2_pharma()
% Run this function WITHIN the experiment directory. 
% i.e. 'C:\matlabroot\G4_Protocols\nested_RF_protocol2\data\2025_05_09_15_46'. 

exp_folder = cd;

PROJECT_ROOT = "C:\matlabroot\G4_Protocols\nested_RF_protocol2";
% PROJECT_ROOT = "/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2";

% Load metadata - contains 'Frame', 'Age', 'Strain' and 'Side'. Side is the
% side of the arena upon which protocol 1 was run. 
load('currentExp.mat', 'metadata');

process_bars_p2_pharma(exp_folder, metadata, PROJECT_ROOT);

% TBD - bar flash analysis 
% process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT)

cd(exp_folder)

end 

