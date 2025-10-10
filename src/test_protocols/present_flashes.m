function present_flashes(peak_frame, screen_hemi, flash_size)

% THIS DOES NOT SAVE THE DATA

% Uses the patterns generated for Protocol 1 
% After running Protocol 1 the EphysGrid plots are generated.
% From these plots, use the numbers above each subplot to enter as
% 'peak_frame' and this function will present 5 flashes at that location. 

%% SET PARAMETERS
n_flashes = 3;
flash_dur = 1;
bkg_dur = 0.5;
n_reps = 1;

%% Get the patterns:

% Assuming using left at the moment:
if screen_hemi == "L"
    pattern_path = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\LHS\protocol1_10kHz_4reps_12px_6px_LHS_2sbkg_200msfl_50msint_12-13-24_14-33-03\Patterns';
elseif screen_hemi == "R"
    % pattern_path = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\RHS\protocol1_10kHz_4reps_12px_6px_RHS_2sbkg_200msfl_50msint_04-08-25_08-08-42\Patterns'; 
    pattern_path = "C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\RHS2\protocol1_10kHz_4reps_12px_6px_RHS2_2sbkg_200msfl_50msint_81_180_05-05-25_16-18-66\Patterns";
end

cd(pattern_path)

% Define which pattern to use (12px or 6 px flashes)
if flash_size == 12
    mat_files = dir('0001_*');
    mat_name = fullfile(mat_files.folder, mat_files.name);
    pat_files = dir('pat0001*');
    pat_name = fullfile(pat_files.folder, pat_files.name);
elseif flash_size == 6
    mat_files = dir('0002_*');
    mat_name = fullfile(mat_files.folder, mat_files.name);
    pat_files = dir('pat0002*');
    pat_name = fullfile(pat_files.folder, pat_files.name);
end 

% HARD CODE THE DIRECTORY WITH THE POSITION FUNCTION IN
exp_folder = "C:\matlabroot\G4_Protocols\nested_RF_stimulus\src\test_protocols\5flash_protocol";

pattern_folder = fullfile(exp_folder, "Patterns");
if ~exist(pattern_folder,'dir')
    mkdir(pattern_folder);
end 

%% Copy over the pattern files (.mat and .pat) to the "5flash_protocol" folder.
mat_name_new = "0001_pattern.mat";
pat_name_new = "pat0001.pat";

copyfile(mat_name, fullfile(pattern_folder, mat_name_new));
copyfile(pat_name, fullfile(pattern_folder, pat_name_new));

%% Generate and save the position function

function_folder = fullfile(exp_folder, "Functions");
if ~exist(function_folder,'dir')
    mkdir(function_folder);
end

total_dur = generate_5flash_function(peak_frame, n_flashes, flash_dur, bkg_dur, function_folder); 

%% Run the protocol

%Load configuration and start G4 Host
userSettings;
num_conditions = 1;
 
%% Open new Panels controller instance
ctlr = PanelsController();
ctlr.open(true);

%% Change to root directory
ctlr.setRootDirectory(exp_folder);

%% start experiment

for r = 1:n_reps
    for c = 1:num_conditions
        %trial portion
        ctlr.setControlMode(1);
        ctlr.setPatternID(1);
        ctlr.setPatternFunctionID(1);
        trial_t = total_dur;
        ctlr.startDisplay(ceil(trial_t*10)-1); %duration expected in 100ms units
    end
end 

ctlr.stopDisplay()
ctlr.close()
disp('finished');

end 