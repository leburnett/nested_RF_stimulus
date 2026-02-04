function run_protocol2(exp_folder, pattern_order, func_order, trial_dur)
% RUN_PROTOCOL2  Execute Protocol 2 on the G4 LED arena.
%
%   RUN_PROTOCOL2(EXP_FOLDER, PATTERN_ORDER, FUNC_ORDER, TRIAL_DUR)
%   connects to the G4 Panels controller and runs the complete Protocol 2
%   experiment, presenting all patterns with their position functions and
%   logging voltage and frame timing data.
%
%   INPUTS:
%     exp_folder    - Full path to experiment directory containing
%                     currentExp.mat, Patterns/, and Functions/
%     pattern_order - 1xN array of pattern IDs in presentation order
%     func_order    - 1xN array of function IDs matching pattern_order
%     trial_dur     - 1xN array of trial durations in seconds
%
%   EXPERIMENT EXECUTION:
%     - Runs 3 repetitions of all conditions (n_reps = 3)
%     - Control mode 1 (position function control)
%     - Logs data to TDMS files in 'Log Files/' subdirectory
%     - Displays progress to command window
%     - Presents grey screen at end of experiment
%     - Converts TDMS files to .mat format after completion
%
%   PREREQUISITES:
%     - G4 Display Tools installed (PanelsController class available)
%     - userSettings.m configured with correct paths
%     - Arena connected and powered on
%     - No existing log files in experiment folder
%
%   OUTPUT:
%     Creates in exp_folder/Log Files/:
%       - G4_TDMS_Log_*.mat files with voltage and frame timing data
%
%   DURATION:
%     Estimated time is printed to console before experiment starts.
%     Typical Protocol 2 duration: ~45-60 minutes for 3 repetitions.
%
%   See also CREATE_PROTOCOL2, GENERATE_PROTOCOL2, G4_TDMS_FOLDER2STRUCT

    %% Experiment metadata from user input:
    n_reps = 3;
    
    %% set up for experiment
    %Load configuration and start G4 Host
    % Check the use of this
    userSettings;
    
    load(fullfile(exp_folder,'currentExp.mat'));
    num_conditions = numel(pattern_order);
    log_folder = fullfile(exp_folder,'Log Files');
    if ~exist(log_folder,'dir')
        mkdir(log_folder);
    end
    
    %% Open new Panels controller instance
    ctlr = PanelsController();
    ctlr.open(true);
    
    %% Change to root directory
    ctlr.setRootDirectory(exp_folder);
    
    %check if log files already present for this experiment
    assert(~exist(fullfile(exp_folder, 'Log Files\*'),'file'),'unsorted log files present in save folder, remove before restarting experiment\n');
    % assert(~exist(fullfile(exp_folder, 'Results\', fly_name),'dir'),'Results folder already exists with that fly name\n');
    
    %finish setting up for experiment
    exp_seconds = currentExp.totalDuration;
    fprintf(['Estimated experiment duration: ' num2str((exp_seconds*n_reps)/60) ' minutes\n']);
    
    %% start experiment
    ctlr.startLog(); %starts logging data in .tdms files

    for r = 1:n_reps
        for c = 1:num_conditions
            %trial portion
            ctlr.setControlMode(1);
            ctlr.setPatternID(pattern_order(1,c));
            ctlr.setPatternFunctionID(func_order(1,c));
            trial_t = trial_dur(1, c);
            fprintf(['Rep ', num2str(r), ' of ', num2str(n_reps), ', cond ' num2str(c) ' of ' num2str(num_conditions) ': ' strjoin(currentExp.pattern.pattNames(pattern_order(1,c))) '\n']);
            ctlr.startDisplay(ceil(trial_t*10)-1); %duration expected in 100ms units
        end
    end 

    % ADD Greyscale at the end - only at the very end - not at the end of
    % every rep. 
    ctlr.setControlMode(1);
    ctlr.setPatternID(pattern_order(1,1));
    ctlr.setPatternFunctionID(func_order(1,1));
    trial_t = trial_dur(1, 1);
    fprintf('Greyscale interval at the end.');
    ctlr.startDisplay(ceil(trial_t*10)-1); %duration expected in 100ms units

    %rename/move results folder
    ctlr.stopLog('showTimeoutMessage', true);
    ctlr.stopDisplay()
    ctlr.close()
    disp('finished');

    % % Convert TDMS files to mat file - current issues.
    G4_TDMS_folder2struct(log_folder)

end 