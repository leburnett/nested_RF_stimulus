function run_protocol2(exp_folder, pattern_order, func_order, trial_dur)
% Run protocol 2 
    %% Experiment metadata from user input:
    metadata = get_input_parameters();

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
    fprintf(['Estimated experiment duration: ' num2str(exp_seconds/60) ' minutes\n']);
    
    %% start experiment
    ctlr.startLog(); %starts logging data in .tdms files
    
    for c = 1:num_conditions
        %trial portion
        ctlr.setControlMode(1);
        ctlr.setPatternID(pattern_order(1,c));
        ctlr.setPatternFunctionID(func_order(1,c));
        trial_t = trial_dur(1, c);
        fprintf(['Cond ' num2str(c) ' of ' num2str(num_conditions) ': ' strjoin(currentExp.pattern.pattNames(pattern_order(1,c))) '\n']);
        ctlr.startDisplay(ceil(trial_t*10)-1); %duration expected in 100ms units
    end

    %rename/move results folder
    ctlr.stopLog('showTimeoutMessage', true);

    % save metadata
    exp_name = exp_folder(end-15:end);
    save(fullfile(exp_folder, strcat('metadata_', exp_name, '.mat')), 'metadata');

    % Convert TDMS files to mat file.
    G4_TDMS_folder2struct(log_folder)

    ctlr.stopDisplay()
    ctlr.close()
    disp('finished');

end 