function run_protocol2(exp_folder, pattern_order, func_order, trial_dur)
% Run protocol 2 
    %% User-defined experiment conditions
    num_reps = 1; %number of repetitions for each stimuli
    fly_name = 'test_1';
    
    %% set up for experiment
    %Load configuration and start G4 Host
    % Check the use of this
    userSettings;
    
    load([exp_folder '\currentExp.mat']);
    num_conditions = numel(pattern_order);
    if ~exist(fullfile(exp_folder,'Log Files'),'dir')
        mkdir(exp_folder,'Log Files');
    end
    
    %% Open new Panels controller instance
    ctlr = PanelsController();
    ctlr.open(true);
    
    %% Change to root directory
    ctlr.setRootDirectory(exp_folder);
    
    %check if log files already present for this experiment
    assert(~exist([exp_folder '\Log Files\*'],'file'),'unsorted log files present in save folder, remove before restarting experiment\n');
    assert(~exist([exp_folder '\Results\' fly_name],'dir'),'Results folder already exists with that fly name\n');
    
    %finish setting up for experiment
    exp_seconds = currentExp.totalDuration;
    fprintf(['Estimated experiment duration: ' num2str(exp_seconds/60) ' minutes\n']);
    
    %% start experiment
    ctlr.startLog(); %starts logging data in .tdms files
    
    %block trial structure
    for r = 1:num_reps
        for c = 1:num_conditions
            %trial portion
            ctlr.setControlMode(1);
            ctlr.setPatternID(pattern_order(1,c));
            ctlr.setPatternFunctionID(func_order(1,c));
            trial_t = trial_dur(1, c);
            fprintf(['Rep ' num2str(r) ' of ' num2str(num_reps) ', cond ' num2str(c) ' of ' num2str(num_conditions) ': ' strjoin(currentExp.pattern.pattNames(pattern_order(1,c))) '\n']);
            ctlr.startDisplay(ceil(trial_t*10)-1); %duration expected in 100ms units
        end
    end
    
    %rename/move results folder
    ctlr.stopLog('showTimeoutMessage', true);
    movefile([exp_folder '\Log Files\*'],fullfile(exp_folder,'Results',fly_name));
    ctlr.close()
    ctrl.stopDisplay()
    
    disp('finished');
    % Need metadata and Log....
    % Didn't close the controller but ran everything. 
end 