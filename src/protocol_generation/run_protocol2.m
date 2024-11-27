function run_protocol2(exp_folder)
% Run protocol 2 
    
    %% User-defined experiment conditions
    num_reps = 1; %number of repetitions for each stimuli
    exp_mode = 1; %0=streaming, 1=position function, 2=constant rate, 3=position change, 4=Closed-loop (CL), 5=CL+bias, 6=CL+OL
    fly_name = 'testfly1';
    
    trial_duration = 4; %duration (in seconds) of each trial
    
    %% set up for experiment
    %Load configuration and start G4 Host
    % Check the use of this
    userSettings;
    
    load([exp_folder '\currentExp.mat']);
    num_conditions = currentExp.pattern.num_patterns;
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
    
    %create .mat file of experiment order
    exp_order = repmat(1:num_conditions,num_reps,1);
    
    %finish setting up for experiment
    exp_seconds = num_reps*num_conditions*trial_duration;
    fprintf(['Estimated experiment duration: ' num2str(exp_seconds/60) ' minutes\n']);
    save([exp_folder '\Log Files\exp_order.mat'],'exp_order')
    
    %% start experiment
    ctlr.startLog(); %starts logging data in .tdms files
    
    %block trial structure
    for r = 1:num_reps
        for c = 1:num_conditions
            %trial portion
            ctlr.setControlMode(exp_mode);
            ctlr.setPatternID(exp_order(r,c));
            ctlr.setPatternFunctionID(exp_order(r,c));
            fprintf(['Rep ' num2str(r) ' of ' num2str(num_reps) ', cond ' num2str(c) ' of ' num2str(num_conditions) ': ' strjoin(currentExp.pattern.pattNames(exp_order(r,c))) '\n']);
            ctlr.startDisplay((trial_duration*10)-1); %duration expected in 100ms units
        end
    end
    
    %rename/move results folder
    ctlr.stopLog('showTimeoutMessage', true);
    movefile([exp_folder '\Log Files\*'],fullfile(exp_folder,'Results',fly_name));
    ctlr.close()
    
    disp('finished');
end 