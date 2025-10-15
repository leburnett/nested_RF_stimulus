function run_protocol2(exp_folder, pattern_order, func_order, trial_dur, n_reps)
% Run protocol 2 
    %% Experiment metadata from user input:
    % n_reps = 3;

    bar_flash_pattern_slow = max(pattern_order)-1; % The bar flash pattern is the last pattern.
    bar_flash_pattern_fast = max(pattern_order);
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

            % Use different position functions for bar flashes each rep.
            if pattern_order(1,c) == bar_flash_pattern_slow
                if r == 2  
                    ctlr.setPatternID(bar_flash_pattern_slow+1);
                elseif r == 3
                    ctlr.setPatternID(bar_flash_pattern_slow+2);
                else
                    ctlr.setPatternID(pattern_order(1,c));
                end 
            else
                ctlr.setPatternID(pattern_order(1,c));
            end 

            % Use different position functions for bar flashes each rep.
            if pattern_order(1,c) == bar_flash_pattern_fast
                if r == 2  
                    ctlr.setPatternID(bar_flash_pattern_fast+1);
                elseif r == 3
                    ctlr.setPatternID(bar_flash_pattern_fast+2);
                else
                    ctlr.setPatternID(pattern_order(1,c));
                end 
            else
                ctlr.setPatternID(pattern_order(1,c));
            end

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