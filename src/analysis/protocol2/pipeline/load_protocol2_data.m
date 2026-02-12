function [date_str, time_str, Log, params, pfnparam] = load_protocol2_data(exp_folder)
% LOAD_PROTOCOL2_DATA  Load recorded data and metadata from Protocol 2 experiment.
%
%   [DATE_STR, TIME_STR, LOG, PARAMS, PFNPARAM] = LOAD_PROTOCOL2_DATA(EXP_FOLDER)
%   loads the TDMS log file, parameter files, and position function data
%   from a completed Protocol 2 experiment folder.
%
%   INPUT:
%     exp_folder - Full path to the experiment directory
%                  (e.g., 'C:\...\nested_RF_protocol2\2025_05_09_15_46')
%
%   OUTPUTS:
%     date_str - Experiment date string (format: 'YYYY_MM_DD')
%     time_str - Experiment time string (format: 'HH_MM')
%     Log      - Structure from G4_TDMS_Log*.mat containing:
%                .Voltage - Recorded voltage trace (10kHz sampling)
%                .Frame   - Frame number at each time point
%     params   - Structure with stimulus parameters (.x, .y, .on_off, etc.)
%     pfnparam - Position function parameters (.func, .dur, .frames, etc.)
%
%   FILE LOCATIONS:
%     - Log Files/G4_TDMS_Log*.mat - Main recording data
%     - params/*.mat - Stimulus generation parameters
%     - Functions/0001*.mat - Position function parameters
%
%   USAGE:
%     Called by PROCESS_BARS_P2 and PROCESS_FLASH_P2 at the start of
%     data analysis pipeline.
%
%   See also PROCESS_PROTOCOL2, PARSE_BAR_DATA, PARSE_FLASH_DATA

    % Input experiment folder with data from protocol 2:
    cd(exp_folder)

    if contains(exp_folder, '/')
        strrs = split(exp_folder, '/');
    elseif contains(exp_folder, '\')
        strrs = split(exp_folder, '\');
    end 
    exp_meta = strrs{end};
    
    date_str = exp_meta(end-15:end-6);
    time_str = exp_meta(end-4:end);
    
    % Find the 'G4_TDMS_Log..mat' file and load it:
    log_folder = fullfile(exp_folder, "Log Files"); 
    cd(log_folder);
    log_file = dir('G4_TDMS*');
    if isempty(log_file)
        warning("G4_TDMS log file does not exist.")
    else
        load(log_file.name, 'Log');
    end 

    % Find out whether the flashes presented were ON or OFF pixels:
    params_folder = fullfile(exp_folder, "params"); 
    cd(params_folder);
    param_file = dir('*.mat');
    if isempty(param_file)
        warning("params file does not exist.")
    else
       % % % % % % % % % % % % THIS CURRENTLY OPENS THE FIRST FILE - 4px
       % flashes - currently this is only used for file saving so shouldn't
       % matter with new 6px stimulus. 
        load(param_file(1).name, 'params');
    end 

    func_folder = fullfile(exp_folder, "Functions");
    cd(func_folder)
    func_file = dir("0001*");
    if isempty(func_file)
        warning("function file does not exist.")
    else
        load(func_file.name, 'pfnparam');
    end 

end 