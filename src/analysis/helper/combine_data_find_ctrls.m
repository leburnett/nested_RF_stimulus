function [data_align_all, ctrl_flies] = combine_data_find_ctrls(data_folder)
% COMBINE_DATA_FIND_CTRLS  Aggregate bar results and identify control flies.
%
%   [DATA_ALIGN_ALL, CTRL_FLIES] = COMBINE_DATA_FIND_CTRLS(DATA_FOLDER)
%   loads all bar results files from a folder, combines the aligned
%   timeseries data, and identifies which recordings are from control flies.
%
%   INPUT:
%     data_folder - Path to folder containing peak_vals*.mat files
%
%   OUTPUTS:
%     data_align_all - 32 x n_exp cell array of mean timeseries
%                      Rows 1-16: slow bar directions, Rows 17-32: fast bars
%                      Each column is one experiment/cell
%     ctrl_flies     - n_exp x 1 logical array
%                      true if filename contains 'CTL', false otherwise
%
%   FILE NAMING:
%     Control flies are identified by 'CTL' in the filename.
%     Results files expected format: peak_vals_<strain>_<type>_<date>_<time>.mat
%
%   USAGE:
%     Used for population analyses comparing experimental vs control groups.
%     Combined data can be used for group averaging and statistical tests.
%
%   See also ANALYSE_BAR_RESPONSE_COMBINE_DATA, PROCESS_BARS_P2

    cd(data_folder)
    file_list = dir('*peak_vals*');
    n_exp = height(file_list);
    
    % Initialise empty arrays for filling in:
    data_align_all = cell(32, n_exp);
    ctrl_flies = zeros(n_exp, 1);
    
    %% Combine the data from all of the cells of this type:
    
    for i = 1:n_exp
    
        % Load needed variables:
        load(file_list(i).name);
    
        % Determine if the data is from a control or test fly.
        if contains(file_list(i).name, 'CTL')
            ctrl_flies(i) = 1;
        end 
    
        % Combine the "mean" aligned timeseries from all of the flies.  
        data_align_all(:, i) = data_aligned(:, 4); % Combine the mean timeseries only.
    end 
    
    % Convert 'ctrl_flies' from double to logical for logical indexing. 
    ctrl_flies = logical(ctrl_flies);

end 