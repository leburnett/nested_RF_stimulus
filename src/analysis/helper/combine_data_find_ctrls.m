function [data_align_all, ctrl_flies] = combine_data_find_ctrls(data_folder)

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