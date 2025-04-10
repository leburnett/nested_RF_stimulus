function [median_vs] = combine_median_v(data_folder)
    % Find the median voltage per recording and store them altogether in a
    % single array.

    cd(data_folder)
    file_list = dir('*peak_vals*');
    n_exp = height(file_list);
   
    median_vs = zeros(1, n_exp);
    
    for i = 1:n_exp
        % Load the data from one cell:
        load(file_list(i).name);
        median_vs(i) = median_voltage;
    end 

end 