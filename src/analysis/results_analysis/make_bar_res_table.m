function make_bar_res_table(results_folder, results_save_folder)
% 
    % results_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/bar_results';
    cd(results_folder)
    
    res_files = dir('peak_vals_*');
    n_files = height(res_files);
    
    bar_res_table = table();
    
    for f = 1:n_files
        fname = res_files(f).name;
        load(fname, 'bar_results');
        bar_results.Type = string(bar_results.Type);
        bar_results.Strain = string(bar_results.Strain);
        bar_res_table = vertcat(bar_res_table, bar_results);
    end
    
    % results_save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results';
    save(fullfile(results_save_folder, 'bar_res_table_250213.mat'), 'bar_res_table');

end 
