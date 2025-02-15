function make_rf_res_table(results_folder, results_save_folder)
% 
    % results_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/rf_results';
    cd(results_folder)
    
    res_files = dir('rf_results_*');
    n_files = height(res_files);
    
    rf_res_table = table();
    
    for f = 1:n_files
        fname = res_files(f).name;
        load(fname, 'rf_results');
        rf_results.Type = string(rf_results.Type);
        rf_results.Strain = string(rf_results.Strain);
        rf_res_table = vertcat(rf_res_table, rf_results);
    end
    
    % results_save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results';
    save(fullfile(results_save_folder, 'rf_res_table_250213.mat'), 'rf_res_table');

end 
