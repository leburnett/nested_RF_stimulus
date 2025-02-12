function make_qual_res_table(results_folder, results_save_folder)

    % results_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/qual_results';
    cd(results_folder)
    
    res_files = dir('qual_res*');
    n_files = height(res_files);
    
    qual_res_table = table();
    
    for f = 1:n_files
        fname = res_files(f).name;
        load(fname, 'qual_table');
        qual_table.Type = string(qual_table.Type);
        qual_table.Strain = string(qual_table.Strain);
        qual_res_table = vertcat(qual_res_table, qual_table);
    end
    
    % results_save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results';
    save(fullfile(results_save_folder, 'qual_res_table_250212.mat'), 'qual_res_table');

end 
