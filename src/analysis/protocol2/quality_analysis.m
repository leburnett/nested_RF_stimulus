%% Quality analysis

close all
clear

% 1 - Load the data: 
% Input experiment folder with data from protocol 2:
date_folder = cd;
strrs = split(date_folder, '/');
% date_str = strrs{end};
date_str = date_folder(end-15:end-6);
time_str = date_folder(end-4:end);
type_str = string(strrs{end-1});
strain_str = string(strrs{end-2});

% Find the 'G4_TDMS_Log..mat' file and load it:
log_folder = fullfile(date_folder, "Log Files"); cd(log_folder);
log_file = dir('G4_TDMS*');
load(log_file.name, 'Log');

% Load the voltage data:
v_data = Log.ADC.Volts(2, :)*10;

%% Save quality metrics as a table:
qual_table = table();
qual_table.Date = date_str;
qual_table.Time = time_str;
qual_table.Strain = strain_str;
qual_table.Type = type_str;
qual_table.median_v = median(v_data);
qual_table.var_v = var(v_data);
qual_table.max_v = prctile(v_data, 98);
qual_table.min_v = prctile(v_data, 2);

qual_save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/qual_results';
save(fullfile(qual_save_folder, strcat('qual_res_', date_str, '_', time_str, '_', strain_str, '_', type_str, '.mat')), 'qual_table');
