function process_flash_p2(exp_folder)

% exp_folder = "/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/control/ON/2024_12_18_15_07";
strain_str = "TmY3"; % eventually will get this from metadata. 

[date_str, time_str, Log, params, pfnparam] = load_protocol2_data(exp_folder);
on_off = params.on_off;

sampling_rate = 10000;

f_data = Log.ADC.Volts(1, :); % frame data
diff_f_data = diff(f_data);
idx = find(diff_f_data == min(diff_f_data)); % where the flash stimuli end.

v_data = Log.ADC.Volts(2, :)*10; % voltage data
median_v = median(v_data);
v2_data = v_data - median_v; % Get the median-subtracted voltage.

% STD and the upper/lower bounds for finding response groups later on. 
std_v = std(v2_data);
upper_bound_med = std_v;
lower_bound_med = -std_v/2;

% % - - Figure - check quality of the recording:
% plot_quality_check_data_full_rec(f_data, v_data)

% % - - Figure - check bound positions
% plot_quality_check_data_bounds(v_data, median_v, upper_bound_med, lower_bound_med)

% Voltage variance - store in metadata 
filtered_voltage_data = movmean(v_data, 20000);
var_filtered_v = var(filtered_voltage_data);

% Flash duration(s)
flash_dur_s = pfnparam.dur;
flash_dur_ms = 976700; % flash_dur_s*sampling_rate
slow_flashes_dur = 5000; % 340ms + 160ms bkg.

% % - - Figure - check flash timing
% plot_quality_check_flash_timing(f_data, flash_dur_ms)

% 'fc' is the function - each number corresponds to the frame that is being
% presented. There is one data point every 2ms. 
fc = pfnparam2.func; 

% Determine if the flashes presented were ON (bright) or OFF (dark)
max_frame = max(fc);
if max_frame < 198
    on_off = "off";
else 
    on_off = "on";
end 

[data_comb,...
    cmap_id,...
    var_across_reps,...
    var_within_reps,...
    diff_mean,...
    max_data,...
    min_data] = parse_flash_data(f_data, flash_dur_ms);

% Rescale the combined data to be between 0 and 1.
data_comb2 = rescale(data_comb, 0, 1);

f_timeseries = plot_rf_estimate_timeseries_line(data_comb2, cmap_id, f_data, v2_data, slow_flashes_dur, idx, flash_dur_ms, on_off);
f_heatmap = plot_heatmap_flash_responses(data_comb2);

% Generate plot with Gausssian RF estimates:
exc_data = data_comb2;
inh_data = data_comb;
inh_data(cmap_id~=2)=0;
[optEx, R_squared, optInh, R_squaredi, ~, ~] = gaussian_RF_estimate(exc_data, inh_data);

% Combine all of the results into a table:
rf_results = table();
rf_results.Date = date_str;
rf_results.Time = time_str;
rf_results.Strain = strain_str;
rf_results.Type = on_off;
rf_results.data_comb = {data_comb};
rf_results.max_data = {max_data};
rf_results.min_data = {min_data};
rf_results.diff_mean = {diff_mean};
rf_results.cmap_id = {cmap_id};
rf_results.max_val = prctile(reshape(max_data, [1, 196]), 98);
rf_results.min_val = prctile(reshape(min_data, [1, 196]), 2);
rf_results.var_within_reps = {var_within_reps};
rf_results.var_across_reps = {var_across_reps};
rf_results.var_filtered_v = var_filtered_v;
rf_results.med_var_X_reps = med_var_X_reps;
rf_results.med_var_W_reps = med_var_W_reps;
rf_results.med_diff_mean = med_diff_mean;
rf_results.R_squared = R_squared;
rf_results.sigma_x_exc = optEx(4);
rf_results.sigma_y_exc = optEx(5);
rf_results.optExc = {optEx};
rf_results.R_squaredi = R_squaredi;
rf_results.optInh = {optInh};
rf_results.sigma_x_inh = optInh(4);
rf_results.sigma_y_inh = optInh(5);

% save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/rf_results';
% save(fullfile(save_folder, strcat('rf_results_', date_str,'_', strain_str, '_', type_str, '.mat')), 'rf_results');

end 




