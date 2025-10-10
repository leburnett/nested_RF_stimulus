function process_flash_p2(exp_folder, metadata, PROJECT_ROOT, resultant_angle)

results_folder = fullfile(PROJECT_ROOT, "results", "flash_results");
if ~isfolder(results_folder)
    mkdir(results_folder);
end

figures_folder = fullfile(PROJECT_ROOT, "figures", "flash_stimuli");
if ~isfolder(figures_folder)
    mkdir(figures_folder);
end

% strain_str = "T4T5"; % eventually will get this from metadata. 

[date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);
on_off = params.on_off;
params.date = date_str;
params.time = time_str;
params.strain = metadata.Strain;
params.resultant_angle = resultant_angle;

% sampling_rate = 10000;

f_data = Log.ADC.Volts(1, :); % frame data
diff_f_data = diff(f_data);

% Find 'idx' the timepoints of the first 4px flash and the first 6 px flash
% of each rep.
n_flashes_4px = 196;
n_flashes_6px = 100;
if on_off == "off"
    % conveniently for the dark flashes, the first flash of both the 4px
    % and 6px flashes is frame 1. 
    idx = find(diff_f_data == 1 & f_data(2:end) == 1); % First flash.
elseif on_off == "on"
    % 4 px flash light flashes:
    idx_4 = find(diff_f_data == 1 + n_flashes_4px & f_data(2:end) == 1 + n_flashes_4px); % First flash.
    idx_4(:, [2,4,6]) = []; % Remove timepoints within the 6px flashes. 
    % 6 px light flashes
    idx_6 = find(diff_f_data == 1 + n_flashes_6px & f_data(2:end) == 1 + n_flashes_6px); % Should be 3 values.
    % combine the indices for the start of the 4px and 6 px flashes - sort
    % to order them.
    idx = sort(horzcat(idx_4, idx_6));
end 

% TEST 
% Plot the values of 'idx'
% figure; plot(f_data);
% hold on;
% for kk = 1:numel(idx)
%     plot([idx(kk), idx(kk)], [0 400], 'r');
% end 

v_data = Log.ADC.Volts(2, :)*10; % voltage data
median_v = median(v_data);
v2_data = v_data - median_v; % Get the median-subtracted voltage.

% % - - Figure - check quality of the recording:
% plot_quality_check_data_full_rec(f_data, v_data, save_fig, PROJECT_ROOT)

% Voltage variance - store in metadata 
filtered_voltage_data = movmean(v_data, 20000);
var_filtered_v = var(filtered_voltage_data);

% Flash duration(s)
% flash_dur_s = pfnparam.dur;
% flash_dur_ms = 976700; % flash_dur_s*sampling_rate
% slow_flashes_dur = 5000; % 340ms + 160ms FLASH.
% fast_flashes_dur = 2500; % 170ms + 80ms FLASH.

% % - - Figure - check flash timing
% plot_quality_check_flash_timing(f_data, flash_dur_ms, save_fig, PROJECT_ROOT)
% rf_results = struct();
% rf_results.Date = date_str;
% rf_results.Time = time_str;
% rf_results.Strain = metadata.Strain;
% rf_results.Type = on_off;

for px_size = [4, 6]

    slow_fast = "slow";
    speed_str = "160ms";

    if px_size == 4
        n_flashes = n_flashes_4px;
    elseif px_size == 6
        n_flashes = n_flashes_6px;
    end 
    
    [data_comb,cmap_id,var_across_reps,var_within_reps,diff_mean,max_data,min_data] = parse_flash_data(f_data, v_data, on_off, slow_fast, px_size, PROJECT_ROOT);
    
    % med_var_X_reps = median(reshape(var_across_reps, [1, n_flashes]));
    % med_var_W_reps = var(reshape(var_within_reps, [1, n_flashes]));
    % med_diff_mean = median(reshape(diff_mean, [1, n_flashes]));
    
    % Rescale the combined data to be between 0 and 1.
    data_comb2 = rescale(data_comb, 0, 1);
    
    % Timeseries plot:
    f_timeseries = plot_rf_estimate_timeseries_line(data_comb2, cmap_id, f_data, v2_data, slow_fast, px_size, idx, params);
    fname = fullfile(figures_folder, strcat('Timeseries_', date_str, '_', time_str, '_', metadata.Strain, '_', on_off, "_", speed_str, ".pdf"));
    exportgraphics(f_timeseries ...
            , fname ...
            , 'ContentType', 'vector' ...
            , 'BackgroundColor', 'none' ...
            ); 
    
    % Simple heat map plot:
    f_heatmap = plot_heatmap_flash_responses(data_comb2);
    fname = fullfile(figures_folder, strcat('Heatmap_', date_str, '_', time_str, '_', metadata.Strain, '_', on_off, "_", speed_str, ".pdf"));
    exportgraphics(f_heatmap ...
            , fname ...
            , 'ContentType', 'vector' ...
            , 'BackgroundColor', 'none' ...
            ); 
    
    % Generate plot with Gausssian RF estimates:
    % exc_data = data_comb2;
    % inh_data = data_comb;
    % inh_data(cmap_id~=2)=0;
    % [optEx, R_squared, optInh, R_squaredi, ~, ~] = gaussian_RF_estimate(exc_data, inh_data);
    
    % Combine all of the results into a table:
    % rf_results.(slow_fast).data_comb = {data_comb};
    % rf_results.(slow_fast).max_data = {max_data};
    % rf_results.(slow_fast).min_data = {min_data};
    % rf_results.(slow_fast).diff_mean = {diff_mean};
    % rf_results.(slow_fast).cmap_id = {cmap_id};
    % rf_results.(slow_fast).max_val = prctile(reshape(max_data, [1, n_flashes]), n_flashes/2);
    % rf_results.(slow_fast).min_val = prctile(reshape(min_data, [1, n_flashes]), 2);
    % rf_results.(slow_fast).var_within_reps = {var_within_reps};
    % rf_results.(slow_fast).var_across_reps = {var_across_reps};
    % rf_results.(slow_fast).var_filtered_v = var_filtered_v;
    % rf_results.(slow_fast).med_var_X_reps = med_var_X_reps;
    % rf_results.(slow_fast).med_var_W_reps = med_var_W_reps;
    % rf_results.(slow_fast).med_diff_mean = med_diff_mean;
    % rf_results.(slow_fast).R_squared = R_squared;
    % rf_results.(slow_fast).sigma_x_exc = optEx(4);
    % rf_results.(slow_fast).sigma_y_exc = optEx(5);
    % rf_results.(slow_fast).optExc = {optEx};
    % rf_results.(slow_fast).R_squaredi = R_squaredi;
    % rf_results.(slow_fast).optInh = {optInh};
    % rf_results.(slow_fast).sigma_x_inh = optInh(4);
    % rf_results.(slow_fast).sigma_y_inh = optInh(5);

end 

% save(fullfile(results_folder, strcat('rf_results_', date_str,'_', time_str, '_', metadata.Strain, '_', on_off, '.mat')), 'rf_results');

end 




