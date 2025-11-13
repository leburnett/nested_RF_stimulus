function process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT)
% Bar flashes occur at 2 speeds currently. 80ms and 140ms. 
% Acquisition rate = 10,000Hz - therefore one sample every 100us. 
% 80 ms = 80,000us = 800 samples. 

% 80ms + 920ms gap = 1s.
% 14ms + 486ms  = 500ms = 0.5s . 


results_folder = fullfile(PROJECT_ROOT, "results", "bar_flash_results");
if ~isfolder(results_folder)
    mkdir(results_folder);
end

figures_folder = fullfile(PROJECT_ROOT, "figures", "bar_flash_stimuli");
if ~isfolder(figures_folder)
    mkdir(figures_folder);
end

% strain_str = "T4T5"; % eventually will get this from metadata. 

[date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);
on_off = params.on_off;
params.date = date_str;
params.time = time_str;
params.strain = metadata.Strain;

f_data = Log.ADC.Volts(1, :); % frame data
    
v_data = Log.ADC.Volts(2, :)*10; % voltage data
% median_v = median(v_data);

prop_int_2_show = 0.75; % Proportion of interval time to keep in timeseries.
[data_slow, data_fast, mean_slow, mean_fast]  = parse_bar_flash_data(f_data, v_data, prop_int_2_show);

save_fig = 1;
interval_between_flashes_slow = 10000; % number of sampling points (10,000 samples * 100us (time per sample) = 1s)
n_samples_int_slow = interval_between_flashes_slow*prop_int_2_show; % How many samples of the interval are in the timeseries.
fig_slow = plot_bar_flash_data(data_slow, mean_slow, n_samples_int_slow);
sgtitle(strcat("Barflashes-80ms-", string(n_samples_int_slow/10),"msInt-", strrep(date_str, '_', '-'),'-', strrep(time_str, '_', '-'), '-', strrep(metadata.Strain, '_', '-')), "FontSize", 25);

if save_fig
    fname_slow = fullfile(figures_folder, strcat("Bar_flashes_80ms_",  strrep(date_str, '_', '-'),'_', strrep(time_str, '_', '-'), '_', strrep(metadata.Strain, '_', '-'), '.pdf')); 
    exportgraphics(fig_slow, fname_slow, 'ContentType', 'vector', 'BackgroundColor', 'none');
    % close
end 

interval_between_flashes_fast = 5000; % time in ms.
n_samples_int_fast = interval_between_flashes_fast*prop_int_2_show; 
fig_fast = plot_bar_flash_data(data_fast, mean_fast, n_samples_int_fast); % time in ms.;
sgtitle(strcat("Barflashes-14ms-", string(n_samples_int_fast/10),"msInt-", strrep(date_str, '_', '-'),'-', strrep(time_str, '_', '-'), '-', strrep(metadata.Strain, '_', '-')), "FontSize", 25);

if save_fig
    fname_fast = fullfile(figures_folder, strcat("Bar_flashes_14ms_",  strrep(date_str, '_', '-'),'_', strrep(time_str, '_', '-'), '_', strrep(metadata.Strain, '_', '-'), '.pdf')); 
    exportgraphics(fig_fast, fname_fast, 'ContentType', 'vector', 'BackgroundColor', 'none');
    % close
end 

save(fullfile(results_folder, strcat('bar_flash_results_', date_str,'_', time_str, '_', metadata.Strain, '_', on_off, '.mat')), 'data_slow', 'data_fast', 'mean_slow', 'mean_fast');

end 