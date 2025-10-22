function process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT)

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
median_v = median(v_data);

[data_slow, data_fast, mean_slow, mean_fast]  = parse_bar_flash_data(f_data, v_data);

save_fig = 1;
fig_slow = plot_bar_flash_data(data_slow, mean_slow, median_v);
fig_fast = plot_bar_flash_data(data_fast, mean_fast, median_v);

if save_fig

    fname_slow = fullfile(figures_folder, strcat("Bar_flashes_80ms_",  strrep(date_str, '_', '-'),'_', strrep(time_str, '_', '-'), '_', strrep(metadata.Strain, '_', '-'), '.pdf')); 
    exportgraphics(fig_slow, fname_slow, 'ContentType', 'vector', 'BackgroundColor', 'none');

    fname_fast = fullfile(figures_folder, strcat("Bar_flashes_14ms_",  strrep(date_str, '_', '-'),'_', strrep(time_str, '_', '-'), '_', strrep(metadata.Strain, '_', '-'), '.pdf')); 
    exportgraphics(fig_fast, fname_fast, 'ContentType', 'vector', 'BackgroundColor', 'none');

end 

save(fullfile(results_folder, strcat('bar_flash_results_', date_str,'_', time_str, '_', metadata.Strain, '_', on_off, '.mat')), 'data_slow', 'data_fast', 'mean_slow', 'mean_fast');

end 