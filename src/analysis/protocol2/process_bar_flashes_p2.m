function process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT)
% PROCESS_BAR_FLASHES_P2  Analyse bar flash responses from Protocol 2.
%
%   PROCESS_BAR_FLASHES_P2(EXP_FOLDER, METADATA, PROJECT_ROOT)
%   loads recorded data from the experiment folder, parses the bar flash
%   stimuli at two flash durations (80ms slow, 14ms fast), generates tiled
%   figures showing individual and mean responses across 8 orientations x
%   11 positions, and saves results.
%
%   INPUTS:
%     exp_folder   - char | string
%                    Full path to the experiment directory (containing
%                    Log Files/, Patterns/, Functions/, currentExp.mat).
%     metadata     - struct
%                    Experiment metadata with fields:
%                      .Strain (char) - fly genotype string
%                      (plus .Frame, .Side, .Age from GET_INPUT_PARAMETERS)
%     PROJECT_ROOT - char | string
%                    Root directory for saving results and figures.
%                    Results go to PROJECT_ROOT/results/bar_flash_results/
%                    Figures go to PROJECT_ROOT/figures/bar_flash_stimuli/
%
%   OUTPUTS (saved to disk):
%     Figures:
%       Bar_flashes_80ms_<date>_<time>_<strain>.pdf  - slow flash tiled plot
%       Bar_flashes_14ms_<date>_<time>_<strain>.pdf  - fast flash tiled plot
%     Results:
%       bar_flash_results_<date>_<time>_<strain>_<on_off>.mat containing:
%         data_slow  - 11x8x3 cell (positions x orientations x reps)
%         data_fast  - 11x8x3 cell
%         mean_slow  - 11x8 cell (mean across reps)
%         mean_fast  - 11x8 cell
%
%   PARAMETERS:
%     prop_int_2_show = 0.75 : fraction of the inter-flash interval kept
%       in the extracted timeseries. This controls how much baseline
%       before/after each flash is visible in the plots.
%
%   FLASH TIMING:
%     Slow : 80ms flash + 920ms gap = 1s per cycle (10000 samples)
%     Fast : 14ms flash + 486ms gap = 0.5s per cycle (5000 samples)
%     Acquisition rate = 10 kHz (1 sample = 100 us)
%
%   See also PARSE_BAR_FLASH_DATA, PLOT_BAR_FLASH_DATA, LOAD_PROTOCOL2_DATA


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