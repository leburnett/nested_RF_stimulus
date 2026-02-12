function process_bar_flashes_p2(exp_folder, metadata, PROJECT_ROOT)
% PROCESS_BAR_FLASHES_P2  Analyze bar flash stimulus responses.
%
%   PROCESS_BAR_FLASHES_P2(EXP_FOLDER, METADATA, PROJECT_ROOT) extracts
%   voltage responses to bar flash stimuli, computes mean across reps,
%   generates visualization figures, and saves results.
%
%   INPUTS:
%     exp_folder   - Full path to the Protocol 2 experiment directory
%     metadata     - Structure containing: .Strain, .Frame, .Side, .Age
%     PROJECT_ROOT - Base directory for saving results and figures
%
%   ANALYSIS PIPELINE:
%     1. Loads voltage and frame data from TDMS log files
%     2. Loads bar flash pattern to determine frame count
%     3. Parses bar flash responses using PARSE_BAR_FLASH_DATA
%     4. Generates visualization figures for slow and fast conditions
%     5. Saves results to .mat file
%
%   OUTPUT FILES:
%     Saves to PROJECT_ROOT/results/bar_flash_data/:
%       bar_flash_results_<date>_<time>_<strain>_<on_off>.mat
%     Saves to PROJECT_ROOT/figures/bar_flash_stimuli/:
%       PDF figures for slow and fast conditions
%
%   See also PARSE_BAR_FLASH_DATA, PLOT_BAR_FLASH_DATA, GET_STIMULUS_METADATA

% Guard against load_protocol2_data changing the working directory
orig_dir = cd;
restoreDir = onCleanup(@() cd(orig_dir));

results_folder = fullfile(PROJECT_ROOT, "results", "bar_flash_data");
if ~isfolder(results_folder)
    mkdir(results_folder);
end

figures_folder = fullfile(PROJECT_ROOT, "figures", "bar_flash_stimuli");
if ~isfolder(figures_folder)
    mkdir(figures_folder);
end

[date_str, time_str, Log, params, ~] = load_protocol2_data(exp_folder);
on_off = params.on_off;
params.date = date_str;
params.time = time_str;
params.strain = metadata.Strain;

f_data = Log.ADC.Volts(1, :); % frame data
v_data = Log.ADC.Volts(2, :)*10; % voltage data
median_v = median(v_data);

% Load bar flash pattern info if available
bar_flash_pattern = [];
pattern_folder = fullfile(exp_folder, 'Patterns');
if isfolder(pattern_folder)
    mat_files = dir(fullfile(pattern_folder, '*FLASHES*.mat'));
    if ~isempty(mat_files)
        patt_data = load(fullfile(pattern_folder, mat_files(1).name), 'pattern');
        if isfield(patt_data, 'pattern')
            bar_flash_pattern = patt_data.pattern;
        end
    end
end

[data_slow, data_fast, mean_slow, mean_fast] = parse_bar_flash_data(f_data, v_data, bar_flash_pattern);

save_fig = 1;
fig_slow = plot_bar_flash_data(data_slow, mean_slow, median_v);
sgtitle(strcat("Barflashes-slow-", strrep(date_str, '_', '-'), '-', strrep(time_str, '_', '-'), '-', strrep(metadata.Strain, '_', '-')));

if save_fig
    fname_slow = fullfile(figures_folder, strcat("Bar_flashes_slow_", strrep(date_str, '_', '-'), '_', strrep(time_str, '_', '-'), '_', strrep(metadata.Strain, '_', '-'), '.pdf'));
    exportgraphics(fig_slow, fname_slow, 'ContentType', 'vector', 'BackgroundColor', 'none');
    close
end

fig_fast = plot_bar_flash_data(data_fast, mean_fast, median_v);
sgtitle(strcat("Barflashes-fast-", strrep(date_str, '_', '-'), '-', strrep(time_str, '_', '-'), '-', strrep(metadata.Strain, '_', '-')));

if save_fig
    fname_fast = fullfile(figures_folder, strcat("Bar_flashes_fast_", strrep(date_str, '_', '-'), '_', strrep(time_str, '_', '-'), '_', strrep(metadata.Strain, '_', '-'), '.pdf'));
    exportgraphics(fig_fast, fname_fast, 'ContentType', 'vector', 'BackgroundColor', 'none');
    close
end

save(fullfile(results_folder, strcat('bar_flash_results_', date_str, '_', time_str, '_', metadata.Strain, '_', on_off, '.mat')), ...
    'data_slow', 'data_fast', 'mean_slow', 'mean_fast');

end
