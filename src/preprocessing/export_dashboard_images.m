function export_dashboard_images()
% EXPORT_DASHBOARD_IMAGES  Generate all PNG images for the ephys dashboard.
%
%   EXPORT_DASHBOARD_IMAGES() reads the experiment log, filters to T4T5
%   Summer/Autumn 2025 cells (excluding test genotype), and generates:
%     - GridPlot PNGs from Protocol 1 .fig files
%     - Bar sweep polar PNGs from saved bar sweep data
%     - Flash RF timeseries and heatmap PNGs from raw data + rf_results
%     - Bar flash PNGs from saved bar flash data
%     - Thumbnails for gallery view
%     - cell_index.json metadata file
%
%   All images are exported with white backgrounds at 300 DPI.
%
%   PREREQUISITES:
%     - MATLAB with access to the analysis code on the path
%     - All Protocol 2 data and results directories accessible
%
%   Run this once before launching the Python dashboard:
%     >> addpath(genpath('/Users/burnettl/Documents/GitHub/nested_RF_stimulus/src'))
%     >> export_dashboard_images

%% Configuration
EXCEL_PATH = '/Users/burnettl/Documents/GitHub/nested_RF_stimulus/reiser-lab-ephys-experiment-log.xlsx';
P1_BASE = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol1/bkg4/protocol1_10kHz_4reps_12px_6px_LHS_2sbkg_200msfl_150msint_07-10-25_15-45-60/';
P2_SUMMER = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/T4T5_Summer2025/';
P2_AUTUMN = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/T4T5_autumn2025/';
RESULTS_1DRF = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/1DRF/';
FLASH_RESULTS = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/flash_results/';
BAR_SWEEP_DATA = fullfile(RESULTS_1DRF, 'bar_sweep_data');
BAR_FLASH_DATA = fullfile(RESULTS_1DRF, 'bar_flash_data');
% Additional search paths — process_bars_p2 and process_bar_flashes_p2 save
% to these directories (used by batch_generate_flash_results), whereas the
% 1DRF paths above were populated by batch_analyze_1DRF (Autumn only).
BAR_RESULTS_ALT = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/bar_results/';
BAR_FLASH_ALT = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/bar_flash_results/';
OUTPUT_DIR = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/dashboard_images/';

EXPORT_OPTS = {'Resolution', 300, 'BackgroundColor', 'white'};
THUMB_WIDTH = 400; % pixels

%% Create output directories
subdirs = {'gridplots', 'bar_polar', 'flash_rf', 'flash_rf_heatmap', 'bar_flash', ...
           'thumbs/gridplots', 'thumbs/bar_polar', 'thumbs/flash_rf', 'thumbs/flash_rf_heatmap', 'thumbs/bar_flash'};
for i = 1:numel(subdirs)
    d = fullfile(OUTPUT_DIR, subdirs{i});
    if ~isfolder(d)
        mkdir(d);
    end
end

%% Read and filter experiment log
T = readtable(EXCEL_PATH, 'Sheet', 'Form responses 1', 'VariableNamingRule', 'preserve');

% Rename columns with special chars for easier access
T.Properties.VariableNames = strrep(T.Properties.VariableNames, ' ', '_');
T.Properties.VariableNames = strrep(T.Properties.VariableNames, '?', '_');
T.Properties.VariableNames = strrep(T.Properties.VariableNames, '/', '_');

% Filter to Summer2025 and Autumn2025 project types
mask = contains(T.project_type, 'Summer2025') | contains(T.project_type, 'Autumn2025');

% Exclude test genotype
mask = mask & ~contains(T.strain, 'test', 'IgnoreCase', true);

cells = T(mask, :);
n_cells = height(cells);
fprintf('Processing %d cells (Summer + Autumn 2025, excluding test)\n', n_cells);

%% Process each cell
for i = 1:n_cells
    date_str = cells.date{i};
    time_str = cells.time{i};
    cell_id = sprintf('%s_%s', date_str, time_str);
    strain = cells.strain{i};
    pc = cells.preferred_contrast{i};

    fprintf('\n[%d/%d] %s (%s, %s)\n', i, n_cells, cell_id, strain, pc);

    % Determine ON/OFF for file naming
    if strcmpi(pc, 'ON')
        on_off = 'on';
    else
        on_off = 'off';
    end

    % Find the P2 experiment folder
    exp_folder = '';
    summer_path = fullfile(P2_SUMMER, cell_id);
    autumn_path = fullfile(P2_AUTUMN, cell_id);
    if isfolder(summer_path)
        exp_folder = summer_path;
    elseif isfolder(autumn_path)
        exp_folder = autumn_path;
    end

    %% 1. GridPlots from Protocol 1 (4 separate figures per .fig file)
    gridplot_1_path = fullfile(OUTPUT_DIR, 'gridplots', [cell_id '_gridplot_1.png']);
    if ~isfile(gridplot_1_path)
        fig_path = find_matching_p1_gridplot(P1_BASE, date_str, time_str);
        if ~isempty(fig_path) && isfile(fig_path)
            try
                fprintf('  Exporting GridPlots (4 figures)...\n');
                figs = openfig(fig_path, 'invisible');
                n_figs = numel(figs);
                fprintf('    openfig returned %d figure handles\n', n_figs);
                for fi = 1:n_figs
                    figs(fi).Color = 'w';
                    gp_out = fullfile(OUTPUT_DIR, 'gridplots', ...
                        sprintf('%s_gridplot_%d.png', cell_id, fi));
                    exportgraphics(figs(fi), gp_out, EXPORT_OPTS{:});
                    make_thumbnail(gp_out, ...
                        fullfile(OUTPUT_DIR, 'thumbs', 'gridplots', ...
                        sprintf('%s_gridplot_%d.png', cell_id, fi)), THUMB_WIDTH);
                    fprintf('    Exported gridplot %d/%d\n', fi, n_figs);
                end
                close(figs);
            catch ME
                warning('  GridPlot failed for %s: %s', cell_id, ME.message);
            end
            close all force
        else
            fprintf('  No matching P1 GridPlot found\n');
        end
    else
        fprintf('  GridPlots already exist, skipping\n');
    end

    %% 2. Bar sweep polar plot — regenerate from saved bar sweep data
    bar_polar_path = fullfile(OUTPUT_DIR, 'bar_polar', [cell_id '_bar_polar.png']);
    if ~isfile(bar_polar_path)
        % Find saved bar sweep results — search multiple directories
        sweep_path = find_peak_vals(BAR_SWEEP_DATA, BAR_RESULTS_ALT, strain, on_off, date_str, time_str);

        if ~isempty(sweep_path)
            try
                fprintf('  Regenerating bar polar from saved data...\n');
                S = load(sweep_path, 'data', 'bar_results');
                params_polar = struct('date', date_str, 'time', time_str, ...
                    'strain', strain, 'on_off', on_off);
                plot_timeseries_polar_bars(S.data, S.bar_results.median_voltage, ...
                    params_polar, false, '');
                fig = gcf;
                fig.Color = 'w';
                exportgraphics(fig, bar_polar_path, EXPORT_OPTS{:});
                make_thumbnail(bar_polar_path, ...
                    fullfile(OUTPUT_DIR, 'thumbs', 'bar_polar', [cell_id '_bar_polar.png']), THUMB_WIDTH);
            catch ME
                warning('  Bar polar failed for %s: %s', cell_id, ME.message);
            end
            close all force
        else
            fprintf('  No bar sweep data found for %s\n', cell_id);
        end
    else
        fprintf('  Bar polar already exists, skipping\n');
    end

    %% 3. Flash RF timeseries + heatmap
    flash_rf_path = fullfile(OUTPUT_DIR, 'flash_rf', [cell_id '_flash_rf.png']);
    flash_rf_heatmap_path = fullfile(OUTPUT_DIR, 'flash_rf_heatmap', [cell_id '_flash_rf_heatmap.png']);
    if (~isfile(flash_rf_path) || ~isfile(flash_rf_heatmap_path)) && ~isempty(exp_folder)
        % Find matching rf_results .mat file
        rf_pattern = fullfile(FLASH_RESULTS, sprintf('rf_results_%s_*_%s.mat', cell_id, on_off));
        rf_files = dir(rf_pattern);

        % Try broader pattern if exact match fails
        if isempty(rf_files)
            rf_pattern2 = fullfile(FLASH_RESULTS, sprintf('rf_results_%s_*%s*.mat', date_str, on_off));
            rf_files = dir(rf_pattern2);
        end

        if ~isempty(rf_files)
            rf_path = fullfile(FLASH_RESULTS, rf_files(1).name);
            fprintf('  Generating flash RF images (timeseries + heatmap)...\n');
            gen_success = generate_flash_rf_images(cell_id, exp_folder, rf_path, ...
                flash_rf_path, flash_rf_heatmap_path, 0);
            if gen_success
                make_thumbnail(flash_rf_path, ...
                    fullfile(OUTPUT_DIR, 'thumbs', 'flash_rf', [cell_id '_flash_rf.png']), THUMB_WIDTH);
                make_thumbnail(flash_rf_heatmap_path, ...
                    fullfile(OUTPUT_DIR, 'thumbs', 'flash_rf_heatmap', [cell_id '_flash_rf_heatmap.png']), THUMB_WIDTH);
            end
            close all force
        else
            fprintf('  WARNING: No rf_results .mat file found for %s (%s). Run batch_generate_flash_results first.\n', cell_id, on_off);
        end
    elseif isfile(flash_rf_path) && isfile(flash_rf_heatmap_path)
        fprintf('  Flash RF images already exist, skipping\n');
    elseif isempty(exp_folder)
        fprintf('  No experiment folder found for flash RF\n');
    end

    %% 4. Bar flash — regenerate from saved bar flash data
    bar_flash_path = fullfile(OUTPUT_DIR, 'bar_flash', [cell_id '_bar_flash.png']);
    if ~isfile(bar_flash_path)
        % Find saved bar flash results — search multiple directories
        bf_path = find_bar_flash(BAR_FLASH_DATA, BAR_FLASH_ALT, cell_id, strain, on_off, date_str, time_str);

        if ~isempty(bf_path)
            try
                fprintf('  Regenerating bar flash from saved data...\n');
                S = load(bf_path, 'data_slow', 'mean_slow');
                % n_int = number of inter-flash interval samples (slow = 0.75 * 10000)
                n_int = 7500;
                fig = plot_bar_flash_data(S.data_slow, S.mean_slow, n_int);
                fig.Color = 'w';
                exportgraphics(fig, bar_flash_path, EXPORT_OPTS{:});
                make_thumbnail(bar_flash_path, ...
                    fullfile(OUTPUT_DIR, 'thumbs', 'bar_flash', [cell_id '_bar_flash.png']), THUMB_WIDTH);
            catch ME
                warning('  Bar flash failed for %s: %s', cell_id, ME.message);
            end
            close all force
        else
            fprintf('  No bar flash data found for %s\n', cell_id);
        end
    else
        fprintf('  Bar flash already exists, skipping\n');
    end
end

%% Build cell index JSON
fprintf('\nBuilding cell index...\n');
build_cell_index(cells, OUTPUT_DIR);

fprintf('\nDone! Dashboard images saved to:\n  %s\n', OUTPUT_DIR);

end


function make_thumbnail(src_path, thumb_path, max_width)
% MAKE_THUMBNAIL  Create a resized thumbnail of a PNG image.
    try
        img = imread(src_path);
        [h, w, ~] = size(img);
        if w > max_width
            scale = max_width / w;
            new_h = round(h * scale);
            thumb = imresize(img, [new_h, max_width]);
            imwrite(thumb, thumb_path);
        else
            copyfile(src_path, thumb_path);
        end
    catch
        % If imresize fails, just copy the original
        copyfile(src_path, thumb_path);
    end
end


function found_path = find_peak_vals(primary_dir, alt_dir, strain, on_off, date_str, time_str)
% FIND_PEAK_VALS  Search for peak_vals .mat file across multiple directories.
%
%   Searches in:
%     1. primary_dir (results/1DRF/bar_sweep_data/) — flat directory
%     2. alt_dir (results/bar_results/) — recursively, including summer2025/
%
%   Returns full path to file, or '' if not found.

    found_path = '';

    % 1. Search primary directory (exact pattern, then broad)
    pattern1 = fullfile(primary_dir, sprintf('peak_vals_%s_%s_%s_%s.mat', strain, on_off, date_str, time_str));
    files = dir(pattern1);
    if isempty(files)
        pattern1b = fullfile(primary_dir, sprintf('peak_vals_*%s*%s*%s*.mat', on_off, date_str, time_str));
        files = dir(pattern1b);
    end
    if ~isempty(files)
        found_path = fullfile(primary_dir, files(1).name);
        return;
    end

    % 2. Search alternative directory recursively
    if isfolder(alt_dir)
        pattern2 = fullfile(alt_dir, '**', sprintf('peak_vals_%s_%s_%s_%s.mat', strain, on_off, date_str, time_str));
        files = dir(pattern2);
        if isempty(files)
            % Broader pattern: match by date/time and on_off (strain may differ)
            pattern2b = fullfile(alt_dir, '**', sprintf('peak_vals_*_%s_%s_%s.mat', on_off, date_str, time_str));
            files = dir(pattern2b);
        end
        if ~isempty(files)
            found_path = fullfile(files(1).folder, files(1).name);
            return;
        end
    end
end


function found_path = find_bar_flash(primary_dir, alt_dir, cell_id, strain, on_off, date_str, time_str)
% FIND_BAR_FLASH  Search for bar_flash_results .mat file across directories.
%
%   Searches in:
%     1. primary_dir (results/1DRF/bar_flash_data/) — flat directory
%     2. alt_dir (results/bar_flash_results/) — including 1DRF/ subdirectory
%
%   Returns full path to file, or '' if not found.

    found_path = '';

    % 1. Search primary directory
    pattern1 = fullfile(primary_dir, sprintf('bar_flash_results_%s_%s_%s.mat', cell_id, strain, on_off));
    files = dir(pattern1);
    if isempty(files)
        pattern1b = fullfile(primary_dir, sprintf('bar_flash_results_%s*%s*.mat', date_str, on_off));
        files = dir(pattern1b);
    end
    if ~isempty(files)
        found_path = fullfile(primary_dir, files(1).name);
        return;
    end

    % 2. Search alternative directory recursively
    if isfolder(alt_dir)
        pattern2 = fullfile(alt_dir, '**', sprintf('bar_flash_results_%s_%s_%s.mat', cell_id, strain, on_off));
        files = dir(pattern2);
        if isempty(files)
            pattern2b = fullfile(alt_dir, '**', sprintf('bar_flash_results_%s*%s*.mat', date_str, on_off));
            files = dir(pattern2b);
        end
        if ~isempty(files)
            found_path = fullfile(files(1).folder, files(1).name);
            return;
        end
    end
end
