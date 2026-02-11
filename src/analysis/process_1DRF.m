function process_1DRF()
    % Batch processing script for all 1DRF experiments.
    %
    % Iterates over all experiment folders in the 1DRF data directory. For
    % each experiment:
    %   1. Runs bar sweep analysis and determines preferred direction
    %   2. Creates a GIF of the preferred direction bar sweep stimulus
    %   3. Plots bar flash spatial profiles along the PD axis
    %
    % Saves per-cell figures in individual folders and a summary .mat file.

    %% Pathsdebug_PD_mapping('/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/1DRF/2025_10_23_10_31')

    REPO_ROOT = '/Users/burnettl/Documents/GitHub/nested_RF_stimulus';
    DATA_ROOT = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/1DRF';
    FIGURES_ROOT = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/1DRF/figures';

    % Add analysis code to path
    addpath(genpath(fullfile(REPO_ROOT, 'src', 'analysis')));

    % Create output directory
    if ~isfolder(FIGURES_ROOT)
        mkdir(FIGURES_ROOT);
    end

    %% Get experiment folders
    exp_dirs = dir(DATA_ROOT);
    exp_dirs = exp_dirs([exp_dirs.isdir]);
    exp_dirs = exp_dirs(~ismember({exp_dirs.name}, {'.', '..'}));

    n_experiments = numel(exp_dirs);
    fprintf('Found %d experiment folders in %s\n\n', n_experiments, DATA_ROOT);

    %% Process each experiment
    results = cell(n_experiments, 1);
    failed = {};

    for i = 1:n_experiments
        exp_folder = fullfile(DATA_ROOT, exp_dirs(i).name);
        fprintf('=== Experiment %d/%d: %s ===\n', i, n_experiments, exp_dirs(i).name);

        try
            results{i} = process_1DRF_single(exp_folder, FIGURES_ROOT);
        catch ME
            warning('Failed on %s: %s', exp_dirs(i).name, ME.message);
            fprintf('  Stack trace:\n');
            for s = 1:numel(ME.stack)
                fprintf('    %s (line %d)\n', ME.stack(s).name, ME.stack(s).line);
            end
            fprintf('\n');
            failed{end+1} = exp_dirs(i).name;
            results{i} = [];
        end
    end

    %% Build summary table
    valid_idx = ~cellfun(@isempty, results);
    valid_results = results(valid_idx);

    if ~isempty(valid_results)
        summary = struct();
        summary.date = cellfun(@(r) r.date, valid_results, 'UniformOutput', false);
        summary.time = cellfun(@(r) r.time, valid_results, 'UniformOutput', false);
        summary.strain = cellfun(@(r) r.strain, valid_results, 'UniformOutput', false);
        summary.on_off = cellfun(@(r) r.on_off, valid_results, 'UniformOutput', false);
        summary.pd_angle_deg = cellfun(@(r) r.pd_angle_deg, valid_results);
        summary.pd_speed_label = cellfun(@(r) r.pd_speed_label, valid_results, 'UniformOutput', false);
        summary.pd_max_v = cellfun(@(r) r.pd_max_v, valid_results);
        summary.median_v = cellfun(@(r) r.median_v, valid_results);
        summary.bf_col = cellfun(@(r) r.bf_col, valid_results);
        summary.vector_sum_magnitude_slow = cellfun(@(r) r.vector_sum_magnitude_slow, valid_results);
        summary.vector_sum_angle_slow = cellfun(@(r) r.vector_sum_angle_slow, valid_results);
        summary.DSI_slow = cellfun(@(r) r.DSI_slow, valid_results);
        summary.DSI_pdnd_slow = cellfun(@(r) r.DSI_pdnd_slow, valid_results);
        summary.sym_ratio_slow = cellfun(@(r) r.sym_ratio_slow, valid_results);
        summary.fwhm_slow = cellfun(@(r) r.fwhm_slow, valid_results);
        summary.cv_slow = cellfun(@(r) r.cv_slow, valid_results);
        summary.n_processed = numel(valid_results);
        summary.n_failed = numel(failed);
        summary.failed_folders = failed;

        % Save summary
        summary_file = fullfile(FIGURES_ROOT, 'summary_1DRF.mat');
        save(summary_file, 'summary', 'results');
        fprintf('\nSummary saved: %s\n', summary_file);
    end

    %% Print summary
    fprintf('\n========================================\n');
    fprintf('Processing complete.\n');
    fprintf('  Successful: %d / %d\n', sum(valid_idx), n_experiments);
    fprintf('  Failed:     %d\n', numel(failed));
    if ~isempty(failed)
        fprintf('  Failed folders:\n');
        for j = 1:numel(failed)
            fprintf('    - %s\n', failed{j});
        end
    end
    fprintf('  Figures saved to: %s\n', FIGURES_ROOT);
    fprintf('========================================\n');

end
