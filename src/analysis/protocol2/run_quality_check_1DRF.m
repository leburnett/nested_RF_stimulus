function summary = run_quality_check_1DRF(data_root, opts)
% RUN_QUALITY_CHECK_1DRF  Batch quality check for all 1DRF experiments.
%
%   SUMMARY = RUN_QUALITY_CHECK_1DRF(DATA_ROOT) recursively finds all
%   experiment folders under DATA_ROOT, runs quality_check_recording on
%   each, and returns a summary table.
%
%   SUMMARY = RUN_QUALITY_CHECK_1DRF(DATA_ROOT, OPTS) passes options
%   through to quality_check_recording and controls batch behavior.
%
%   INPUTS:
%     data_root - Path to the 1DRF data directory
%                 (e.g., '/path/to/protocol2/data/1DRF')
%     opts      - (Optional) structure with fields:
%                   .save_summary  - Save summary table to .mat
%                                    (default: true)
%                   .summary_dir   - Where to save summary
%                                    (default: data_root)
%                   Plus all opts from quality_check_recording
%
%   OUTPUT:
%     summary - Table with one row per experiment containing: folder,
%               date_str, time_str, on_off, group, FLAG, mean_v, std_v,
%               median_v, mad_v, drift_v, range_v, p2_v, p98_v, var_v,
%               fail_reasons
%
%   EXAMPLE:
%     data_root = '/path/to/protocol2/data/1DRF';
%     summary = run_quality_check_1DRF(data_root);
%
%   See also QUALITY_CHECK_RECORDING, PROCESS_PROTOCOL2

    if nargin < 2, opts = struct(); end
    opts = set_batch_quality_defaults(opts, data_root);

    %% Discover experiment folders
    file_list = dir(fullfile(data_root, '**', 'currentExp.mat'));
    all_exp_folders = {file_list.folder};
    n_exp = numel(all_exp_folders);
    fprintf('Found %d experiment folders under %s\n', n_exp, data_root);

    if n_exp == 0
        summary = table();
        return;
    end

    %% Process each experiment
    results = cell(n_exp, 1);

    for exp_idx = 1:n_exp
        exp_folder = all_exp_folders{exp_idx};
        fprintf('\n[%d/%d] %s\n', exp_idx, n_exp, exp_folder);

        try
            quality = quality_check_recording(exp_folder, opts);
            r = build_result_row(exp_folder, quality);
            results{exp_idx} = r;
        catch ME
            fprintf('  ERROR: %s\n', ME.message);
            r = build_error_row(exp_folder, ME.message);
            results{exp_idx} = r;
        end
    end

    %% Build summary table
    results = vertcat(results{:});
    summary = struct2table(results, 'AsArray', true);

    %% Print summary
    valid_mask = ~isnan(summary.FLAG);
    n_pass = sum(summary.FLAG == 0 & valid_mask);
    n_fail = sum(summary.FLAG == 1 & valid_mask);
    n_error = sum(~valid_mask);

    fprintf('\n=== Quality Check Summary ===\n');
    fprintf('  Pass:  %d / %d\n', n_pass, n_exp);
    fprintf('  Fail:  %d / %d\n', n_fail, n_exp);
    if n_error > 0
        fprintf('  Error: %d / %d\n', n_error, n_exp);
    end

    %% Save summary
    if opts.save_summary
        if ~isfolder(opts.summary_dir)
            mkdir(opts.summary_dir);
        end
        save_path = fullfile(opts.summary_dir, 'quality_check_summary.mat');
        save(save_path, 'summary');
        fprintf('\nSummary saved to: %s\n', save_path);
    end

end


%% ========================= Helper Functions ==============================

function r = build_result_row(exp_folder, quality)
% BUILD_RESULT_ROW  Create a result struct from quality check output.

    [date_str, time_str, on_off, group] = parse_path_metadata(exp_folder);

    r.folder       = string(exp_folder);
    r.date_str     = string(date_str);
    r.time_str     = string(time_str);
    r.on_off       = string(on_off);
    r.group        = string(group);
    r.FLAG         = quality.FLAG;
    r.mean_v       = quality.mean_v;
    r.std_v        = quality.std_v;
    r.median_v     = quality.median_v;
    r.mad_v        = quality.mad_v;
    r.drift_v      = quality.drift_v;
    r.range_v      = quality.range_v;
    r.p2_v         = quality.p2_v;
    r.p98_v        = quality.p98_v;
    r.var_v        = quality.var_v;
    r.fail_reasons = {quality.fail_reasons};

end


function r = build_error_row(exp_folder, err_msg)
% BUILD_ERROR_ROW  Create a NaN-filled result struct for failed experiments.

    [date_str, time_str, on_off, group] = parse_path_metadata(exp_folder);

    r.folder       = string(exp_folder);
    r.date_str     = string(date_str);
    r.time_str     = string(time_str);
    r.on_off       = string(on_off);
    r.group        = string(group);
    r.FLAG         = NaN;
    r.mean_v       = NaN;
    r.std_v        = NaN;
    r.median_v     = NaN;
    r.mad_v        = NaN;
    r.drift_v      = NaN;
    r.range_v      = NaN;
    r.p2_v         = NaN;
    r.p98_v        = NaN;
    r.var_v        = NaN;
    r.fail_reasons = {{err_msg}};

end


function [date_str, time_str, on_off, group] = parse_path_metadata(exp_folder)
% PARSE_PATH_METADATA  Extract metadata from the nested folder path.
%
%   Path structure: .../1DRF/<protocol>/<group>/<ON|OFF>/YYYY_MM_DD_HH_MM

    parts = split(exp_folder, filesep);
    exp_name = parts{end};

    % Date and time from folder name
    date_str = exp_name(1:10);
    time_str = exp_name(12:16);

    % ON/OFF and control/ttl from parent directories
    on_off = parts{end-1};
    group  = parts{end-2};

end


%% ========================= Default Options ==============================

function opts = set_batch_quality_defaults(opts, data_root)
% SET_BATCH_QUALITY_DEFAULTS  Fill in default values for batch options.

    if ~isfield(opts, 'save_summary'), opts.save_summary = true;      end
    if ~isfield(opts, 'summary_dir'),  opts.summary_dir = data_root;  end

end
