function batch_generate_flash_results()
% BATCH_GENERATE_FLASH_RESULTS  Generate missing analysis results.
%
%   BATCH_GENERATE_FLASH_RESULTS() scans T4T5 Summer and Autumn 2025
%   experiment folders and runs analysis pipelines for experiments that
%   are missing results:
%
%   Phase 1 — experiments missing rf_results:
%     1. Runs process_bars_p2 to get the resultant_angle
%     2. Runs process_flash_p2 to generate rf_results and figures
%     3. Runs process_bar_flashes_p2 to generate bar flash results
%
%   Phase 2 — experiments that have rf_results but missing bar_flash:
%     1. Runs process_bar_flashes_p2 to generate bar flash results
%
%   PREREQUISITES:
%     >> addpath(genpath('/Users/burnettl/Documents/GitHub/nested_RF_stimulus/src'))
%
%   See also PROCESS_FLASH_P2, PROCESS_BARS_P2, PROCESS_BAR_FLASHES_P2

%% Configuration
PROJECT_ROOT = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2';
P2_DIRS = {
    fullfile(PROJECT_ROOT, 'data', 'T4T5_Summer2025')
    fullfile(PROJECT_ROOT, 'data', 'T4T5_autumn2025')
};
FLASH_RESULTS = fullfile(PROJECT_ROOT, 'results', 'flash_results');

%% Collect all experiment folders
exp_folders = {};
for d = 1:numel(P2_DIRS)
    contents = dir(P2_DIRS{d});
    contents = contents([contents.isdir]);
    contents = contents(~ismember({contents.name}, {'.', '..'}));
    for k = 1:numel(contents)
        exp_folders{end+1} = fullfile(P2_DIRS{d}, contents(k).name); %#ok<AGROW>
    end
end

fprintf('Found %d experiment folders\n', numel(exp_folders));

%% Identify experiments missing rf_results
existing_rf = dir(fullfile(FLASH_RESULTS, 'rf_results_*.mat'));
existing_rf_names = {existing_rf.name};

missing = {};
for i = 1:numel(exp_folders)
    folder = exp_folders{i};
    [~, folder_name] = fileparts(folder);

    % Check if any rf_results file matches this experiment date/time
    has_results = false;
    for j = 1:numel(existing_rf_names)
        if contains(existing_rf_names{j}, folder_name)
            has_results = true;
            break;
        end
    end

    if ~has_results
        missing{end+1} = folder; %#ok<AGROW>
    end
end

fprintf('%d experiments missing rf_results (out of %d total)\n', numel(missing), numel(exp_folders));

%% Process each missing experiment
orig_dir = pwd;
n_success = 0;
n_fail = 0;

if isempty(missing)
    fprintf('Phase 1: Nothing to do — all experiments have rf_results files.\n');
else

for i = 1:numel(missing)
    folder = missing{i};
    [~, folder_name] = fileparts(folder);
    fprintf('\n[%d/%d] Processing %s\n', i, numel(missing), folder_name);

    % Verify currentExp.mat exists
    meta_path = fullfile(folder, 'currentExp.mat');
    if ~isfile(meta_path)
        fprintf('  SKIP: no currentExp.mat found\n');
        n_fail = n_fail + 1;
        continue;
    end

    try
        % Restore working directory on exit (load_protocol2_data uses cd)
        cleanup = onCleanup(@() cd(orig_dir));

        % Load metadata
        S = load(meta_path, 'metadata');
        metadata = S.metadata;
        fprintf('  Strain: %s\n', metadata.Strain);

        % Step 1: Run bar analysis to get preferred direction
        fprintf('  Running bar analysis...\n');
        resultant_angle = process_bars_p2(folder, metadata, PROJECT_ROOT);
        close all force

        % Step 2: Run flash analysis
        fprintf('  Running flash analysis...\n');
        process_flash_p2(folder, metadata, PROJECT_ROOT, resultant_angle);
        close all force

        % Step 3: Run bar flash analysis
        fprintf('  Running bar flash analysis...\n');
        process_bar_flashes_p2(folder, metadata, PROJECT_ROOT);
        close all force

        fprintf('  Done.\n');
        n_success = n_success + 1;

    catch ME
        fprintf('  FAILED: %s\n', ME.message);
        n_fail = n_fail + 1;
        cd(orig_dir);
        close all force
    end
end

fprintf('\n=== Phase 1 Summary ===\n');
fprintf('  Processed: %d\n', n_success);
fprintf('  Failed:    %d\n', n_fail);
fprintf('  Skipped:   %d\n', numel(missing) - n_success - n_fail);
end  % if ~isempty(missing)

%% Phase 2: Generate missing bar flash results for experiments that already
%  have rf_results but are missing bar_flash_results (e.g., Summer cells
%  that were processed before bar flash analysis was added).
BAR_FLASH_RESULTS = fullfile(PROJECT_ROOT, 'results', 'bar_flash_results');
if ~isfolder(BAR_FLASH_RESULTS)
    mkdir(BAR_FLASH_RESULTS);
end

existing_bf = dir(fullfile(BAR_FLASH_RESULTS, 'bar_flash_results_*.mat'));
% Also check 1DRF subdirectory
existing_bf_1drf = dir(fullfile(BAR_FLASH_RESULTS, '1DRF', 'bar_flash_results_*.mat'));
existing_bf_names = [reshape({existing_bf.name}, 1, []), reshape({existing_bf_1drf.name}, 1, [])];

missing_bf = {};
for i = 1:numel(exp_folders)
    folder = exp_folders{i};
    [~, folder_name] = fileparts(folder);

    has_bf = false;
    for j = 1:numel(existing_bf_names)
        if contains(existing_bf_names{j}, folder_name)
            has_bf = true;
            break;
        end
    end

    if ~has_bf
        % Also check if we just processed this in Phase 1 (skip if so)
        was_in_phase1 = false;
        for m = 1:numel(missing)
            if strcmp(folder, missing{m})
                was_in_phase1 = true;
                break;
            end
        end
        if ~was_in_phase1
            missing_bf{end+1} = folder; %#ok<AGROW>
        end
    end
end

fprintf('\n%d experiments missing bar_flash_results (excluding Phase 1)\n', numel(missing_bf));

n_bf_success = 0;
n_bf_fail = 0;

for i = 1:numel(missing_bf)
    folder = missing_bf{i};
    [~, folder_name] = fileparts(folder);
    fprintf('\n[Phase 2: %d/%d] Bar flash for %s\n', i, numel(missing_bf), folder_name);

    meta_path = fullfile(folder, 'currentExp.mat');
    if ~isfile(meta_path)
        fprintf('  SKIP: no currentExp.mat found\n');
        n_bf_fail = n_bf_fail + 1;
        continue;
    end

    try
        cleanup = onCleanup(@() cd(orig_dir));
        S = load(meta_path, 'metadata');
        metadata = S.metadata;
        fprintf('  Strain: %s\n', metadata.Strain);

        process_bar_flashes_p2(folder, metadata, PROJECT_ROOT);
        close all force

        fprintf('  Done.\n');
        n_bf_success = n_bf_success + 1;
    catch ME
        fprintf('  FAILED: %s\n', ME.message);
        n_bf_fail = n_bf_fail + 1;
        cd(orig_dir);
        close all force
    end
end

fprintf('\n=== Phase 2 Summary (Bar Flash) ===\n');
fprintf('  Processed: %d\n', n_bf_success);
fprintf('  Failed:    %d\n', n_bf_fail);

end
