function n_deleted = delete_tdms_files(data_root, opts)
% DELETE_TDMS_FILES  Remove .tdms and .tdms_index files from experiment folders.
%
%   N_DELETED = DELETE_TDMS_FILES(DATA_ROOT) recursively finds all .tdms
%   and .tdms_index files under DATA_ROOT and deletes them. Defaults to
%   dry-run mode (lists files without deleting).
%
%   N_DELETED = DELETE_TDMS_FILES(DATA_ROOT, OPTS) uses options to control
%   behavior:
%     opts.dry_run - If true (default), list files without deleting.
%                    Set to false to actually delete files.
%
%   INPUTS:
%     data_root - Path to the data directory containing experiment folders
%                 (e.g., '/path/to/protocol2/data/1DRF')
%     opts      - (Optional) structure with fields:
%                   .dry_run - Logical, default true (safe mode)
%
%   OUTPUT:
%     n_deleted - Number of files deleted (or that would be deleted in
%                 dry-run mode)
%
%   EXAMPLE:
%     data_root = '/path/to/protocol2/data/1DRF';
%
%     % Preview what would be deleted:
%     delete_tdms_files(data_root);
%
%     % Actually delete:
%     opts.dry_run = false;
%     n = delete_tdms_files(data_root, opts);
%
%   See also RUN_QUALITY_CHECK_1DRF, QUALITY_CHECK_RECORDING

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'dry_run'), opts.dry_run = true; end

    %% Find all TDMS files recursively
    tdms_files = dir(fullfile(data_root, '**', '*.tdms'));
    tdms_idx_files = dir(fullfile(data_root, '**', '*.tdms_index'));
    all_files = [tdms_files; tdms_idx_files];
    n_deleted = numel(all_files);

    if n_deleted == 0
        fprintf('No .tdms or .tdms_index files found under %s\n', data_root);
        return;
    end

    %% Calculate total size
    total_bytes = sum([all_files.bytes]);

    if opts.dry_run
        fprintf('[DRY RUN] Found %d TDMS files (%.1f GB) under %s\n', ...
            n_deleted, total_bytes / 1e9, data_root);
        fprintf('[DRY RUN] Set opts.dry_run = false to delete.\n\n');

        for f = 1:numel(all_files)
            fprintf('  %s (%.1f MB)\n', ...
                fullfile(all_files(f).folder, all_files(f).name), ...
                all_files(f).bytes / 1e6);
        end
    else
        fprintf('Deleting %d TDMS files (%.1f GB) under %s\n', ...
            n_deleted, total_bytes / 1e9, data_root);

        n_success = 0;
        n_fail = 0;
        for f = 1:numel(all_files)
            fpath = fullfile(all_files(f).folder, all_files(f).name);
            try
                delete(fpath);
                n_success = n_success + 1;
            catch ME
                fprintf('  WARNING: Could not delete %s: %s\n', fpath, ME.message);
                n_fail = n_fail + 1;
            end
        end

        fprintf('Deleted %d files. ', n_success);
        if n_fail > 0
            fprintf('%d files could not be deleted.\n', n_fail);
        else
            fprintf('%.1f GB freed.\n', total_bytes / 1e9);
        end
        n_deleted = n_success;
    end

end
