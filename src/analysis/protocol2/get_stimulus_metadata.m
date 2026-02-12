function stim_meta = get_stimulus_metadata(exp_folder)
% GET_STIMULUS_METADATA  Decode stimulus parameters from pattern files.
%
%   STIM_META = GET_STIMULUS_METADATA(EXP_FOLDER) loads the pattern files
%   from an experiment folder, empirically determines bar orientations and
%   motion directions from the pixel data, and returns a comprehensive
%   metadata struct.
%
%   INPUT:
%     exp_folder - Full path to an experiment directory containing:
%                  Patterns/, Functions/, and currentExp.mat
%
%   OUTPUT:
%     stim_meta - Struct with fields:
%       .bar_sweep      - Bar sweep metadata (see below)
%       .bar_flash      - Bar flash metadata (see below)
%       .cross_reference - Mapping between sweep and flash orientations
%       .convention     - String describing the direction convention
%
%   BAR SWEEP FIELDS (stim_meta.bar_sweep):
%     .direction_deg    - 32x1 motion direction for each raw data row
%     .direction_rad    - 32x1 same in radians
%     .orientation_deg  - 32x1 bar orientation [0, 180)
%     .orientation_rad  - 32x1 same in radians
%     .speed_label      - 32x1 string: '28dps' or '56dps'
%     .speed_dps        - 32x1 numeric speed
%     .pattern_file     - 32x1 string: pattern filename
%     .is_flip          - 32x1 logical: flip function used
%     .plot_order       - 16x1 reordering for sequential angles
%     .summary_table    - MATLAB table with all above
%
%   BAR FLASH FIELDS (stim_meta.bar_flash):
%     .frame_to_orientation_idx - 88x1 orientation index per frame
%     .frame_to_position_idx   - 88x1 position index per frame
%     .orientation_deg         - 8x1 bar orientation per index
%     .orientation_rad         - 8x1 same in radians
%     .raw_cell_map            - 11x8 struct documenting actual mapping
%
%   DIRECTION CONVENTION:
%     0° = East (vertical bar moving left to right)
%     90° = North (horizontal bar moving bottom to top)
%     Angles increase counterclockwise.
%     Orientations are in [0, 180), perpendicular to motion direction.
%
%   NOTE: MATLAB row 1 = bottom of LED screen (vertically flipped).
%         This function accounts for the flip when computing angles.
%
%   See also VERIFY_STIMULUS_METADATA, REINDEX_BAR_FLASH_DATA

    %% Load experiment info
    exp_mat = fullfile(exp_folder, 'currentExp.mat');
    assert(isfile(exp_mat), 'currentExp.mat not found in %s', exp_folder);
    loaded = load(exp_mat, 'currentExp', 'pattern_order', 'func_order');
    currentExp = loaded.currentExp;
    pattern_order = loaded.pattern_order;
    func_order = loaded.func_order;

    pattern_folder = fullfile(exp_folder, 'Patterns');
    function_folder = fullfile(exp_folder, 'Functions');

    %% Identify pattern and function files
    pat_mat_files = dir(fullfile(pattern_folder, '*.mat'));
    func_mat_files = dir(fullfile(function_folder, '*.mat'));

    %% Classify patterns: flash patterns vs bar patterns vs bar flash patterns
    n_patterns = numel(pat_mat_files);
    pat_names = {pat_mat_files.name};

    % Bar flash patterns contain 'FLASHES' in the name
    is_bar_flash_pat = contains(pat_names, 'FLASHES', 'IgnoreCase', true);
    % Flash patterns contain 'square_RF' or 'px_square' in the name
    is_flash_pat = contains(pat_names, 'square_RF', 'IgnoreCase', true) | ...
                   contains(pat_names, 'px_square', 'IgnoreCase', true);
    % Bar sweep patterns are the rest (excluding flash and bar flash)
    is_bar_pat = ~is_flash_pat & ~is_bar_flash_pat;

    bar_pat_indices = find(is_bar_pat);
    bar_flash_pat_indices = find(is_bar_flash_pat);

    %% Classify functions: flash, bar (forward/flip), bar flash, static
    func_names = {func_mat_files.name};
    is_flash_fn = contains(func_names, 'flashes', 'IgnoreCase', true) & ...
                  ~contains(func_names, 'bar_flash', 'IgnoreCase', true);
    is_bar_flash_fn = contains(func_names, 'bar_flash', 'IgnoreCase', true);
    is_static_fn = contains(func_names, 'static', 'IgnoreCase', true);
    is_bar_fn = ~is_flash_fn & ~is_bar_flash_fn & ~is_static_fn;

    bar_fn_indices = find(is_bar_fn);

    % Determine which bar functions are forward vs flip
    is_flip_fn = contains(func_names, 'FLIP', 'IgnoreCase', true);

    % Determine speed from function filenames
    func_speed_dps = zeros(1, numel(func_names));
    for f = 1:numel(func_names)
        fn = func_names{f};
        % Extract speed from pattern like "28dps" or "56dps"
        speed_match = regexp(fn, '(\d+)dps', 'tokens');
        if ~isempty(speed_match)
            func_speed_dps(f) = str2double(speed_match{1}{1});
        end
    end

    %% Analyze each bar sweep pattern: determine orientation and direction
    n_bar_pats = numel(bar_pat_indices);
    bar_orientations_deg = zeros(n_bar_pats, 1);
    bar_forward_dirs_deg = zeros(n_bar_pats, 1);
    bar_pat_names = cell(n_bar_pats, 1);

    for bp = 1:n_bar_pats
        idx = bar_pat_indices(bp);
        patt_file = fullfile(pattern_folder, pat_mat_files(idx).name);
        bar_pat_names{bp} = pat_mat_files(idx).name;

        patt_data = load(patt_file, 'pattern');
        Pats = patt_data.pattern.Pats;

        [orient_deg, fwd_dir_deg] = analyze_bar_pattern(Pats);
        bar_orientations_deg(bp) = orient_deg;
        bar_forward_dirs_deg(bp) = fwd_dir_deg;
    end

    %% Map data rows to stimuli using pattern_order and func_order
    % The bar portion of pattern_order starts after the flash entries.
    % From create_protocol2.m:
    %   pattern_order = [1, flash_patt, repmat(repelem(bar_patt, n_dir), [1, n_speeds])]
    %   func_order = [fn_static, fns_flash, repmat(fns_28dps, [1,n_bar_patt]), repmat(fns_56dps, [1,n_bar_patt])]
    %
    % We need to find which entries in pattern_order correspond to bar patterns.

    n_trials = numel(pattern_order);
    is_bar_trial = false(1, n_trials);
    for t = 1:n_trials
        pat_idx = pattern_order(t);
        if pat_idx >= 1 && pat_idx <= n_patterns
            is_bar_trial(t) = is_bar_pat(pat_idx);
        end
    end

    bar_trial_indices = find(is_bar_trial);
    n_bar_trials = numel(bar_trial_indices);

    % Build mapping for each bar trial (which becomes a data row)
    direction_deg = zeros(n_bar_trials, 1);
    direction_rad = zeros(n_bar_trials, 1);
    orientation_deg = zeros(n_bar_trials, 1);
    orientation_rad = zeros(n_bar_trials, 1);
    speed_dps = zeros(n_bar_trials, 1);
    speed_label = strings(n_bar_trials, 1);
    pattern_file = strings(n_bar_trials, 1);
    is_flip = false(n_bar_trials, 1);

    for row = 1:n_bar_trials
        trial_idx = bar_trial_indices(row);
        pat_idx = pattern_order(trial_idx);
        fn_idx = func_order(trial_idx);

        % Find which bar pattern this is (index into bar_pat_indices)
        bp = find(bar_pat_indices == pat_idx, 1);
        if isempty(bp)
            warning('Could not find bar pattern index %d in bar patterns list', pat_idx);
            continue;
        end

        orient = bar_orientations_deg(bp);
        fwd_dir = bar_forward_dirs_deg(bp);

        % Determine if flip
        flip = is_flip_fn(fn_idx);
        is_flip(row) = flip;

        if flip
            dir_deg = mod(fwd_dir + 180, 360);
        else
            dir_deg = fwd_dir;
        end

        direction_deg(row) = dir_deg;
        direction_rad(row) = deg2rad(dir_deg);
        orientation_deg(row) = orient;
        orientation_rad(row) = deg2rad(orient);
        speed_dps(row) = func_speed_dps(fn_idx);
        speed_label(row) = sprintf('%ddps', func_speed_dps(fn_idx));
        pattern_file(row) = bar_pat_names{bp};
    end

    %% Build bar_sweep struct
    stim_meta.bar_sweep.direction_deg = direction_deg;
    stim_meta.bar_sweep.direction_rad = direction_rad;
    stim_meta.bar_sweep.orientation_deg = orientation_deg;
    stim_meta.bar_sweep.orientation_rad = orientation_rad;
    stim_meta.bar_sweep.speed_dps = speed_dps;
    stim_meta.bar_sweep.speed_label = speed_label;
    stim_meta.bar_sweep.pattern_file = pattern_file;
    stim_meta.bar_sweep.is_flip = is_flip;

    % Build summary table
    Row = (1:n_bar_trials)';
    stim_meta.bar_sweep.summary_table = table(Row, direction_deg, direction_rad, ...
        orientation_deg, orientation_rad, speed_dps, speed_label, pattern_file, is_flip);

    % Determine plot_order: reorder first 16 rows (slow bars) to sequential angles
    if n_bar_trials >= 16
        slow_dirs = direction_deg(1:16);
        [~, plot_order] = sort(slow_dirs);
        stim_meta.bar_sweep.plot_order = plot_order;
    end

    %% Analyze bar flash pattern
    if ~isempty(bar_flash_pat_indices)
        bf_idx = bar_flash_pat_indices(1);
        bf_patt_file = fullfile(pattern_folder, pat_mat_files(bf_idx).name);
        bf_data = load(bf_patt_file, 'pattern');
        bf_Pats = bf_data.pattern.Pats;

        % Number of orientations in bar flash pattern = 8 (always)
        n_bar_flash_orientations = 8;
        [bf_orient_deg, bf_orient_rad, frame_to_orient, frame_to_pos] = ...
            analyze_bar_flash_pattern(bf_Pats, n_bar_flash_orientations);

        stim_meta.bar_flash.orientation_deg = bf_orient_deg;
        stim_meta.bar_flash.orientation_rad = bf_orient_rad;
        stim_meta.bar_flash.frame_to_orientation_idx = frame_to_orient;
        stim_meta.bar_flash.frame_to_position_idx = frame_to_pos;

        % Document the raw cell mapping (off-by-one due to frame 1 = grey)
        stim_meta.bar_flash.raw_cell_map = build_raw_cell_map(frame_to_orient, frame_to_pos, bf_orient_deg);
    end

    %% Cross-reference bar flash orientations with bar sweep orientations
    if isfield(stim_meta, 'bar_flash') && ~isempty(stim_meta.bar_flash.orientation_deg)
        [flash_to_sweep, sweep_to_flash] = cross_reference_orientations( ...
            stim_meta.bar_flash.orientation_deg, ...
            stim_meta.bar_sweep.orientation_deg);

        stim_meta.cross_reference.flash_orient_to_sweep_rows = flash_to_sweep;
        stim_meta.cross_reference.sweep_row_to_flash_orient = sweep_to_flash;
    end

    %% Convention string
    stim_meta.convention = sprintf([...
        'Direction convention:\n' ...
        '  0° = East (vertical bar, left to right)\n' ...
        '  90° = North (horizontal bar, bottom to top)\n' ...
        '  Angles increase counterclockwise.\n' ...
        '  Orientation = bar long axis angle in [0, 180).\n' ...
        '  Direction = perpendicular to orientation, indicating motion.\n' ...
        'Coordinate system: as seen on LED arena screen\n' ...
        '  (MATLAB row 1 = bottom of screen, flipped from imagesc default).']);

end


%% ========== LOCAL FUNCTIONS ==========

function [orient_deg, fwd_dir_deg] = analyze_bar_pattern(Pats)
% ANALYZE_BAR_PATTERN  Determine orientation and forward direction of a bar pattern.
%
%   Uses PCA on non-background pixels to find bar orientation, and
%   centroid tracking across frames to find direction of motion.

    % Handle compressed patterns (192 x 3 x N → stretch to 48 x 192 x N)
    Pats = stretch_pattern_if_needed(Pats);

    n_rows = size(Pats, 1);
    n_frames = size(Pats, 3);

    % Background value: mode of frame 1
    bkg = mode(double(Pats(:, :, 1)), 'all');

    % Choose two frames for centroid comparison
    frame_early = min(20, floor(n_frames * 0.15));
    frame_late = min(50, floor(n_frames * 0.4));
    frame_early = max(frame_early, 2);
    frame_late = max(frame_late, frame_early + 5);

    % Find non-background pixels in each frame
    [orient_deg, ~] = compute_orientation_from_frame(Pats(:, :, round((frame_early + frame_late) / 2)), bkg, n_rows);

    % Compute centroids for direction
    [cy_early, cx_early] = compute_centroid(Pats(:, :, frame_early), bkg, n_rows);
    [cy_late, cx_late] = compute_centroid(Pats(:, :, frame_late), bkg, n_rows);

    % Direction vector (in screen coordinates)
    dx = cx_late - cx_early;
    dy = cy_late - cy_early; % Already in screen-y (flipped)

    fwd_dir_deg = mod(atan2d(dy, dx), 360);
end


function [orient_deg, orient_rad] = compute_orientation_from_frame(frame, bkg, n_rows)
% COMPUTE_ORIENTATION_FROM_FRAME  Use PCA on non-background pixels.

    mask = double(frame) ~= bkg;
    [row_idx, col_idx] = find(mask);

    if numel(row_idx) < 3
        orient_deg = 0;
        orient_rad = 0;
        return;
    end

    % Convert to screen coordinates (flip y-axis)
    screen_y = n_rows - row_idx + 1;
    screen_x = col_idx;

    % PCA
    coords = [screen_x, screen_y];
    coords_centered = coords - mean(coords);
    [~, ~, V] = svd(coords_centered, 'econ');

    % First principal component = bar long axis
    pc1 = V(:, 1);
    angle = atan2d(pc1(2), pc1(1));

    % Normalize to [0, 180) since orientation is undirected
    orient_deg = mod(angle, 180);
    orient_rad = deg2rad(orient_deg);
end


function [cy, cx] = compute_centroid(frame, bkg, n_rows)
% COMPUTE_CENTROID  Centroid of non-background pixels in screen coordinates.

    mask = double(frame) ~= bkg;
    [row_idx, col_idx] = find(mask);

    if isempty(row_idx)
        cy = 0;
        cx = 0;
        return;
    end

    % Convert to screen coordinates
    screen_y = n_rows - row_idx + 1;
    cx = mean(col_idx);
    cy = mean(screen_y);
end


function [orient_deg, orient_rad, frame_to_orient, frame_to_pos] = ...
        analyze_bar_flash_pattern(Pats, n_orientations)
% ANALYZE_BAR_FLASH_PATTERN  Determine orientation for each frame group.
%
%   The bar flash pattern has frame 1 = grey, then groups of 11 frames
%   (one per position) for each of 8 orientations.

    Pats = stretch_pattern_if_needed(Pats);
    n_rows = size(Pats, 1);
    n_frames = size(Pats, 3);
    n_positions = 11;

    bkg = mode(double(Pats(:, :, 1)), 'all');

    orient_deg = zeros(n_orientations, 1);
    orient_rad = zeros(n_orientations, 1);
    frame_to_orient = zeros(n_frames - 1, 1); % Frames 2 to n_frames
    frame_to_pos = zeros(n_frames - 1, 1);

    for o = 1:n_orientations
        % Frames for this orientation (11 frames per orientation, starting at frame 2)
        frame_start = 2 + (o - 1) * n_positions;
        frame_end = frame_start + n_positions - 1;

        if frame_end > n_frames
            break;
        end

        % Use the middle frame of this orientation group for analysis
        mid_frame = round((frame_start + frame_end) / 2);
        [od, or] = compute_orientation_from_frame(Pats(:, :, mid_frame), bkg, n_rows);
        orient_deg(o) = od;
        orient_rad(o) = or;

        % Map frame numbers to orientation and position indices
        for p = 1:n_positions
            frame_idx = frame_start + p - 1;
            if frame_idx <= n_frames
                frame_to_orient(frame_idx - 1) = o;  % -1 because we skip frame 1
                frame_to_pos(frame_idx - 1) = p;
            end
        end
    end
end


function raw_cell_map = build_raw_cell_map(frame_to_orient, frame_to_pos, orient_deg)
% BUILD_RAW_CELL_MAP  Document what each (row,col) in the stored 11x8 cell represents.
%
%   Due to the off-by-one (frame 1 = grey, flash frames start at 2),
%   data_rep{flash_frame_num} stores data at linear index = flash_frame_num
%   in an 11x8 cell. This means the mapping is shifted.

    n_positions = 11;
    n_orientations = 8;

    raw_cell_map = struct('row', {}, 'col', {}, 'frame_num', {}, ...
        'actual_orientation_idx', {}, 'actual_position_idx', {}, ...
        'actual_orientation_deg', {});

    for col = 1:n_orientations
        for row = 1:n_positions
            linear_idx = (col - 1) * n_positions + row; % MATLAB column-major
            entry = struct();
            entry.row = row;
            entry.col = col;
            entry.linear_idx = linear_idx;

            if linear_idx == 1
                % Frame 1 = grey, never stored
                entry.frame_num = 1;
                entry.actual_orientation_idx = NaN;
                entry.actual_position_idx = NaN;
                entry.actual_orientation_deg = NaN;
                entry.status = 'empty (grey frame)';
            elseif linear_idx <= numel(frame_to_orient) + 1
                frame_num = linear_idx; % frame number = linear index
                fi = linear_idx - 1; % index into frame_to_orient (which skips frame 1)
                entry.frame_num = frame_num;
                if fi >= 1 && fi <= numel(frame_to_orient)
                    entry.actual_orientation_idx = frame_to_orient(fi);
                    entry.actual_position_idx = frame_to_pos(fi);
                    if frame_to_orient(fi) >= 1 && frame_to_orient(fi) <= numel(orient_deg)
                        entry.actual_orientation_deg = orient_deg(frame_to_orient(fi));
                    else
                        entry.actual_orientation_deg = NaN;
                    end
                    entry.status = 'valid';
                else
                    entry.actual_orientation_idx = NaN;
                    entry.actual_position_idx = NaN;
                    entry.actual_orientation_deg = NaN;
                    entry.status = 'out of range';
                end
            else
                entry.frame_num = linear_idx;
                entry.actual_orientation_idx = NaN;
                entry.actual_position_idx = NaN;
                entry.actual_orientation_deg = NaN;
                entry.status = 'out of range';
            end

            raw_cell_map(row, col) = entry;
        end
    end
end


function [flash_to_sweep, sweep_to_flash] = cross_reference_orientations(flash_orient_deg, sweep_orient_deg)
% CROSS_REFERENCE_ORIENTATIONS  Match bar flash and bar sweep orientations.

    n_flash = numel(flash_orient_deg);
    n_sweep = numel(sweep_orient_deg);

    flash_to_sweep = cell(n_flash, 1);
    sweep_to_flash = zeros(n_sweep, 1);

    for fo = 1:n_flash
        % Find sweep rows with matching orientation (within tolerance)
        angle_diff = abs(sweep_orient_deg - flash_orient_deg(fo));
        % Handle wrapping at 180
        angle_diff = min(angle_diff, 180 - angle_diff);
        matching_rows = find(angle_diff < 10); % 10 degree tolerance
        flash_to_sweep{fo} = matching_rows;
    end

    for sr = 1:n_sweep
        angle_diff = abs(flash_orient_deg - sweep_orient_deg(sr));
        angle_diff = min(angle_diff, 180 - angle_diff);
        [~, best_match] = min(angle_diff);
        if angle_diff(best_match) < 10
            sweep_to_flash(sr) = best_match;
        end
    end
end


function Pats = stretch_pattern_if_needed(Pats)
% STRETCH_PATTERN_IF_NEEDED  Handle compressed 192x3xN patterns.
%
%   If the pattern is stored as 192x3xN (compressed), stretch each of
%   the 3 rows to 16 rows to get a 48x192xN array.

    [d1, d2, d3] = size(Pats);

    if d1 == 192 && d2 == 3
        % Compressed format: 192 x 3 x N → need to transpose and stretch
        % Actually this means rows=192 (width), cols=3 (compressed height)
        stretched = zeros(48, 192, d3);
        for f = 1:d3
            frame = Pats(:, :, f); % 192 x 3
            frame_t = frame'; % 3 x 192
            % Stretch 3 rows to 48 rows (each row → 16 rows)
            for r = 1:3
                row_start = (r - 1) * 16 + 1;
                row_end = r * 16;
                stretched(row_start:row_end, :, f) = repmat(frame_t(r, :), 16, 1);
            end
        end
        Pats = stretched;
    elseif d2 == 192 && d1 == 3
        % Also compressed but transposed: 3 x 192 x N
        stretched = zeros(48, 192, d3);
        for f = 1:d3
            frame = Pats(:, :, f); % 3 x 192
            for r = 1:3
                row_start = (r - 1) * 16 + 1;
                row_end = r * 16;
                stretched(row_start:row_end, :, f) = repmat(frame(r, :), 16, 1);
            end
        end
        Pats = stretched;
    end
    % If already 48x192xN, no change needed
end
