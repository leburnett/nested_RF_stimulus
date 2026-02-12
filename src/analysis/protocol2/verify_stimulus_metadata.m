function verify_stimulus_metadata(exp_folder)
% VERIFY_STIMULUS_METADATA  Generate diagnostic plots to verify stimulus mappings.
%
%   VERIFY_STIMULUS_METADATA(EXP_FOLDER) runs GET_STIMULUS_METADATA on the
%   given experiment folder and produces 7 diagnostic figures plus console
%   output to verify the correctness of all stimulus-to-data mappings.
%
%   INPUT:
%     exp_folder - Full path to a Protocol 2 experiment directory
%
%   FIGURES GENERATED:
%     Figure 1 - Bar pattern frames with annotated orientation/direction
%     Figure 2 - Compass diagram of bar sweep directions
%     Figure 3 - Bar flash orientation verification
%     Figure 4 - Cross-reference: flash vs sweep pattern comparison
%     Figure 5 - Bar flash epoch detection overlay on frame data
%     Figure 6 - Bar flash frame number distribution
%     Figure 7 - Bar flash timeseries sanity check (8x11 grid)
%
%   See also GET_STIMULUS_METADATA, REINDEX_BAR_FLASH_DATA, PARSE_BAR_FLASH_DATA

    fprintf('=== Stimulus Metadata Verification ===\n');
    fprintf('Experiment: %s\n\n', exp_folder);

    %% Run the metadata decoder
    stim_meta = get_stimulus_metadata(exp_folder);

    %% Print summary table
    fprintf('--- Bar Sweep Summary ---\n');
    disp(stim_meta.bar_sweep.summary_table);

    fprintf('\n--- Direction Convention ---\n');
    fprintf('%s\n\n', stim_meta.convention);

    %% Check forward/flip pairs
    fprintf('--- Forward/Flip Pair Check ---\n');
    n_rows = numel(stim_meta.bar_sweep.direction_deg);
    for row = 1:2:min(n_rows, 16)
        if row + 1 <= n_rows
            dir1 = stim_meta.bar_sweep.direction_deg(row);
            dir2 = stim_meta.bar_sweep.direction_deg(row + 1);
            diff_angle = abs(mod(dir2 - dir1 + 180, 360) - 180);
            ok_str = '';
            if abs(diff_angle - 180) > 5
                ok_str = ' *** WARNING: not 180° apart!';
            end
            fprintf('  Rows %2d-%2d: %.1f° and %.1f° (diff=%.1f°)%s\n', ...
                row, row+1, dir1, dir2, diff_angle, ok_str);
        end
    end
    fprintf('\n');

    %% Cross-reference table
    if isfield(stim_meta, 'cross_reference')
        fprintf('--- Cross-Reference: Flash Orientation → Sweep Rows ---\n');
        for fo = 1:numel(stim_meta.bar_flash.orientation_deg)
            matching = stim_meta.cross_reference.flash_orient_to_sweep_rows{fo};
            fprintf('  Flash orient %d (%.1f°) → Sweep rows: [%s]\n', ...
                fo, stim_meta.bar_flash.orientation_deg(fo), ...
                strjoin(string(matching'), ', '));
        end
        fprintf('\n');
    end

    %% Build green LED colormap
    green_cmap = zeros(256, 3);
    green_cmap(:, 2) = linspace(0, 1, 256); % Black to green

    %% FIGURE 1: Bar pattern frames with annotated orientation/direction
    pattern_folder = fullfile(exp_folder, 'Patterns');
    pat_mat_files = dir(fullfile(pattern_folder, '*.mat'));
    pat_names = {pat_mat_files.name};

    % Find bar patterns (not flash, not bar flash)
    is_bar_flash_pat = contains(pat_names, 'FLASHES', 'IgnoreCase', true);
    is_flash_pat = contains(pat_names, 'square_RF', 'IgnoreCase', true) | ...
                   contains(pat_names, 'px_square', 'IgnoreCase', true);
    is_bar_pat = ~is_flash_pat & ~is_bar_flash_pat;
    bar_indices = find(is_bar_pat);
    n_bars = numel(bar_indices);

    if n_bars > 0
        n_cols_fig = min(4, n_bars);
        n_rows_fig = ceil(n_bars / n_cols_fig);

        figure('Name', 'Fig 1: Bar Pattern Orientation & Direction');
        for bp = 1:n_bars
            subplot(n_rows_fig, n_cols_fig, bp);

            patt_data = load(fullfile(pattern_folder, pat_mat_files(bar_indices(bp)).name), 'pattern');
            Pats = stretch_pattern_if_needed(patt_data.pattern.Pats);

            n_rows_pat = size(Pats, 1);
            bkg = mode(double(Pats(:, :, 1)), 'all');

            % Show middle frame, flipped for screen coordinates
            mid_frame = round(size(Pats, 3) * 0.25);
            frame = double(Pats(:, :, mid_frame));
            frame_flipped = flipud(frame);

            imagesc(frame_flipped);
            colormap(gca, green_cmap);
            clim([0 15]);
            hold on;
            axis equal tight;

            % Find this pattern's index in stim_meta
            pat_name = pat_mat_files(bar_indices(bp)).name;
            pat_match = find(stim_meta.bar_sweep.pattern_file == string(pat_name) & ...
                            ~stim_meta.bar_sweep.is_flip, 1);

            if ~isempty(pat_match)
                dir_deg = stim_meta.bar_sweep.direction_deg(pat_match);
                orient_deg = stim_meta.bar_sweep.orientation_deg(pat_match);

                % Draw direction arrow (in image coordinates: y is flipped)
                cx = size(Pats, 2) / 2;
                cy = n_rows_pat / 2;
                arrow_len = 15;
                dx = arrow_len * cosd(dir_deg);
                dy = -arrow_len * sind(dir_deg); % Negative because image y is down

                quiver(cx, cy, dx, dy, 0, 'w', 'LineWidth', 2, 'MaxHeadSize', 2);

                % Draw orientation line
                orient_dx = 20 * cosd(orient_deg);
                orient_dy = -20 * sind(orient_deg);
                plot([cx - orient_dx, cx + orient_dx], [cy - orient_dy, cy + orient_dy], ...
                    'y--', 'LineWidth', 1.5);

                title(sprintf('Dir=%.0f° Ori=%.0f°\n%s', dir_deg, orient_deg, ...
                    strrep(pat_name(1:min(30,end)), '_', '\_')), 'FontSize', 8);
            else
                title(strrep(pat_name(1:min(30,end)), '_', '\_'), 'FontSize', 8);
            end

            set(gca, 'XTick', [], 'YTick', []);
        end
        set(gcf, 'Position', [50 400 1200 600]);
    end

    %% FIGURE 2: Compass diagram of bar sweep directions
    figure('Name', 'Fig 2: Bar Sweep Direction Compass');

    % Slow bars
    n_total = numel(stim_meta.bar_sweep.direction_deg);
    n_per_speed = min(16, n_total);

    polaraxes; hold on;

    colors_speed = {[0.2 0.4 0.7], [0.4 0.8 1]};
    speed_labels = {'28dps', '56dps'};

    for row = 1:n_total
        dir_rad = stim_meta.bar_sweep.direction_rad(row);
        spd = stim_meta.bar_sweep.speed_dps(row);

        if spd <= 30
            col = colors_speed{1};
            speed_idx = 1;
        else
            col = colors_speed{2};
            speed_idx = 2;
        end

        if stim_meta.bar_sweep.is_flip(row)
            line_style = '--';
        else
            line_style = '-';
        end

        polarplot([0 dir_rad], [0 1], line_style, 'Color', col, 'LineWidth', 1.5);

        % Label
        text_offset = 1.1;
        [tx, ty] = pol2cart(dir_rad, text_offset);
        text(dir_rad, text_offset, sprintf('%d', row), 'FontSize', 7, ...
            'HorizontalAlignment', 'center', 'Color', col);
    end

    title('Bar Sweep Directions (row numbers)');
    legend(speed_labels, 'Location', 'best');

    %% FIGURE 3: Bar flash orientation verification
    bar_flash_pats = find(is_bar_flash_pat);
    if ~isempty(bar_flash_pats)
        bf_patt_data = load(fullfile(pattern_folder, pat_mat_files(bar_flash_pats(1)).name), 'pattern');
        bf_Pats = stretch_pattern_if_needed(bf_patt_data.pattern.Pats);

        n_orient = numel(stim_meta.bar_flash.orientation_deg);

        figure('Name', 'Fig 3: Bar Flash Orientations');
        for o = 1:n_orient
            subplot(1, n_orient, o);

            % Show the central position frame for this orientation
            % Central position = position 6 (middle of 11)
            frame_idx = 2 + (o - 1) * 11 + 5; % Frame 1=grey, 11 per orient, position 6
            if frame_idx <= size(bf_Pats, 3)
                frame_flipped = flipud(double(bf_Pats(:, :, frame_idx)));
                imagesc(frame_flipped);
                colormap(gca, green_cmap);
                clim([0 15]);
            end
            axis equal tight;
            title(sprintf('Orient %d\n%.1f°', o, stim_meta.bar_flash.orientation_deg(o)), 'FontSize', 9);
            set(gca, 'XTick', [], 'YTick', []);
        end
        set(gcf, 'Position', [50 100 1600 200]);
    end

    %% FIGURE 4: Cross-reference — flash vs sweep pattern comparison
    if isfield(stim_meta, 'cross_reference') && ~isempty(bar_flash_pats)
        n_orient = numel(stim_meta.bar_flash.orientation_deg);

        figure('Name', 'Fig 4: Cross-Reference Flash vs Sweep');
        for o = 1:n_orient
            % Flash frame
            subplot(2, n_orient, o);
            frame_idx = 2 + (o - 1) * 11 + 5;
            if frame_idx <= size(bf_Pats, 3)
                imagesc(flipud(double(bf_Pats(:, :, frame_idx))));
                colormap(gca, green_cmap);
                clim([0 15]);
            end
            axis equal tight;
            title(sprintf('Flash %d (%.0f°)', o, stim_meta.bar_flash.orientation_deg(o)), 'FontSize', 8);
            set(gca, 'XTick', [], 'YTick', []);
            if o == 1; ylabel('Bar Flash'); end

            % Matching sweep frame
            subplot(2, n_orient, n_orient + o);
            matching_rows = stim_meta.cross_reference.flash_orient_to_sweep_rows{o};
            if ~isempty(matching_rows)
                % Use the first matching row's pattern
                match_row = matching_rows(1);
                sweep_pat_name = stim_meta.bar_sweep.pattern_file(match_row);
                sweep_patt = load(fullfile(pattern_folder, sweep_pat_name), 'pattern');
                sweep_Pats = stretch_pattern_if_needed(sweep_patt.pattern.Pats);
                mid = round(size(sweep_Pats, 3) * 0.25);
                imagesc(flipud(double(sweep_Pats(:, :, mid))));
                colormap(gca, green_cmap);
                clim([0 15]);
                title(sprintf('Sweep row %d (%.0f°)', match_row, ...
                    stim_meta.bar_sweep.orientation_deg(match_row)), 'FontSize', 8);
            end
            axis equal tight;
            set(gca, 'XTick', [], 'YTick', []);
            if o == 1; ylabel('Bar Sweep'); end
        end
        set(gcf, 'Position', [50 300 1600 400]);
    end

    %% FIGURE 5: Bar flash epoch detection on frame data
    % Load raw data (guard against cd change in load_protocol2_data)
    orig_dir = cd;
    restoreDir = onCleanup(@() cd(orig_dir));
    [~, ~, Log, ~, ~] = load_protocol2_data(exp_folder);
    f_data = Log.ADC.Volts(1, :);
    v_data = Log.ADC.Volts(2, :) * 10;

    % Get bar flash pattern for passing to parse
    bar_flash_pattern = [];
    if ~isempty(bar_flash_pats)
        bar_flash_pattern = bf_patt_data.pattern;
    end

    [~, ~, ~, ~, debug_info] = parse_bar_flash_data(f_data, v_data, bar_flash_pattern);

    figure('Name', 'Fig 5: Bar Flash Epoch Detection');
    plot(f_data, 'Color', [0.5 0.5 0.5]);
    hold on;

    epoch_colors = {[1 0.3 0.3], [0.3 0.3 1], [1 0.6 0.3], [0.3 0.8 0.3], [0.8 0.3 0.8], [0.3 0.8 0.8]};
    bf_epochs = debug_info.bar_flash_epochs;
    for e = 1:numel(bf_epochs)
        col = epoch_colors{mod(e-1, numel(epoch_colors)) + 1};
        ep_range = bf_epochs(e).start:bf_epochs(e).stop;

        % Shade the epoch region
        patch([ep_range(1) ep_range(end) ep_range(end) ep_range(1)], ...
              [0 0 max(f_data)*0.5 max(f_data)*0.5], col, ...
              'FaceAlpha', 0.2, 'EdgeColor', 'none');

        speed_str = 'unknown';
        if e <= numel(debug_info.epoch_speeds)
            if debug_info.epoch_speeds(e) == 1
                speed_str = 'slow';
            elseif debug_info.epoch_speeds(e) == 2
                speed_str = 'fast';
            end
        end

        text(double(bf_epochs(e).start), double(max(f_data)*0.4), ...
            sprintf('E%d (%s)\n%d trans', e, speed_str, bf_epochs(e).n_transitions), ...
            'FontSize', 7, 'Color', col);
    end

    xlabel('Sample');
    ylabel('Frame #');
    title(sprintf('Bar Flash Epoch Detection (%d epochs found, %d expected)', ...
        debug_info.n_bar_flash_epochs, 6));
    set(gcf, 'Position', [50 50 1800 300]);

    fprintf('--- Epoch Detection ---\n');
    fprintf('Total 3s+ gaps: %d\n', debug_info.n_gaps);
    fprintf('Bar flash epochs detected: %d (expected 6)\n', debug_info.n_bar_flash_epochs);
    for e = 1:numel(bf_epochs)
        fprintf('  Epoch %d: samples %d-%d, max_frame=%d, transitions=%d, speed=%d\n', ...
            e, bf_epochs(e).start, bf_epochs(e).stop, bf_epochs(e).max_frame, ...
            bf_epochs(e).n_transitions, debug_info.epoch_speeds(e));
    end
    fprintf('\n');

    %% FIGURE 6: Frame number distribution
    % Re-parse one rep to get frame number distribution
    figure('Name', 'Fig 6: Bar Flash Frame Number Distribution');

    if ~isempty(bf_epochs)
        ep = bf_epochs(1);
        ep_frames = f_data(ep.start:ep.stop);
        above_zero = ep_frames > 0;
        flash_ups = find(diff(above_zero) == 1);
        flash_downs = find(diff(above_zero) == -1);

        n_flashes = min(numel(flash_ups), numel(flash_downs));
        frame_nums = zeros(1, n_flashes);

        for f = 1:n_flashes
            st = flash_ups(f) + ep.start - 1;
            nd = flash_downs(f) + ep.start - 1;
            frame_nums(f) = max(f_data(st:nd));
        end

        histogram(frame_nums, 'BinMethod', 'integers');
        xlabel('Frame Number');
        ylabel('Count');
        title(sprintf('Frame Numbers in Epoch 1 (%d flashes, range %d-%d)', ...
            n_flashes, min(frame_nums), max(frame_nums)));

        fprintf('--- Frame Number Check (Epoch 1) ---\n');
        fprintf('  Flashes detected: %d (expected 88)\n', n_flashes);
        fprintf('  Frame range: %d to %d\n', min(frame_nums), max(frame_nums));
        unique_frames = unique(frame_nums);
        fprintf('  Unique frames: %d (expected 88)\n', numel(unique_frames));
        missing = setdiff(2:89, frame_nums);
        if ~isempty(missing)
            fprintf('  Missing frame numbers: [%s]\n', strjoin(string(missing), ', '));
        else
            fprintf('  No missing frame numbers.\n');
        end
        fprintf('\n');
    end

    %% FIGURE 7: Bar flash timeseries sanity check
    [data_slow, ~, mean_slow, ~] = parse_bar_flash_data(f_data, v_data, bar_flash_pattern);

    figure('Name', 'Fig 7: Bar Flash Timeseries (slow, raw indexing)');
    tiledlayout(8, 11, 'TileSpacing', 'compact', 'Padding', 'compact');
    median_v = median(v_data);

    for col = 1:8
        for row = 1:11
            nexttile;
            d = mean_slow{row, col};
            if ~isempty(d)
                plot(d, 'k', 'LineWidth', 0.8);
                hold on;
                plot([1 numel(d)], [median_v median_v], 'Color', [0.7 0.7 0.7]);
                ylim([median_v - 20, median_v + 20]);
                xlim([1 numel(d)]);
            end
            set(gca, 'XTick', [], 'YTick', []);
            if row == 1
                title(sprintf('C%d', col), 'FontSize', 7);
            end
            if col == 1 && row == 1
                ylabel('R1', 'FontSize', 7);
            end
        end
    end
    sgtitle('Bar Flash Responses (raw indexing) - rows=positions, cols=orientations');
    set(gcf, 'Position', [50 50 1600 900]);

    % Check for empty cells
    n_empty = sum(cellfun(@isempty, mean_slow), 'all');
    fprintf('--- Bar Flash Data Check ---\n');
    fprintf('  Empty cells in mean_slow: %d out of %d\n', n_empty, numel(mean_slow));
    if n_empty > 1
        fprintf('  WARNING: More than 1 empty cell (only (1,1) should be empty due to off-by-one).\n');
    end

    fprintf('\n=== Verification Complete ===\n');
    fprintf('Check the 7 figures for visual confirmation.\n');

end


%% ========== LOCAL FUNCTION ==========

function Pats = stretch_pattern_if_needed(Pats)
% Duplicated from get_stimulus_metadata for standalone use.
    [d1, d2, d3] = size(Pats);
    if d1 == 192 && d2 == 3
        stretched = zeros(48, 192, d3);
        for f = 1:d3
            frame = Pats(:, :, f);
            frame_t = frame';
            for r = 1:3
                row_start = (r - 1) * 16 + 1;
                row_end = r * 16;
                stretched(row_start:row_end, :, f) = repmat(frame_t(r, :), 16, 1);
            end
        end
        Pats = stretched;
    elseif d2 == 192 && d1 == 3
        stretched = zeros(48, 192, d3);
        for f = 1:d3
            frame = Pats(:, :, f);
            for r = 1:3
                row_start = (r - 1) * 16 + 1;
                row_end = r * 16;
                stretched(row_start:row_end, :, f) = repmat(frame(r, :), 16, 1);
            end
        end
        Pats = stretched;
    end
end
