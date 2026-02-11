function plot_bar_flash_PD_positions(data_slow, mean_slow, bf_col, positions_order, median_v, params, save_folder, pd_angle_rad, exp_folder)
    % Plot voltage timeseries for bar flashes at 11 positions along the
    % preferred direction axis. Generates three figures:
    %   1. Raw voltage
    %   2. Baseline-subtracted (per-rep baseline from samples 2500:5000)
    %   3. Diagnostic: shows per-rep baselines on full-resolution data
    %
    % Each main figure has 2 rows x 11 columns:
    %   Row 1 (~25% height): Stimulus stills from the bar flash pattern
    %   Row 2 (~75% height): Voltage timeseries (downsampled 100x)
    %
    % Baseline subtraction: each rep has its own baseline computed from
    % samples 2500:5000, then the baseline-subtracted mean is computed
    % from the corrected reps. Peak position (column 6) is determined
    % from the baseline-subtracted mean traces.
    %
    % Inputs
    % ------
    %   data_slow : cell [11 x 8 x 3]
    %       Bar flash timeseries per position, orientation, repetition.
    %
    %   mean_slow : cell [11 x 8]
    %       Mean timeseries across repetitions (used for raw plot only).
    %
    %   bf_col : int
    %       Bar flash orientation column (1-8) aligned with PD axis.
    %
    %   positions_order : array [1 x 11]
    %       Position indices ordered leading to trailing.
    %
    %   median_v : float
    %       Median voltage across the recording.
    %
    %   params : struct
    %       Contains date, time, strain, on_off fields.
    %
    %   save_folder : str
    %       Directory to save the figures.
    %
    %   pd_angle_rad : float
    %       Preferred direction angle in radians.
    %
    %   exp_folder : str
    %       Path to experiment data folder (for loading pattern file).

    n_positions = 11;
    n_reps = 3;
    ds_factor = 100; % Downsample factor
    baseline_range = [2500, 5000]; % Full-resolution datapoint indices for baseline

    bold_color = [0.2, 0.4, 0.7];
    gray_color = [0.75, 0.75, 0.75];

    %% Check for data
    has_data = false;
    for p = 1:n_positions
        if ~isempty(mean_slow{p, bf_col})
            has_data = true;
            break;
        end
    end
    if ~has_data
        warning('No bar flash data found for column %d. Skipping figure.', bf_col);
        return;
    end

    %% Per-rep baseline subtraction
    % For each position and rep, compute baseline from that rep's samples
    % 2500:5000, subtract it, then compute the baseline-subtracted mean
    % by averaging the corrected reps.
    bl_sub_means = cell(1, n_positions);    % baseline-subtracted mean traces
    bl_sub_reps = cell(n_positions, n_reps); % baseline-subtracted individual reps
    rep_baselines = NaN(n_positions, n_reps); % per-rep baselines for diagnostics

    for p = 1:n_positions
        sub_reps = {};
        for r = 1:n_reps
            d_rep = data_slow{p, bf_col, r};
            if ~isempty(d_rep)
                d_rep = d_rep(:);
                bl_end = min(baseline_range(2), numel(d_rep));
                if baseline_range(1) <= numel(d_rep)
                    bl = mean(d_rep(baseline_range(1):bl_end), 'omitnan');
                else
                    bl = 0;
                end
                rep_baselines(p, r) = bl;
                bl_sub_reps{p, r} = d_rep - bl;
                sub_reps{end+1} = d_rep - bl; %#ok<AGROW>
            end
        end
        % Compute mean of baseline-subtracted reps (NaN-pad to equal length)
        if ~isempty(sub_reps)
            max_len = max(cellfun(@numel, sub_reps));
            padded = NaN(max_len, numel(sub_reps));
            for k = 1:numel(sub_reps)
                padded(1:numel(sub_reps{k}), k) = sub_reps{k};
            end
            bl_sub_means{p} = nanmean(padded, 2);
        end
    end

    %% Determine peak position from baseline-subtracted mean traces
    peak_vals = -Inf(1, n_positions);
    for p = 1:n_positions
        if ~isempty(bl_sub_means{p})
            peak_vals(p) = max(bl_sub_means{p});
        end
    end
    [~, best_pos] = max(peak_vals);

    %% Reorder positions so peak response is in column 6
    idx_in_order = find(positions_order == best_pos, 1);
    rot = mod(6 - idx_in_order, n_positions);
    display_order = circshift(positions_order, rot);

    fprintf('  Bar flash peak position: %d (bl-sub max = %.2f mV), placed in column 6\n', ...
        best_pos, peak_vals(best_pos));

    %% Load bar flash pattern for stimulus stills
    Pats = [];
    try
        pattern_folder = fullfile(exp_folder, 'Patterns');
        pattern_files = dir(fullfile(pattern_folder, '*.mat'));
        [~, si] = sort({pattern_files.name});
        pattern_files = pattern_files(si);
        bf_pattern_file = fullfile(pattern_folder, pattern_files(end).name);
        loaded = load(bf_pattern_file, 'pattern');
        Pats = loaded.pattern.Pats;
    catch
        warning('Could not load bar flash pattern file. Stimulus stills will be blank.');
    end

    % Green colormap for stimulus display
    n_levels = 16;
    green_cmap = zeros(n_levels, 3);
    for k = 1:n_levels
        green_cmap(k, :) = [0, (k-1)/(n_levels-1), 0];
    end

    %% Diagnostic figure: show per-rep baselines on full-resolution data
    diag_fig = figure;
    diag_fig.Position = [50, 50, 1700, 400];
    tiledlayout(1, n_positions, 'TileSpacing', 'compact', 'Padding', 'compact');

    rep_colors = lines(n_reps); % distinct colors for each rep

    for c = 1:n_positions
        pos = display_order(c);
        nexttile; hold on;

        % Plot individual reps (downsampled) with their baselines
        has_any_rep = false;
        for r = 1:n_reps
            d_rep = data_slow{pos, bf_col, r};
            if ~isempty(d_rep)
                has_any_rep = true;
                d_rep = d_rep(:);
                n_pts = numel(d_rep);
                d_ds = d_rep(1:ds_factor:end);
                x_ds = (0:numel(d_ds)-1) * ds_factor;
                plot(x_ds, d_ds, 'Color', [rep_colors(r,:), 0.4], 'LineWidth', 0.7);

                % Show per-rep baseline as dashed line
                if ~isnan(rep_baselines(pos, r))
                    plot([0, n_pts], [rep_baselines(pos, r), rep_baselines(pos, r)], ...
                        '--', 'Color', rep_colors(r,:), 'LineWidth', 1);
                end
            end
        end

        if has_any_rep
            % Plot baseline-subtracted mean (shifted back up by mean of baselines for display)
            % Actually plot the raw mean_slow for context
            d_mean = mean_slow{pos, bf_col};
            if ~isempty(d_mean)
                d_mean = d_mean(:);
                n_pts = numel(d_mean);
                d_ds = d_mean(1:ds_factor:end);
                x_ds = (0:numel(d_ds)-1) * ds_factor;
                plot(x_ds, d_ds, 'Color', bold_color, 'LineWidth', 1.5);
            end

            % Shade the baseline region in red
            yl_data = [];
            for r = 1:n_reps
                d_rep = data_slow{pos, bf_col, r};
                if ~isempty(d_rep)
                    d_rep = d_rep(:);
                    yl_data = [yl_data; d_rep(1:ds_factor:end)]; %#ok<AGROW>
                end
            end
            yl = [min(yl_data)*1.02, max(yl_data)*0.98];
            if yl(1) >= yl(2)
                yl = [yl(1)-1, yl(2)+1];
            end
            max_n = 0;
            for r = 1:n_reps
                d_rep = data_slow{pos, bf_col, r};
                if ~isempty(d_rep)
                    max_n = max(max_n, numel(d_rep));
                end
            end
            bl_end = min(baseline_range(2), max_n);
            patch([baseline_range(1), bl_end, bl_end, baseline_range(1)], ...
                  [yl(1), yl(1), yl(2), yl(2)], ...
                  'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

            ylim(yl);

            % Title with per-rep baselines
            bl_strs = {};
            for r = 1:n_reps
                if ~isnan(rep_baselines(pos, r))
                    bl_strs{end+1} = sprintf('%.1f', rep_baselines(pos, r)); %#ok<AGROW>
                end
            end
            title(sprintf('Pos %d\nbl=[%s]', pos, strjoin(bl_strs, ', ')), 'FontSize', 7);
        else
            title(sprintf('Pos %d\nempty', pos), 'FontSize', 8);
        end
        box off;
        if c > 1
            set(gca, 'YTickLabel', []);
        else
            ylabel('Voltage (mV)');
        end
    end

    annotation('textbox', [0 0.93 1 0.07], 'String', ...
        sprintf('DIAGNOSTIC: Per-rep baselines (samples %d:%d), dashed=rep baselines, bold=raw mean - %s %s', ...
        baseline_range(1), baseline_range(2), params.date, params.time), ...
        'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
        'FontSize', 10, 'FontWeight', 'bold', 'FitBoxToText', 'off');

    diag_fname = fullfile(save_folder, sprintf('%s_%s_bar_flash_baseline_diagnostic.pdf', ...
        params.date, params.time));
    exportgraphics(diag_fig, diag_fname, 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('  Baseline diagnostic saved: %s\n', diag_fname);
    close(diag_fig);

    %% Generate two versions: raw and baseline-subtracted
    versions = {'raw', 'baseline_sub'};

    for v = 1:numel(versions)
        version = versions{v};
        is_baseline_sub = strcmp(version, 'baseline_sub');

        % Compute y-limits for this version
        all_vals = [];
        for p = 1:n_positions
            if is_baseline_sub
                d = bl_sub_means{p};
            else
                d = mean_slow{p, bf_col};
                if ~isempty(d)
                    d = d(:);
                end
            end
            if ~isempty(d)
                d_ds = d(1:ds_factor:end);
                all_vals = [all_vals; d_ds]; %#ok<AGROW>
            end
        end

        y_range = [min(all_vals), max(all_vals)];
        y_pad = (y_range(2) - y_range(1)) * 0.05;
        y_lims = [y_range(1) - y_pad, y_range(2) + y_pad];

        % Create figure with unequal row heights
        fig = figure;
        fig.Position = [50, 200, 1700, 500];
        t = tiledlayout(2, n_positions, 'TileSpacing', 'compact', 'Padding', 'loose');

        % Set row heights: stills ~25%, timeseries ~75%
        try
            t.RowHeight = {'1x', '3x'};
        catch
            % RowHeight not available in older MATLAB versions
        end

        % --- Row 1: Stimulus stills ---
        for c = 1:n_positions
            pos = display_order(c);
            nexttile(c);

            if ~isempty(Pats)
                frame_num = (bf_col - 1) * 11 + pos;
                still = Pats(:, :, frame_num + 1);
                imagesc(still);
                colormap(gca, green_cmap);
                clim([0 15]);
            end

            axis image;
            set(gca, 'XTick', [], 'YTick', []);
            title(sprintf('Pos %d', pos), 'FontSize', 9);
        end

        % --- Row 2: Voltage timeseries ---
        for c = 1:n_positions
            pos = display_order(c);
            nexttile(n_positions + c);
            hold on;

            if is_baseline_sub
                ref_line = 0;

                % Plot baseline-subtracted individual reps
                for r = 1:n_reps
                    d_rep_sub = bl_sub_reps{pos, r};
                    if ~isempty(d_rep_sub)
                        d_rep_ds = d_rep_sub(1:ds_factor:end);
                        plot(1:numel(d_rep_ds), d_rep_ds, 'Color', gray_color, 'LineWidth', 0.7);
                    end
                end

                % Plot baseline-subtracted mean
                d_mean_sub = bl_sub_means{pos};
                if isempty(d_mean_sub)
                    continue;
                end
                d_mean_ds = d_mean_sub(1:ds_factor:end);
                x_ds = 1:numel(d_mean_ds);

                plot([1, numel(d_mean_ds)], [ref_line, ref_line], ...
                    'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
                plot(x_ds, d_mean_ds, 'Color', bold_color, 'LineWidth', 1.5);

            else
                d_mean = mean_slow{pos, bf_col};
                if isempty(d_mean)
                    continue;
                end
                d_mean = d_mean(:);
                d_mean_ds = d_mean(1:ds_factor:end);
                x_ds = 1:numel(d_mean_ds);

                % Individual reps (raw, no subtraction)
                for r = 1:n_reps
                    d_rep = data_slow{pos, bf_col, r};
                    if ~isempty(d_rep)
                        d_rep = d_rep(:);
                        d_rep_ds = d_rep(1:ds_factor:end);
                        plot(1:numel(d_rep_ds), d_rep_ds, 'Color', gray_color, 'LineWidth', 0.7);
                    end
                end

                plot(x_ds, d_mean_ds, 'Color', bold_color, 'LineWidth', 1.5);
            end

            ylim(y_lims);
            xlim([1, numel(d_mean_ds)]);

            if c == 1
                ylabel('Voltage (mV)');
            else
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'XTick', []);
            box off;
        end

        % Figure title
        if is_baseline_sub
            version_str = ' [baseline subtracted]';
        else
            version_str = '';
        end

        title_str = sprintf('Bar Flash - PD axis (%.0f deg) - %s - %s - %s - %s%s', ...
            rad2deg(pd_angle_rad), ...
            strrep(params.date, '_', '-'), ...
            strrep(params.time, '_', '-'), ...
            strrep(params.strain, '_', '-'), ...
            params.on_off, ...
            version_str);

        annotation('textbox', [0 0.95 1 0.05], 'String', title_str, ...
            'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
            'FontSize', 11, 'FontWeight', 'bold', 'FitBoxToText', 'off');

        % Save
        if is_baseline_sub
            fname = fullfile(save_folder, sprintf('%s_%s_bar_flash_PD_positions_baseline_sub.pdf', ...
                params.date, params.time));
        else
            fname = fullfile(save_folder, sprintf('%s_%s_bar_flash_PD_positions.pdf', ...
                params.date, params.time));
        end
        exportgraphics(fig, fname, 'ContentType', 'vector', 'BackgroundColor', 'none');
        fprintf('  Bar flash PD figure saved: %s\n', fname);
        close(fig);
    end

end
