function plot_bar_flash_PD_positions(data_slow, mean_slow, bf_col, positions_order, median_v, params, save_folder, pd_angle_rad, exp_folder)
    % Plot voltage timeseries for bar flashes at 11 positions along the
    % preferred direction axis. Generates two figures:
    %   1. Raw voltage
    %   2. Baseline-subtracted (mean of datapoints 5000:7500 set to 0)
    %
    % Each figure has 2 rows x 11 columns:
    %   Row 1: Stimulus stills from the bar flash pattern
    %   Row 2: Voltage timeseries (downsampled 100x)
    %
    % Positions are ordered leading-to-trailing but rotated so that the
    % peak-responding position sits in column 6.
    %
    % Inputs
    % ------
    %   data_slow : cell [11 x 8 x 3]
    %       Bar flash timeseries per position, orientation, repetition.
    %
    %   mean_slow : cell [11 x 8]
    %       Mean timeseries across repetitions.
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

    %% Reorder positions so peak response is in column 6
    % Find the position with the strongest peak response
    peak_vals = -Inf(1, n_positions);
    for p = 1:n_positions
        d = mean_slow{p, bf_col};
        if ~isempty(d)
            peak_vals(p) = max(d);
        end
    end
    [~, best_pos] = max(peak_vals);

    % Rotate positions_order so best_pos lands at column 6
    idx_in_order = find(positions_order == best_pos, 1);
    rot = mod(6 - idx_in_order, n_positions);
    display_order = circshift(positions_order, rot);

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

    % Green colormap for stimulus display (gs_val=4, values 0-15)
    n_levels = 16;
    green_cmap = zeros(n_levels, 3);
    for k = 1:n_levels
        green_cmap(k, :) = [0, (k-1)/(n_levels-1), 0];
    end

    %% Generate two versions: raw and baseline-subtracted
    versions = {'raw', 'baseline_sub'};

    for v = 1:numel(versions)
        version = versions{v};

        % Compute y-limits for this version
        all_vals = [];
        for p = 1:n_positions
            d = mean_slow{p, bf_col};
            if ~isempty(d)
                if strcmp(version, 'baseline_sub')
                    baseline = mean(d(5000:min(7500, numel(d))));
                    d = d - baseline;
                end
                d_ds = d(1:ds_factor:end);
                all_vals = [all_vals; d_ds(:)]; %#ok<AGROW>
            end
        end

        y_range = [min(all_vals), max(all_vals)];
        y_pad = (y_range(2) - y_range(1)) * 0.05;
        y_lims = [y_range(1) - y_pad, y_range(2) + y_pad];

        % Create figure
        fig = figure;
        t = tiledlayout(2, n_positions, 'TileSpacing', 'compact', 'Padding', 'compact');

        % --- Row 1: Stimulus stills ---
        for c = 1:n_positions
            pos = display_order(c);
            nexttile(c);

            if ~isempty(Pats)
                frame_num = (bf_col - 1) * 11 + pos;
                still = Pats(:, :, frame_num + 1); % +1 because frame 1 is grey background
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

            d_mean = mean_slow{pos, bf_col};
            if isempty(d_mean)
                continue;
            end

            % Compute baseline for subtraction version
            if strcmp(version, 'baseline_sub')
                baseline = mean(d_mean(5000:min(7500, numel(d_mean))));
            else
                baseline = 0;
            end

            % Reference line (median_v for raw, 0 for baseline-subtracted)
            if strcmp(version, 'baseline_sub')
                ref_line = 0;
            else
                ref_line = median_v;
            end

            % Downsample and plot
            d_mean_sub = d_mean - baseline;
            d_mean_ds = d_mean_sub(1:ds_factor:end);
            x_ds = 1:numel(d_mean_ds);

            plot([1, numel(d_mean_ds)], [ref_line, ref_line], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);

            % Plot individual reps in grey
            for r = 1:n_reps
                d_rep = data_slow{pos, bf_col, r};
                if ~isempty(d_rep)
                    if strcmp(version, 'baseline_sub')
                        rep_baseline = mean(d_rep(5000:min(7500, numel(d_rep))));
                        d_rep = d_rep - rep_baseline;
                    end
                    d_rep_ds = d_rep(1:ds_factor:end);
                    plot(1:numel(d_rep_ds), d_rep_ds, 'Color', gray_color, 'LineWidth', 0.7);
                end
            end

            % Plot mean in bold colour
            plot(x_ds, d_mean_ds, 'Color', bold_color, 'LineWidth', 1.5);

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

        % Supertitle
        if strcmp(version, 'baseline_sub')
            version_str = ' [baseline subtracted]';
        else
            version_str = '';
        end

        sgtitle(sprintf('Bar Flash - PD axis (%.0f deg) - %s - %s - %s - %s%s', ...
            rad2deg(pd_angle_rad), ...
            strrep(params.date, '_', '-'), ...
            strrep(params.time, '_', '-'), ...
            strrep(params.strain, '_', '-'), ...
            params.on_off, ...
            version_str));

        fig.Position = [50, 200, 1700, 500];

        % Save
        if strcmp(version, 'baseline_sub')
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
