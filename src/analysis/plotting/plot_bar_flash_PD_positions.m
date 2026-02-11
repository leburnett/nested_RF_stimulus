function f = plot_bar_flash_PD_positions(data_slow, mean_slow, bf_col, positions_order, median_v, params, save_folder, pd_angle_rad)
    % Plot voltage timeseries for bar flashes at 11 positions along the
    % preferred direction axis, ordered from leading to trailing.
    %
    % Inputs
    % ------
    %   data_slow : cell [11 x 8 x 3]
    %       Bar flash timeseries per position (rows), orientation (cols),
    %       and repetition (layers). From parse_bar_flash_data.
    %
    %   mean_slow : cell [11 x 8]
    %       Mean timeseries across repetitions.
    %
    %   bf_col : int
    %       Bar flash orientation column (1-8) aligned with PD axis.
    %
    %   positions_order : array [1 x 11]
    %       Position indices ordered from leading to trailing.
    %
    %   median_v : float
    %       Median voltage across the recording.
    %
    %   params : struct
    %       Contains date, time, strain, on_off fields.
    %
    %   save_folder : str
    %       Directory to save the figure.
    %
    %   pd_angle_rad : float
    %       Preferred direction angle in radians.

    n_positions = 11;
    n_reps = 3;

    bold_color = [0.2, 0.4, 0.7];
    gray_color = [0.75, 0.75, 0.75];

    % Determine consistent y-limits across all 11 positions
    all_vals = [];
    for p = 1:n_positions
        d = mean_slow{p, bf_col};
        if ~isempty(d)
            all_vals = [all_vals; d(:)];
        end
    end

    if isempty(all_vals)
        warning('No bar flash data found for column %d. Skipping figure.', bf_col);
        f = [];
        return;
    end

    y_range = [min(all_vals), max(all_vals)];
    y_pad = (y_range(2) - y_range(1)) * 0.05;
    y_lims = [y_range(1) - y_pad, y_range(2) + y_pad];

    % Create figure
    figure;
    tiledlayout(1, n_positions, 'TileSpacing', 'compact', 'Padding', 'compact');

    for c = 1:n_positions
        pos = positions_order(c); % original position index
        nexttile;
        hold on;

        d_mean = mean_slow{pos, bf_col};

        % Plot median voltage reference line
        if ~isempty(d_mean)
            plot([1, numel(d_mean)], [median_v, median_v], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
        end

        % Plot individual reps in grey
        for r = 1:n_reps
            d_rep = data_slow{pos, bf_col, r};
            if ~isempty(d_rep)
                plot(1:numel(d_rep), d_rep, 'Color', gray_color, 'LineWidth', 0.7);
            end
        end

        % Plot mean in bold colour
        if ~isempty(d_mean)
            plot(1:numel(d_mean), d_mean, 'Color', bold_color, 'LineWidth', 1.5);
        end

        ylim(y_lims);
        if ~isempty(d_mean)
            xlim([1, numel(d_mean)]);
        end

        title(sprintf('Pos %d', pos), 'FontSize', 9);

        if c == 1
            ylabel('Voltage (mV)');
        else
            set(gca, 'YTickLabel', []);
        end
        set(gca, 'XTick', []);
        box off;
    end

    sgtitle(sprintf('Bar Flash - PD axis (%.0f deg) - %s - %s - %s - %s', ...
        rad2deg(pd_angle_rad), ...
        strrep(params.date, '_', '-'), ...
        strrep(params.time, '_', '-'), ...
        strrep(params.strain, '_', '-'), ...
        params.on_off));

    f = gcf;
    f.Position = [50, 300, 1700, 300];

    % Save
    fname = fullfile(save_folder, sprintf('%s_%s_bar_flash_PD_positions.pdf', ...
        params.date, params.time));
    exportgraphics(f, fname, 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('  Bar flash PD figure saved: %s\n', fname);

end
