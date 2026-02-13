function fig = plot_bar_flash_heatmap(data_slow, mean_slow, median_v, ...
    bar_flash_col, pos_order, orient_labels, title_str)
% PLOT_BAR_FLASH_HEATMAP  Full 8x11 bar flash response grid with heatmap background.
%
%   FIG = PLOT_BAR_FLASH_HEATMAP(DATA_SLOW, MEAN_SLOW, MEDIAN_V, ...
%       BAR_FLASH_COL, POS_ORDER, ORIENT_LABELS, TITLE_STR)
%   creates an 8x11 tiled layout showing all bar flash responses across
%   8 orientations and 11 spatial positions. Background colour intensity
%   encodes response magnitude relative to baseline. The row corresponding
%   to the preferred direction orientation is highlighted in red.
%
%   INPUTS:
%     data_slow      - 11x8x3 cell array of individual rep timeseries from
%                      parse_bar_flash_data (positions x orientations x reps)
%     mean_slow      - 11x8 cell array of mean timeseries
%     median_v       - Median voltage for heatmap colour normalization
%     bar_flash_col  - Column index (1-8) for the PD orientation, which is
%                      highlighted with red traces and red border
%     pos_order      - 1x11 spatial ordering array for the PD row (ND to PD).
%                      Other rows use default 1:11 ordering.
%     orient_labels  - 8x1 cell array of orientation label strings for the
%                      y-axis (e.g., {'Orient: 90°', 'Orient: 68°', ...})
%     title_str      - Figure title string
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     8 rows (orientations) x 11 columns (positions). Each tile shows
%     individual reps (gray lines) overlaid with the mean (black, or red
%     for the PD row). Background colour ranges from white (no response) to
%     red (maximum depolarization). The PD row has a red border and bold
%     red y-axis label marked with *PD*.
%
%   See also PARSE_BAR_FLASH_DATA, PLOT_BAR_FLASH_1X11, ANALYZE_SINGLE_EXPERIMENT
% ________________________________________________________________________

    n_orient = 8;
    n_pos    = 11;

    %% Compute normalised heatmap background values
    % Find overall max for consistent y-limits
    max_overall = -Inf;
    for i = 1:numel(mean_slow)
        d = mean_slow{i};
        if ~isempty(d)
            n_pts = numel(d);
            max_d = prctile(d(ceil(n_pts * 0.5):ceil(n_pts * 0.75)), 98);
            if max_d > max_overall
                max_overall = max_d;
            end
        end
    end

    % Per-tile peak values for heatmap colour
    max_vals = zeros(n_pos, n_orient);
    for orient_idx = 1:n_orient
        for pos_idx = 1:n_pos
            d = mean_slow{pos_idx, orient_idx};
            if ~isempty(d)
                n_pts = numel(d);
                max_vals(pos_idx, orient_idx) = ...
                    prctile(d(ceil(n_pts * 0.5):ceil(n_pts * 0.75)), 98);
            end
        end
    end

    % Normalise relative to median voltage
    max_vals_med = max_vals - median_v;
    max_vals_med(max_vals_med < 0) = 0;
    if max(max_vals_med(:)) > 0
        normalizedArray = 1 - max_vals_med / max(max_vals_med(:));
    else
        normalizedArray = ones(size(max_vals_med));
    end

    %% Create figure
    fig = figure('Name', 'Bar Flashes — All Orientations (80ms)');
    tiledlayout(n_orient, n_pos, 'TileSpacing', 'compact', 'Padding', 'compact');

    for orient_idx = 1:n_orient
        for pos_idx_raw = 1:n_pos
            % For the PD row, use spatial ordering (ND to PD)
            if orient_idx == bar_flash_col
                pos_idx = pos_order(pos_idx_raw);
            else
                pos_idx = pos_idx_raw;
            end

            nexttile;
            hold on;

            % Heatmap background rectangle
            d = mean_slow{pos_idx, orient_idx};
            if ~isempty(d)
                bg_val = normalizedArray(pos_idx, orient_idx);
                rectangle('Position', [0, -70, numel(d), 40], ...
                    'FaceColor', [1, bg_val, bg_val] * 0.9);
            end

            % Individual reps (gray)
            for r = 1:3
                ts = data_slow{pos_idx, orient_idx, r};
                if ~isempty(ts)
                    plot(1:numel(ts), ts, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7);
                end
            end

            % Mean trace (red for PD row, black otherwise)
            if ~isempty(d)
                if orient_idx == bar_flash_col
                    plot(1:numel(d), d, 'Color', [0.8 0 0], 'LineWidth', 1.2);
                else
                    plot(1:numel(d), d, 'k', 'LineWidth', 1);
                end
            end

            ylim([-70 max_overall * 0.9]);
            if ~isempty(d)
                xlim([0 numel(d)]);
            end
            set(gca, 'XTick', [], 'YTick', []);

            % Row labels (leftmost column only)
            if pos_idx_raw == 1
                if orient_idx == bar_flash_col
                    ylabel(sprintf('%s *PD*', orient_labels{orient_idx}), ...
                        'FontSize', 7, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                else
                    ylabel(orient_labels{orient_idx}, 'FontSize', 7);
                end
            end

            % Column labels (top row only)
            if orient_idx == 1
                title(sprintf('Pos %d', pos_idx_raw), 'FontSize', 6);
            end

            % Highlight PD row border
            if orient_idx == bar_flash_col
                set(gca, 'XColor', [0.8 0 0], 'YColor', [0.8 0 0], 'LineWidth', 1.5);
                box on;
            end
        end
    end

    sgtitle(title_str);
    set(fig, 'Position', [45 128 1675 902]);

end
