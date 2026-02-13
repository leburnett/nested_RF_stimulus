function fig = plot_bar_flash_1x11(flash_data, flash_mean, pos_order, ...
    title_str, opts)
% PLOT_BAR_FLASH_1X11  Single-row bar flash plot along one orientation axis.
%
%   FIG = PLOT_BAR_FLASH_1X11(FLASH_DATA, FLASH_MEAN, POS_ORDER, TITLE_STR)
%   creates a 1x11 tiled layout showing baseline-subtracted bar flash
%   responses across 11 spatial positions for a single bar orientation.
%
%   FIG = PLOT_BAR_FLASH_1X11(..., OPTS) uses the options structure to
%   override default plotting parameters.
%
%   This function is used for both the PD-ND axis plot and the orthogonal
%   orientation plot — the only difference is the data and title passed in.
%
%   INPUTS:
%     flash_data - 11x1x3 cell array of individual rep timeseries for one
%                  orientation (i.e., one column of the data_slow output
%                  from parse_bar_flash_data)
%     flash_mean - 11x1 cell array of mean timeseries for the same
%                  orientation (one column of mean_slow)
%     pos_order  - 1x11 array specifying the spatial ordering of positions
%                  (ND side to PD side, left to right)
%     title_str  - Figure title string
%     opts       - (Optional) structure with fields:
%                    .baseline_samples - Sample indices for computing the
%                                        pre-flash baseline mean that is
%                                        subtracted from each trace
%                                        (default: 1:5000, i.e. 500ms at 10kHz)
%                    .y_limits         - [ymin ymax] for all subplots in mV
%                                        (default: [-15 25])
%                    .fig_position     - [x y width height] figure position
%                                        in pixels (default: [50 400 1800 300])
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     1x11 tiled layout. Each tile shows 3 individual reps (light gray)
%     and the mean (thick black), all baseline-subtracted. The first tile
%     is labelled 'ND', the sixth 'Center', and the eleventh 'PD'. The
%     y-axis shows voltage change relative to baseline (DeltamV).
%
%   BASELINE SUBTRACTION:
%     For each trace, the mean voltage over the baseline_samples window
%     is computed and subtracted. This normalises responses so that the
%     pre-flash period sits at 0 mV, making depolarisation magnitudes
%     directly comparable across positions and experiments.
%
%   See also PLOT_BAR_FLASH_HEATMAP, PARSE_BAR_FLASH_DATA,
%            ANALYZE_SINGLE_EXPERIMENT
% ________________________________________________________________________

    % Set defaults
    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'baseline_samples'), opts.baseline_samples = 1:5000; end
    if ~isfield(opts, 'y_limits'),         opts.y_limits         = [-15 25]; end
    if ~isfield(opts, 'fig_position'),     opts.fig_position     = [50 400 1800 300]; end

    n_pos = 11;
    bl_samples = opts.baseline_samples;

    fig = figure('Name', title_str);
    tiledlayout(1, n_pos, 'TileSpacing', 'compact', 'Padding', 'compact');

    for pos_idx = 1:n_pos
        flash_pos = pos_order(pos_idx);
        nexttile;
        hold on;

        % Individual reps (light gray, baseline-subtracted)
        for r = 1:3
            ts = flash_data{flash_pos, 1, r};
            if ~isempty(ts)
                bl = mean(ts(bl_samples));
                plot(1:numel(ts), ts - bl, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            end
        end

        % Mean trace (thick black, baseline-subtracted)
        ts_mean = flash_mean{flash_pos};
        if ~isempty(ts_mean)
            bl_mean = mean(ts_mean(bl_samples));
            plot(1:numel(ts_mean), ts_mean - bl_mean, 'k', 'LineWidth', 1.5);
        end

        ylim(opts.y_limits);
        set(gca, 'XTick', []);

        % Y-axis label on first tile only
        if pos_idx == 1
            ylabel('\DeltamV');
        else
            set(gca, 'YTickLabel', []);
        end

        % Position labels
        if pos_idx == 1
            title('ND');
        elseif pos_idx == 6
            title('Center');
        elseif pos_idx == 11
            title('PD');
        else
            title(sprintf('%d', pos_idx));
        end

        box on;
    end

    sgtitle(title_str);
    set(fig, 'Position', opts.fig_position);

end
