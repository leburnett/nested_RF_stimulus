function fig = plot_flash_1x11_population(traces_ctrl, traces_ttl, title_str, opts)
% PLOT_FLASH_1X11_POPULATION  Population 1x11 bar flash plot with shaded spread.
%
%   FIG = PLOT_FLASH_1X11_POPULATION(TRACES_CTRL, TRACES_TTL, TITLE_STR)
%   creates a 1x11 tiled layout showing the population center bar flash
%   response at each spatial position, with a shaded spread band. Control
%   and ttl groups are overlaid on the same axes.
%
%   FIG = PLOT_FLASH_1X11_POPULATION(..., OPTS) uses the options structure
%   to override default plotting parameters.
%
%   INPUTS:
%     traces_ctrl - Cell array (one per control cell). Each element is an
%                   11xN numeric matrix of baseline-subtracted mean flash
%                   traces (11 positions x N timepoints). Positions are
%                   ordered ND (row 1) to PD (row 11).
%     traces_ttl  - Cell array (one per ttl cell), same format.
%     title_str   - Figure title string
%     opts        - (Optional) structure with fields:
%                     .y_limits     - [ymin ymax] in mV (default: [-15 35])
%                     .fig_position - [x y w h] in pixels
%                                     (default: [50 400 1800 300])
%                     .plot_type    - 'pd_nd' or 'orthogonal'
%                                     (default: 'pd_nd'). Controls whether
%                                     tiles 1/11 are labelled PD/ND.
%                     .stat_method  - 'median_mad' (default) or 'mean_sem'
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     1x11 tiled layout. Each tile shows the center trace for control
%     (black) and ttl (red) with shaded spread (gray and light red
%     respectively). The shaded band represents MAD (default) or SEM
%     depending on stat_method. For 'pd_nd' plots, tile 1 is labelled
%     'ND', tile 6 'Center', and tile 11 'PD'. For 'orthogonal' plots,
%     tile 1 is '1', tile 6 'Center', and tile 11 '11'.
%
%   See also PLOT_BAR_FLASH_1X11, BATCH_ANALYZE_1DRF, ANALYZE_SINGLE_EXPERIMENT
% ________________________________________________________________________

    % Set defaults
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'y_limits'),     opts.y_limits     = [-15 35]; end
    if ~isfield(opts, 'fig_position'), opts.fig_position = [50 400 1800 300]; end
    if ~isfield(opts, 'plot_type'),    opts.plot_type    = 'pd_nd'; end
    if ~isfield(opts, 'stat_method'),  opts.stat_method  = 'median_mad'; end

    n_pos = 11;

    % Colours
    col_ctrl_line = [0 0 0];            % black
    col_ctrl_fill = [0.80 0.80 0.80];   % light gray
    col_ttl_line  = [1 0 0];            % red
    col_ttl_fill  = [1 0.70 0.70];      % light red
    alpha_val     = 0.40;

    fig = figure('Name', title_str);
    tiledlayout(1, n_pos, 'TileSpacing', 'compact', 'Padding', 'compact');

    for pos_idx = 1:n_pos
        nexttile;
        hold on;

        % Extract traces for this position from all cells, compute stats
        [ctr_ctrl, spr_ctrl, n_ctrl] = compute_trace_stats(traces_ctrl, pos_idx, opts.stat_method);
        [ctr_ttl,  spr_ttl,  n_ttl]  = compute_trace_stats(traces_ttl, pos_idx, opts.stat_method);

        % Plot control: shaded spread then center line
        if ~isempty(ctr_ctrl)
            x = 1:numel(ctr_ctrl);
            plot_shaded_trace(x, ctr_ctrl, spr_ctrl, ...
                col_ctrl_line, col_ctrl_fill, alpha_val, 1.5);
        end

        % Plot ttl: shaded spread then center line
        if ~isempty(ctr_ttl)
            x = 1:numel(ctr_ttl);
            plot_shaded_trace(x, ctr_ttl, spr_ttl, ...
                col_ttl_line, col_ttl_fill, alpha_val, 1.5);
        end

        ylim(opts.y_limits);
        set(gca, 'XTick', []);

        % Y-axis label on first tile only
        if pos_idx == 1
            ylabel('\DeltamV');
        else
            set(gca, 'YTickLabel', []);
        end

        % Position labels — show PD/ND only for pd_nd plots
        if strcmpi(opts.plot_type, 'pd_nd')
            if pos_idx == 1
                title('ND');
            elseif pos_idx == 6
                title('Center');
            elseif pos_idx == 11
                title('PD');
            else
                title(sprintf('%d', pos_idx));
            end
        else
            % Orthogonal plot: no PD/ND labels
            if pos_idx == 6
                title('Center');
            else
                title(sprintf('%d', pos_idx));
            end
        end

        box off;
        ax = gca;
        ax.LineWidth = 1.2;
        ax.TickDir = 'out';
        ax.TickLength = [0.015 0.015];
        ax.FontSize = 14;

        % Add legend to first tile only
        if pos_idx == 1
            h = [];
            labs = {};
            if ~isempty(ctr_ctrl)
                h(end+1) = plot(NaN, NaN, '-', 'Color', col_ctrl_line, 'LineWidth', 1.5);
                labs{end+1} = sprintf('control (n=%d)', n_ctrl);
            end
            if ~isempty(ctr_ttl)
                h(end+1) = plot(NaN, NaN, '-', 'Color', col_ttl_line, 'LineWidth', 1.5);
                labs{end+1} = sprintf('ttl (n=%d)', n_ttl);
            end
            if ~isempty(h)
                legend(h, labs, 'Location', 'northwest', 'FontSize', 8);
            end
        end
    end

    sgtitle(title_str, 'FontSize', 16);
    set(fig, 'Position', opts.fig_position);

end


%% ========================= Local Functions ============================

function [center_trace, spread_trace, n] = compute_trace_stats(traces_cell, pos_idx, stat_method)
% COMPUTE_TRACE_STATS  Stack traces at one position across cells, compute center/spread.
%   stat_method: 'median_mad' (default) or 'mean_sem'.

    center_trace = [];
    spread_trace = [];
    n            = 0;

    if isempty(traces_cell)
        return;
    end

    % First pass: find minimum trace length across cells
    min_len = Inf;
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            min_len = min(min_len, size(mat, 2));
        end
    end
    if isinf(min_len)
        return;
    end

    % Second pass: collect traces truncated to common length
    all_traces = [];
    for k = 1:numel(traces_cell)
        mat = traces_cell{k};
        if ~isempty(mat) && pos_idx <= size(mat, 1)
            all_traces = [all_traces; mat(pos_idx, 1:min_len)]; %#ok<AGROW>
        end
    end

    if isempty(all_traces)
        return;
    end

    n = size(all_traces, 1);

    if strcmpi(stat_method, 'mean_sem')
        center_trace = mean(all_traces, 1, 'omitnan');
        sd           = std(all_traces, 0, 1, 'omitnan');
        spread_trace = sd ./ sqrt(n);
    else  % 'median_mad'
        center_trace = median(all_traces, 1, 'omitnan');
        spread_trace = mad(all_traces, 1, 1);  % median absolute deviation (flag=1)
    end

end


function plot_shaded_trace(x, center_vals, spread_vals, line_col, fill_col, alpha_val, lw)
% PLOT_SHADED_TRACE  Plot center line with shaded spread band on current axes.

    x = x(:)';
    m = center_vals(:)';
    s = spread_vals(:)';

    upper = m + s;
    lower = m - s;

    % Shaded spread band
    fill([x, fliplr(x)], [upper, fliplr(lower)], fill_col, ...
        'FaceAlpha', alpha_val, 'EdgeColor', 'none');

    % Center line on top
    plot(x, m, '-', 'Color', line_col, 'LineWidth', lw);

end
