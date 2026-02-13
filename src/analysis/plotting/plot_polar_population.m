function fig = plot_polar_population(aligned_data_ctrl, aligned_data_ttl, title_str)
% PLOT_POLAR_POPULATION  PD-aligned polar tuning curve with mean +/- SEM.
%
%   FIG = PLOT_POLAR_POPULATION(ALIGNED_DATA_CTRL, ALIGNED_DATA_TTL, TITLE_STR)
%   creates a polar plot showing the population-averaged directional tuning
%   curve for two groups (control and ttl). Each cell's tuning data has been
%   rotated so that the preferred direction aligns to pi/2 (90 degrees,
%   pointing up). The mean is shown as a solid line and the SEM as a
%   semi-transparent shaded band.
%
%   INPUTS:
%     aligned_data_ctrl - Cell array of 16x2 matrices, one per control cell.
%                         Each matrix has columns [aligned_angle_rad, response].
%                         PD is at pi/2 after alignment by find_PD_and_order_idx.
%     aligned_data_ttl  - Cell array of 16x2 matrices, one per ttl cell.
%                         Same format as aligned_data_ctrl.
%     title_str         - Figure title string
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     Polar axes showing directional tuning. Control group in black with
%     gray SEM shading. TTL group in red with light red SEM shading.
%     Legend shows group names and sample sizes.
%
%   APPROACH:
%     Uses a dual-axis technique: polar axes for the mean lines and grid,
%     with a transparent Cartesian axes overlay for the SEM patches
%     (MATLAB's polaraxes do not natively support patch objects).
%
%   See also FIND_PD_AND_ORDER_IDX, POLAR_MEAN_BY_STRAIN,
%            BATCH_ANALYZE_1DRF
% ________________________________________________________________________

    % Compute group statistics
    stats_ctrl = compute_group_stats(aligned_data_ctrl);
    stats_ttl  = compute_group_stats(aligned_data_ttl);

    % Get canonical theta from whichever group has data
    if ~isempty(stats_ctrl.theta)
        theta = stats_ctrl.theta;
    elseif ~isempty(stats_ttl.theta)
        theta = stats_ttl.theta;
    else
        error('No valid aligned data in either group.');
    end

    % Colours
    col_ctrl_line = [0 0 0];            % black
    col_ctrl_fill = [0.80 0.80 0.80];   % light gray
    col_ttl_line  = [1 0 0];            % red
    col_ttl_fill  = [1 0.70 0.70];      % light red

    % Create figure with polar axes
    fig = figure('Name', title_str);
    axPolar = polaraxes;
    hold(axPolar, 'on');

    % Set r-limits based on data
    rmax = estimate_rmax(stats_ctrl, stats_ttl);
    if ~isnan(rmax) && rmax > 0
        rlim(axPolar, [0 rmax * 1.1]);
    end

    % Create overlay Cartesian axes for SEM patches
    axFill = axes('Position', axPolar.Position, 'Color', 'none', ...
                  'XColor', 'none', 'YColor', 'none', 'HitTest', 'off');
    axis(axFill, 'equal');
    hold(axFill, 'on');

    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);

    % Plot control group
    hCtrl = gobjects(1, 1);
    if ~isempty(stats_ctrl.mean_resp)
        hCtrl = plot_with_shade(axPolar, axFill, theta, ...
            stats_ctrl.mean_resp, stats_ctrl.sem, ...
            col_ctrl_line, col_ctrl_fill, 0.30, 2);
    end

    % Plot ttl group
    hTtl = gobjects(1, 1);
    if ~isempty(stats_ttl.mean_resp)
        hTtl = plot_with_shade(axPolar, axFill, theta, ...
            stats_ttl.mean_resp, stats_ttl.sem, ...
            col_ttl_line, col_ttl_fill, 0.30, 2);
    end

    % Keep polar axes on top
    uistack(axPolar, 'top');

    % Legend
    hs = []; labs = {};
    if ~isempty(stats_ctrl.mean_resp)
        hs(end+1)   = hCtrl;
        labs{end+1}  = sprintf('control (n=%d)', stats_ctrl.n);
    end
    if ~isempty(stats_ttl.mean_resp)
        hs(end+1)   = hTtl;
        labs{end+1}  = sprintf('ttl (n=%d)', stats_ttl.n);
    end
    if ~isempty(hs)
        legend(axPolar, hs, labs, 'Location', 'best');
    end

    axPolar.FontSize = 12;
    title(axPolar, title_str);

end


%% ========================= Local Functions ============================

function stats = compute_group_stats(aligned_data)
% COMPUTE_GROUP_STATS  Stack aligned responses and compute mean/SEM.

    stats.theta     = [];
    stats.mean_resp = [];
    stats.sem       = [];
    stats.n         = 0;

    if isempty(aligned_data)
        return;
    end

    vals = [];
    for k = 1:numel(aligned_data)
        d = aligned_data{k};
        if isnumeric(d) && isequal(size(d), [16 2])
            if isempty(stats.theta)
                stats.theta = d(:, 1);
            end
            vals = [vals, d(:, 2)]; %#ok<AGROW>
        end
    end

    if isempty(vals)
        return;
    end

    stats.n         = size(vals, 2);
    stats.mean_resp = mean(vals, 2, 'omitnan');
    sd              = std(vals, 0, 2, 'omitnan');
    stats.sem       = sd ./ sqrt(stats.n);

end


function rmax = estimate_rmax(sc, st)
% ESTIMATE_RMAX  Find maximum r-value across both groups for axis scaling.

    allVals = [];
    if ~isempty(sc.mean_resp)
        allVals = [allVals; sc.mean_resp(:) + sc.sem(:)];
    end
    if ~isempty(st.mean_resp)
        allVals = [allVals; st.mean_resp(:) + st.sem(:)];
    end
    if isempty(allVals)
        rmax = NaN;
    else
        rmax = max(allVals, [], 'omitnan');
        if ~isfinite(rmax) || rmax <= 0
            rmax = NaN;
        end
    end

end


function hLine = plot_with_shade(axPolar, axFill, theta, meanVals, bandVals, ...
    lineColor, fillColor, alphaVal, lw)
% PLOT_WITH_SHADE  Draw mean line on polar axes and SEM patch on Cartesian overlay.

    th = theta(:);
    m  = meanVals(:);
    b  = bandVals(:);

    % Ensure non-negative radii
    upper = max(m + b, 0);
    lower = max(m - b, 0);

    % Build closed polygon in polar, convert to Cartesian for patch
    th_poly = [th; th(1); flipud(th); th(1)];
    r_poly  = [upper; upper(1); flipud(lower); lower(1)];
    [x_poly, y_poly] = pol2cart(th_poly, r_poly);

    patch('XData', x_poly, 'YData', y_poly, ...
          'FaceColor', fillColor, 'FaceAlpha', alphaVal, ...
          'EdgeColor', 'none', 'Parent', axFill, 'HitTest', 'off');

    % Draw the mean line on polar axes (closed loop)
    th_line = [th; th(1)];
    r_line  = [m; m(1)];
    hLine   = polarplot(axPolar, th_line, r_line, '-', ...
        'Color', lineColor, 'LineWidth', lw);

end
