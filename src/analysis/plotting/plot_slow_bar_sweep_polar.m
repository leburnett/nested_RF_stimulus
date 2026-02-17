function fig = plot_slow_bar_sweep_polar(bar_data, max_v, lut_directions, ...
    plot_order, median_v, title_str)
% PLOT_SLOW_BAR_SWEEP_POLAR  Radial timeseries figure with central polar plot.
%
%   FIG = PLOT_SLOW_BAR_SWEEP_POLAR(BAR_DATA, MAX_V, LUT_DIRECTIONS, ...
%       PLOT_ORDER, MEDIAN_V, TITLE_STR)
%   creates a publication-quality figure with 16 timeseries subplots
%   arranged in a circle at their LUT-verified directions, and a central
%   polar plot showing peak response vs direction with a vector sum arrow.
%
%   INPUTS:
%     bar_data       - Nx(R+1) cell array from parse_bar_data. Rows 1-16
%                      are slow bar directions. Last column is the mean.
%     max_v          - 16x1 depolarization amplitudes ordered by subplot
%                      position (from compute_bar_sweep_responses)
%     lut_directions - 16x1 LUT motion directions in degrees, indexed by
%                      data row
%     plot_order     - 1x16 subplot-to-data-row mapping
%     median_v       - Median voltage for baseline reference line
%     title_str      - Figure title string
%
%   OUTPUT:
%     fig - Figure handle
%
%   FIGURE LAYOUT:
%     16 small axes positioned radially at each LUT direction. Each subplot
%     shows individual repetitions (gray) and the mean trace (blue). A
%     horizontal gray line marks the median voltage. The central polaraxes
%     displays the directional tuning curve with a vector sum arrow
%     indicating the preferred direction.
%
%   See also COMPUTE_BAR_SWEEP_RESPONSES, VECTOR_SUM_POLAR,
%            ADD_ARROW_TO_POLARPLOT, PLOT_TIMESERIES_POLAR_BARS
% ________________________________________________________________________

    numPlots = numel(plot_order);
    col = [0.2 0.4 0.7];

    % Layout parameters
    centerX = 0.5;
    centerY = 0.5;
    radius  = 0.35;
    subW    = 0.15;
    subH    = 0.15;

    % LUT directions for each subplot position (degrees and radians)
    lut_dirs_ordered     = lut_directions(plot_order);
    lut_dirs_ordered_rad = deg2rad(lut_dirs_ordered);

    fig = figure('Name', 'Slow Bar Sweep Polar Timeseries (28 dps)');

    %% Radial timeseries subplots
    for subplot_idx = 1:numPlots
        data_row  = plot_order(subplot_idx);
        angle_rad = lut_dirs_ordered_rad(subplot_idx);

        % Position subplot at the LUT direction
        x_pos = centerX + radius * cos(angle_rad);
        y_pos = centerY + radius * sin(angle_rad);
        ax = axes('Position', [x_pos - subW/2, y_pos - subH/2, subW, subH]);
        hold on;

        n_reps = size(bar_data, 2) - 1;

        for r = 1:n_reps + 1
            d2plot = bar_data{data_row, r};
            x_vals = 1:numel(d2plot);

            if r == 1
                % Median voltage reference line
                plot([1 x_vals(end)], [median_v, median_v], ...
                    'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            end

            if r < n_reps + 1
                plot(x_vals, d2plot, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
            else
                plot(x_vals, d2plot, 'Color', col, 'LineWidth', 1.2);
            end
        end

        ylim([-80 -10]);
        axis(ax, 'off');
    end

    %% Central polar plot
    centralSize = (2 * radius) * 0.65;
    centralPosition = [centerX - centralSize/2, centerY - centralSize/2, ...
                       centralSize, centralSize];

    % Sort by direction for a smooth curve
    dir_and_val = sortrows([lut_dirs_ordered_rad, max_v], 1);
    polar_theta = [dir_and_val(:,1); dir_and_val(1,1)];  % close the loop
    polar_rho   = [dir_and_val(:,2); dir_and_val(1,2)];

    axCentral = polaraxes('Position', centralPosition);
    hold on;
    polarplot(polar_theta, polar_rho, 'Color', col, 'LineWidth', 2, ...
        'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w');

    % Vector sum arrow
    rho_for_vecsum   = max_v';
    theta_for_vecsum = lut_dirs_ordered_rad';
    [~, resultant_angle] = vector_sum_polar(rho_for_vecsum, theta_for_vecsum);

    arrow_magnitude = max(polar_rho) * 0.9;
    add_arrow_to_polarplot(arrow_magnitude, resultant_angle, [0.3 0.3 0.3]);

    axCentral.LineWidth = 1.2;
    axCentral.ThetaTick = 0:22.5:337.5;
    axCentral.ThetaTickLabel = {};

    set(fig, 'Position', [303 78 961 969]);
    sgtitle(title_str);

end
