function resultant_angle = plot_polar_with_arrow(max_v, median_voltage, params, save_fig, fig_save_folder)
% PLOT_POLAR_WITH_ARROW  Create polar plot with vector sum direction arrow.
%
%   RESULTANT_ANGLE = PLOT_POLAR_WITH_ARROW(MAX_V, MEDIAN_VOLTAGE, ...
%       PARAMS, SAVE_FIG, FIG_SAVE_FOLDER)
%   generates a polar plot showing directional tuning with an arrow
%   indicating the preferred direction computed via vector sum.
%
%   INPUTS:
%     max_v          - 16x2 array of max response per direction
%                      Column 1: slow bars, Column 2: fast bars
%     median_voltage - Baseline voltage for normalization
%     params         - Structure with: .date, .time, .strain, .on_off
%     save_fig       - Boolean, true to save figure
%     fig_save_folder - Directory for saving figures
%
%   OUTPUT:
%     resultant_angle - Preferred direction in radians from vector sum
%
%   FIGURE ELEMENTS:
%     - Polar line plot showing response magnitude at each direction
%     - Dark blue line: slow (28 dps) bars with filled markers
%     - Light blue line: fast (56 dps) bars with filled markers
%     - Gray arrow: vector sum direction (preferred direction)
%
%   VECTOR SUM:
%     The arrow direction is computed as the angle of the weighted
%     sum of unit vectors, where weights are the response magnitudes.
%     Arrow magnitude is fixed for visualization purposes.
%
%   OUTPUT FILES:
%     <date>_<time>_<strain>_polarplot_bar_DS_wArrow.png
%     <date>_<time>_<strain>_polarplot_bar_DS_wArrow.pdf
%
%   See also VECTOR_SUM_POLAR, ADD_ARROW_TO_POLARPLOT, PLOT_TIMESERIES_POLAR_BARS

    % Convert max values for both conditions into polar format
    max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % slow bars
    max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % fast bars

    colors = {[0.2 0.4 0.7], [0.4 0.8 1]};
 
    theta = linspace(0, 2*pi, 17);
    rho = [max_v_polar1 - median_voltage]';
    [~, resultant_angle] = vector_sum_polar(rho, theta);
    
    resultant_magnitude = 30; 
    col = [0.3 0.3 0.3];
    
    figure;
    polarplot(theta, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    hold on
    polarplot(theta, max_v_polar2 - median_voltage, 'Color', colors{2}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    ax = gca;
    ax.LineWidth = 1.2;
    ax.FontSize = 15;
    ax.ThetaTick = rad2deg(theta);
    ax.ThetaTickLabel = {};
    hold on
    add_arrow_to_polarplot(resultant_magnitude, resultant_angle, col)

    if save_fig

        fig_folder = fullfile(fig_save_folder, params.on_off);
        if ~isfolder(fig_folder)
            mkdir(fig_folder)
        end 

        % Save as PNG
        exportgraphics(gca ...
                , strcat(params.date, "_", params.time, "_", params.strain, "_polarplot_bar_DS_wArrow.png") ...
                );
        % Save as PDF
        exportgraphics(gca ...
                , strcat(params.date, "_", params.time, "_", params.strain, "_polarplot_bar_DS_wArrow.pdf") ...
                , 'ContentType', 'vector' ...
                , 'BackgroundColor', 'none' ...
                ); 
    end 


end 