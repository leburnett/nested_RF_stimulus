function resultant_angle = plot_polar_with_arrow_flex(max_v, median_voltage, params, save_fig, fig_save_folder)
% PLOT_POLAR_WITH_ARROW_FLEX  Speed-flexible polar plot with vector sum direction arrows.
%
%   RESULTANT_ANGLE = PLOT_POLAR_WITH_ARROW_FLEX(MAX_V, MEDIAN_VOLTAGE, ...
%       PARAMS, SAVE_FIG, FIG_SAVE_FOLDER)
%   generates a polar plot showing directional tuning with arrows
%   indicating the preferred direction computed via vector sum.
%   Automatically detects the number of bar sweep speeds from the
%   max_v columns and adjusts plotting accordingly.
%
%   INPUTS:
%     max_v          - 16 x n_speeds array of max response per direction
%                      Column 1: slow bars (28 dps)
%                      Column 2: fast bars (56 dps)
%                      Column 3: very fast bars (168 dps, if present)
%                      Values are raw 98th percentile voltages (negative).
%     median_voltage - Baseline voltage for normalization (subtracted for plotting)
%     params         - Structure with: .date, .time, .strain, .on_off
%     save_fig       - Boolean, true to save figure
%     fig_save_folder - Directory for saving figures
%
%   OUTPUT:
%     resultant_angle - Preferred direction in radians from vector sum
%                       (computed from slow-speed data, 28 dps)
%
%   FIGURE ELEMENTS:
%     - Polar line plot showing response magnitude at each direction
%     - One line per speed, colored by speed condition
%     - Arrows showing vector sum direction for each speed
%     - Dark blue = slow (28 dps), light blue = fast (56 dps),
%       purple = very fast (168 dps, if present)
%
%   SPEED DETECTION:
%     n_speeds = size(max_v, 2)
%     Supports 2 speeds (Summer 2025) and 3 speeds (Autumn 2025)
%
%   See also VECTOR_SUM_POLAR, ADD_ARROW_TO_POLARPLOT, PLOT_TIMESERIES_POLAR_BARS_FLEX

    % Detect number of speeds
    n_speeds = size(max_v, 2);

    % Colors for each speed
    all_colors = {[0.2 0.4 0.7], [0.4 0.8 1], [0.45, 0.0, 0.55]};
    all_arrow_colors = {[0.3 0.3 0.3], [0.6 0.6 0.6], [0.85 0.85 0.85]};
    colors = all_colors(1:n_speeds);
    arrow_colors = all_arrow_colors(1:n_speeds);

    theta = linspace(0, 2*pi, 17);
    resultant_magnitude = 30;

    figure;

    % Plot polar traces for each speed — subtract median for positive values
    for sp = 1:n_speeds
        max_v_polar = vertcat(max_v(:, sp), max_v(1, sp));
        polarplot(theta, max_v_polar - median_voltage, 'Color', colors{sp}, ...
            'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
        hold on
    end

    ax = gca;
    ax.LineWidth = 1.2;
    ax.FontSize = 15;
    ax.ThetaTick = rad2deg(theta);
    ax.ThetaTickLabel = {};
    hold on

    % Add arrows for each speed
    for sp = 1:n_speeds
        max_v_polar = vertcat(max_v(:, sp), max_v(1, sp));
        rho = [max_v_polar(1:end-1, 1) - median_voltage]';
        [~, ra] = vector_sum_polar(rho, theta(1:end-1));

        % Store the slow-speed (sp=1) resultant angle as the return value
        if sp == 1
            resultant_angle = ra;
        end

        add_arrow_to_polarplot(resultant_magnitude, ra, arrow_colors{sp});
    end

    if save_fig
        fig_folder = fullfile(fig_save_folder, params.on_off);
        if ~isfolder(fig_folder)
            mkdir(fig_folder)
        end

        % Save as PNG
        exportgraphics(gca, fullfile(fig_folder, ...
            strcat(params.date, "_", params.time, "_", params.strain, "_polarplot_bar_DS_wArrow.png")));
        % Save as PDF
        exportgraphics(gca, fullfile(fig_folder, ...
            strcat(params.date, "_", params.time, "_", params.strain, "_polarplot_bar_DS_wArrow.pdf")), ...
            'ContentType', 'vector', 'BackgroundColor', 'none');
    end

end
