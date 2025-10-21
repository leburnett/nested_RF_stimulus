function resultant_angle = plot_polar_with_arrow(max_v, median_voltage, params, save_fig, fig_save_folder)

    % Convert max values for both conditions into polar format
    max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % slow bars
    max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % fast bars
    max_v_polar3 = vertcat(max_v(:, 3), max_v(1, 3)); % very fast bars

    
 
    theta  = linspace(0, 2*pi, 17);
    resultant_magnitude = 30; 
    colors = {[0.2 0.4 0.7], [0.4 0.8 1], [0.45, 0.0, 0.55]};
    arrows_colors = {[0.3 0.3 0.3], [0.6 0.6 0.6], [0.85 0.85 0.85]};

    
    figure;
    polarplot(theta, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    hold on
    polarplot(theta, max_v_polar2 - median_voltage, 'Color', colors{2}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    polarplot(theta, max_v_polar3 - median_voltage, 'Color', colors{3}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    ax = gca;
    ax.LineWidth = 1.2;
    ax.FontSize = 15;
    ax.ThetaTick = rad2deg(theta);
    ax.ThetaTickLabel = {};
    hold on

    for sp = 1:3

        if sp == 1
            array = max_v_polar1;
        elseif sp == 2
            array = max_v_polar2;
        elseif sp == 3
            array = max_v_polar3;
        end 
        rho = [array(1:end-1, 1) - median_voltage]';
        [~, resultant_angle] = vector_sum_polar(rho, theta(1:end-1));

        add_arrow_to_polarplot(resultant_magnitude, resultant_angle, arrows_colors{sp})
    end 

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