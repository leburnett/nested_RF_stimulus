function resultant_angle = plot_polar_with_arrow(max_v, median_voltage, params, save_fig)

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
    polarplot(theta, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 3);
    hold on
    polarplot(theta, max_v_polar2 - median_voltage, 'Color', colors{2}, 'LineWidth', 3);
    ax = gca;
    ax.LineWidth = 1.2;
    ax.FontSize = 15;
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
                , strcat(params.date, "_", params.time, "_", params.strain, "_polarplot_20dps_bar_DS_wArrow.png") ...
                );
        % Save as PDF
        exportgraphics(gca ...
                , strcat(params.date, "_", params.time, "_", params.strain, "_polarplot_20dps_bar_DS_wArrow.pdf") ...
                , 'ContentType', 'vector' ...
                , 'BackgroundColor', 'none' ...
                ); 
    end 


end 