function resultant_angle = plot_polar_with_arrow_pharma(max_v, median_voltage, params, save_fig, fig_save_folder)

    % Convert max values for both conditions into polar format
    max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % 28 dps
    max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % 56 dps
    max_v_polar3 = vertcat(max_v(:, 3), max_v(1, 3)); % 168 dps
    max_v_polar4 = vertcat(max_v(:, 4), max_v(1, 4)); % 250 dps
    max_v_polar5 = vertcat(max_v(:, 5), max_v(1, 5)); % 500 dps

    theta  = linspace(0, 2*pi, 17);
    resultant_magnitude = 30; 
    colors = {[0.2 0.4 0.7], [0.4 0.8 1], [0.45, 0.0, 0.55], [0.99, 0.78, 0.8], [0.85, 0 , 0.7]}; 
    % arrows_colors = {[0.3 0.3 0.3], [0.6 0.6 0.6], [0.85 0.85 0.85]};
    arrows_colors = {[0.2 0.2 0.2], [0.45 0.45 0.45], [0.6 0.6 0.6], [0.75 0.75 0.75], [0.9 0.9 0.9]};

    
    figure;
    % polarplot(theta, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    % hold on
    % polarplot(theta, max_v_polar2 - median_voltage, 'Color', colors{2}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    % polarplot(theta, max_v_polar3 - median_voltage, 'Color', colors{3}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
   
    polarplot(theta, max_v_polar1, 'Color', colors{1}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    hold on
    polarplot(theta, max_v_polar2, 'Color', colors{2}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    polarplot(theta, max_v_polar3, 'Color', colors{3}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    polarplot(theta, max_v_polar4, 'Color', colors{4}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    polarplot(theta, max_v_polar5, 'Color', colors{5}, 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    
    ax = gca;
    ax.LineWidth = 1.2;
    ax.FontSize = 15;
    ax.ThetaTick = rad2deg(theta);
    ax.ThetaTickLabel = {};
    hold on

    for sp = 1:5

        if sp == 1
            array = max_v_polar1;
        elseif sp == 2
            array = max_v_polar2;
        elseif sp == 3
            array = max_v_polar3;
        elseif sp == 4
            array = max_v_polar4;
        elseif sp == 5
            array = max_v_polar5;
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

        % % Save as PNG
        % exportgraphics(gca ...
        %         , strcat(params.date, "_", params.time, "_", params.strain, '_', params.application, '_', params.drug, "_polarplot_bar_DS_wArrow.png") ...
        %         );
        % Save as PDF
        exportgraphics(gca ...
                , fullfile(fig_folder, strcat(params.date, "_", params.time, "_", params.strain, '_', params.application, '_', params.drug, "_polarplot_bar_DS_wArrow.pdf")) ...
                , 'ContentType', 'vector' ...
                , 'BackgroundColor', 'none' ...
                ); 
    end 


end 