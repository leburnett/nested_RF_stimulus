function plot_each_rep_28dps_pharma(data, params, save_fig, fig_save_folder)
    
    % This is from looking at the movement of the bar stimuli in real life. 
    % The order should always be the same, this should not change. 
    % Bars move left to right, then right to left, then in a counter clockwise
    % fashion. Always moving in one direction then the opposite.
    
    % For plotting - in MATLAB, the first plot is positioned in the 'E'
    % position of a compass and moves counter clockwise.
    
    % PD and ND pairs are ordered sequentially in 'data'. e.g. data{1, 4} is
    % left to right and data{2, 4} is right to left.
    
    % Dark bars - two directions are mixed up. Occurred when making patterns with
    % bkg4... 
    plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

    angls = linspace(0, 2*pi, 17); % 17 points to include 0 and 2π
   
    % Define colors for the three repetitions
    colors = {[0.53, 0.81, 0.98], [0.0, 0.0, 0.8], [0.10, 0.10, 0.44]};

    interval_len = [9000, 9000, 5000, 4000, 4000];

     % Preallocate max/min voltage arrays
    max_v = zeros(16, 3);

    %% Create the figure
    figure
    
    tiledlayout(1, 4, 'TileSpacing', 'compact');

    for rep = 1:3
  
        for i = 1:16 % directions 
            
            % Get the data index - which row in 'data' contains the data for
            % the desired direction.
            % d_idx = plot_order(i) + n_conditions*(sp-1);
            d_idx = plot_order(i); % only 28 dps. 
            % d_idx = plot_order(i)+16; % only 56 dps. 

            % Store max/min values from the repetition
            d = data{d_idx, rep}; 

            int_l = interval_len(rep);

            % Find the mean voltage in the 900ms before the flash. 
            d_before_flash = d(1:int_l);
            mean_before = mean(d_before_flash);

            % Remove the 900ms before and last 700ms after stimulus to look
            % for min / max response.
            d = d(int_l:end-(int_l*0.75));

            % max_v(i, sp) = prctile(d, 98);
            max_v(i, rep) = abs(diff([prctile(d, 98), mean_before])); % ignore the first 8th of the condition
        end
        nexttile
        max_v_polar = vertcat(max_v(:, rep), max_v(1, rep)); %
        polarplot(angls, max_v_polar, 'Color', colors{rep}, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
        title(strcat("Rep ", string(rep)), 'FontSize', 15);
        rlim([0 30])
    end
    
    nexttile
    for rr = 1:3
        max_v_polar = vertcat(max_v(:, rr), max_v(1, rr)); %
        polarplot(angls, max_v_polar, 'Color', colors{rr}, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
        hold on
    end 
    title("All reps");
    rlim([0 30])
    
    % Add title
    sgtitle(sprintf("28 dps - 4 pixel bar stimuli - 30 pix square - %s - %s - %s - %s - %s %s", ...
                    strrep(params.date, '_', '-'), strrep(params.time, '_', '-'), strrep(params.strain, '_', '-'), params.on_off, params.application, params.drug), 'FontSize', 16);
    
    % Set figure size
    set(gcf, 'Position', [34  710  1156  319]);

    if save_fig
        fig_folder = fullfile(fig_save_folder, params.on_off);
        if ~isfolder(fig_folder)
            mkdir(fig_folder)
        end 
        f = gcf; 
        fname = fullfile(fig_folder, strcat(params.strain, '_', params.on_off, '_', params.date, '_', params.time, '_', params.application, '_', params.drug, '_bar_polar_per_rep.pdf'));
        exportgraphics(f ...
                , fname ...
                , 'ContentType', 'vector' ...
                , 'BackgroundColor', 'none' ...
                ); 
    end 

end 