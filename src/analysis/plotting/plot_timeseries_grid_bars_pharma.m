function plot_timeseries_grid_bars_pharma(data, median_v, params, save_fig, fig_save_folder)

    % Inputs
    % ------
    %       data: cell array [n_conditions * 2, n_reps +1]
    %           Timeseries voltage data for each condition (rows). The
    %           first (n_reps) columns are the repetition data and the
    %           (n_reps +1) column is the mean data across the reps. 

    %       median_voltage: float
    %           Median voltage across the entire recording.

    %       params: struct
    %           Contains meta data about the experiment. Used for saving. 

    %       save_fig : bool
    %           Boolean value for whether to save the figures of not. If
    %           True, then the figures are saved to `fig_sve_folder`, if
    %           False then the figures are created but not saved. 

    %       fig_save_folder : str
    %           Path where to save the figures.           
    % ________________________________________________________________________________________
    
    % Number of subplots
    numDir = 16;
    numSpeeds = 5;
    
    %% The order in which the data for the different directions are stored:
    
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
    
    % Define colors for the two conditions
    colors = {[0.2 0.4 0.7], [0.4 0.8 1], [0.45, 0.0, 0.55], [0.99, 0.78, 0.8], [0.85, 0 , 0.7]};  % Dark blue (28 dps) and Light blue (56 dps), purple = 168, light pink = 350, magenta = 500

    interval_len = [9000, 9000, 5000, 3000, 3000];

    n_conditions = numel(plot_order);
    
    %% Create the figure
    figure
    idx = 1;

    % Loop to create subplots
    for i = 1:numDir

        for sp = 1:numSpeeds

            col = colors{sp}; % Get color for current condition
            int_l = interval_len(sp);
    
            subplot(numDir, numSpeeds, idx);
            hold on;
    
            % Get the data index - which row in 'data' contains the data for
            % the desired direction.
            d_idx = plot_order(i) + n_conditions*(sp-1);
    
            n_reps = size(data, 2)-1;
    
            % Plot the timeseries data for each repetition and the mean
            for r = 1:n_reps+1
    
                % Extract data
                d2plot = data{d_idx, r};
    
                % x values
                x_vals = 1:numel(d2plot);
    
                if r  == 1
                    % Plot median voltage background
                    plot([1 x_vals(end)], [median_v, median_v], 'Color', [0.7 0.7 0.7], 'LineWidth', 1); 
                end 
    
                % Plot repetitions in gray, last rep (mean) in condition color
                if r < n_reps+1
                    plot(x_vals, d2plot, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
                else 
                    plot(x_vals, d2plot, 'Color', col, 'LineWidth', 1.2);
                    plot([int_l int_l],[-80 -10], 'Color', 'k', 'LineWidth', 0.1) % Start of bar stimulus
                    % plot([numel(x_vals)-9000 numel(x_vals)-9000],[-80 -10], 'Color', 'k', 'LineWidth', 0.1) % End of bar stimulus.
                end
    
             end 
    
            % Set consistent Y-limits
            ylim([-70 -30])
            xlim([0 max(x_vals)])

            if idx > 1
                ax = gca;
                ax.YAxis.Visible = "off";
                if idx > 5
                    ax.XAxis.Visible = "off";
                end 
            end 

            if idx == 1
                title("28 dps");
            elseif idx == 2
                title("56 dps");
            elseif idx == 3
                title("168 dps");
            elseif idx == 4
                title("250 dps");
            elseif idx == 5
                title("500 dps");
            end 

            idx = idx+1;
        end
    end
   
    % Add title
    sgtitle(sprintf("28 / 56 / 168 / 250 / 500 dps - %s - %s - %s - %s - %s %s", ...
                    strrep(params.date, '_', '-'), strrep(params.time, '_', '-'), strrep(params.strain, '_', '-'), params.on_off, params.application, params.drug), 'FontSize', 16);
    
    % Set figure size
    set(gcf, 'Position', [2023  -522  747  1588]);

    if save_fig
        fig_folder = fullfile(fig_save_folder, params.on_off);
        if ~isfolder(fig_folder)
            mkdir(fig_folder)
        end 
        f = gcf; 
        fname = fullfile(fig_folder, strcat(params.strain, '_', params.on_off, '_', params.date, '_', params.time, '_', params.application, '_', params.drug, '_bar_grid.pdf'));
        exportgraphics(f ...
                , fname ...
                , 'ContentType', 'vector' ...
                , 'BackgroundColor', 'none' ...
                ); 
    end 

end 























