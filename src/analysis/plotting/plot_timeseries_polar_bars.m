function [max_v, min_v] = plot_timeseries_polar_bars(data, median_voltage, params, save_fig, fig_save_folder)
    % Central polar plot of maximum voltage for each direction, with timeseries of voltage data for each moving bar stimulus.
    % Responses to fast bars are in light blue and responses to slow bars in dark blue. 

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

    % Outputs
    % -------

    %       max_v: arrray [n_conditions, n_speeds]
    %           98th percentile value of the mean voltage during each bar condition. Range of data includes the time 
    %           when bar stimulus is shown + 200ms after the bar stimulus ends. Excludes the interval time before 
    %           the bar stimulus starts.

    %       min_v: arrray [n_conditions, n_speeds]
    %           2nd percentile value of the mean voltage during the second half of each bar condition. Range of data 
    %           includes the time when bar stimulus is shown + 200ms after the bar stimulus ends. Excludes the interval
    %           time before the bar stimulus starts.   

    % ________________________________________________________________________________________
    
    % Number of subplots
    numPlots = 16;
    theta = linspace(0, 2*pi, numPlots+1); % Angles (add 2*pi to complete the circle)
    theta = theta(1:end-1); % Remove redundant last point
    
    % Center and radius of the circle
    centerX = 0.5;
    centerY = 0.5;
    radius = 0.35;
    
    % Define central polar plot position
    centralSize = (2 * radius) * 0.65; 
    centralPosition = [centerX - centralSize/2, centerY - centralSize/2, centralSize, centralSize];
    
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

    angls = linspace(0, 2*pi, 17); % 17 points to include 0 and 2π
    
    % Preallocate max/min voltage arrays
    max_v = zeros(numPlots, 2);
    min_v = zeros(numPlots, 2);
    
    % Define colors for the two conditions
    colors = {[0.2 0.4 0.7], [0.4 0.8 1], [0.45, 0.0, 0.55]};  % Dark blue (28 dps) and Light blue (56 dps)

    n_conditions = numel(plot_order);
    
    %% Create the figure
    figure
    
    for sp = 1:1:3
        col = colors{sp}; % Get color for current condition
    
        % Loop to create subplots
        for i = 1:numPlots

            % Compute subplot position
            x = centerX + radius * cos(theta(i)); 
            y = centerY + radius * sin(theta(i));
            subplotWidth = 0.15; 
            subplotHeight = 0.15;
            subplotPosition = [x - subplotWidth/2, y - subplotHeight/2, subplotWidth, subplotHeight];
    
            % Create subplot axes
            ax = axes('Position', subplotPosition);
            hold on
    
            % Get the data index - which row in 'data' contains the data for
            % the desired direction.
            d_idx = plot_order(i) + n_conditions*(sp-1);
    
            n_reps = size(data, 2)-1;

            % Plot the timeseries data for each repetition and the mean
            for r = 1:n_reps+1
                d2plot = data{d_idx, r};
                x_vals = 1:numel(d2plot);
    
                if r == 1 
                    % Plot median voltage background
                    plot([1 x_vals(end)], [median_voltage, median_voltage], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);  
                end 
    
                % Plot repetitions in gray, last rep (mean) in condition color
                if r < n_reps+1
                    plot(x_vals, d2plot, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
                else 
                    plot(x_vals, d2plot, 'Color', col, 'LineWidth', 1.2);
                    % plot([9000 9000],[-80 -10], 'Color', 'k', 'LineWidth', 0.1) % Start of bar stimulus
                    % plot([numel(x_vals)-9000 numel(x_vals)-9000],[-80 -10], 'Color', 'k', 'LineWidth', 0.1) % End of bar stimulus.
                end
            end 
    
            % Set consistent Y-limits
            ylim([-80 -10])

            % Store max/min values from the mean voltage per condition.
            d = data{d_idx, n_reps+1}; 

            % Find the mean voltage in the 900ms before the flash. 
            d_before_flash = d(1000:9000);
            mean_before = mean(d_before_flash);

            % Remove the 900ms before and last 700ms after stimulus to look
            % for min / max response.
            d = d(9000:end-7000);

            % Number of data points now:
            n_vals_d = numel(d);

            % max_v(i, sp) = prctile(d, 98);
            max_v(i, sp) = abs(diff([prctile(d, 98), mean_before])); % ignore the first 8th of the condition

            min_v(i, sp) = prctile(d(floor(n_vals_d/2):end), 2); % find the min during the 2nd half of the condition.
    
            % Turn off axes for better visualization
            axis(ax, 'off');

            % % % TODO - - - Add lines when the stimulus starts / stops? - 9000
            % datapoints added before and after. 
        end
    end
    
    % Add polar plot in the center
    axCentral = polaraxes('Position', centralPosition, 'ThetaTick', rad2deg(angls));
    hold on
    
    % Convert max values for both conditions into polar format
    max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % slow bars
    max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % fast bars
    max_v_polar3 = vertcat(max_v(:, 3), max_v(1, 3)); % fast bars
    
    % % Plot the polar data
    % polarplot(angls, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 2);
    % polarplot(angls, max_v_polar2 - median_voltage, 'Color', colors{2}, 'LineWidth', 2);
    % polarplot(angls, max_v_polar3 - median_voltage, 'Color', colors{3}, 'LineWidth', 2);
    
    % Plot the polar data
    polarplot(angls, max_v_polar1, 'Color', colors{1}, 'LineWidth', 2);
    polarplot(angls, max_v_polar2, 'Color', colors{2}, 'LineWidth', 2);
    polarplot(angls, max_v_polar3, 'Color', colors{3}, 'LineWidth', 2);
    
    % Add title
    sgtitle(sprintf("28 / 56 / 168 dps - 4 pixel bar stimuli - 30 pix square - %s - %s - %s - %s", ...
                    strrep(params.date, '_', '-'), strrep(params.time, '_', '-'), strrep(params.strain, '_', '-'), params.on_off));
    
    % Set figure size
    set(gcf, 'Position', [303 78 961 969]);

    if save_fig
        fig_folder = fullfile(fig_save_folder, params.on_off);
        if ~isfolder(fig_folder)
            mkdir(fig_folder)
        end 
        f = gcf; 
        fname = fullfile(fig_folder, strcat(params.strain, '_', params.on_off, '_', params.date, '_', params.time,'_bar_polar.pdf'));
        exportgraphics(f ...
                , fname ...
                , 'ContentType', 'vector' ...
                , 'BackgroundColor', 'none' ...
                ); 
    end 

end 