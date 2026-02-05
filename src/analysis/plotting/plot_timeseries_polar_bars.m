function [max_v, min_v] = plot_timeseries_polar_bars(data, median_voltage, params, save_fig, fig_save_folder)
% PLOT_TIMESERIES_POLAR_BARS  Create circular figure with polar plot and timeseries.
%
%   [MAX_V, MIN_V] = PLOT_TIMESERIES_POLAR_BARS(DATA, MEDIAN_VOLTAGE, ...
%       PARAMS, SAVE_FIG, FIG_SAVE_FOLDER)
%   generates a publication-quality figure showing bar stimulus responses
%   arranged radially with a central polar plot summarizing direction tuning.
%
%   INPUTS:
%     data           - 32x4 cell array from PARSE_BAR_DATA:
%                      Rows 1-16: slow (28 dps) bar directions
%                      Rows 17-32: fast (56 dps) bar directions
%                      Columns 1-3: individual repetitions
%                      Column 4: mean across repetitions
%     median_voltage - Baseline voltage for normalization
%     params         - Structure with: .date, .time, .strain, .on_off
%     save_fig       - Boolean, true to save figure as PDF
%     fig_save_folder - Directory for saving figures
%
%   OUTPUTS:
%     max_v - 16x2 array of 98th percentile response per direction/speed
%             Column 1: slow bars, Column 2: fast bars
%     min_v - 16x2 array of 2nd percentile response per direction/speed
%
%   FIGURE LAYOUT:
%     - 16 small subplots arranged in a circle, one per direction
%     - Each subplot shows all 3 reps (gray) and mean (colored)
%     - Central polar plot shows max response vs direction
%     - Dark blue = slow (28 dps), light blue = fast (56 dps)
%
%   DIRECTION ORDER:
%     Subplots start at East (0 degrees) and proceed counterclockwise.
%     Data is reordered from forward/backward pairs to sequential angles.
%
%   OUTPUT FILE:
%     <strain>_<on_off>_<date>_<time>_bar_polar.pdf
%
%   See also PARSE_BAR_DATA, PLOT_POLAR_WITH_ARROW, PROCESS_BARS_P2
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
    if params.on_off == "on"
        % dark bars: 
        plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    elseif params.on_off == "off"
        % light bars: 
        plot_order = [1,3,5,7,9,11,14,16,2,4,6,8,10,12,13,15];
    end 
    
    angls = linspace(0, 2*pi, 17); % 17 points to include 0 and 2π
    
    % Preallocate max/min voltage arrays
    max_v = zeros(numPlots, 2);
    min_v = zeros(numPlots, 2);
    
    % Define colors for the two conditions
    colors = {[0.2 0.4 0.7], [0.4 0.8 1]};  % Dark blue (28 dps) and Light blue (56 dps)

    n_conditions = size(data, 1)/2;
    
    %% Create the figure
    figure
    
    for sp = 1:2
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
                    plot([9000 9000],[-80 -10], 'Color', 'k', 'LineWidth', 0.1) % Start of bar stimulus
                    % plot([numel(x_vals)-9000 numel(x_vals)-9000],[-80 -10], 'Color', 'k', 'LineWidth', 0.1) % End of bar stimulus.
                end
            end 
    
            % Set consistent Y-limits
            ylim([-80 -10])
    
            % % % % TODO -- -- Update this to make the range of time over
            % which we are looking at - 9000 start / end. 

            % Store max/min values from the mean voltage per condition.
            d = data{d_idx, n_reps+1}; 

            % Remove the 900ms before and last 700ms after stimulus to look
            % for min / max response.
            d = d(9000:end-7000);

            % Number of data points now:
            n_vals_d = numel(d);

            max_v(i, sp) = prctile(d, 98); % ignore the first 8th of the condition
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
    
    % Plot the polar data
    polarplot(angls, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 2);
    polarplot(angls, max_v_polar2 - median_voltage, 'Color', colors{2}, 'LineWidth', 2);
    
    % Add title
    sgtitle(sprintf("28 / 56 dps - 4 pixel bar stimuli - 30 pix square - %s - %s - %s - %s", ...
                    strrep(params.date, '_', '-'), strrep(params.time, '_', '-'), params.strain, params.on_off));
    
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