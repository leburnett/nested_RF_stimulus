function [max_v, min_v] = plot_timeseries_polar_bars_flex(data, median_voltage, params, save_fig, fig_save_folder)
% PLOT_TIMESERIES_POLAR_BARS_FLEX  Speed-flexible circular figure with polar plot and timeseries.
%
%   [MAX_V, MIN_V] = PLOT_TIMESERIES_POLAR_BARS_FLEX(DATA, MEDIAN_VOLTAGE, ...
%       PARAMS, SAVE_FIG, FIG_SAVE_FOLDER)
%   generates a publication-quality figure showing bar stimulus responses
%   arranged radially with a central polar plot summarizing direction tuning.
%   Automatically detects the number of bar sweep speeds from the data
%   dimensions and adjusts plotting accordingly.
%
%   INPUTS:
%     data           - Nx4 cell array from PARSE_BAR_DATA:
%                      N = 16 * n_speeds (32 for 2 speeds, 48 for 3)
%                      Rows 1-16: slow (28 dps) bar directions
%                      Rows 17-32: fast (56 dps) bar directions
%                      Rows 33-48: very fast (168 dps) bar directions (if present)
%                      Columns 1-3: individual repetitions
%                      Column 4: mean across repetitions
%     median_voltage - Baseline voltage for normalization
%     params         - Structure with: .date, .time, .strain, .on_off
%     save_fig       - Boolean, true to save figure as PDF
%     fig_save_folder - Directory for saving figures
%
%   OUTPUTS:
%     max_v - 16 x n_speeds array of 98th percentile response per direction/speed
%     min_v - 16 x n_speeds array of 2nd percentile response per direction/speed
%
%   FIGURE LAYOUT:
%     - 16 small subplots arranged in a circle, one per direction
%     - Each subplot shows the mean trace colored by speed
%     - Central polar plot shows max response vs direction
%     - Dark blue = slow (28 dps), light blue = fast (56 dps),
%       purple = very fast (168 dps, if present)
%
%   SPEED DETECTION:
%     n_speeds = size(data, 1) / 16
%     Supports 2 speeds (Summer 2025) and 3 speeds (Autumn 2025)
%
%   See also PARSE_BAR_DATA, PLOT_POLAR_WITH_ARROW_FLEX, PROCESS_BARS_P2
% ________________________________________________________________________________________

    % Detect number of speeds from data dimensions
    n_data_rows = size(data, 1);
    n_speeds = n_data_rows / 16;

    % Number of subplots (directions)
    numPlots = 16;
    theta = linspace(0, 2*pi, numPlots+1);
    theta = theta(1:end-1);

    % Center and radius of the circle
    centerX = 0.5;
    centerY = 0.5;
    radius = 0.35;

    % Define central polar plot position
    centralSize = (2 * radius) * 0.65;
    centralPosition = [centerX - centralSize/2, centerY - centralSize/2, centralSize, centralSize];

    %% Direction ordering
    % ON/OFF-aware plot_order: OFF cells have two directions swapped due to
    % a pattern generation issue with dark bars (bkg4).
    if params.on_off == "on"
        plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    elseif params.on_off == "off"
        plot_order = [1,3,5,7,9,11,14,16,2,4,6,8,10,12,13,15];
    else
        plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];
    end

    angls = linspace(0, 2*pi, 17); % 17 points to include 0 and 2pi

    % Preallocate max/min voltage arrays — dynamic columns
    max_v = zeros(numPlots, n_speeds);
    min_v = zeros(numPlots, n_speeds);

    % Colors for each speed condition
    all_colors = {[0.2 0.4 0.7], [0.4 0.8 1], [0.45, 0.0, 0.55]};
    colors = all_colors(1:n_speeds);

    % Speed labels for title
    all_speed_labels = {'28', '56', '168', '250', '500'};
    speed_labels = all_speed_labels(1:min(n_speeds, numel(all_speed_labels)));

    n_conditions = numPlots; % 16 directions per speed

    %% Create the figure
    figure

    for sp = 1:n_speeds
        col = colors{sp};

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

            % Get the data index
            d_idx = plot_order(i) + n_conditions*(sp-1);

            n_reps = size(data, 2) - 1;

            % Plot the timeseries data — mean trace in color
            for r = 1:n_reps+1
                d2plot = data{d_idx, r};
                x_vals = 1:numel(d2plot);

                if r == 1
                    % Plot median voltage background
                    plot([1 x_vals(end)], [median_voltage, median_voltage], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
                end

                % Plot mean trace in condition color (skip individual reps for clarity)
                if r == n_reps + 1
                    plot(x_vals, d2plot, 'Color', col, 'LineWidth', 1.2);
                end
            end

            % Set consistent Y-limits
            ylim([-80 -10])

            % Compute max/min values from the mean voltage per condition
            d = data{d_idx, n_reps+1};

            % Remove the 900ms before and last 700ms after stimulus
            d = d(9000:end-7000);
            n_vals_d = numel(d);

            % 98th percentile of response (raw value, Summer convention)
            max_v(i, sp) = prctile(d, 98);
            % 2nd percentile during second half
            min_v(i, sp) = prctile(d(floor(n_vals_d/2):end), 2);

            % Turn off axes for better visualization
            axis(ax, 'off');
        end
    end

    % Add polar plot in the center
    axCentral = polaraxes('Position', centralPosition, 'ThetaTick', rad2deg(angls));
    hold on

    % Plot polar traces for each speed — subtract median for positive values
    for sp = 1:n_speeds
        max_v_polar = vertcat(max_v(:, sp), max_v(1, sp));
        polarplot(angls, max_v_polar - median_voltage, 'Color', colors{sp}, 'LineWidth', 2);
    end

    % Add title with speed labels
    speed_str = strjoin(speed_labels, ' / ');
    sgtitle(sprintf("%s dps - 4 pixel bar stimuli - 30 pix square - %s - %s - %s - %s", ...
        speed_str, strrep(params.date, '_', '-'), strrep(params.time, '_', '-'), ...
        strrep(params.strain, '_', '-'), params.on_off));

    % Set figure size
    set(gcf, 'Position', [303 78 961 969]);

    if save_fig
        fig_folder = fullfile(fig_save_folder, params.on_off);
        if ~isfolder(fig_folder)
            mkdir(fig_folder)
        end
        f = gcf;
        fname = fullfile(fig_folder, strcat(params.strain, '_', params.on_off, '_', params.date, '_', params.time, '_bar_polar.pdf'));
        exportgraphics(f, fname, 'ContentType', 'vector', 'BackgroundColor', 'none');
    end

end
