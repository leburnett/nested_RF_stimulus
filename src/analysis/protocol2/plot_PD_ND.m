%% Plot PD and ND 

%% Figure with the PD and ND plotted on top of each other. 

angls_deg = [0, 180, 22.5, 202.5, 45, 225, 67.5, 247.5, 90, 270, 112.5, 292.5, 135, 315, 157.5, 337.5];

 % 0    22.5   45   67.5   90   112.5   135   157.5   
 % 180  202.5  225  247.5  270  292.5   315   337.5

figure;
numSubplots = 8; % Number of subplots (8 rows)
rowsPerPlot = 2; % Each subplot contains 2 rows of data

for i = 1:numSubplots
    subplot(2, numSubplots, i); % Create subplot
    
    % Get the row indices for this subplot
    row1 = (i-1) * rowsPerPlot + 1; % First row in the pair
    row2 = row1 + 1; % Second row in the pair

    ang1 = angls_deg(row1);
    ang2 = angls_deg(row2);

    for ii = 1:3
    plot(data{row1, ii}, 'Color', [0.8 0.8 0.8]); 
    plot(data{row2, ii}, 'Color', [1 0.8 0.8]); 
    hold on
    end 

    % Mean 
    plot(data{row1, 4}, 'k', 'LineWidth', 1.5); % Plot first row (blue)
    plot(data{row2, 4}, 'r', 'LineWidth', 1.5); % Plot second row (red)
    
    title(sprintf('Dir %.1f (k) & %.1f (r)', ang1, ang2));
    hold off;

    ylim([-70 -30])
end

for i = 1:numSubplots
    subplot(2, numSubplots, i+8); % Create subplot
    
    % Get the row indices for this subplot
    row1 = (i-1) * rowsPerPlot + 17; % First row in the pair
    row2 = row1 + 1; % Second row in the pair

    ang1 = angls_deg(row1-16);
    ang2 = angls_deg(row2-16);

    for ii = 1:3
        plot(data{row1, ii}, 'Color', [0.8 0.8 0.8]); 
        plot(data{row2, ii}, 'Color', [1 0.8 0.8]); 
        hold on
    end 

    plot(data{row1, 4}, 'k',  'LineWidth', 1.5); 
    hold on;
    plot(data{row2, 4}, 'r',  'LineWidth', 1.5); 
    
    title(sprintf('Dir %.1f (k) & %.1f (r)', ang1, ang2));
    hold off;

    ylim([-65 -40])
end

f = gcf;
f.Position = [2  621  1799  350];


%% Plot the PD and ND - linear timeseries plot

angls_deg = [0, 180, 22.5, 202.5, 45, 225, 67.5, 247.5, 90, 270, 112.5, 292.5, 135, 315, 157.5, 337.5];

figure;
numSubplots = 8; % Number of subplots (8 rows)
rowsPerPlot = 2; % Each subplot contains 2 rows of data

% Function to plot data for a given row pair
plot_data_pair = @(row1, row2, ang1, ang2, yLimits) ...
    arrayfun(@(ii) plot(data{row1, ii}, 'Color', [0.8 0.8 0.8]), 1:3, 'UniformOutput', false); % Plot the first 3 columns

% Combined loop for both sets of subplots
for i = 1:numSubplots
    subplot(2, numSubplots, i); % Create subplot for the first 8 subplots

    % Calculate row indices for the first set (rows 1 to 16)
    row1 = (i-1) * rowsPerPlot + 1; % First row in the pair (1-16)
    row2 = row1 + 1; % Second row in the pair (1-16)

    ang1 = angls_deg(row1);
    ang2 = angls_deg(row2);

    % Plot the first 3 columns for both rows
    hold on
    plot_data_pair(row1, row2, ang1, ang2, [-70 -30]);

    % Plot the mean (4th column) for both rows
    plot(data{row1, 4}, 'k', 'LineWidth', 1.5); % Plot first row (black)
    plot(data{row2, 4}, 'r', 'LineWidth', 1.5); % Plot second row (red)

    % Title and formatting
    title(sprintf('Dir %.1f (k) & %.1f (r)', ang1, ang2));
    ylim([-70 -30]);
    hold off;
    
    % Create subplot for the second set (rows 17 to 32)
    subplot(2, numSubplots, i + numSubplots); % Create subplot for second set

    % Calculate row indices for the second set (rows 17-32)
    row1 = (i-1) * rowsPerPlot + 17; % First row in the pair (17-32)
    row2 = row1 + 1; % Second row in the pair (17-32)

    ang1 = angls_deg(row1 - 16);  % Adjust index for the second set
    ang2 = angls_deg(row2 - 16);  % Adjust index for the second set

    % Plot the first 3 columns for both rows
    hold on
    plot_data_pair(row1, row2, ang1, ang2, [-65 -40]);

    % Plot the mean (4th column) for both rows
    plot(data{row1, 4}, 'k', 'LineWidth', 1.5); % Plot first row (black)
    plot(data{row2, 4}, 'r', 'LineWidth', 1.5); % Plot second row (red)

    % Title and formatting
    title(sprintf('Dir %.1f (k) & %.1f (r)', ang1, ang2));
    ylim([-65 -40]);
    hold off;
end

% Set figure size
f = gcf;
f.Position = [2 621 1799 350];


%% Plot a circle showing the bar direction in degrees

figure; 
% Generate angles for circle points
theta = 0:pi/8:2*pi;
% Plot the circle
centers = [0, 0]; % Center point
radii = 2; 
h = viscircles(centers, radii, 'EdgeColor', [0.8 0.8 0.8]); 
hold on 
for ang = theta(1:end-1)
    if rad2deg(ang) <180
        col = 'k-';
    else
        col = 'r-';
    end 
    plot([x_center, x_center + radius*cos(ang)], [y_center, y_center + radius*sin(ang)], col); 
    text( x_center + radius*cos(ang), y_center + radius*sin(ang), string(rad2deg(ang)), 'FontSize', 20, 'HorizontalAlignment','center')
end 
axis equal; 
xlim([-2.5 2.5])
ylim([-2.5 2.5])
box off
axis off

%%