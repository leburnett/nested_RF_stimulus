%% Plot PD and ND 

results_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/control_RNAi/ON';
cd(results_folder)
files = dir('*_peak_vals*');

% Load the data: 
load(files(1).name)

%% Plot the PD and ND - linear timeseries plot

angls_deg = [0, 180, 22.5, 202.5, 45, 225, 67.5, 247.5, 90, 270, 112.5, 292.5, 135, 315, 157.5, 337.5];

figure;
numSubplots = 8; % Number of subplots (8 rows)
rowsPerPlot = 2; % Each subplot contains 2 rows of data

plot_data_pair = @(row1, col1, col2, col3) ...
    arrayfun(@(ii) plot(data{row1, ii}, 'Color', [col1, col2, col3]), 1:3, 'UniformOutput', false);

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
    plot_data_pair(row1, 0.8, 0.8, 0.8);
    plot_data_pair(row2, 1, 0.8, 0.8);

    % Plot the mean (4th column) for both rows
    plot(data{row1, 4}, 'k', 'LineWidth', 1.5); % Plot first row (black)
    plot(data{row2, 4}, 'r', 'LineWidth', 1.5); % Plot second row (red)

    % Plot the median voltage for the recording:
    plot([0 numel(data{row1, 4})], [median_voltage, median_voltage], 'k')

    % Title and formatting
    title(sprintf('Dir %.1f (k) & %.1f (r)', ang1, ang2));
    if i == 1
        ylabel('Voltage (mV)')
        xlabel('Time (s)')
    end 
    xticks([0 10000 20000])
    xticklabels({'0', '1', '2'})
    ylim([-70 -30]);
    hold off;
    
    % Create subplot for the second set (rows 17 to 32)
    subplot(2, numSubplots, i + numSubplots); % Create subplot for second set

    % Calculate row indices for the second set (rows 17-32)
    row3 = (i-1) * rowsPerPlot + 17; % First row in the pair (17-32)
    row4 = row3 + 1; % Second row in the pair (17-32)

    ang1 = angls_deg(row3 - 16);  % Adjust index for the second set
    ang2 = angls_deg(row4 - 16);  % Adjust index for the second set

    % Plot the first 3 columns for both rows
    hold on
    plot_data_pair(row3, 0.8, 0.8, 0.8);
    plot_data_pair(row4, 1, 0.8, 0.8);

    % Plot the mean (4th column) for both rows
    plot(data{row3, 4}, 'k', 'LineWidth', 1.5); % Plot first row (black)
    plot(data{row4, 4}, 'r', 'LineWidth', 1.5); % Plot second row (red)

    % Plot the median voltage for the recording:
    plot([0 numel(data{row3, 4})], [median_voltage, median_voltage], 'k')

    % Title and formatting
    ylim([-65 -40]);
    xticks([0 5000 10000])
    xticklabels({'0', '0.5', '1'})
    if i == 1
        ylabel('Voltage (mV)')
        xlabel('Time (s)')
    end 
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

data2 = cell(size(data));

for j = 1:height(data)
    for k = 1:4
        % Subtract the first element of each array in the cell array
        data2{j, k} = data{j, k} - data{j, k}(1);
    end
end

    
%% For initial value subtracted data:

figure;
numSubplots = 8; % Number of subplots (8 rows)
rowsPerPlot = 2; % Each subplot contains 2 rows of data

plot_data_pair = @(row1, col1, col2, col3) ...
    arrayfun(@(ii) plot(data2{row1, ii}, 'Color', [col1, col2, col3]), 1:3, 'UniformOutput', false);

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
    plot_data_pair(row1, 0.8, 0.8, 0.8);
    plot_data_pair(row2, 1, 0.8, 0.8);

    % Plot the mean (4th column) for both rows
    plot(data2{row1, 4}, 'k', 'LineWidth', 1.5); % Plot first row (black)
    plot(data2{row2, 4}, 'r', 'LineWidth', 1.5); % Plot second row (red)

    % Plot initial voltage
    plot([0 numel(data2{row1, 4})], [0, 0], 'k')

    % Title and formatting
    title(sprintf('Dir %.1f (k) & %.1f (r)', ang1, ang2));
    if i == 1
        ylabel('Voltage (mV)')
        xlabel('Time (s)')
    end 
    xticks([0 10000 20000])
    xticklabels({'0', '1', '2'})
    ylim([-5 35]);
    hold off;
    
    % Create subplot for the second set (rows 17 to 32)
    subplot(2, numSubplots, i + numSubplots); % Create subplot for second set

    % Calculate row indices for the second set (rows 17-32)
    row3 = (i-1) * rowsPerPlot + 17; % First row in the pair (17-32)
    row4 = row3 + 1; % Second row in the pair (17-32)

    ang1 = angls_deg(row3 - 16);  % Adjust index for the second set
    ang2 = angls_deg(row4 - 16);  % Adjust index for the second set

    % Plot the first 3 columns for both rows
    hold on
    plot_data_pair(row3, 0.8, 0.8, 0.8);
    plot_data_pair(row4, 1, 0.8, 0.8);

    % Plot the mean (4th column) for both rows
    plot(data2{row3, 4}, 'k', 'LineWidth', 1.5); % Plot first row (black)
    plot(data2{row4, 4}, 'r', 'LineWidth', 1.5); % Plot second row (red)

    % Plot initial voltage
    plot([0 numel(data2{row3, 4})], [0, 0], 'k')

    % Title and formatting
    ylim([-5 25]);
    xticks([0 5000 10000])
    xticklabels({'0', '0.5', '1'})
    if i == 1
        ylabel('Voltage (mV)')
        xlabel('Time (s)')
    end 
    hold off;
end

% Set figure size
f = gcf;
f.Position = [2 621 1799 350];






















