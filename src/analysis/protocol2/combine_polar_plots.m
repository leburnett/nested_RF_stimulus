
% Load in the data from the 'results' folder
% data_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/control_RNAi';

file_list = dir('*_peak_vals*');
n_exp = height(file_list);
col = [0.8 0.8 0.8];

% For polar plots - max vals. 
dd = [];

data_aligned = cell(32, n_exp);
order_idx = nan(16, n_exp);

plot_order= [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

% % % data structure to save aligned, centred and re-ordered data across
% experiments, 

figure
for i = 1:n_exp

    % Load the data from one cell:
    load(file_list(i).name);

    %% 1 - rearrange data to be in sequential order for each angle. 
    % Before, it was in the order in which the stimuli were presented. Same
    % orientations, in the opposite directions sequentially. Whereas now, they
    % move around a circle in 1/16*pi sections. 
   
    data_ordered = cell(32, 1);
    
    for k = 1:32
    
        if k >16
            data_idx = plot_order(k-16);
            data_idx = data_idx+16;
        else 
            data_idx = plot_order(k);
        end 
        data_ordered(k, 1) = data(data_idx, 4);
    end 
    
    angls = angls(1:end-1)'; 


    %% 2 - plot the timeseries responses of the individual cell to the 16 different directions
    % Plots the mean of 3 reps per each direction on an individual subplot.
    % slow - 28 dps = top row. 
    % fast - 56 dps - bottom row. 
    
    % f1 = figure;
    % for i2 = 1:32
    % 
    %     subplot(2,16,i2)
    % 
    %     plot(data_ordered{i2, 1}*10, 'k');
    %     n_xvals = numel(data_ordered{i2, 1});
    %     hold on 
    %     plot([0 n_xvals], [median_voltage median_voltage], 'r', 'LineWidth', 0.3)
    % 
    %     box off
    %     ax = gca;
    %     % ax.XAxis.Visible = 'off';
    %     ax.TickLength = [0.04 0.04];
    %     ax.TickDir = 'out';
    % 
    %     if i2>16
    %         ylim([-75 -40])
    %         ax.XTickLabel = {'0', '0.5', '1'};
    %     else
    %         ylim([-75 -30])
    %         ax.XTickLabel = {'0', '1', '2'};
    %     end 
    % 
    %     if i2 ==1 || i2 == 17
    %         ylabel('Voltage (mV)');
    %         xlabel('Time (s)')
    %     else
    %         ax.YAxis.Visible = 'off';
    %     end 
    % 
    %     if i2 == 8 
    %         title('28 dps')
    %     elseif i2 == 24
    %         title('56 dps')
    %     end 
    % 
    % end 
    % 
    % f1.Position = [5 743 1795 283];



    %% 3 - polar plot - re-orientated so that the direction nearest to the direction of the vector sum points upwards - only 28 dps

    responses = max_v_polar(1:end-1)-median_voltage;  % peak neural responses
    
    % Compute vector components
    x_component = responses.*cos(angls);
    y_component = responses.*sin(angls);

    % Compute vector sum
    vector_sum_x = sum(x_component)';
    vector_sum_y = sum(y_component)';

    % Compute magnitude and direction
    magnitude = sqrt(vector_sum_x.^2 + vector_sum_y.^2);
    % This is the angle in which the vector sum points.
    angle_rad = atan2(vector_sum_y, vector_sum_x);

    % Find the closest direction in 'directions' to the vector sum angle. 
    if angle_rad < 0 
        angle_rad = angle_rad+2*pi;
    end

    diff_angl = abs(angls - angle_rad);
    loc_min = find(diff_angl == min(diff_angl));

    % Find order of directions to plot - if the 'peak' direction is in
    % subplot (position) 8. 
    for j = 1:16

        id = (j - loc_min)+8;
        if id<1
            id = id+16;
        elseif id>16
            id = id-16;
        end 
        order_idx(j, i) = id;
    end 

    target_direction = pi/2;  % Target direction (π/2 radians)
    current_peak_direction = angls(loc_min);
    rotation_offset = target_direction - current_peak_direction;
    
    aligned_directions = mod(angls + rotation_offset, 2*pi);
    d = [aligned_directions, responses];
    d = sortrows(d, 1);
    
    % if i ==1 
    %     f_polar = figure;
    % else
    %     figure(f_polar)
    % end 

    % PLOT 
    data_to_plot = vertcat(d, d(1, :)); % repeat the first value to close the circle. 
    data_to_plot(:,2) = data_to_plot(:, 2)/max(data_to_plot(:, 2)); % normalise max to 1 - 
    polarplot(data_to_plot(:, 1), data_to_plot(:, 2), 'Color', [1 0.8 0.8], 'LineWidth', 1);
    hold on
    rlim([0 1.1])

    if i == 1
        dd = d;
    else
        dd(:, i+1) = d(:, 2);
    end 

    %% 4 - reorder timeseries data, so that PD data is all in the same column. 

    for kk = 1:32
        if kk<17
            ord_id = order_idx(kk, i);
        else 
            ord_id = order_idx(kk-16, i)+16;
        end 
        data_aligned(ord_id, i) = data_ordered(kk, 1);
    end 


    %% 5 - Plot the individual timeseries - now with the PD plotted in the middle subplot

    % f2 = figure;
    % for i3 = 1:32
    % 
    %     subplot(2,16, i3)
    % 
    %     plot(data_aligned{i3, i}*10, 'k');
    %     n_xvals = numel(data_aligned{i3, i});
    %     hold on 
    %     plot([0 n_xvals], [median_voltage median_voltage], 'r', 'LineWidth', 0.3)
    % 
    %     box off
    %     ax = gca;
    %     % ax.XAxis.Visible = 'off';
    %     ax.TickLength = [0.04 0.04];
    %     ax.TickDir = 'out';
    % 
    %     if i3>16
    %         ylim([-80 -40])
    %         ax.XTickLabel = {'0', '0.5', '1'};
    %     else
    %         ylim([-80 -30])
    %         ax.XTickLabel = {'0', '1', '2'};
    %     end 
    % 
    %     if i3 ==1 || i3 == 17
    %         ylabel('Voltage (mV)');
    %         xlabel('Time (s)')
    %     else
    %         ax.YAxis.Visible = 'off';
    %     end 
    % 
    %     if i3 == 8 
    %         title('28 dps')
    %     elseif i3 == 24
    %         title('56 dps')
    %     end 
    % 
    % end 
    % 
    % f2.Position = [5 743 1795 283];

end 

%% Plot the mean onto the polarplot

% Plot the mean on top
ang = dd(:, 1);
mag = mean(dd(:, 2:end), 2);

ang = vertcat(ang, ang(1));
mag = vertcat(mag, mag(1));

mag = mag/max(mag);

% polarplot(ang, mag, 'Color', [1 0.8 0.8], 'LineWidth', 2);
% figure(f_polar)
polarplot(ang, mag, 'Color', 'r', 'LineWidth', 2);
title('RNAi ttl - ON')




%% FInd the mean across all of the cells - plot the timeseries - PEAKS ARE NOT ALIGNED. 

f3 = figure;
for i4 = 1:32

    subplot(2,16,i4)

    if i4<17
        data = nan(n_exp, 24000);
    else
        data = nan(n_exp, 12000);
    end 

    for jj = 1:n_exp
        d_cell = data_aligned{i4, jj}*10;
        plot(d_cell, 'Color', col);
        hold on
        n_cell = numel(d_cell);
        data(jj, 1:n_cell) = d_cell;
    end 

    mean_dir = nanmean(data);
    mean_dir = mean_dir(~isnan(mean_dir));

    plot(mean_dir, 'k', 'LineWidth', 2);
    n_xvals = numel(mean_dir);
    % hold on 
    % plot([0 n_xvals], [median_voltage median_voltage], 'r', 'LineWidth', 0.3)

    box off
    ax = gca;
    ax.TickLength = [0.04 0.04];
    ax.TickDir = 'out';

    if i4>16
        ylim([-75 -40])
        ax.XTickLabel = {'0', '0.5', '1'};
    else
        ylim([-75 -30])
        ax.XTickLabel = {'0', '1', '2'};
    end 

    if i4 ==1 || i4 == 17
        ylabel('Voltage (mV)');
        xlabel('Time (s)')
    else
        ax.YAxis.Visible = 'off';
    end 

    if i4 == 8 
        title('28 dps')
    elseif i4 == 24
        title('56 dps')
    end 

end 

f3.Position = [5 743 1795 283];


%% With the peaks of the data aligned. 


figure

for d = 1:32
    
    subplot(2, 16, d)

    arrays = data_aligned(d, :);
    num_arrays = length(arrays);
    
    % Find peak indices
    peak_indices = zeros(1, num_arrays);
    array_lengths = zeros(1, num_arrays);
    
    for q = 1:num_arrays
        [~, peak_indices(q)] = max(arrays{q});
        array_lengths(q) = length(arrays{q});
    end
    
    % Determine the maximum peak index
    max_peak = max(peak_indices);
    
    % Compute shifts needed
    shifts = max_peak - peak_indices;
    
    % Determine the new length (maximum shifted length)
    new_length = max(array_lengths + shifts);
    
    % Create shifted arrays
    aligned_arrays = nan(num_arrays, new_length);  % Use NaN padding
    
    for i2 = 1:num_arrays
        if i2<=6
            col = [0.8 0.8 0.8];
        else 
            col = [1 0.8 0.8];
        end 

        ddd = arrays{i2}; % data 
        x = [1:1:numel(ddd)]+shifts(i2);
        xq = 1:1:new_length;
        vq = interp1(x, ddd, xq, 'nearest', 'extrap'); % extrapolate data 
        aligned_arrays(i2, :) = vq;
        plot(aligned_arrays(i2, :)*10, 'Color', col)
        hold on
    end

    plot(nanmean(aligned_arrays(1:6, :))*10, 'Color', 'k', 'LineWidth', 1.5)
    plot(nanmean(aligned_arrays(7:12, :))*10, 'Color', 'r', 'LineWidth', 1.5)

    box off
    ax = gca;
    ax.TickLength = [0.04 0.04];
    ax.TickDir = 'out';

    if d>16
        ylim([-73 -45])
        ax.XTickLabel = {'0', '0.5', '1'};
    else
        ylim([-73 -35])
        ax.XTickLabel = {'0', '1', '2'};
    end 

    if d ==1 || d == 17
        ylabel('Voltage (mV)');
        xlabel('Time (s)')
    else
        ax.YAxis.Visible = 'off';
    end 

    if d == 8 
        title('28 dps')
    elseif d == 24
        title('56 dps')
    end 
end 

f = gcf;
f.Position = [18  688  1762 316];

%% Median voltage per recording

median_vs = zeros(1, n_exp);

for i = 1:n_exp

    % Load the data from one cell:
    load(file_list(i).name);
    median_vs(i) = median_voltage;
end 
%% Plot each individual cell 

% Only 28 dps at the moment

figure
for d = 1:16  % Run through the 16 directions. Columns. 

    arrays = data_aligned(d, :);
    num_arrays = length(arrays);
    
    
    % Find peak indices
    peak_indices = zeros(1, num_arrays);
    array_lengths = zeros(1, num_arrays);
    xv = cell(1, num_arrays);
    
    for q = 1:num_arrays
        [~, peak_indices(q)] = max(arrays{q});
        array_lengths(q) = length(arrays{q});
    end
    
    % Determine the maximum peak index
    max_peak = max(peak_indices);
    
    % Compute shifts needed
    shifts = max_peak - peak_indices;

    % Shift timing 
    for q2 = 1:num_arrays
        data_to_plot = arrays{q2}*10;
        maxval = max(data_to_plot)*0.95;
        minval = min(data_to_plot);
        n_dp = numel(data_to_plot);
        x_vals = 1:1:n_dp;
        x_vals_shifted = x_vals-peak_indices(q2);
        xv{q2} = x_vals_shifted;

        subplot(num_arrays, 16, 16*(q2-1) + d)
        % Plot median
        mv = nanmean(data_to_plot(1:3000));
        hold on

        plot([-20000 20000], [median_vs(q2), median_vs(q2)], '-', 'Color', [1 0.8 0.8], 'LineWidth', 1)
        rectangle('Position', [-20000 minval 40000 abs(median_vs(q2)-minval)], 'FaceColor', 'k', 'FaceAlpha', 0.2)

        % plot([-20000 20000], [mv, mv], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1) 
        % rectangle('Position', [-20000 minval 40000 abs(mv-minval)], 'FaceColor', 'k', 'FaceAlpha', 0.2)

        plot([0 0], [-90, -30], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1)
        % Plot the trace for the individ dir and indiv cell 
        plot(x_vals_shifted, data_to_plot, 'Color', 'k', 'LineWidth', 1.2); 
        
        ylim([minval maxval])
        xlim([-13000 15000])
        box off
        ax = gca;
        ax.TickLength = [0.04 0.04];
        ax.TickDir = 'out';
        ax.XAxis.Visible = 'off';

        if q2 ==1 
            if d == 1 || d == 16
                title('ND')
            elseif d == 8
                title('PD')
            end 
        end 
    end 

end 

f = gcf;
f.Position = [1  75  1800 972];
sgtitle('RNAi control - ON')

    % 
    % 
    % if d>16
    %     ylim([-75 -40])
    %     ax.XTickLabel = {'0', '0.5', '1'};
    % else
    %     ylim([-75 -30])
    %     ax.XTickLabel = {'0', '1', '2'};
    % end 
    % 
    % if d ==1 || d == 17
    %     ylabel('Voltage (mV)');
    %     xlabel('Time (s)')
    % else
    %     ax.YAxis.Visible = 'off';
    % end 
    % 
    % if d == 8 
    %     title('28 dps')
    % elseif d == 24
    %     title('56 dps')
    % end 
























































%% 3 - Plot polar plots - re-orientated so that the direction nearest to the direction of the vector sum points upwards - only 28 dps

dd = [];

col = [0.8 0.8 0.8];

order_idx = nan(16, n_exp);

figure 

for i = 1:n_exp

    d_timeseries = data(1:16, 4);
    
    %% Vector sum
    
    responses = max_v_polar(1:end-1)-median_voltage;  % peak neural responses
    directions = angls(1:end-1)'; 
    
    % Compute vector components
    x_component = responses.*cos(directions);
    y_component = responses.*sin(directions);

    % Compute vector sum
    vector_sum_x = sum(x_component)';
    vector_sum_y = sum(y_component)';

    % Compute magnitude and direction
    magnitude = sqrt(vector_sum_x.^2 + vector_sum_y.^2);
    % This is the angle in which the vector sum points.
    angle_rad = atan2(vector_sum_y, vector_sum_x);

    % Find the closest direction in 'directions' to the vector sum angle. 
    if angle_rad < 0 
        angle_rad = angle_rad+2*pi;
    end

    diff_angl = abs(directions - angle_rad);
    loc_min = find(diff_angl == min(diff_angl));

    for j = 1:16
        id = j - loc_min;
        if id<1
            id = id+16;
        end 
        order_idx(j, i) = id;
    end 

    target_direction = pi/2;  % Target direction (π/2 radians)
    current_peak_direction = directions(loc_min);
    rotation_offset = target_direction - current_peak_direction;
    
    aligned_directions = mod(directions + rotation_offset, 2*pi);
    d = [aligned_directions, responses];
    d = sortrows(d, 1);
    
    % PLOT 
    data_to_plot = vertcat(d, d(1, :)); % repeat the first value to close the circle. 
    polarplot(data_to_plot(:, 1), data_to_plot(:, 2), 'Color', col, 'LineWidth', 1);
    hold on

    if i == 1
        dd = d;
    else
        dd(:, i+1) = d(:, 2);
    end 
    
end 

% Plot the mean on top
ang = dd(:, 1);
mag = mean(dd(:, 2:end), 2);

ang = vertcat(ang, ang(1));
mag = vertcat(mag, mag(1));

% polarplot(ang, mag, 'Color', [1 0.8 0.8], 'LineWidth', 2);
polarplot(ang, mag, 'Color', 'k', 'LineWidth', 2);

sgtitle('control RNAi ON')

%% 









% Plot with peaks of each trace in the middle: 

figure

for i = 1:16
    
    subplot(1,16,i)

    arrays = d_time_all(i, :);
    num_arrays = length(arrays);
    
    % Find peak indices
    peak_indices = zeros(1, num_arrays);
    array_lengths = zeros(1, num_arrays);
    
    for q = 1:num_arrays
        [~, peak_indices(q)] = max(arrays{q});
        array_lengths(q) = length(arrays{q});
    end
    
    % Determine the maximum peak index
    max_peak = max(peak_indices);
    
    % Compute shifts needed
    shifts = max_peak - peak_indices;
    
    % Determine the new length (maximum shifted length)
    new_length = max(array_lengths + shifts);
    
    % Create shifted arrays
    aligned_arrays = nan(num_arrays, new_length);  % Use NaN padding
    
    for i2 = 1:num_arrays
        start_idx = shifts(i2) + 1;
        end_idx = start_idx + length(arrays{i2}) - 1;
        aligned_arrays(i2, start_idx:end_idx) = arrays{i2};
    end

    for j = 1:n_exp
        plot(aligned_arrays(j, :)*10, 'Color', [0.8 0.8 0.8], 'LineWidth',1)
        hold on
    end 
    plot(nanmean(aligned_arrays(:, 3000:24000))*10, 'Color', 'k', 'LineWidth', 2)
    box off
    ylim([-80 -30])

    if i == 8 
        title('PD')
    end 

end 

f = gcf;
f.Position = [18  829  1771  175];



%% Plot every time series on it's own plot.

figure
for j = 1:n_exp
    arrays = d_time_all(j, :);
    for i = 1:16
        subplot(n_exp,16, 16*(j-1)+i)
        plot(aligned_arrays(i, :)*10, 'Color', [0.8 0.8 0.8], 'LineWidth',1)
        hold on
    end 
end 





%% 














