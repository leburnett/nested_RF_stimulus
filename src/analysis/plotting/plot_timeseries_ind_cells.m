function plot_timeseries_ind_cells(data_folder)

[data_align_all, ctrl_flies] = combine_data_find_ctrls(data_folder);

[median_vs] = combine_median_v(data_folder);

%% Plot linear timeseries plot - each cell individually 

figure
for d = 1:16  % Run through the 16 directions. 28 dps only. Columns. 

    arrays = data_align_all(d, :);
    num_arrays = length(arrays); % Number of cells.
    
    % Find peak indices
    peak_indices = zeros(1, num_arrays);
    array_lengths = zeros(1, num_arrays);
    xv = cell(1, num_arrays);
    
    for q = 1:num_arrays
        [~, peak_indices(q)] = max(arrays{q});
        array_lengths(q) = length(arrays{q});
    end
    
    % Determine the maximum peak index
    % max_peak = max(peak_indices);
    
    % Compute shifts needed
    % shifts = max_peak - peak_indices;

    % Shift timing 
    for q2 = 1:num_arrays

        if ctrl_flies(q2)==0
            col1 = [0.8 0.8 0.8];
            col2 = 'k';
        else 
            col1 = [1 0.8 0.8];
            col2 = 'r';
        end 

        data_to_plot = arrays{q2};
        maxval = prctile(max(data_to_plot), 95);
        minval = min(data_to_plot);
        n_dp = numel(data_to_plot);
        x_vals = 1:1:n_dp;
        x_vals_shifted = x_vals-peak_indices(q2);
        xv{q2} = x_vals_shifted;

        subplot(num_arrays, 16, 16*(q2-1) + d)
        % Plot median
        mv = nanmean(data_to_plot(1:3000));
        hold on

        plot([-20000 20000], [median_vs(q2), median_vs(q2)], '-', 'Color', col1, 'LineWidth', 1)
        rectangle('Position', [-20000 minval 40000 abs(median_vs(q2)-minval)], 'FaceColor', col2, 'FaceAlpha', 0.2)

        % Plot a vertical line through 0. 
        plot([0 0], [-90, -30], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1)
        % Plot the trace for the individ dir and indiv cell 
        plot(x_vals_shifted, data_to_plot, 'Color', col2, 'LineWidth', 1.2); 
        
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









%% Plot timeseries data vertically:
% figure
% offset = 0;
% celll = 7;
% for pp = 1:16
%     da = data_aligned{pp, celll}*10;
%     plot(da-offset, 'LineWidth', 0.2);
%     hold on
%     offset = offset+5;
% end 
% box off
% title(string(celll))
% 
% f = gcf;
% f.Position = [1   543   254   504];
% 




