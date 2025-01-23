
% Load in the data from the 'results' folder
data_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/control_RNAi';

file_list = dir('peak_vals*');
n_exp = height(file_list);

% ang_vals = zeros(n_exp, 17);
% vals = zeros(n_exp, 17);
dd = [];

col = [0.8 0.8 0.8]; 

%% 

figure 

for i = 1:n_exp

    load(file_list(i).name);
   
    % Find the max voltage for each direction to the 28 dps stimulus. 
    % max_v_polar = vertcat(max_v(:, 1), max_v(1, 1));
    
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
    angle_rad = atan2(vector_sum_y, vector_sum_x);

    [~, peak_index] = max(responses);

    target_direction = 1.5708;  % Target direction (π/2 radians)
    current_peak_direction = directions(peak_index);
    rotation_offset = target_direction - current_peak_direction;
    
    aligned_directions = mod(directions + rotation_offset, 2*pi);
    d = [aligned_directions, responses];
    d = sortrows(d, 1);
    
    % PLOT 
  

    % subplot(1,3,1)
    % polarplot(angls, max_v_polar-median_voltage, 'Color', col, 'LineWidth', 1.5);
    % hold on
    % 
    % subplot(1,3,2)
    % polarplot([0 angle_rad], [0 magnitude], 'r-')
    % hold on
    % 
    % % % Find how much to shift the angles to in order to be facing vertically.
    % % shift_angl = 90*pi/180 - angle_rad;
    % % % Shift the current angles by this much:
    % % new_angls = angls + shift_angl;
    % 
    % subplot(1,3,3)
    % data_to_plot = vertcat(d, d(1, :));
    % polarplot(data_to_plot(:, 1), data_to_plot(:, 2), 'Color', col, 'LineWidth', 1);
    % hold on
    % polarplot(new_angls, max_v_polar-median_voltage, 'Color', col, 'LineWidth', 1);
    % hold on

    % ang_vals(i, :) = d;
    % vals(i, :) = max_v_polar-median_voltage;
    if i == 1
        dd = d;
    else
        dd(:, i+1) = d(:, 2);
    end 
    
end 

f = gcf;
f.Position = [11 642  1778   390];

% Plot the mean on top

ang = dd(:, 1);
mag = mean(dd(:, 2:end), 2);

ang = vertcat(ang, ang(1));
mag = vertcat(mag, mag(1));

subplot(1,3,3)
polarplot(ang, mag, 'Color', 'r', 'LineWidth', 2);

sgtitle('control RNAi')


figure;
polarplot(ang, mag, 'Color', 'r', 'LineWidth', 2);