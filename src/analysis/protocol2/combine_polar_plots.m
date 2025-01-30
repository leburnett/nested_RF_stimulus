
% Load in the data from the 'results' folder
data_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/control_RNAi/ON';
cd(data_folder)
strrs = split(data_folder, '/');

if strrs{end-1} == "control_RNAi"
    strain = "control";
else 
    strain = "ttl";
end 

cell_type = strrs{end};

file_list = dir('*_peak_vals*');
n_exp = height(file_list);

% Initialise empty arrays for filling in:
dd = [];
data_aligned = cell(32, n_exp);
order_idx = nan(16, n_exp);

%% Align data and generate aligned polar plot for all cells of the type and strain:

figure

for i = 1:n_exp

    % Load the data from one cell:
    load(file_list(i).name);

    %% 1 - rearrange data to be in sequential order for each angle. 
    % Before, it was in the order in which the stimuli were presented. Same
    % orientations, in the opposite directions sequentially. Whereas now, they
    % move around a circle in 1/16*pi sections. 
   
    % data_ordered = align_data_by_seq_angles(data);

    %% 2 - plot the timeseries responses of the individual cell to the 16 different directions
    % Plots the mean of 3 reps per each direction on an individual subplot.
    % slow - 28 dps = top row. 
    % fast - 56 dps - bottom row. 
    
    % plot_linear_timeseries_16D(data_ordered)

    %% 3 - polar plot - re-orientated so that the direction nearest to the direction of the vector sum points upwards - only 28 dps

    % [d, ord] = find_PD_and_order_idx(max_v_polar, median_voltage);

    % Add this cell's order to the combined array.
    order_idx(:, i) = ord;

    % Repeat the first value to close the circle.
    data_to_plot = vertcat(d, d(1, :));  

    % Normalise responses so that max = 1.
    data_to_plot(:,2) = data_to_plot(:, 2)/max(data_to_plot(:, 2)); 

    % Plot the polar plot for this cell - data aligned with PD at 90
    % degrees (up). 

    if strain == "control"
        col = [0.8 0.8 0.8];
    elseif strain == "ttl"
        col = [1 0.8 0.8];
    end 

    polarplot(data_to_plot(:, 1), data_to_plot(:, 2), 'Color', col, 'LineWidth', 1);
    hold on
    rlim([0 1.1])

    % Combine the peak response arrays across cells.
    % All cells have their PD to pi/2 here.
    if i == 1
        dd = d;
    else
        dd(:, i+1) = d(:, 2); % second column contains the peak voltage responses.
    end 

    %% 4 - reorder the timeseries data, so that PD data is all in the same row. 
    % 
    % for kk = 1:32 % for each direction
    % 
    %     if kk<17
    %         ord_id = ord(kk);
    %     else 
    %         ord_id = ord(kk-16)+16;
    %     end 
    % 
    %     data_aligned(ord_id, i) = data_ordered(kk, 4); % Combine the mean timeseries only.
    % end 

end 

%% Plot the mean onto the polarplot

% Plot the mean on top
ang = dd(:, 1);
mag = mean(dd(:, 2:end), 2);

ang = vertcat(ang, ang(1));
mag = vertcat(mag, mag(1));

% Normalise max to 1.
mag = mag/max(mag);

if strain == "control"
    col = 'k';
elseif strain == "ttl"
    col = "r";
end 

polarplot(ang, mag, 'Color', col, 'LineWidth', 2);
title(strcat(strain, '-', cell_type))


% %% Median voltage per recording
% 
% median_vs = zeros(1, n_exp);
% 
% for i = 1:n_exp
% 
%     % Load the data from one cell:
%     load(file_list(i).name);
%     median_vs(i) = median_voltage;
% end 


