function plot_timeseries_init0_alignedpeaks(data_folder)

% 'comb_ind' - set to 'comb' if data_folder has files from both control
% and turtle cells. Set to 'ind' if it's only one.

%% Plot timeseries data of voltage responses to bar stimuli.
% With the peaks of the data aligned AND the initial voltage subtracted so 
% that they all start at 0mV.

% data_folder = "/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/OFF_all"
cd(data_folder)
file_list = dir('*peak_vals*');
n_exp = height(file_list);

% Initialise empty arrays for filling in:
data_align_all = cell(32, n_exp);
ctrl_flies = zeros(n_exp, 1);

%% Combine the data from all of the cells of this type:

for i = 1:n_exp

    % Load needed variables:
    load(file_list(i).name);

    % Determine if the data is from a control or test fly.
    if contains(file_list(i).name, 'CTL')
        ctrl_flies(i) = 1;
    end 

    % Combine the "mean" aligned timeseries from all of the flies.  
    data_align_all(:, i) = data_aligned(:, 4); % Combine the mean timeseries only.
end 

% Convert 'ctrl_flies' from double to logical for logical indexing. 
ctrl_flies = logical(ctrl_flies);

%% Plot both control and turtle flies together. 
% Individual subplots for each direction. 
% Faint lines are individual cells and the bold lines are the mean across
% the cells. 

% Peaks are aligned in time. 
% Voltage traces have the initial value subtracted so that all of the
% traces start from zero. 
% Timeseries are filled with NaNs to make the arrays all the same length
% and then extrapolated. 

figure
for d = 1:32
    
    subplot(2, 16, d)

    arrays = data_align_all(d, :);
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

        if ctrl_flies(i2)==0
            col = [0.8 0.8 0.8];
        else 
            col = [1 0.8 0.8];
        end 

        ddd = arrays{i2}; % data 
        x = [1:1:numel(ddd)]+shifts(i2);
        xq = 1:1:new_length;
        vq = interp1(x, ddd, xq, 'nearest', 'extrap'); % extrapolate data 
        vq = vq-vq(1);
        aligned_arrays(i2, :) = vq;
        plot(aligned_arrays(i2, :), 'Color', col)
        hold on
    end

    plot(nanmean(aligned_arrays(ctrl_flies, :)), 'Color', 'k', 'LineWidth', 1.5)
    plot(nanmean(aligned_arrays(~ctrl_flies, :)), 'Color', 'r', 'LineWidth', 1.5)

    box off
    ax = gca;
    ax.TickLength = [0.04 0.04];
    ax.TickDir = 'out';
    xlim([0 new_length])

    if d>16
        % ylim([-73 -45])
        ax.XTickLabel = {'0', '0.5', '1'};
    else
        % ylim([-73 -35])
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

    ylim([-10 35])
end 

f = gcf;
f.Position = [18  688  1762 316];
