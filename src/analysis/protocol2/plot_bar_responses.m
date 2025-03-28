%% Analysing responses to bar stimuli - Jin Yong ephys data from T4T5 cells
% Dec 2024. 

% This processes the responses to the bar stimuli and plots a polar plot
% with the timeseries around the polar plot for the 2 speeds. 
% It saves this figure and saves the data in a results folder. 

close all
clear

date_folder = cd;
date_str = date_folder(end-15:end-6);
time_str = date_folder(end-4:end);

strrs = split(date_folder, '/');
type_str = string(strrs(end-1));
strain_str = string(strrs{end-2});

strrs = split(date_folder, '/');
cell_type = strrs{end-1};
if strrs{end-2}=="control"
    strain = "CTL";
else
    strain = "TTL";
end 

fig_save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/figures/bar_polar';
results_save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results';

%% Go to experiment folder with data from protocol 2:

% Find the 'G4_TDMS_Log..mat' file and load it:
log_folder = fullfile(date_folder, "Log Files"); cd(log_folder);
log_file = dir('G4_TDMS*');
load(log_file.name, 'Log');

f_data = Log.ADC.Volts(1, :);
v_data = Log.ADC.Volts(2, :)*10;

%% % % % % % %  % % % %  % % CHECK DATA:
% Plot the frame position data and the raw voltage trace. 

% figure; subplot(3, 1, 1); plot(f_data);
% subplot(3,1,2:3); plot(v_data)
% f2 = gcf;
% f2.Position =  [59  651  1683 375]; 

%%
% Any date after 12/12/24 - only 2 speeds and 30 pixel square area. 
% 12/12/24 or before = 3 speeds and 45 pixel square area. 
disp(date_str)

% 1 - find the difference between frame positions:
diff_f_data = diff(f_data);
idx = find(diff_f_data == min(diff_f_data));

%% % % % % % %  % % % %  % % CHECK DATA:
% Plot the difference in frame data with the position of each time 
% the difference is the minimum (-392).

% figure; 
% subplot(2,1,1);
% plot(f_data); hold on
% for ii = 1:numel(idx)
%     plot([idx(ii), idx(ii)], [0 400], 'r')
% end 
% subplot(2,1,2);
% plot(diff_f_data); hold on
% for ii = 1:numel(idx)
%     plot([idx(ii), idx(ii)], [-400 400], 'r')
% end 

%% Find the duration of each single bar movement in : 

% After 12_12_24
dur_t = (2.273+1.155)*10000*16; % (dur_bar_slow + dur_bar_fast) * acq_speed * n_directions. 

% Find the timings of when the bar stimuli start and end: 
% 1120 is added because idx(2) etc is the end of the last flash and then we
% add 1120 because that is the gap between the last flash and the beginning
% of the bar stimuli. 
rep1_rng = [idx(2)+1120, idx(2)+dur_t];
rep2_rng = [idx(4)+1120, idx(4)+dur_t];
rep3_rng = [idx(6)+1120, numel(f_data)]; % til the end of the recording. 

% For experiments pre 12_12_24
% rep1_rng = [1474710, 2022040];
% rep2_rng = [3496300, 4043600];
% rep3_rng = [5517920, 6065030];

%% % % % % % %  % % % %  % % CHECK DATA:
% Plot vertical lines at the beginning and end of the bar stimuli. 
% Check that the frame position data matches. 

% figure; plot(f_data);
% hold on;
% plot([rep1_rng(1), rep1_rng(1)], [0 400], 'm');
% plot([rep1_rng(2), rep1_rng(2)], [0 400], 'm');
% plot([rep2_rng(1), rep2_rng(1)], [0 400], 'm');
% plot([rep2_rng(2), rep2_rng(2)], [0 400], 'm');
% plot([rep3_rng(1), rep3_rng(1)], [0 400], 'm');
% plot([rep3_rng(2), rep3_rng(2)], [0 400], 'm');

%% Find the median voltage across the entire recording:
median_voltage = median(v_data);

%% Combine the voltage timeseries data from the three repetitions into 
% one data structure. 

rep_ranges = {rep1_rng, rep2_rng, rep3_rng};
idxs_all = cell(1, 3); 

for i = 1:3
    st_val = rep_ranges{i}(1);
    end_val = rep_ranges{i}(2);
    
    frames_rep = f_data(st_val:end_val);
    
    r_min = islocalmin(frames_rep);
    r_min_vals = find(r_min);
    r_st = [st_val, r_min_vals + st_val, end_val];

    r_max = islocalmax(frames_rep);
    r_max_vals = find(r_max);
    r_nd = r_max_vals + st_val;
    
    all_idxs = [r_st, r_nd];
    idxs_all{i} = sort(all_idxs); % Store sorted indices
end

% Extract data segments using the computed indices
num_segments = numel(idxs_all{1}) - 1; % Assuming all idxs_all{i} have the same size
data = cell(num_segments, 3); % Preallocate data cell array

for i = 1:num_segments
    for j = 1:3
        data{i, j} = v_data(idxs_all{j}(i):idxs_all{j}(i+1)-1);
    end
end

%% % % % % % %  % % % %  % % CHECK DATA:
% Test the chopping up of the bar stimuli:

% figure; 
% plot(f_data, 'k') % Plot the data in black for better visibility
% hold on
% 
% % Define colors for different repetitions (optional)
% colors = {'r', 'g', 'b'}; % Red, Green, Blue for each repetition
% 
% % Loop through repetitions
% for rep = 1:3
%     idxs = idxs_all{rep}; % Get indices for current repetition
% 
%     % Plot vertical lines for segmentation
%     for i = 1:numel(idxs)
%         plot([idxs(i), idxs(i)], [min(f_data), max(f_data)], colors{rep}, 'LineWidth', 1.5)
%     end
% end
% 
% hold off
% xlabel('Frame Number')
% ylabel('Signal')
% title('Segment Boundaries for Different Repetitions')
% legend('f\_data', 'Rep 1 Boundaries', 'Rep 2 Boundaries', 'Rep 3 Boundaries')

%% Add the mean of the three repetitions as the 4th column:

for j = 1:height(data)
    % Extract data for all repetitions and ensure column vectors
    d = cellfun(@(x) x(:), data(j, 1:3), 'UniformOutput', false); 
    n_per_col = cellfun(@numel, d); % Get the number of elements in each repetition

    % Find minimum length and trim all series to match
    min_val = min(n_per_col);
    d_trimmed = cellfun(@(x) x(1:min_val), d, 'UniformOutput', false);

    % Convert to matrix and compute mean across rows (time series)
    d_matrix = horzcat(d_trimmed{:}); % Convert cell array to matrix (columns = repetitions)
    mean_resp = nanmean(d_matrix, 2); % Compute mean across repetitions (columns)

    % Store mean time series in the 4th column
    data{j, 4} = mean_resp';
end

%% GENERATE PLOT 
% Central polar plot with timeseries around in location 
% representative of bar direction. Responses to fast bars are in light 
% blue and responses to slow bars in dark blue. 

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

plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

angls = linspace(0, 2*pi, 17); % 17 points to include 0 and 2π

% Preallocate max/min voltage arrays
max_v = zeros(numPlots, 2);
min_v = zeros(numPlots, 2);

% Define colors for the two conditions
colors = {[0.2 0.4 0.7], [0.4 0.8 1]};  % Dark blue (28 dps) and Light blue (56 dps)

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
        d_idx = plot_order(i) + 16*(sp-1);

        % Plot data for each repetition
        for r = 1:4
            d2plot = data{d_idx, r};
            x_vals = 1:numel(d2plot);

            if r == 1 
                % Plot median voltage background
                plot([1 x_vals(end)], [median_voltage, median_voltage], 'Color', [0.7 0.7 0.7], 'LineWidth', 1);  
            end 

            % Plot repetitions in gray, last rep in condition color
            if r < 4
                plot(x_vals, d2plot, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
            else 
                plot(x_vals, d2plot, 'Color', col, 'LineWidth', 1.2);
            end
        end 

        % Set consistent Y-limits
        ylim([-80 -10])

        % Store max/min values from last column of data = mean over the 3
        % repetitions:
        d = data{d_idx, 4}; 
        n_vals_d = numel(d);
        max_v(i, sp) = prctile(d(floor(n_vals_d/8):end), 98); % ignore the first 8th of the condition
        min_v(i, sp) = prctile(d(floor(n_vals_d/2):end), 2); % find the min during the 2nd half of the condition.

        % Turn off axes for better visualization
        axis(ax, 'off');
    end
end

% Add polar plot in the center
axCentral = polaraxes('Position', centralPosition);
hold on

% Convert max values for both conditions into polar format
max_v_polar1 = vertcat(max_v(:, 1), max_v(1, 1)); % slow bars
max_v_polar2 = vertcat(max_v(:, 2), max_v(1, 2)); % fast bars

% Plot the polar data
polarplot(angls, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 2);
polarplot(angls, max_v_polar2 - median_voltage, 'Color', colors{2}, 'LineWidth', 2);

% Add title
sgtitle(sprintf("28 / 56 dps - 4 pixel bar stimuli - 30 pix square - %s - %s - %s - %s", ...
                strrep(date_str, '_', '-'), strrep(time_str, '_', '-'), strain, cell_type));

% Set figure size
set(gcf, 'Position', [303 78 961 969]);

fig_folder = fullfile(fig_save_folder, cell_type);
if ~isfolder(fig_folder)
    mkdir(fig_folder)
end 

% f = gcf; 
% fname = fullfile(fig_folder, strcat(strain, '_', cell_type, '_', date_str, '_', time_str,'_bar_polar.pdf'));
% exportgraphics(f ...
%         , fname ...
%         , 'ContentType', 'vector' ...
%         , 'BackgroundColor', 'none' ...
%         ); 


% Add an arrow on top in the direction of the vector sum

theta = angls;
rho = [max_v_polar1 - median_voltage]';

[aaa, resultant_angle] = vector_sum_polar(rho, theta);

resultant_magnitude = 30; 
col = [0.3 0.3 0.3];

figure;
polarplot(angls, max_v_polar1 - median_voltage, 'Color', colors{1}, 'LineWidth', 3);
ax = gca;
ax.LineWidth = 1.2;
ax.FontSize = 15;
ax.ThetaTickLabel = {};
hold on
add_arrow_to_polarplot(resultant_magnitude, resultant_angle, col)


% Save as PNG
exportgraphics(gca ...
        , "2024_12_18_15_07_polarplot_20dps_bar_DS_wArrow.png" ...
        );

exportgraphics(gca ...
        , "2024_12_18_15_07_polarplot_20dps_bar_DS_wArrow.pdf" ...
        , 'ContentType', 'vector' ...
        , 'BackgroundColor', 'none' ...
        ); 


%% Plot heat map of the max voltage reached during each rep. 

figure; 
imagesc(max_v); 

hcb = colorbar;
cm_inferno=inferno(1000);
colormap(cm_inferno)
ax_c= gca;
ax_c.TickDir = 'out';
ax_c.LineWidth = 1;
ax_c.FontSize = 12; 
box off

yticks([1,5,9, 13])
yticklabels({'0', '90', '180', '270'});
ylabel('Direction - deg')

xticks([1,2])
xticklabels({'28', '56'})
xlabel('Speed - dps')

colorTitleHandle = get(hcb,'Title');
set(colorTitleHandle ,'String','Max voltage (mV)');

f2 = gcf;
f2.Position = [620   386   190   581];

%% Align the data
data_ordered = align_data_by_seq_angles(data);

%% Find PR and the order to re-order the data to align PD to pi/2.
% max_v_polar = mean([max_v_polar1, max_v_polar2], 2);
% max_v_polar1_norm = max_v_polar1/max(max_v_polar1);
[d, ord, magnitude_slow, angle_rad_slow, fwhm_slow, cv_slow, thetahat_slow, kappa_slow] = find_PD_and_order_idx(max_v_polar1, median_voltage); % use the max v polar for the slower bars.

data_aligned = cell(size(data_ordered));
for kk = 1:32 
    if kk<17
        ord_id = ord(kk);
    else 
        ord_id = ord(kk-16)+16;
    end 

    data_aligned(ord_id, :) = data_ordered(kk, :); % Combine the mean timeseries only.
end 

[d_fast, ord_fast, magnitude_fast, angle_rad_fast, fwhm_fast, cv_fast, thetahat_fast, kappa_fast] = find_PD_and_order_idx(max_v_polar2, median_voltage);

%% Symmetry index:
% 'd' - column 1 = angle in radians.
% column 2 = max_v for that angle. Adjusted so that the peak is 90 deg.
% (1.5708 rad). 
% Row 5 = the row corresponding to 90 deg. 
% 4-6
% 3-7
% 2-8
% 1-9
% 16-10
% 15-11
% 14-12

vals1 = d([4,3,2,1,16,15,14], 2); % Increases in voltage from median per angle. 
vals2 = d([6,7,8,9,10,11,12], 2);
diff_vals = abs(vals1-vals2); % Find the diff voltage between the peak responses from the 2 halves.
sym_val_slow = sum(diff_vals)/sum(d(:, 2)); % Sum these differences and normalise by the sum of all of the responses.
sym_ratio_slow = 1 - sym_val_slow; % Closer to 1 = more symmetric.

% FAST - 56 dps
vals1 = d_fast([4,3,2,1,16,15,14], 2);
vals2 = d_fast([6,7,8,9,10,11,12], 2);
diff_vals = abs(vals1-vals2); 
sym_val_fast = sum(diff_vals)/sum(d_fast(:, 2));
sym_ratio_fast = 1 - sym_val_fast;

%% DSI - vector sum method

% SLOW 

% Compute the vector sum
vector_sum_slow = sum(d(:,2) .* exp(1i * d(:,1)));
% Compute the DSI
DSI_slow = abs(vector_sum_slow) / sum(d(:, 2));
% DSI: - PD - ND
DSI_pdnd_slow = (d(5,2)-d(13,2))/(d(5,2)+d(13,2));

% FAST 

% Compute the vector sum
vector_sum_fast = sum(d_fast(:,2) .* exp(1i * d_fast(:,1)));
% Compute the DSI
DSI_fast = abs(vector_sum_fast) / sum(d_fast(:, 2));
% DSI: - PD - ND
DSI_pdnd_fast = (d_fast(5,2)-d_fast(13,2))/(d_fast(5,2)+d_fast(13,2));

%% SAVE THE DATA

bar_results = table();
bar_results.Date = date_str;
bar_results.Time = time_str;
bar_results.Strain = strain_str;
bar_results.Type = type_str;
bar_results.max_v_polar1 = {max_v_polar1};
bar_results.max_v_polar2 = {max_v_polar2};
bar_results.min_v = {min_v};
bar_results.max_v = {max_v};
bar_results.median_voltage = median_voltage;

% output of vector sum:
bar_results.magnitude_slow = magnitude_slow;
bar_results.angle_rad_slow = angle_rad_slow;
bar_results.fwhm_slow = fwhm_slow;
bar_results.cv_slow = cv_slow;
bar_results.thetahat_slow = thetahat_slow;
bar_results.kappa_slow = kappa_slow;

bar_results.magnitude_fast = magnitude_fast;
bar_results.angle_rad_fast = angle_rad_fast;
bar_results.fwhm_fast = fwhm_fast;
bar_results.cv_fast = cv_fast;
bar_results.thetahat_fast = thetahat_fast;
bar_results.kappa_fast = kappa_fast;

% symmetry index
bar_results.sym_ratio_slow = sym_ratio_slow;
bar_results.sym_ratio_fast = sym_ratio_fast;

% DSI - vector sum 
bar_results.vector_sum_slow = vector_sum_slow;
bar_results.DSI_vector_slow = DSI_slow;
bar_results.DSI_pdnd_slow = DSI_pdnd_slow;
bar_results.vector_sum_fast = vector_sum_fast;
bar_results.DSI_vector_fast = DSI_fast;
bar_results.DSI_pdnd_fast = DSI_pdnd_fast;

% res_folder = strcat(results_save_folder, '/', cell_type);
res_folder = strcat(results_save_folder, '/bar_results');

if ~isfolder(res_folder)
    mkdir(res_folder)
end 

save(fullfile(res_folder, strcat('peak_vals_', strain, '_', cell_type, '_', date_str, '_', time_str, '.mat'))...
    , "bar_results"...
    , 'data' ...
    , 'data_ordered'...
    , 'data_aligned'...
    , 'ord'...
    , 'd'...
    );
