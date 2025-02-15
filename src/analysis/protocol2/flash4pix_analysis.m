%% Generate plots of the receptive field of T4T5 cells: 
% First generate a receptive field estimate for the slow flashes and then
% for the fast flashes. 
close all
clear

% 1 - Load the data: 

% Input experiment folder with data from protocol 2:
date_folder = cd;
strrs = split(date_folder, '/');
date_str = strrs{end};

date_str2 = date_folder(end-15:end-6);
time_str2 = date_folder(end-4:end);

type_str = string(strrs{end-1});
strain_str = string(strrs{end-2});

% Find the 'G4_TDMS_Log..mat' file and load it:
log_folder = fullfile(date_folder, "Log Files"); cd(log_folder);
log_file = dir('G4_TDMS*');
load(log_file.name, 'Log');

sampling_rate = 10000;

% % % %  Load the frame data:
f_data = Log.ADC.Volts(1, :);

diff_f_data = diff(f_data);
% Find where the flash stimuli end. 
idx = find(diff_f_data == min(diff_f_data));

%% % % % % % % % % % % % % % % % % Check data:
% figure; plot(f_data);
% hold on;
% plot([idx(1), idx(1)], [0 200], 'm');
% plot([idx(2), idx(2)], [0 200], 'm');
% plot([idx(3), idx(3)], [0 200], 'm');
% plot([idx(1), idx(1)], [0 200], 'm');

%% % % %  Load the voltage data:
v_data = Log.ADC.Volts(2, :)*10;  % transform the data in mV
median_v = median(v_data); % Find the median voltage across the entire recording.

%% % % % % % % % % % % % % % % % % Check data:
qual_fig_folder= '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/figures/RF_estimate/quality';
% 
figure; 
subplot(5,1,1)
plot(f_data)
ax = gca;
ax.XAxis.Visible = 'off';
title('f-data')
xlim([0 numel(f_data)])
subplot(5,1,2:3)
plot(v_data)
xlim([0 numel(v_data)])
ylabel('v-data')
ax = gca;
ax.XAxis.Visible = 'off';
title(strcat("var - v-data = ", string(var(v_data))))
subplot(5,1,4:5)
plot(movmean(v_data, 20000))
xlim([0 numel(v_data)])
ylabel('movmean(v-data, 20000)')
hold on 
% plot(movmean(v_data, 40000))
% plot(movmean(v_data, 60000))
plot(movmean(v_data, 100000))
title(strcat("var - movmean - 100000 = ", string(var(movmean(v_data, 100000)))))
f = gcf;
f.Position = [25 629 1121 392];
% 
% savefig(gcf, fullfile(qual_fig_folder, strcat('Full_rec_', date_str, '_', strain_str, '_', type_str)))
% 
filtered_voltage_data = movmean(v_data, 20000);
var_filtered_v = var(filtered_voltage_data);

%%
v2_data = v_data - median_v; % Get the median-subtracted voltage.

% Find the standard deviation and the upper / lower bounds for finding
% response groups later on. 
std_v = std(v2_data);
upper_bound_med = std_v;
lower_bound_med = -std_v/2;

%% % % % % % % % % % % % % % % Check data:
% TEST - plot with median and upper/lower bounds plotted for understanding
% how groups are found. 

% % figure; 
% plot(v_data);
% hold on;
% plot([1, numel(v_data)],[median_v, median_v], 'r');
% % mm = movmean(v_data, 10000*30);
% % plot(mm, 'c')
% plot([1, numel(v_data)],[median_v+upper_bound_med, median_v+upper_bound_med], 'm');
% plot([1, numel(v_data)],[median_v+lower_bound_med, median_v+lower_bound_med], 'm');

%% Function for slow flashes 
func_folder = fullfile(date_folder, "Functions");
cd(func_folder)
func_file = dir("0001*");
load(func_file.name, 'pfnparam');
pfnparam2 = pfnparam;
dur_slowflashes = pfnparam2.dur; % seconds. 
% dur_ms = dur_slowflashes*sampling_rate; % 
dur_ms = 976700;

%% % % % % % % % % % % % % % % Check data:
% TEST - check timing
% figure; plot(f_data);
% hold on;
% plot([dur_ms , dur_ms], [0 400], 'r', 'LineWidth', 1.2)
% % plot([dur_ms_fast , dur_ms_fast], [0 400], 'm', 'LineWidth', 1.2)
% 

%%
% 'fc' is the function - each number corresponds to the frame that is being
% presented. There is one data point every 2ms. 

fc = pfnparam2.func; 

% Find the max frame number to determine if dark or light squares were
% presented. 
max_frame = max(fc);
if max_frame < 198
    on_off = "off";
else 
    on_off = "on";
end 

%% RECEPTIVE FIELD PLOT WITH TIMESERIES 
% Plot the timeseries responses to each flash at its corresponding position - colour the background of the plot
% Look at the difference between the peak and the trough in the timeseries.

% - Check for the variance across repetitions - if the variance is high - set as grey. 

%%  SLOW FLASHES

% % 340 ms OFF - bkg - 160 ms FLASH. 500 ms = 0.5s. 

slow_flashes_dur = 5000; % 0.5s * sampling rate. 

% Empty arrays to fill - 196 flashes/values. 
data_comb = zeros(14, 14);
cmap_id = zeros(14, 14);
var_across_reps = zeros(14, 14);
var_within_reps = zeros(14,14);
diff_mean = zeros(14,14);
max_data = zeros(14,14);
min_data = zeros(14,14);
i_num = zeros(14,14);

% Run through the responses to each flash.
for i = 1:196

    % Array to collect the data for each flash over the three repetitions.
    data_flash = ones(3, slow_flashes_dur); 
    data_flash_raw = ones(3, slow_flashes_dur); 

    for r = 1:3

        if r == 1 % rep 1
            end_t = idx(1); %3865; 
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        elseif r == 2 % rep 2 
            end_t = idx(3);%start_t = 2025440; %2025520;
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        elseif r == 3 % rep3 
            end_t = idx(5); %start_t = 4047040; %4047120;
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        end 

        if i < 196
            d = f_data(edge_vals(i):edge_vals(i+1)-1); % frame data during flash. 
            v = v2_data(edge_vals(i):edge_vals(i+1)-1); % median-subtracted voltage data
            v_raw = v_data(edge_vals(i):edge_vals(i+1)-1);
        elseif i == 196
            d = f_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
            v = v2_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
            v_raw = v_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
        end 

        data_flash(r, :) = v;
        data_flash_raw(r, :) = v_raw;
    end 
    
    % Check for the variance across reps:
    % var_X_rep = var(var(data_flash));

    % Standard deviation within each repetition (across time points)
    std_within_trial = std(data_flash_raw, 0, 2);  % std along rows (across time points)
    
    % Mean within each repetition (across time points)
    mean_within_trial = mean(data_flash_raw, 2);  % mean along rows (across time points)
    
    % Coefficient of Variation (CV) within each repetition
    cv_within_trial = std_within_trial ./ mean_within_trial;
    
    % Mean CV across all repetitions (mean within trial CV)
    mean_cv_within_trial = abs(mean(cv_within_trial));

    % Standard deviation across time points (across trials)
    std_across_trials = std(data_flash_raw, 0, 1);  % std along columns (across trials)
    
    % Mean across time points (across trials)
    mean_across_trials = mean(data_flash_raw, 1);  % mean along columns (across trials)
    
    % Coefficient of Variation (CV) across repetitions (across trials)
    cv_across_trials = std_across_trials ./ mean_across_trials;

    % max_cv_across_trials = prctile(cv_across_trials, 98);
    mean_cv_across_trials = abs(mean(cv_across_trials));

    mean_data_flash = mean(data_flash);
    n_vals = size(mean_data_flash, 2);

    % Max and min of the mean flash response:
    max_val_flash = prctile(mean_data_flash(500:end), 98); % Max 
    % min_val_flash = prctile(mean_data_flash(1250:end), 2); % Min in the second half - ignore if min is early. 
    min_val_flash = prctile(mean_data_flash(2500:end), 2); % Min in the second half - ignore if min is early. 

    diff_resp = max_val_flash - min_val_flash; 
    diff_med = median_v - min_val_flash;

    if abs(max_val_flash)>=abs(min_val_flash) % larger excitatory peak
        
        if diff_resp>3
                val = max_val_flash;
                cm = 1;
        else 
            val = mean(mean_data_flash(n_vals*0.75:end));
            cm = 3;
        end 
    elseif abs(max_val_flash)<abs(min_val_flash) % larger inhibitory peak 
        if diff_resp>2.8
            val = min_val_flash;
            cm = 2;
        else 
            val = mean(mean_data_flash(n_vals*0.75:end));
            cm = 3;
        end 
    end 

    flash_frame_num = max(d)-1;
    % fnum(i) = flash_frame_num;
    % disp(flash_frame_num)

    if on_off == "on" % from 196
        rows = 14 - mod((flash_frame_num - 196), 14);   % Rows decrease from 14 to 1
        cols = floor((flash_frame_num - 196) / 14) + 1; % Columns increase normally
    elseif on_off == "off" % 1- 196
        rows = 14 - mod(flash_frame_num, 14);   % Rows decrease from 14 to 1
        cols = floor(flash_frame_num / 14) + 1; % Columns increase normally
    end

    data_comb(rows, cols) = val;
    cmap_id(rows, cols) = cm;
    var_across_reps(rows, cols) = mean_cv_across_trials;
    var_within_reps(rows, cols) = mean_cv_within_trial;
    diff_mean(rows, cols) = diff_resp;
    max_data(rows, cols) = max_val_flash;
    min_data(rows, cols) = min_val_flash;
    i_num(rows, cols) = i;
end 

%% Generate a plot with some 'quality' check plots:

% With median value across all flashes in title.
figure; 
subplot(2,3,1) % 1 - variance across reps - consistency
imagesc(var_across_reps); colorbar
med_var_X_reps = median(reshape(var_across_reps, [1, 196]));
title(strcat("CoV across reps - ", string(med_var_X_reps)))

subplot(2,3,4)
histogram(var_across_reps);
xlabel('Coeff. of var. across reps - 98% val')

subplot(2,3,2) % 2 - variance within reps - strength of response
imagesc(var_within_reps); colorbar
med_var_W_reps = var(reshape(var_within_reps, [1, 196]));
title(strcat("var within reps - ", string(med_var_W_reps)))

subplot(2,3,5)
histogram(var_within_reps);
xlabel('Var. within reps - mean')

subplot(2,3,3) % 3 - difference between max and min per flash for mean.
imagesc(diff_mean); colorbar
med_diff_mean = median(reshape(diff_mean, [1, 196]));
title(strcat("diff mean - ", string(med_diff_mean)))

subplot(2,3,6)
histogram(diff_mean);
xlabel('Diff between max and min of mean')

f = gcf;
f.Position = [1   684   827   363];
savefig(gcf, fullfile(qual_fig_folder, strcat('Var_', date_str, '_', strain_str, '_', type_str)))

figure; 
subplot(2,1,1)
imagesc(max_data); colorbar; title('max - 98th prctile')
subplot(2,1,2)
histogram(max_data)
f=gcf;
f.Position = [620   501   275   466];
savefig(gcf, fullfile(qual_fig_folder, strcat('Max_', date_str, '_', strain_str, '_', type_str)))

figure; 
subplot(2,1,1)
imagesc(min_data*-1); colorbar; title('min - 98th prctile')
subplot(2,1,2)
histogram(min_data)
f=gcf;
f.Position = [620   501   275   466];
savefig(gcf, fullfile(qual_fig_folder, strcat('Min_', date_str, '_', strain_str, '_', type_str)))

%% Plot the timeseries responses to the flashes at each location. 

data_comb2 = rescale(data_comb, 0, 1);

% f = plot_rf_estimate_timeseries(data_comb2, cmap_id, f_data, v2_data, slow_flashes_dur, idx, dur_ms, on_off);
f = plot_rf_estimate_timeseries_line(data_comb2, cmap_id, f_data, v2_data, slow_flashes_dur, idx, dur_ms, on_off);

fig_save_folder1 = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/figures/RF_estimate/timeseries_line';
% fname = fullfile(fig_save_folder1, strcat(date_str, '_', strain_str, '_', type_str,  ".pdf"));
% exportgraphics(f ...
%         , fname ...
%         , 'ContentType', 'vector' ...
%         , 'BackgroundColor', 'none' ...
%         ); 

%% Heatmap without traces: 

figure; imagesc(data_comb2)
cmap = redblue();
colormap(cmap)

med_val = median(data_comb2(:));
clim([med_val-0.5 med_val+0.5])
colorbar
axis square

% fig_save_folder2 = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/figures/RF_estimate/heatmap';
% savefig(gcf, fullfile(fig_save_folder2, strcat(date_str, '_', strain_str, '_', type_str)))

%% Gaussian fits:
exc_data = data_comb2;

inh_data = data_comb;
inh_data(cmap_id~=2)=0;
% inh_data = sign(inh_data) .* log(1 + abs(inh_data)); 
[optEx, R_squared, optInh, R_squaredi, f1, f2] = gaussian_RF_estimate(exc_data, inh_data);

% Save the figures:

% fig_save_folder4 = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/figures/RF_estimate/gaussfits';
% fname2 = fullfile(fig_save_folder4, strcat(date_str, '_', strain_str, '_', type_str, ".pdf"));
% exportgraphics(f2 ...
%         , fname2 ...
%         , 'ContentType', 'vector' ...
%         , 'BackgroundColor', 'none' ...
%         ); 

%% Combine results into a table. 
rf_results = table();
rf_results.Date = date_str2;
rf_results.Time = time_str2;
rf_results.Strain = strain_str;
rf_results.Type = type_str;
rf_results.data_comb = {data_comb};
rf_results.max_data = {max_data};
rf_results.min_data = {min_data};
rf_results.diff_mean = {diff_mean};
rf_results.cmap_id = {cmap_id};
rf_results.max_val = prctile(reshape(max_data, [1, 196]), 98);
rf_results.min_val = prctile(reshape(min_data, [1, 196]), 2);
rf_results.var_within_reps = {var_within_reps};
rf_results.var_across_reps = {var_across_reps};
rf_results.var_filtered_v = var_filtered_v;
rf_results.med_var_X_reps = med_var_X_reps;
rf_results.med_var_W_reps = med_var_W_reps;
rf_results.med_diff_mean = med_diff_mean;
rf_results.R_squared = R_squared;
rf_results.sigma_x_exc = optEx(4);
rf_results.sigma_y_exc = optEx(5);
rf_results.optExc = {optEx};
rf_results.R_squaredi = R_squaredi;
rf_results.optInh = {optInh};
rf_results.sigma_x_inh = optInh(4);
rf_results.sigma_y_inh = optInh(5);

% Save data:
save_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/results/rf_results';
save(fullfile(save_folder, strcat('rf_results_', date_str,'_', strain_str, '_', type_str, '.mat')), 'rf_results');

% , "data_comb"...
    % , "cmap_id"...
    % , "var_within_reps"...
    % , "var_across_reps"...
    % , "min_data"...
    % , "max_data"...
    % , "diff_mean"...
    % , "R_squaredi"...
    % , "optInh"...
    % , "R_squared"...
    % , "optEx"...
    % )

% % % % % Test
% figure; plot(f_data)
% hold on 
% plot([idx(3), idx(3)], [0 400], 'g')
% plot([idx(5), idx(5)], [0 400], 'g')
% plot([idx(3)-dur_ms, idx(3)-dur_ms], [0 400], 'm')
% plot([idx(5)-dur_ms, idx(5)-dur_ms], [0 400], 'm')
% 


%% FAST FLASHES
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% FAST FLASHES
% fast_flashes_dur = 2500; % 0.25s * sampling rate. 
% 
% % 196 flashes / values. 
% data_comb = zeros(14, 14);
% cmap_id = zeros(14, 14);
% fnum = [];
% 
% for i = 1:196
% 
%     data_flash = ones(3, fast_flashes_dur); 
% 
%     for r = 1:3
% 
%         if r == 1
%             edge_vals = idx(1)+4860:fast_flashes_dur:idx(2);
%         elseif r == 2
%             edge_vals = idx(3)+4860:fast_flashes_dur:idx(4);
%         elseif r == 3
%             edge_vals = idx(5)+4860:fast_flashes_dur:idx(6);
%         end 
% 
%         if i < 196
%             d = f_data(edge_vals(i):edge_vals(i+1)-1);
%             v = v2_data(edge_vals(i):edge_vals(i+1)-1);
%         elseif i == 196
%             d = f_data(edge_vals(196):edge_vals(196)+fast_flashes_dur-1);
%             v = v2_data(edge_vals(196):edge_vals(196)+fast_flashes_dur-1);
%         end 
% 
%         data_flash(r, :) = v;
%     end 
% 
%     % In the future - add a check for outlier reps. 
% 
%     mean_data_flash = mean(data_flash);
% 
%     max_val_flash = max(mean_data_flash(500:2000));
%     min_val_flash = min(mean_data_flash(1250:2500));
% 
%     % choose colour of background. 
% 
%     if max_val_flash > upper_bound_med % excitation - red
%         % % % % % % different for fast and slow 
%         val = mean(mean_data_flash(400:2000));
%         cm = 1;
%     elseif min_val_flash < lower_bound_med % inhibition - blue 
%         val = mean(mean_data_flash(1250:2500));
%         cm = 2;
%     elseif max_val_flash <= upper_bound_med && min_val_flash >= lower_bound_med  % neither. 
%         val = mean(mean_data_flash);
%         cm = 3;
%     else 
%         disp('error')
%         val = mean(mean_data_flash);
%         cm = 3;
%     end 
% 
%     flash_frame_num = max(d)-1;
%     fnum(i) = flash_frame_num;
% 
%     if on_off == "on" % from 196
%         rows = 14 - floor((flash_frame_num - 196) / 14);
%         cols = mod((flash_frame_num - 196), 14) + 1;
%     elseif on_off == "off" % 1- 196
%         rows = 14 - floor(flash_frame_num/14);
%         cols = mod(flash_frame_num, 14) + 1;
%     end 
% 
%     data_comb(rows, cols) = val;
%     cmap_id(rows, cols) = cm;
% 
% end 
% 
% data_comb2 = rescale(data_comb, 0, 1);
% % figure; imagesc(data_comb); title('mean')
% % caxis([-68 -51])
% 
% 
% % Plot figure:
% figure
% for i = 1:196
% 
%     data_flash = ones(3, fast_flashes_dur); 
% 
%     for r = 1:3 
% 
%          if r == 1
%             edge_vals = idx(1)+4860:fast_flashes_dur:idx(2);
%         elseif r == 2
%             edge_vals = idx(3)+4860:fast_flashes_dur:idx(4);
%         elseif r == 3
%             edge_vals = idx(5)+4860:fast_flashes_dur:idx(6);
%         end 
% 
%         if i < 196
%             d = f_data(edge_vals(i):edge_vals(i+1)-1);
%             v = v2_data(edge_vals(i):edge_vals(i+1)-1);
%         elseif i == 196
%             d = f_data(edge_vals(196):edge_vals(196)+fast_flashes_dur-1);
%             v = v2_data(edge_vals(196):edge_vals(196)+fast_flashes_dur-1);
%         end 
% 
%         data_flash(r, :) = v;
%     end 
% 
%     mean_data_flash = mean(data_flash);
% 
%     flash_frame_num = max(d)-1;
% 
%     if on_off == "on" % from 196
%         rows = 14 - floor((flash_frame_num - 196) / 14);
%         cols = mod((flash_frame_num - 196), 14) + 1;
%     elseif on_off == "off" % 1- 196
%         rows = 14 - floor(flash_frame_num/14);
%         cols = mod(flash_frame_num, 14) + 1;
%     end 
% 
%     val  = data_comb2(rows, cols);
%     cm = cmap_id(rows, cols);
% 
%     X = (rows - 1) * 14 + cols;
%     subplot(14, 14, X)
% 
%     if cm == 1 % RED 
%         c = [1, val, val];
%     elseif cm == 2 % blue
%         c = [0, 0, 1-val];
%     elseif cm == 3 % grey 
%         c = [1-val, 1-val, 1-val];
%     end 
% 
%     rectangle('Position', [0 -25 fast_flashes_dur, 50], "FaceColor", c, "EdgeColor", 'none', "FaceAlpha", 0.5);
%     hold on
%     plot(mean_data_flash, 'Color', 'k', 'LineWidth', 2)
%     hold on 
%     xmax = numel(v);
%     plot([1 xmax], [0, 0], 'Color', [0.7 0.7 0.7])
%     ylim([-10 25])
%     axis off
%     box off
%     axis square
% 
% end 
% 
% sgtitle('80ms flashes - 170ms interval')
% f = gcf;
% f.Position = [77  173  1057  874]; %[77  76   1379   971];
% 
% 
% figure; imagesc(data_comb2)
% cmap = redblue();
% colormap(cmap)
% clim([-0.15 1.15])

% % % % Test
% figure; plot(f_data)
% hold on 
% plot([idx(1), idx(1)], [0 400], 'g')
% plot([idx(3), idx(3)], [0 400], 'g')
% plot([idx(5), idx(5)], [0 400], 'g')
% plot([idx(2), idx(2)], [0 400], 'g')
% plot([idx(4), idx(4)], [0 400], 'g')
% plot([idx(6), idx(6)], [0 400], 'g')


% % % % % Important timing test - how the flashes are bieng chopped up. 
% figure;
% plot(f_data)
% hold on 
% for i = 1:numel(edge_vals)
%     plot([edge_vals(i), edge_vals(i)], [0 200], 'c');
% end 
% 














%% Plot the time series responses to individual flashes, grouped by which group they are in. 

% TEST - - - - - Plot the traces from the 3 different response types: 

% figure
% for i = 1:196
% 
%     data_flash = ones(3, slow_flashes_dur); 
% 
%     for r = 1:3 
% 
%         if r == 1
%             start_t = 3865; 
%             edge_vals = start_t:slow_flashes_dur:dur_ms;
%         elseif r == 2
%             start_t = 2025520;
%             edge_vals = start_t:slow_flashes_dur:start_t+dur_ms;
%         elseif r == 3
%             start_t = 4047120;
%             edge_vals = start_t:slow_flashes_dur:start_t+dur_ms;
%         end 
% 
%         if i < 196
%             d = f_data(edge_vals(i):edge_vals(i+1)-1);
%             v = v2_data(edge_vals(i):edge_vals(i+1)-1);
%         elseif i == 196
%             d = f_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
%             v = v2_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
%         end 
% 
%         data_flash(r, :) = v;
%     end 
% 
%     mean_data_flash = mean(data_flash);
%     max_val_flash = max(mean_data_flash);
%     min_val_flash = min(mean_data_flash);
% 
%     if max_val_flash > upper_bound_med % excitation - red
%         val = max_val_flash; 
%         cm = 1;
%     elseif min_val_flash < lower_bound_med % inhibition - blue 
%         val = min_val_flash;
%         cm = 2;
%     elseif max_val_flash <= upper_bound_med && min_val_flash >= lower_bound_med  % neither. 
%         val = mean(mean_data_flash);
%         cm = 3;
%     else 
%         disp('error')
%         val = mean(mean_data_flash);
%         cm = 3;
%     end 
% 
%     subplot(1, 3, cm)
%     hold on;
%     plot(mean_data_flash)
% 
% end 









