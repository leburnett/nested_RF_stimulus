%% Generate plots of the receptive field of T4T5 cells: 
% First generate a receptive field estimate for the slow flashes and then
% for the fast flashes. 
close all
clear

% 1 - Load the data: 

% Input experiment folder with data from protocol 2:
date_folder = cd;

% Find the 'G4_TDMS_Log..mat' file and load it:
log_folder = fullfile(date_folder, "Log Files"); cd(log_folder);
log_file = dir('G4_TDMS*');
load(log_file.name, 'Log');

sampling_rate = 10000;

% % % %  Load the frame data:
f_data = Log.ADC.Volts(1, :);

diff_f_data = diff(f_data);
idx = find(diff_f_data == min(diff_f_data));

% figure; plot(f_data);
% hold on;
% plot([idx(1), idx(1)], [0 200], 'm');
% plot([idx(2), idx(2)], [0 200], 'm');
% plot([idx(3), idx(3)], [0 200], 'm');
% plot([idx(1), idx(1)], [0 200], 'm');

% % % %  Load the voltage data:
v_data = Log.ADC.Volts(2, :)*10;
median_v = median(v_data);
v2_data = v_data - median_v; 

% Find the standard deviation and the upper / lower bounds for finding
% response groups later on. 
std_v = std(v2_data);
upper_bound_med = std_v;
lower_bound_med = -std_v/2;

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

% Function slow 
func_folder = fullfile(date_folder, "Functions");
cd(func_folder)
func_file = dir("0001*");
load(func_file.name, 'pfnparam');
pfnparam2 = pfnparam;
dur_slowflashes = pfnparam2.dur; % seconds. 
% dur_ms = dur_slowflashes*sampling_rate; % 
dur_ms = 976700;

% TEST - check timing
% figure; plot(f_data);
% hold on;
% plot([dur_ms , dur_ms], [0 400], 'r', 'LineWidth', 1.2)
% plot([dur_ms_fast , dur_ms_fast], [0 400], 'm', 'LineWidth', 1.2)
% 

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

%%  SLOW FLASHES

% % 340 ms OFF - bkg - 160 ms FLASH. 
% % 500 ms = 0.5s. 
slow_flashes_dur = 5000; % 0.5s * sampling rate. 

% 196 flashes / values. 
data_comb = zeros(14, 14);
cmap_id = zeros(14, 14);

for i = 1:196

    data_flash = ones(3, slow_flashes_dur); 

    for r = 1:3

        if r == 1
            end_t = idx(1); %3865; 
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        elseif r == 2
            end_t = idx(3);%start_t = 2025440; %2025520;
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        elseif r == 3
            end_t = idx(5); %start_t = 4047040; %4047120;
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        end 

        if i < 196
            d = f_data(edge_vals(i):edge_vals(i+1)-1);
            v = v2_data(edge_vals(i):edge_vals(i+1)-1);
        elseif i == 196
            d = f_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
            v = v2_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
        end 

        data_flash(r, :) = v;
    end 
    
    % In the future - add a check for outlier reps. 

    mean_data_flash = mean(data_flash);

    max_val_flash = max(mean_data_flash(350:2500));
    min_val_flash = min(mean_data_flash(1250:4250));
    
    % choose colour of background. 

    if max_val_flash > upper_bound_med % excitation - red
        val = mean(mean_data_flash(350:2500));
        cm = 1;
    elseif min_val_flash < lower_bound_med % inhibition - blue 
        val = mean(mean_data_flash(1250:4250));
        cm = 2;
    elseif max_val_flash <= upper_bound_med && min_val_flash >= lower_bound_med  % neither. 
        val = mean(mean_data_flash);
        cm = 3;
    else 
        disp('error')
        val = mean(mean_data_flash);
        cm = 3;
    end 

    flash_frame_num = max(d)-1;
    fnum(i) = flash_frame_num;
    % disp(flash_frame_num)

    if on_off == "on" % from 196
        rows = 14 - floor((flash_frame_num - 196) / 14);
        cols = mod((flash_frame_num - 196), 14) + 1;
    elseif on_off == "off" % 1- 196
        rows = 14 - floor(flash_frame_num/14);
        cols = mod(flash_frame_num, 14) + 1;
    end 

    data_comb(rows, cols) = val;
    cmap_id(rows, cols) = cm;

end 

data_comb2 = rescale(data_comb, 0, 1);
% figure; imagesc(data_comb); title('mean')
% caxis([-68 -51])

figure
for i = 1:196

    data_flash = ones(3, slow_flashes_dur); 

    for r = 1:3 

         if r == 1
            end_t = idx(1); %3865; 
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        elseif r == 2
            end_t = idx(3);%start_t = 2025440; %2025520;
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        elseif r == 3
            end_t = idx(5); %start_t = 4047040; %4047120;
            edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;
        end 

        if i < 196
            d = f_data(edge_vals(i):edge_vals(i+1)-1);
            v = v2_data(edge_vals(i):edge_vals(i+1)-1);
        elseif i == 196
            d = f_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
            v = v2_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
        end 

        data_flash(r, :) = v;
    end 

    mean_data_flash = mean(data_flash);

    flash_frame_num = max(d)-1;

    if on_off == "on" % from 196
        rows = 14 - floor((flash_frame_num - 196) / 14);
        cols = mod((flash_frame_num - 196), 14) + 1;
    elseif on_off == "off" % 1- 196
        rows = 14 - floor(flash_frame_num/14);
        cols = mod(flash_frame_num, 14) + 1;
    end 

    val  = data_comb2(rows, cols);
    cm = cmap_id(rows, cols);

    X = (rows - 1) * 14 + cols;
    subplot(14, 14, X)

    if cm == 1 % RED 
        c = [1, val, val];
    elseif cm == 2 % blue
        c = [0, 0, 1-val];
    elseif cm == 3 % grey 
        c = [1-val, 1-val, 1-val];
    end 

    rectangle('Position', [0 -25 5000, 50], "FaceColor", c, "EdgeColor", 'none', "FaceAlpha", 0.5);
    hold on
    plot(mean_data_flash, 'Color', 'k', 'LineWidth', 2)
    hold on 
    xmax = numel(v);
    plot([1 xmax], [0, 0], 'Color', [0.7 0.7 0.7])
    ylim([-10 25])
    axis off
    box off
    axis square

end 

sgtitle('160ms flashes - 340ms interval')
f = gcf;
f.Position = [77  173  1057  874]; %[77  76   1379   971];


%% Heatmap without traces: 

figure; imagesc(data_comb2)
cmap = redblue();
colormap(cmap)
clim([-0.15 1.15])


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













