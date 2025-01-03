%% Generate plots of the receptive field of T4T5 cells: 

% 1 - Load the data: 

% Input experiment folder with data from protocol 2:
date_folder = cd;

% Find the 'G4_TDMS_Log..mat' file and load it:
log_folder = fullfile(date_folder, "Log Files"); cd(log_folder);
log_file = dir('G4_TDMS*');
load(log_file.name, 'Log');

% Load the frame and voltage data:
f_data = Log.ADC.Volts(1, :);

diff_f_data = diff(f_data);
idx = find(diff_f_data == min(diff_f_data));

% figure; plot(f_data);
% hold on;
% plot([idx(1), idx(1)], [0 200], 'm');
% plot([idx(2), idx(2)], [0 200], 'm');
% plot([idx(3), idx(3)], [0 200], 'm');
% plot([idx(1), idx(1)], [0 200], 'm');


v_data = Log.ADC.Volts(2, :)*10;
median_v = median(v_data);
v2_data = v_data - median_v; 

std_v = std(v2_data);

upper_bound_med = std_v;
lower_bound_med =  -std_v;

% TEST 
% figure; 
% plot(v_data);
% hold on;
% plot([1, numel(v_data)],[median_v, median_v], 'r');
% % mm = movmean(v_data, 10000*30);
% % plot(mm, 'c')
% plot([1, numel(v_data)],[median_v+std_v/2, median_v+std_v/2], 'g');
% plot([1, numel(v_data)],[median_v-std_v/2, median_v-std_v/2], 'g');

% pattern with 4 pix squares. 
pat_folder = fullfile(date_folder, "Patterns");
cd(pat_folder)
pat_file = dir("0001*");
load(pat_file.name, 'pattern');

% Function slow 
func_folder = fullfile(date_folder, "Functions");
cd(func_folder)
func_file = dir("0001*");
load(func_file.name, 'pfnparam');
pfnparam2 = pfnparam;
func_file2 = dir("0002*");
load(func_file2.name, 'pfnparam');

dur_slowflashes = pfnparam2.dur; % seconds. 
dur_fastflashes = pfnparam.dur; % seconds. 
sampling_rate = 10000;
% End of slow flashes - rep1 
% dur_ms = dur_slowflashes*sampling_rate; % 
dur_ms = 976700;
% End of fast flashes - rep1 
dur_ms_fast = dur_ms + dur_fastflashes*sampling_rate;

% TEST - check timing
% figure; plot(f_data);
% hold on;
% plot([dur_ms , dur_ms], [0 400], 'r', 'LineWidth', 1.2)
% plot([dur_ms_fast , dur_ms_fast], [0 400], 'm', 'LineWidth', 1.2)
% 
fc = pfnparam2.func; 
min_frame = min(fc);
max_frame = max(fc);

if max_frame < 198
    on_off = "off";
else 
    on_off = "on";
end 

%%  PLOT % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% fast_flashes_dur = 2500; % 0.25s * sampling rate. 
% edge_vals2 = dur_ms+1:fast_flashes_dur:dur_ms_fast;

% SLOW FLASHES: 

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

    max_val_flash = max(mean_data_flash);
    min_val_flash = min(mean_data_flash);
    
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

f = gcf;
f.Position = [77  173  1057  874]; %[77  76   1379   971];

figure; imagesc(data_comb2)
cmap = redblue();
colormap(cmap)
caxis([-0.15 1.15])


% Test
hold on 
plot([idx(3), idx(3)], [0 400], 'g')
plot([idx(5), idx(5)], [0 400], 'g')
plot([idx(3)-dur_ms, idx(3)-dur_ms], [0 400], 'm')
plot([idx(5)-dur_ms, idx(5)-dur_ms], [0 400], 'm')















%% OLD EXPERIMENTS - 20kHZ - 200ms flash - 50 ms wait



% 4 pixel flashes - 2 speeds. 

% f_data = Log.ADC.Volts(1, :);
% v_data = Log.ADC.Volts(2, :)*10;
% median_v = median(v_data);
% 
% allf = pattern.Pats;
% 
% dur_slowflashes = pfnparam.dur; % seconds. 
% sampling_rate = 20000;
% dur_ms = dur_slowflashes*sampling_rate; 
% 
% % rep1 
% 
% edge_vals = 6719:5000:dur_ms;
% % 
% figure; plot(f_data); hold on;
% for i = 1:floor(numel(edge_vals)/2)
%     plot([edge_vals(i) , edge_vals(i)], [0 400], 'g', 'LineWidth', 1)
% end 

% 196 flashes / values. 
% 
% figure
% 
% % Rep 1
% for i = 1:ceil(numel(edge_vals)/2)
% 
%     if i < ceil(numel(edge_vals)/2)
%         d = f_data(edge_vals(i):edge_vals(i+1)-1);
%         v = v_data(edge_vals(i):edge_vals(i+1)-1);
%     elseif i == ceil(numel(edge_vals)/2)
%         d = f_data(edge_vals(ceil(numel(edge_vals)/2)):edge_vals(ceil(numel(edge_vals)/2))+5000);
%         v = v_data(edge_vals(ceil(numel(edge_vals)/2)):edge_vals(ceil(numel(edge_vals)/2))+5000);
%     end 
% 
%     flash_frame_num = max(d)-1;
% 
%     rows = 14 - floor((flash_frame_num) / 14);
%     cols = mod((flash_frame_num), 14) + 1;
% 
%     X = (rows - 1) * 14 + cols;
% 
%     subplot(14, 14, X)
%     plot(v, 'Color', [0.4 0.4 0.4])
%     hold on 
%     plot([4500 4500], [-70 -40], 'r')
%     xmax = numel(v);
%     plot([1 xmax], [median_v, median_v], 'Color', [0.7 0.7 0.7])
%     ylim([-70 -40])
%     box off
%     axis off
% end 
% 
% %% Rep 2 
% 
% start_t = 2025520;
% edge_vals = start_t:5000:start_t+dur_ms;
% 
% % Test
% % for i = 1:197
% %     plot([edge_vals(i) , edge_vals(i)], [0 400], 'g', 'LineWidth', 1)
% % end 
% 
% for i = 1:196
% 
%     if i < 196
%         d = f_data(edge_vals(i):edge_vals(i+1)-1);
%         v = v_data(edge_vals(i):edge_vals(i+1)-1);
%     elseif i == 196
%         d = f_data(edge_vals(196):edge_vals(196)+5000);
%         v = v_data(edge_vals(196):edge_vals(196)+5000);
%     end 
% 
%     flash_frame_num = max(d)-1;
% 
%     rows = 14 - floor((flash_frame_num - 196) / 14);
%     cols = mod((flash_frame_num - 196), 14) + 1;
% 
%     X = (rows - 1) * 14 + cols;
% 
%     subplot(14, 14, X)
%     plot(v, 'Color', [0.4 0.4 0.4])
%     hold on 
%     plot([1600 1600], [-70 -40], 'r')
%     xmax = numel(v);
%     plot([1 xmax], [median_v, median_v], 'Color', [0.7 0.7 0.7])
%     ylim([-70 -40])
%     box off
%     axis off
% end 
% 
% 
% 
% %% Rep 3
% 
% start_t = 4047120;
% edge_vals = start_t:5000:start_t+dur_ms;
% 
% % % Test
% % for i = 1:197
% %     plot([edge_vals(i) , edge_vals(i)], [0 400], 'y', 'LineWidth', 1)
% % end 
% 
% for i = 1:196
% 
%     if i < 196
%         d = f_data(edge_vals(i):edge_vals(i+1)-1);
%         v = v_data(edge_vals(i):edge_vals(i+1)-1);
%     elseif i == 196
%         d = f_data(edge_vals(196):edge_vals(196)+5000);
%         v = v_data(edge_vals(196):edge_vals(196)+5000);
%     end 
% 
%     flash_frame_num = max(d)-1;
% 
%     rows = 14 - floor((flash_frame_num - 196) / 14);
%     cols = mod((flash_frame_num - 196), 14) + 1;
% 
%     X = (rows - 1) * 14 + cols;
% 
%     subplot(14, 14, X)
%     plot(v, 'Color', [0.2 0.2 0.2])
%     hold on 
%     % plot([1600 1600], [-70 -40], 'r')
%     xmax = numel(v);
%     plot([1 xmax], [median_v, median_v], 'Color', [0.7 0.7 0.7])
%     ylim([-70 -40])
%     box off
%     axis off
% end 
% 





% % REP 1 - 
% slow_flashes_dur = 5000; % 500 ms. 
% fast_flashes_dur = 2500; % 250 ms.
% 
% edge_vals = 3865:slow_flashes_dur:dur_ms;
% edge_vals2 = dur_ms+1:fast_flashes_dur:dur_ms_fast;
% 
% % 196 flashes / values. 
% data_comb = ones(14, 14);
% 
% % v_data_ord = reshape(permute(data_comb,[2 1]), 1, 196);
% % min_v_ord = min(v_data_ord);
% % max_v_ord = max(v_data_ord);
% % v_ord_norm = (v_data_ord- min_v_ord) / (max_v_ord - min_v_ord);
% 
% figure
% 
% % for idx = 1:196
% %     subplot(14, 14, idx)
% %     c = v_ord_norm(idx);
% %     rectangle('Position', [0 -70 5000, 40], "FaceColor", [c c c], "EdgeColor", 'none', "FaceAlpha", 0.5);
% %     axis off
% %     box off
% %     axis square
% %     hold on
% % end 
% 
% % Rep 1
% for i = 1:196
% 
%     if i < 196
%         d = f_data(edge_vals(i):edge_vals(i+1)-1);
%         v = v_data(edge_vals(i):edge_vals(i+1)-1);
%     elseif i == 196
%         d = f_data(edge_vals(196):edge_vals(196)+5000);
%         v = v_data(edge_vals(196):edge_vals(196)+5000);
%     end 
% 
%     flash_frame_num = max(d)-1;
% 
%     rows = 14 - floor((flash_frame_num - 196) / 14);
%     cols = mod((flash_frame_num - 196), 14) + 1;
% 
%     X = (rows - 1) * 14 + cols;
% 
%     subplot(14, 14, X)
%     plot(v, 'Color', 'k')
%     hold on 
%     % plot([1600 1600], [-70 -40], 'k')
%     xmax = numel(v);
%     plot([1 xmax], [median_v, median_v], 'Color', [1 0.7 0.7])
%     ylim([-70 -30])
%     box off
%     axis off
% 
%     data_comb(rows, cols) = data_comb(rows, cols) + mean(v); 
% end 
% 
% 
% %% Rep 2 
% 
% start_t = 2025520;
% edge_vals = start_t:5000:start_t+dur_ms;
% 
% % Test
% % for i = 1:197
% %     plot([edge_vals(i) , edge_vals(i)], [0 400], 'g', 'LineWidth', 1)
% % end 
% 
% for i = 1:196
% 
%     if i < 196
%         d = f_data(edge_vals(i):edge_vals(i+1)-1);
%         v = v_data(edge_vals(i):edge_vals(i+1)-1);
%     elseif i == 196
%         d = f_data(edge_vals(196):edge_vals(196)+5000);
%         v = v_data(edge_vals(196):edge_vals(196)+5000);
%     end 
% 
%     flash_frame_num = max(d)-1;
% 
%     rows = 14 - floor((flash_frame_num - 196) / 14);
%     cols = mod((flash_frame_num - 196), 14) + 1;
% 
%     X = (rows - 1) * 14 + cols;
% 
%     subplot(14, 14, X)
%     plot(v, 'Color', 'k')
%     hold on 
%     % plot([1600 1600], [-70 -40], 'r')
%     xmax = numel(v);
%     % plot([1 xmax], [median_v, median_v], 'Color', [1 0.7 0.7])
%     ylim([-70 -30])
%     box off
%     axis off
%     data_comb(rows, cols) = data_comb(rows, cols) + mean(v); 
% end 
% 
% 
% 
% %% Rep 3
% 
% start_t = 4047120;
% edge_vals = start_t:5000:start_t+dur_ms;
% 
% % % Test
% % for i = 1:197
% %     plot([edge_vals(i) , edge_vals(i)], [0 400], 'y', 'LineWidth', 1)
% % end 
% 
% for i = 1:196
% 
%     if i < 196
%         d = f_data(edge_vals(i):edge_vals(i+1)-1);
%         v = v_data(edge_vals(i):edge_vals(i+1)-1);
%     elseif i == 196
%         d = f_data(edge_vals(196):edge_vals(196)+5000);
%         v = v_data(edge_vals(196):edge_vals(196)+5000);
%     end 
% 
%     flash_frame_num = max(d)-1;
% 
%     rows = 14 - floor((flash_frame_num - 196) / 14);
%     cols = mod((flash_frame_num - 196), 14) + 1;
%     % 
%     X = (rows - 1) * 14 + cols;
% 
%     subplot(14, 14, X)
%     plot(v, 'Color', 'k')
%     hold on 
%     % plot([1600 1600], [-70 -40], 'r')
%     xmax = numel(v);
%     % plot([1 xmax], [median_v, median_v], 'Color', [1 0.7 0.7])
%     ylim([-70 -30])
%     box off
%     axis off
% 
%     data_comb(rows, cols) = data_comb(rows, cols) + mean(v); 
% end 
% 
% 
% data_comb = data_comb / 3;
% %%
% 
% % % Initialize the new 7x7 non-overlapping grid
% % new_grid = zeros(7, 7);
% % 
% % % Loop through the 14x14 array in blocks of 2x2
% % for i = 1:7
% %     for j = 1:7
% %         % Compute the mean of the corresponding 2x2 block
% %         new_grid(i, j) = mean(mean(data_comb(2*i-1:2*i, 2*j-1:2*j)));
% %     end
% % end
% % 
% % % Display the resulting grid
% % figure; imagesc(new_grid); title('mean')
% 
% figure; imagesc(data_comb); title('mean')








% Plot the traces from the 3 different response types: 

figure
for i = 1:196

    data_flash = ones(3, slow_flashes_dur); 

    for r = 1:3 

        if r == 1
            start_t = 3865; 
            edge_vals = start_t:slow_flashes_dur:dur_ms;
        elseif r == 2
            start_t = 2025520;
            edge_vals = start_t:slow_flashes_dur:start_t+dur_ms;
        elseif r == 3
            start_t = 4047120;
            edge_vals = start_t:slow_flashes_dur:start_t+dur_ms;
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
    max_val_flash = max(mean_data_flash);
    min_val_flash = min(mean_data_flash);

    if max_val_flash > upper_bound_med % excitation - red
        val = max_val_flash; 
        cm = 1;
    elseif min_val_flash < lower_bound_med % inhibition - blue 
        val = min_val_flash;
        cm = 2;
    elseif max_val_flash <= upper_bound_med && min_val_flash >= lower_bound_med  % neither. 
        val = mean(mean_data_flash);
        cm = 3;
    else 
        disp('error')
        val = mean(mean_data_flash);
        cm = 3;
    end 

    subplot(1, 3, cm)
    hold on;
    plot(mean_data_flash)

end 


 %% Create figure with coloured squares in the background: 

for idx = 1:196
    subplot(14, 14, idx)
    c = v_ord_norm(idx);
    rectangle('Position', [0 -70 5000, 40], "FaceColor", [c c c], "EdgeColor", 'none', "FaceAlpha", 0.5);
    axis off
    box off
    axis square
    hold on
end 

% Then plot the traces on top: 











