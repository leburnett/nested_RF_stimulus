
% Input experiment folder with data from protocol 2:
date_folder = '/Users/burnettl/Documents/Projects/nested_RF_stimulus/protocol2/data/2024_12_12_12_12';

% Find the 'G4_TDMS_Log..mat' file and load it:
log_folder = fullfile(date_folder, "Log Files"); cd(log_folder);
log_file = dir('G4_TDMS*');
load(log_file.name, 'Log');

% Load the frame and voltage data:
f_data = Log.ADC.Volts(1, :);
v_data = Log.ADC.Volts(2, :)*10;
median_v = median(v_data);

% pattern with 4 pix squares. 
pat_folder = fullfile(date_folder, "Patterns");
cd(pat_folder)
pat_file = dir("0001*");
load(pat_file.name, 'pattern');

% Function slow 
func_folder = fullfile(date_folder, "Functions");
cd(pat_folder)
func_file = dir("0001*");
load(func_file.name, 'pfnparam');

dur_slowflashes = pfnparam.dur; % seconds. 
sampling_rate = 10000;
dur_ms = dur_slowflashes*sampling_rate; 

% 
% figure; plot(f_data);
% hold on;
% plot([dur_ms , dur_ms], [0 400], 'r', 'LineWidth', 1.2)
% 
% fc = pfnparam.func; 
% 
% % 340 ms OFF - bkg - 160 ms FLASH. 
% % 500 ms = 0.5s. 
% dur_cond = 0.5*10000;
% 
% plot([dur_cond , dur_cond], [0 400], 'r', 'LineWidth', 1)
% 

%%  PLOT % % % % % % % % % % % % % % % % % % % % % % % % 


%% UPDATE TO RUN THROUGH ALL THREE REPS FOR EACH SUBPLOT AT THE SAME TIME. 

% FIRST - run through to get overall data and find mean. 

% THEN - determine whether E / I - use different time windows to find min / max then make background colour plot. 

% THEN - run through a second time and plot the mean traces on top. 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% REP 1 - 

edge_vals = 3865:5000:dur_ms;

% 196 flashes / values. 
data_comb = ones(14, 14);

% v_data_ord = reshape(data_comb, 1, 196);
v_data_ord = reshape(permute(data_comb,[2 1]), 1, 196);
min_v_ord = min(v_data_ord);
max_v_ord = max(v_data_ord);
v_ord_norm = (v_data_ord- min_v_ord) / (max_v_ord - min_v_ord);

figure

for idx = 1:196
    subplot(14, 14, idx)
    c = v_ord_norm(idx);
    rectangle('Position', [0 -70 5000, 40], "FaceColor", [c c c], "EdgeColor", 'none', "FaceAlpha", 0.5);
    axis off
    box off
    axis square
    hold on
end 

% Rep 1
for i = 1:196

    if i < 196
        d = f_data(edge_vals(i):edge_vals(i+1)-1);
        v = v_data(edge_vals(i):edge_vals(i+1)-1);
    elseif i == 196
        d = f_data(edge_vals(196):edge_vals(196)+5000);
        v = v_data(edge_vals(196):edge_vals(196)+5000);
    end 

    flash_frame_num = max(d)-1;

    rows = 14 - floor((flash_frame_num - 196) / 14);
    cols = mod((flash_frame_num - 196), 14) + 1;

    X = (rows - 1) * 14 + cols;

    subplot(14, 14, X)
    plot(v, 'Color', 'k')
    hold on 
    % plot([1600 1600], [-70 -40], 'k')
    xmax = numel(v);
    plot([1 xmax], [median_v, median_v], 'Color', [1 0.7 0.7])
    ylim([-70 -30])
    box off
    axis off

    % data_comb(rows, cols) = data_comb(rows, cols) + mean(v); 
end 

%% Rep 2 

start_t = 2025520;
edge_vals = start_t:5000:start_t+dur_ms;

% Test
% for i = 1:197
%     plot([edge_vals(i) , edge_vals(i)], [0 400], 'g', 'LineWidth', 1)
% end 

for i = 1:196

    if i < 196
        d = f_data(edge_vals(i):edge_vals(i+1)-1);
        v = v_data(edge_vals(i):edge_vals(i+1)-1);
    elseif i == 196
        d = f_data(edge_vals(196):edge_vals(196)+5000);
        v = v_data(edge_vals(196):edge_vals(196)+5000);
    end 

    flash_frame_num = max(d)-1;

    rows = 14 - floor((flash_frame_num - 196) / 14);
    cols = mod((flash_frame_num - 196), 14) + 1;

    X = (rows - 1) * 14 + cols;

    subplot(14, 14, X)
    plot(v, 'Color', 'k')
    hold on 
    % plot([1600 1600], [-70 -40], 'r')
    xmax = numel(v);
    % plot([1 xmax], [median_v, median_v], 'Color', [1 0.7 0.7])
    ylim([-70 -30])
    box off
    axis off
    % data_comb(rows, cols) = data_comb(rows, cols) + mean(v); 
end 



%% Rep 3

start_t = 4047120;
edge_vals = start_t:5000:start_t+dur_ms;

% % Test
% for i = 1:197
%     plot([edge_vals(i) , edge_vals(i)], [0 400], 'y', 'LineWidth', 1)
% end 

for i = 1:196

    if i < 196
        d = f_data(edge_vals(i):edge_vals(i+1)-1);
        v = v_data(edge_vals(i):edge_vals(i+1)-1);
    elseif i == 196
        d = f_data(edge_vals(196):edge_vals(196)+5000);
        v = v_data(edge_vals(196):edge_vals(196)+5000);
    end 

    flash_frame_num = max(d)-1;

    rows = 14 - floor((flash_frame_num - 196) / 14);
    cols = mod((flash_frame_num - 196), 14) + 1;
    % 
    X = (rows - 1) * 14 + cols;

    subplot(14, 14, X)
    plot(v, 'Color', 'k')
    hold on 
    % plot([1600 1600], [-70 -40], 'r')
    xmax = numel(v);
    % plot([1 xmax], [median_v, median_v], 'Color', [1 0.7 0.7])
    ylim([-70 -30])
    box off
    axis off

    % data_comb(rows, cols) = data_comb(rows, cols) + mean(v); 
end 


data_comb = data_comb / 3;
%%

% Initialize the new 7x7 non-overlapping grid
new_grid = zeros(7, 7);

% Loop through the 14x14 array in blocks of 2x2
for i = 1:7
    for j = 1:7
        % Compute the mean of the corresponding 2x2 block
        new_grid(i, j) = mean(mean(data_comb(2*i-1:2*i, 2*j-1:2*j)));
    end
end

% Display the resulting grid
figure; imagesc(new_grid); title('mean')

figure; imagesc(data_comb); title('mean')

caxis([-68 -51])











%% OLD EXPERIMENTS - 20kHZ - 200ms flash - 50 ms wait



% 4 pixel flashes - 2 speeds. 

f_data = Log.ADC.Volts(1, :);
v_data = Log.ADC.Volts(2, :)*10;
median_v = median(v_data);

allf = pattern.Pats;

dur_slowflashes = pfnparam.dur; % seconds. 
sampling_rate = 20000;
dur_ms = dur_slowflashes*sampling_rate; 

% rep1 

edge_vals = 6719:5000:dur_ms;
% 
% figure; plot(f_data); hold on;
% for i = 1:floor(numel(edge_vals)/2)
%     plot([edge_vals(i) , edge_vals(i)], [0 400], 'g', 'LineWidth', 1)
% end 

% 196 flashes / values. 

figure

% Rep 1
for i = 1:ceil(numel(edge_vals)/2)

    if i < ceil(numel(edge_vals)/2)
        d = f_data(edge_vals(i):edge_vals(i+1)-1);
        v = v_data(edge_vals(i):edge_vals(i+1)-1);
    elseif i == ceil(numel(edge_vals)/2)
        d = f_data(edge_vals(ceil(numel(edge_vals)/2)):edge_vals(ceil(numel(edge_vals)/2))+5000);
        v = v_data(edge_vals(ceil(numel(edge_vals)/2)):edge_vals(ceil(numel(edge_vals)/2))+5000);
    end 

    flash_frame_num = max(d)-1;

    rows = 14 - floor((flash_frame_num) / 14);
    cols = mod((flash_frame_num), 14) + 1;

    X = (rows - 1) * 14 + cols;

    subplot(14, 14, X)
    plot(v, 'Color', [0.4 0.4 0.4])
    hold on 
    plot([4500 4500], [-70 -40], 'r')
    xmax = numel(v);
    plot([1 xmax], [median_v, median_v], 'Color', [0.7 0.7 0.7])
    ylim([-70 -40])
    box off
    axis off
end 

%% Rep 2 

start_t = 2025520;
edge_vals = start_t:5000:start_t+dur_ms;

% Test
% for i = 1:197
%     plot([edge_vals(i) , edge_vals(i)], [0 400], 'g', 'LineWidth', 1)
% end 

for i = 1:196

    if i < 196
        d = f_data(edge_vals(i):edge_vals(i+1)-1);
        v = v_data(edge_vals(i):edge_vals(i+1)-1);
    elseif i == 196
        d = f_data(edge_vals(196):edge_vals(196)+5000);
        v = v_data(edge_vals(196):edge_vals(196)+5000);
    end 

    flash_frame_num = max(d)-1;

    rows = 14 - floor((flash_frame_num - 196) / 14);
    cols = mod((flash_frame_num - 196), 14) + 1;

    X = (rows - 1) * 14 + cols;

    subplot(14, 14, X)
    plot(v, 'Color', [0.4 0.4 0.4])
    hold on 
    plot([1600 1600], [-70 -40], 'r')
    xmax = numel(v);
    plot([1 xmax], [median_v, median_v], 'Color', [0.7 0.7 0.7])
    ylim([-70 -40])
    box off
    axis off
end 



%% Rep 3

start_t = 4047120;
edge_vals = start_t:5000:start_t+dur_ms;

% % Test
% for i = 1:197
%     plot([edge_vals(i) , edge_vals(i)], [0 400], 'y', 'LineWidth', 1)
% end 

for i = 1:196

    if i < 196
        d = f_data(edge_vals(i):edge_vals(i+1)-1);
        v = v_data(edge_vals(i):edge_vals(i+1)-1);
    elseif i == 196
        d = f_data(edge_vals(196):edge_vals(196)+5000);
        v = v_data(edge_vals(196):edge_vals(196)+5000);
    end 

    flash_frame_num = max(d)-1;

    rows = 14 - floor((flash_frame_num - 196) / 14);
    cols = mod((flash_frame_num - 196), 14) + 1;

    X = (rows - 1) * 14 + cols;

    subplot(14, 14, X)
    plot(v, 'Color', [0.2 0.2 0.2])
    hold on 
    % plot([1600 1600], [-70 -40], 'r')
    xmax = numel(v);
    plot([1 xmax], [median_v, median_v], 'Color', [0.7 0.7 0.7])
    ylim([-70 -40])
    box off
    axis off
end 






