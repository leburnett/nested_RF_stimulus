function f = plot_bar_flash_data(data, meanData, median_v)
% 'data' contains the timeseries for each individual flash, organised per
% rep. 
% 'meanData' is the mean across the reps.

% Find background colour:
max_vals = zeros([8, 11]);
for i = 1:88
    d = meanData{i};
    n_points = numel(d);
    max_d = prctile(d(ceil(n_points*0.5):ceil(n_points*0.75)), 98);
    max_vals(i) = max_d;
end
max_overall = max(max_vals(:));
max_vals_med = max_vals - median_v;
max_vals_med(max_vals_med<0) = 0;
normalizedArray = 1- abs(max_vals_med - 0) / (max(max_vals_med(:)) - 0);

% Plot the data
figure; 
tiledlayout(8, 11);
for i = 1:88
    nexttile
    rectangle('Position', [0, -70, numel(d), 40], 'FaceColor', [1, normalizedArray(i), normalizedArray(i)]*0.9)
    hold on;
    d = meanData{i};
    for r = 1:3
        dd = data(:, :, r);
        d2 = dd{i};
        plot(d2, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7)
    end 
    plot(d, 'Color', 'k', 'LineWidth', 1)
    ylim([-70 max_overall*0.9])
    xlim([0 numel(d)])
    title(string(i))
end 
f = gcf;
f.Position = [45  128  1675 902];

% % Plot the stimulus
% figure; 
% tiledlayout(8, 11);
% for i = 1:88
%     nexttile
%     imagesc(pattern.Pats(:, :, i+1))
%     title(string(i))
%     hold on
%     plot([60, 60], [0 48], 'm')
%     plot([90, 90], [0 48], 'm')
%     plot([0, 192], [14 14], 'm')
%     plot([0, 192], [44 44], 'm')
% 
%     plot([75, 75], [0 48], 'c')
%     plot([0, 192], [29 29], 'c')
% end 


end 