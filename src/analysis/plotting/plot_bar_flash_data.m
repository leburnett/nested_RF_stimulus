function f = plot_bar_flash_data(data, meanData, median_v)
% PLOT_BAR_FLASH_DATA  Plot bar flash responses as an 8x11 grid.
%
%   F = PLOT_BAR_FLASH_DATA(DATA, MEANDATA, MEDIAN_V) creates a figure
%   with 88 subplots (8 orientations x 11 positions) showing the voltage
%   response to each bar flash stimulus.
%
%   INPUTS:
%     data      - 11x8x3 cell array of individual rep timeseries
%     meanData  - 11x8 cell array of mean timeseries across reps
%     median_v  - Median voltage for background coloring
%
%   OUTPUT:
%     f - Figure handle
%
%   See also PROCESS_BAR_FLASHES_P2, PARSE_BAR_FLASH_DATA

% Find background colour based on response amplitude:
max_vals = zeros([11, 8]); % Same shape as meanData (positions x orientations)
for i = 1:88
    d = meanData{i};
    if isempty(d)
        max_vals(i) = median_v;
        continue;
    end
    n_points = numel(d);
    max_d = prctile(d(ceil(n_points*0.5):ceil(n_points*0.75)), 98);
    max_vals(i) = max_d;
end
max_overall = max(max_vals(:));
max_vals_med = max_vals - median_v;
max_vals_med(max_vals_med<0) = 0;
max_range = max(max_vals_med(:));
if max_range == 0
    normalizedArray = ones(size(max_vals_med));
else
    normalizedArray = 1 - abs(max_vals_med - 0) / max_range;
end

% Plot the data
f = figure;
tiledlayout(8, 11);
for i = 1:88
    nexttile
    d = meanData{i};
    if isempty(d)
        title(string(i))
        continue;
    end
    rectangle('Position', [0, -70, numel(d), 40], 'FaceColor', [1, normalizedArray(i), normalizedArray(i)]*0.9)
    hold on;
    for r = 1:3
        dd = data(:, :, r);
        d2 = dd{i};
        if ~isempty(d2)
            plot(d2, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7)
        end
    end
    plot(d, 'Color', 'k', 'LineWidth', 1)
    ylim([-70 max_overall*0.9])
    xlim([0 numel(d)])
    title(string(i))
end
f.Position = [45  128  1675 902];

end
