function f = plot_bar_flash_data(data, meanData, n_int)
% 'data'      : [11 x 8 x 3] cell (columns x rows x reps)
% 'meanData'  : [11 x 8]     cell (columns x rows)
% This version circularly shifts all rows so the row with the global max
% of max_vals_base lands at row 5, preserving relative order.
% 'int_ms' : double - time in ms of the gap between flashes. Use this to 


    % -------------------- Compute PD orientation and shift "data" to have this row in the middle --------------------
    [meanData_out, data_out, ~, ~] = shiftMaxColumnTo5(meanData, data);

    % -------------------- Compute background intensities --------------------
    max_vals = zeros(8,11);
    baseline_vals = zeros(8,11);

    for i = 1:88
        d = meanData_out{i};
        n_points = numel(d);
        max_d = prctile(d(n_int:ceil(n_points*0.8)), 98); % max response per tile
        baseline_v = nanmean(d(ceil(n_int*0.5):n_int));               % baseline per tile
        max_vals(i) = max_d;
        baseline_vals(i) = baseline_v;
    end

    % max_overall = max(max_vals(:));
    max_vals_base = max_vals - baseline_vals;
    max_vals_base(max_vals_base < 0) = 0;

    if all(max_vals_base(:) == 0)
        normalizedArray = ones(8,11); % avoid divide-by-zero; all same bg
    else
        % normalizedArray = 1 - max_vals_base ./ max(max_vals_base(:));
        normalizedArray = 1- abs(max_vals_base - 0) / (max(max_vals_base(:)) - 0);
    end

    % -------------------- Plot --------------------

   % Plot the data
    figure; 
    tiledlayout(8, 11);
    for i = 1:88

        nexttile
        rectangle('Position', [0, -75, numel(d), 45], 'FaceColor', [1, normalizedArray(i), normalizedArray(i)]*0.9)
        hold on;
        d = meanData_out{i};
        % Plot the reps
        for r = 1:3
            dd = data_out(:, :, r);
            d2 = dd{i};
            plot(d2, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7)
        end 
        % Plot the mean
        plot(d, 'Color', 'k', 'LineWidth', 1)

        % Plot the beginning and end of flash
        plot([n_int, n_int], [-75 -30], 'k', 'LineWidth', 0.5);
        plot([numel(d)-n_int, numel(d)-n_int], [-75 -30], 'k', 'LineWidth', 0.5);

        % ylim([-70 max_overall*0.9])
        ylim([-75 -30])
        xlim([0 numel(d)])
        box off
        xticks([])

        title(string(i))

        if i == 45
            ylabel("PD", 'FontSize', 20)
            yticks([-70, -60, -50, -40, -30]);
        elseif i == 1
            ylabel("OD", 'FontSize', 20)
            yticks([-70, -60, -50, -40, -30]);
        else
            yticks([]);
        end 
    end 
    f = gcf;
    f.Position = [45  128  1675 902];

end
