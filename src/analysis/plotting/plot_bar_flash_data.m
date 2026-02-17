function f = plot_bar_flash_data(data, meanData, n_int)
% PLOT_BAR_FLASH_DATA  Tiled figure of bar flash responses across orientations and positions.
%
%   F = PLOT_BAR_FLASH_DATA(DATA, MEANDATA, N_INT)
%   creates an 8x11 tiled figure showing voltage timeseries for each
%   bar flash condition (8 orientations x 11 positions). Orientations are
%   circularly shifted so the preferred direction (PD) is centred at row 5.
%   Individual repetitions are shown in grey, the mean in black, and
%   background tile colour indicates response magnitude (red = strong).
%
%   INPUTS:
%     data     - 11x8x3 cell array
%                Voltage timeseries from PARSE_BAR_FLASH_DATA.
%                Dimensions: positions x orientations x repetitions.
%                Each cell contains a 1xT double timeseries (mV).
%     meanData - 11x8 cell array
%                Mean timeseries across the 3 reps.
%                Each cell contains a Tx1 double vector (mV).
%     n_int    - double
%                Number of samples of inter-flash interval included in
%                each timeseries (before and after the flash). Used to
%                draw vertical lines marking flash onset and offset, and
%                to compute per-tile baseline voltage.
%
%   OUTPUT:
%     f - figure handle
%         Handle to the generated figure (1675x902 pixels).
%
%   FIGURE LAYOUT:
%     8 rows (orientations, PD-aligned) x 11 columns (bar positions).
%     Row 5 = preferred direction (PD), labelled on y-axis.
%     Row 1 = orthogonal direction (OD), labelled on y-axis.
%     Background colour: white = no response, red = strong response.
%     Black vertical lines mark flash onset (at n_int) and offset
%     (at numel(d)-n_int).
%     Y-axis: [-75, -30] mV. Y-ticks shown only for PD and OD rows.
%
%   See also PARSE_BAR_FLASH_DATA, PROCESS_BAR_FLASHES_P2, SHIFTMAXCOLUMNTO5


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
