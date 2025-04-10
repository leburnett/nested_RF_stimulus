function f = plot_rf_estimate_timeseries(data_comb2, cmap_id, f_data, v2_data, slow_flashes_dur, idx, dur_ms, on_off)
% Spatial plot of the timeseries responses to the flashes in each position.
% Colour coded background. 

f = figure;
for i = 1:196

    data_flash = ones(3, slow_flashes_dur); 

    for r = 1:3 

        if r == 1
            end_t = idx(1); %3865; 
        elseif r == 2
            end_t = idx(3);%start_t = 2025440; %2025520;
        elseif r == 3
            end_t = idx(5); %start_t = 4047040; %4047120;
        end 

        edge_vals = end_t-dur_ms:slow_flashes_dur:end_t;

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
        rows = 14 - mod((flash_frame_num - 196), 14);   % Rows decrease from 14 to 1
        cols = floor((flash_frame_num - 196) / 14) + 1; % Columns increase normally
    elseif on_off == "off" % 1- 196
        rows = 14 - mod(flash_frame_num, 14);   % Rows decrease from 14 to 1
        cols = floor(flash_frame_num / 14) + 1; % Columns increase normally
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
    % title(string(i))

end 

sgtitle('160ms flashes - 340ms interval')
f.Position = [77  173  1057  874]; %[77  76   1379   971];


end