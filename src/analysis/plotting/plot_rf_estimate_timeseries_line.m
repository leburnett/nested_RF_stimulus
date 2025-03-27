function f = plot_rf_estimate_timeseries_line(data_comb2, cmap_id, f_data, v2_data, slow_flashes_dur, idx, dur_ms, on_off)
% Spatial plot of the timeseries responses to the flashes in each position.
% Colour coded lines. Red above baseline, blue below. 

f = figure;

% Define subplot grid dimensions
nRows = 14;
nCols = 14;

% Define spacing reduction factor
spacingFactor = 0.1; % Reducing spacing by 75%

% Compute original subplot size
origWidth = 1 / nCols;
origHeight = 1 / nRows;

% Compute new subplot size with reduced spacing
newWidth = origWidth * (1 - spacingFactor);
newHeight = origHeight * (1 - spacingFactor);

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
    mean_data_flash = downsample(mean_data_flash, 10);

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

    % Compute subplot position
    left = (cols - 1) * origWidth + (origWidth - newWidth) / 2;
    bottom = (nRows - rows) * origHeight + (origHeight - newHeight) / 2;

    % Create axes with reduced spacing
    ax = axes('Position', [left, bottom, newWidth, newHeight]);

    xmax = numel(mean_data_flash);
    x = 1:xmax;
    y = mean_data_flash;

    hold on

    if cm ~= 3
        plot(ax, x,y,'Color', [0.3, 0.3, 0.8], 'LineWidth', 3);

        % if cm == 1
            col_r = [1, (1-val), (1-val)];
        % end 

        valsPos = find(y>0);
        if ~isempty(valsPos)
            transitions = find(diff(valsPos)>1);
    
            st_idx = [valsPos(1), valsPos(transitions+1)];
            end_idx = [valsPos(transitions), valsPos(end)];
    
            for kk = 1:numel(st_idx)
                plot(ax, x(st_idx(kk):end_idx(kk)), y(st_idx(kk):end_idx(kk)), 'Color', col_r, 'LineWidth', 3);
            end 
        end 

    elseif cm == 3 % grey 
        col = [0.8 0.8 0.8];
        plot(ax, y, 'Color', col, 'LineWidth', 2)
    end 

    plot([1 xmax], [0, 0], 'Color', [0.7 0.7 0.7]) % Plot '0' = median. 
    ylim([-10 25])

    if i ~= 1
        axis off % keep axes for the first subplot.
    end 

    box off
    axis square

end 

% sgtitle('160ms flashes - 340ms interval')
f.Position = [77   265   862   782]; %[77  76   1379   971];

end