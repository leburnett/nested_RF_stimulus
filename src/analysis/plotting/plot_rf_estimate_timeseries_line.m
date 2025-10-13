function f = plot_rf_estimate_timeseries_line(data_comb2, cmap_id, f_data, v2_data, slow_fast, px_size, idx, params)
% Spatial plot of the timeseries responses to the flashes in each position.
% Colour coded lines. Red above baseline, blue below. 

% Inputs
% ------

% 'data_comb2' - contains the min/max val across the flash trial. Used to
% set the intensity of the coloured background behind each timeseries.

% 'cmap_id' - [n_grid_rows, n_grid_cols] size array containing the values
% either 1,2 or 3. Sets the colour group that should be used for the 
% background of each flash subplot. 

% 'f_data' - frame data across the entire experiment. 10kHz acquisition
% rate. 

% 'v2_data' mediaun subtracted voltage data across the entire experiment. 
% 10kHz acquisition rate. 

% 'slow_fast' - whether the flash speed was the slower speed or the faster
% speed. Currently there are only these two hardcoded options. 

% 'px_size' - size of the flash stimuli in pixels. Either 4 or 6 pixels for
% protocol 2. 

% 'idx' - indices of when the first 4 pixel flash starts for each rep. 

% 'params' - struct containing metadata about the experiment.

% _______________________________________________________________________________________________

f = figure;

% Define subplot grid dimensions
if px_size == 4
    nRows = 14;
    nCols = 14;
    n_flashes = 196;
elseif px_size == 6
    nRows = 10;
    nCols = 10;
    n_flashes = 100;
end 

if params.on_off == "on"
    % For finding the end of the 6 pixel flashes
    drop_at_end = -200; % ON flashes. the last 6 pixel flash is frame 200.
elseif params.on_off == "off"
    % For finding the end of the 6 pixel flashes
    drop_at_end = -100; % ON flashes. the last 6 pixel flash is frame 200.
end 

diff_f_data = diff(f_data);
idx_6px = find(diff_f_data == drop_at_end); % where the flash stimuli end.
idx_6px([1,3,5]) = []; 

% Define spacing reduction factor
spacingFactor = 0.1; % Reducing spacing by 75%

% Compute original subplot size
origWidth = 1 / nCols;
origHeight = 1 / nRows;

% Compute new subplot size with reduced spacing
newWidth = origWidth * (1 - spacingFactor);
newHeight = origHeight * (1 - spacingFactor);

if slow_fast == "slow"
    flashes_dur = 7000; % 0.5s * sampling rate.
    % dur_ms = 976700;
    speed_str = "160ms_flash";
% elseif slow_fast == "fast"
%     flashes_dur = 2500;
%     dur_ms = 976700/2;
%     speed_str = "80ms_flash";
end

for i = 1:n_flashes

    data_frame = ones(3, flashes_dur);
    data_flash = ones(3, flashes_dur); 

    for r = 1:3 

        if slow_fast == "slow"
            if r == 1 % rep 1
                if px_size == 4
                    rng_rep1 = idx(1):idx(2);
                else 
                    rng_rep1 = (idx(2):idx_6px(1));
                end 
                start_idx = rng_rep1(1);
                start_flash_idxs = find(diff(f_data(rng_rep1))>0)+start_idx-1;
            elseif r == 2 % rep 2 
                if px_size == 4
                    rng_rep2 = idx(3):idx(4);
                else
                    rng_rep2 = (idx(4):idx_6px(2));
                end 
                start_idx = rng_rep2(1);
                start_flash_idxs = find(diff(f_data(rng_rep1))>0)+start_idx-1;
            elseif r == 3 % rep3 
                if px_size == 4
                    rng_rep3 = idx(5):idx(6);
                else
                    rng_rep3 = (idx(6):idx_6px(3));
                end
                start_idx = rng_rep3(1);
                start_flash_idxs = find(diff(f_data(rng_rep1))>0)+start_idx-1;
            end
        end 

         % Extract data 1000 timepoints before the flash starts
        % til the end of the interval before the next flash. 
        d = f_data(start_flash_idxs(i)-1000:start_flash_idxs(i)+6000-1); % frame data during flash. 
        v = v2_data(start_flash_idxs(i)-1000:start_flash_idxs(i)+6000-1); % median-subtracted voltage data
               
        data_frame(r, :) = d;
        data_flash(r, :) = v;
    end 

    mean_data_flash = mean(data_flash);
    mean_data_flash = downsample(mean_data_flash, 10);

    flash_frame_num = max(d)-1;

     if params.on_off == "on" % from 196
        rows = nRows - mod((flash_frame_num - n_flashes), nRows);   % Rows decrease from 14 to 1
        cols = floor((flash_frame_num - n_flashes) / nRows) + 1; % Columns increase normally
    elseif params.on_off == "off" % 1- 196
        rows = nRows - mod(flash_frame_num, nRows);   % Rows decrease from 14 to 1
        cols = floor(flash_frame_num / nRows) + 1; % Columns increase normally
    end

    % Get information about how intense the background colour should be:
    val  = data_comb2(rows, cols);
    % Get information about what colour map to use:
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
        plot(ax, x(1:150), y(1:150),'Color', [0.8, 0.8, 0.8], 'LineWidth', 3); hold on
        plot(ax, x(151:end), y(151:end),'Color', [0.3, 0.3, 0.8], 'LineWidth', 3);

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

    plot([1 xmax], [0, 0], 'Color', [0.2 0.2 0.2]) % Plot '0' = median. 
    ylim([-10 25])

    % if px_size == 4
    %     first_subpl = 171;
    % elseif px_size == 6
    %     first_subpl = 60; %% % % Check that this is correct. 
    % end 

    if rows == 1 && cols == 1
        xticks([0, 160, 500])
        yticks([-10, 0, 25])
        ax = gca;
        ax.TickDir = 'out';
        ax.TickLength = [0.05 0.05];
        ax.LineWidth = 1.5;
        xticklabels('')
        yticklabels('')
    else 
        axis off
    end 
    % if i ~= 171
    % axis off % keep axes for the first subplot.
    % end 

    box off
    axis square

end 

% Add an arrow on top representing the vector sum of responses to bar
% stimuli. 
magnitude_arrow = 0.85;
angl = params.resultant_angle;
%  'angl' = direction of vector sum of bar responses in radians.
center_x = 0.5;
center_y = 0.5;
arrow_x = magnitude_arrow * cos(angl);
arrow_y = magnitude_arrow * sin(angl);
annotation('arrow', [center_x - arrow_x/2, center_x + arrow_x/2], ...
                    [center_y - arrow_y/2, center_y + arrow_y/2], ...
           'LineWidth', 6, 'Color', [0.2 0.2 0.2], "HeadLength", 40, "HeadWidth", 40);


% sgtitle('160ms flashes - 340ms interval')
f.Position = [77   265   862   782]; %[77  76   1379   971];

exportgraphics(f ...
        , strcat(params.date, "_", params.time, "_", params.strain, "_spatialRF_timeseries_", speed_str, "_px-size", string(px_size), "_wArrow.pdf") ...
        , 'ContentType', 'vector' ...
        , 'BackgroundColor', 'none' ...
        ); 

exportgraphics(f ...
        , strcat(params.date, "_", params.time, "_", params.strain, "_spatialRF_timeseries_", speed_str, "_px-size", string(px_size), "_wArrow.png") ...
        ); 

end