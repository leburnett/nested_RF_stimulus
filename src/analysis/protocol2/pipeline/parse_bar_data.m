function data = parse_bar_data(f_data, v_data, on_off)
% PARSE_BAR_DATA  Extract individual bar stimulus responses from recording.
%
%   DATA = PARSE_BAR_DATA(F_DATA, V_DATA, ON_OFF) identifies the timing of
%   each moving bar stimulus presentation and extracts the corresponding
%   voltage responses for all 3 repetitions.
%
%   INPUTS:
%     f_data - 1xN array of frame numbers at each time point
%     v_data - 1xN array of voltage values at each time point (10kHz)
%     on_off - 'on' or 'off' indicating which contrast was presented
%
%   OUTPUT:
%     data - Mx4 cell array where M = number of bar stimuli (32 total):
%            data{i,1} - Rep 1 voltage trace for bar stimulus i
%            data{i,2} - Rep 2 voltage trace for bar stimulus i
%            data{i,3} - Rep 3 voltage trace for bar stimulus i
%            data{i,4} - Mean voltage across all 3 reps
%            Each trace includes 9s pre-stimulus and 9s post-stimulus
%
%   BAR STIMULUS STRUCTURE:
%     Protocol 2 presents 16 bar directions at 2 speeds = 32 bar stimuli.
%     Each bar direction includes forward and backward movement.
%     Order: 8 orientations x 2 directions x 2 speeds
%
%   TIMING DETECTION:
%     1. Finds end of flash stimuli (frame drops by -100 or -200)
%     2. Locates start of bar stimuli after gap_between_flash_and_bars
%     3. Identifies individual bar boundaries by frame differences > 9
%     4. Extracts voltage with 9s padding before and after each bar
%
%   NOTE:
%     Assumes 10kHz sampling rate (10000 samples = 1 second)
%     Bar stimuli are preceded by 1s static gray screen
%
%   See also PROCESS_BARS_P2, LOAD_PROTOCOL2_DATA

    diff_f_data = diff(f_data);

    if on_off == "on"
       drop_at_end = -200; % ON flashes. the last 6 pixel flash is frame 200.
    elseif on_off == "off"
        drop_at_end = -100; % OFF flashes. the last 6 pixel flash is frame 100.
    end 

    % This finds when the difference in frame number == drop_at_end. 
    % This should find 6 values. Values 1, 3 and 5 are wthin the 4 pixel
    % flashes and values 2,4 and 6 correspond to the timing of the very end
    % of the last flash.
    idx = find(diff_f_data == drop_at_end); % where the flash stimuli end.
    idx([1,3,5]) = []; % Remove timepoints within flashes - only want the timepoints that correspond to the end of the 6px flashes. 

    %% Check - timepoints of end of 6 pixel flash stimuli.

    % figure; 
    % plot(f_data);
    % hold on;
    % for i = 1:numel(idx)
    %     plot([idx(i), idx(i)], [0 abs(drop_at_end)], 'r');
    % end 

    %% Find the range of timepoints during which all of the bar stimuli are presented - per rep. 

    % 1s of grey is now added before every bar movement. It seems that this is
    % ~9980 datapoints - almost 10000.
    gap_between_flash_and_bars = 4210; % Number of timepoints from end of last 6 pixel flash and ~1000ms before the first bar stimulus. 
    array2find = zeros(1, 20000);
    
    % Find the timings of when the bar stimuli start and end: 
    start_f1 = find(f_data(idx(1)+gap_between_flash_and_bars:end) > 0, 1, 'first') + (idx(1)+gap_between_flash_and_bars);
    end_f1 = strfind(f_data(start_f1:end), array2find) + start_f1; % % % This is very slow. 
    end_f1 = end_f1(1);
    rep1_rng = [start_f1, end_f1];

    start_f2 = find(f_data(idx(2)+gap_between_flash_and_bars:end) > 0, 1, 'first') + (idx(2)+gap_between_flash_and_bars);
    end_f2 = strfind(f_data(start_f2:end), array2find) + start_f2;
    end_f2 = end_f2(1);
    rep2_rng = [start_f2, end_f2];

    start_f3 = find(f_data(idx(3)+gap_between_flash_and_bars:end) > 0, 1, 'first') + (idx(3)+gap_between_flash_and_bars);
    rep3_rng = [start_f3, numel(f_data)]; % til the end of the recording.

    %% Check - timepoints corresponding to the start and end of each repetition of bar stimuli (slow and fast)
    
    % figure; 
    % plot(f_data);
    % hold on;
    % for i = 1:numel(rep1_rng)
    %     plot([rep1_rng(i), rep1_rng(i)], [0 100], 'm');
    % end 
    % for i = 1:numel(rep2_rng)
    %     plot([rep2_rng(i), rep2_rng(i)], [0 100], 'r');
    % end 
    % for i = 1:numel(rep3_rng)
    %     plot([rep3_rng(i), rep3_rng(i)], [0 100], 'g');
    % end 
    
    %% Create cell array of the indices of the beginning and end of each individual bar stimulus in each rep.
    
    rep_ranges = {rep1_rng, rep2_rng, rep3_rng};
    idxs_all = cell(1, 3); 
    
    for i = 1:3
        st_val = rep_ranges{i}(1);
        end_val = rep_ranges{i}(2);
        
        frames_rep = f_data(st_val:end_val);
        diff_vals = diff(frames_rep);
        dd = find(abs(diff_vals)>9);
        all_idxs = [st_val, dd + st_val];
        idxs_all{i} = sort(all_idxs); % Store sorted indices
    end
    
    %% Check - timing of start and end of each individual bar stimuli 
    
    % figure; plot(f_data); hold on;
    % for iii = 1:numel(idxs_all{1,1})
    %     all_rep1 = idxs_all{1, 1};
    %     x_val = all_rep1(iii);
    %     plot([x_val, x_val], [0 75], 'r');
    % end 

    %% Extract the frame ranges for each bar stimulus - including 9000ms before the bar starts and 9000ms after the bar stops.
    
    num_segments = numel(idxs_all{1}) -1; % Assuming all idxs_all{i} have the same size
    data = cell(num_segments, 3); % Preallocate data cell array
    
    interval_t_ms = 9000; % 10,000Hz acquisition - 10000 = 1s.

    for i = 1:2:num_segments % Every other idx = JUST the moving bar segments not the interval. 
        for j = 1:3
            data{i, j} = v_data(idxs_all{j}(i)- interval_t_ms:(idxs_all{j}(i+1)-1)+interval_t_ms);
        end
    end

    %% Combine the data across repetitions into one structure and calculate the mean response across repetitions.

    % Remove empty rows
    emptyRows = all(cellfun(@isempty, data), 2);
    data(emptyRows, :) = [];

    for j = 1:height(data)
        % Extract data for all repetitions and ensure column vectors
        d = cellfun(@(x) x(:), data(j, 1:3), 'UniformOutput', false); 
        n_per_col = cellfun(@numel, d); % Get the number of elements in each repetition
    
        % Find minimum length and trim all series to match
        min_val = min(n_per_col);
        d_trimmed = cellfun(@(x) x(1:min_val), d, 'UniformOutput', false);
    
        % Convert to matrix and compute mean across rows (time series)
        d_matrix = horzcat(d_trimmed{:}); % Convert cell array to matrix (columns = repetitions)
        mean_resp = nanmean(d_matrix, 2); % Compute mean across repetitions (columns)
    
        % Store mean time series in the 4th column
        data{j, 4} = mean_resp';
    end


end 