function data = parse_bar_data(f_data, v_data)


    diff_f_data = diff(f_data);
    idx = find(diff_f_data == min(diff_f_data)); % where the flash stimuli end.

    % BAR stimulus duration 
    dur_t = (2.273+1.155)*10000*16; % (dur_bar_slow + dur_bar_fast) * acq_speed * n_directions. 
    
    % Find the timings of when the bar stimuli start and end: 
    % 1120 is added because idx(2) etc is the end of the last flash and then we
    % add 1120 because that is the gap between the last flash and the beginning
    % of the bar stimuli. 
    rep1_rng = [idx(2)+1120, idx(2)+dur_t];
    rep2_rng = [idx(4)+1120, idx(4)+dur_t];
    rep3_rng = [idx(6)+1120, numel(f_data)]; % til the end of the recording.
    
    %% Combine the voltage timeseries data from the three repetitions into 
    % one data structure. 
    
    rep_ranges = {rep1_rng, rep2_rng, rep3_rng};
    idxs_all = cell(1, 3); 
    
    for i = 1:3
        st_val = rep_ranges{i}(1);
        end_val = rep_ranges{i}(2);
        
        frames_rep = f_data(st_val:end_val);
        
        r_min = islocalmin(frames_rep);
        r_min_vals = find(r_min);
        r_st = [st_val, r_min_vals + st_val, end_val];
    
        r_max = islocalmax(frames_rep);
        r_max_vals = find(r_max);
        r_nd = r_max_vals + st_val;
        
        all_idxs = [r_st, r_nd];
        idxs_all{i} = sort(all_idxs); % Store sorted indices
    end
    
    % Extract data segments using the computed indices
    num_segments = numel(idxs_all{1}) - 1; % Assuming all idxs_all{i} have the same size
    data = cell(num_segments, 3); % Preallocate data cell array
    
    for i = 1:num_segments
        for j = 1:3
            data{i, j} = v_data(idxs_all{j}(i):idxs_all{j}(i+1)-1);
        end
    end

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