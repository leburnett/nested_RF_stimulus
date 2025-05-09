function data = parse_bar_data(f_data, v_data)


    diff_f_data = diff(f_data);
    idx = find(diff_f_data == min(diff_f_data)); % where the flash stimuli end.

    % BAR stimulus duration 
    % dur_t = (2.273+1.155)*10000*16; % (dur_bar_slow + dur_bar_fast) * acq_speed * n_directions. 
    
    % Find the timings of when the bar stimuli start and end: 
    % 1120 is added because idx(2) etc is the end of the last flash and then we
    % add 1120 because that is the gap between the last flash and the beginning
    % of the bar stimuli.
    start_f1 = find(f_data(idx(2)+1121:end) > 0, 1, 'first') + (idx(2)+1120);
    end_f1 = find(f_data(start_f1:end) == 0, 1, 'first') + start_f1;
    rep1_rng = [start_f1, end_f1];

    start_f2 = find(f_data(idx(4)+1121:end) > 0, 1, 'first') + (idx(4)+1120);
    end_f2 = find(f_data(start_f2:end) == 0, 1, 'first') + start_f2;
    rep2_rng = [start_f2, end_f2];

    start_f3 = find(f_data(idx(6)+1121:end) > 0, 1, 'first') + (idx(6)+1120);
    rep3_rng = [start_f3, numel(f_data)]; % til the end of the recording.
    
    %% Combine the voltage timeseries data from the three repetitions into 
    % one data structure. 
    
    rep_ranges = {rep1_rng, rep2_rng, rep3_rng};
    idxs_all = cell(1, 3); 
    
    for i = 1:3
        st_val = rep_ranges{i}(1);
        end_val = rep_ranges{i}(2);
        
        frames_rep = f_data(st_val:end_val);
        
        % Look for min / max 
        r_min = islocalmin(frames_rep);
        r_min_vals = find(r_min);
        r_st = [st_val, r_min_vals + st_val, end_val];
    
        r_max = islocalmax(frames_rep);
        % r_max = find(frames_rep(1:end-1) == 60 & frames_rep(2:end) == 61);
        r_max_vals = find(r_max);
        r_nd = r_max_vals + st_val;
        
        all_idxs = [r_st, r_nd];
        idxs_all{i} = sort(all_idxs); % Store sorted indices
    end
    
% TEST 
% Plot the parsing between stimuli

% figure; plot(f_data); hold on;
% for iii = 1:numel(idxs_all{1,1})
%     all_rep1 = idxs_all{1,1};
%     x_val = all_rep1(iii);
%     plot([x_val, x_val], [0 75], 'r');
% end 
% 
% for kk = 1:8
%     plot([r_max(kk)+ rep_ranges{i}(1), r_max(kk)+ rep_ranges{i}(1)], [0 75], 'g'); hold on
% end 

    % Extract data segments using the computed indices
    num_segments = numel(idxs_all{1}) - 1; % Assuming all idxs_all{i} have the same size
    data = cell(num_segments, 3); % Preallocate data cell array
    
    for i = 1:num_segments
        for j = 1:3
            data{i, j} = v_data(idxs_all{j}(i):idxs_all{j}(i+1)-1);
        end
    end

    % stf = 1828960;
    % endf = 1852243;
    % data{i, j} = v_data(stf:endf);

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