function data = parse_bar_data_pharma(f_data, v_data)
    % The voltage and frame data are sampled at 10,000Hz. Therefore, each
    % value corresponds to (1/10000)s, 0.0001s or 100us.

     % The order of the stimuli (within one repetition)
    % - 10s grey screen 
    % - bar sweeps - 28 dps
    % - 3s grey screen
    % - bar sweeps - 56 dps
    % - 3s grey screen
    % - bar sweeps - 168 dps
    % - 3s grey screen
    % - bar sweeps - 250 dps
    % - 3s grey screen
    % - bar sweeps - 500 dps
    % - 3s grey screen
    % - bar flashes - 80 ms (~28 dps)
    % - 3s grey screen
    % - bar flashes - 14 ms (~ 168dps)

    %% 1 - Find the beginning of each repetition. This will follow a 10s
    % grey screen interval. (This seems long winded, but it's a lot faster 
    % than using "strcomp")

    zero_mask = f_data == 0;
    d = diff([0 zero_mask 0]);
    start_idx = find(d == 1);
    end_idx = find(d == -1) - 1;
    len = end_idx - start_idx + 1;
    idx_3 = start_idx(len >= 30000);

    % Plot all of the 3s gaps
    % figure;
    % plot(f_data);
    % for i = 1:numel(idx_3)
    %     hold on;
    %     plot([idx_3(i), idx_3(i)], [0 200], 'm')
    % end

    %% Frame range of all 3 bar sweep stimuli per rep.

    rep1_rng = [idx_3(1), idx_3(6)];
    rep2_rng = [idx_3(8), idx_3(13)];
    rep3_rng = [idx_3(15), idx_3(20)];

    %% Check - timepoints corresponding to the start and end of each repetition of bar stimuli (slow and fast)
    
    % % y-range for shading
    % y_low = 0; %150;
    % y_high = 100; %200;
    % 
    % % Plot your line (example)
    % figure;
    % hold on;
    % 
    % % Define the ranges to shade
    % ranges = [rep1_rng; rep2_rng; rep3_rng];
    % 
    % % Loop through each range and add a semi-transparent grey patch
    % for i = 1:size(ranges,1)
    %     x1 = ranges(i,1);
    %     x2 = ranges(i,2);
    %     patch([x1 x2 x2 x1], [y_low y_low y_high y_high], ...
    %           [1 0.7 0.7], ...     % grey color (RGB)
    %           'FaceAlpha', 0.3, ...  % transparency (0 = fully transparent, 1 = opaque)
    %           'EdgeColor', 'none');  % remove border
    % end
    % 
    % plot(f_data);
    % ylabel('Frame')
    % f = gcf;
    % f.Position = [18  714  1749  244];

    
    %% Create cell array of the indices of the beginning and end of each individual bar stimulus in each rep.
    
    rep_ranges = {rep1_rng, rep2_rng, rep3_rng};
    idxs_all = cell(1, 3); 
    
    for i = 1:3
        st_val = rep_ranges{i}(1);
        end_val = rep_ranges{i}(2);
        
        frames_rep = f_data(st_val:end_val);
        diff_vals = diff(frames_rep);
        dd = find(abs(diff_vals)>9);
        all_idxs = dd + st_val;
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
    
    for i = 1:2:num_segments % Every other idx = JUST the moving bar segments not the interval. Line at the beginning of the sweep. 
        
        if i < 64
            interval_t_ms = 9000; % 10,000Hz acquisition - 10000 = 1s gap.
        elseif i >= 64 && i < 128 
            interval_t_ms = 5000; % 750ms gap. 
        elseif i >= 128
            interval_t_ms = 4000; % 500ms gap
        end 

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