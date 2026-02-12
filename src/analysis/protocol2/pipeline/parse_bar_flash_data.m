function [data_slow, data_fast, mean_slow, mean_fast] = parse_bar_flash_data(f_data, v_data)
    
    % The voltage and frame data are sampled at 10,000Hz. Therefore, each
    % value corresponds to (1/10000)s, 0.0001s or 100us.
    
     % The order of the stimuli (within one repetition)
    % - 10s grey screen 
    % - 4 pixel flashes
    % - 6 pixel flashes
    % - 3s grey screen
    % - bar sweeps - 28 dps
    % - 3s grey screen
    % - bar sweeps - 56 dps
    % - 3s grey screen
    % - bar sweeps - 168 dps
    % - 3s grey screen
    % - bar flashes - 80 ms (~28 dps)
    % - 3s grey screen
    % - bar flashes - 14 ms (~ 168dps)
    
     %% 1 - Find the beginning of each repetition. This will follow a 10s
    % grey screen interval. (This seems long winded, but it's a lot faster 
    % than using "strcomp")
    
    num_reps = 3;
    num_speeds = 2;
    num_orient = 8;
    num_positions = 11;

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

    for sp = 1:num_speeds

        % Currently, rows = positions and cols = orient.
        % Will flip at end.
        data = cell(num_positions, num_orient, num_reps); 

        % These timings are from the end of the last flash - before the
        % gap.
        if sp == 1 % slow
            rep1_rng = idx_3(5:6);
            rep2_rng = idx_3(11:12);
            rep3_rng = idx_3(17:18);
            gap_between_flashes = 5000; % Gap between flashes is 1s (10000) but only clip half that for the timeseries.
        elseif sp == 2 % fast
            rep1_rng = idx_3(6:7);
            rep2_rng = idx_3(12:13);
            rep3_rng = idx_3(18:19);
            gap_between_flashes = 2500;
        end 

        % Find the timeseries per bar flash position for the first rep. 

        for r = 1:num_reps

            data_rep = cell(num_positions, num_orient); 

            if r == 1
                rng = rep1_rng;
            elseif r == 2
                rng = rep2_rng;
            elseif r == 3
                rng = rep3_rng;
            end 
            
            % Find all of the times that the frame # >1 within this range
            flash_idxs = find(diff(f_data(rng(1):rng(2))>0))+rng(1)-1;
            start_flash_idxs = flash_idxs(1:2:end);
            end_flash_idxs = flash_idxs(2:2:end);

            n_flashes = numel(start_flash_idxs);

            for f = 1:n_flashes
                st = start_flash_idxs(f);
                nd = end_flash_idxs(f);
                flash_frame_num = max(f_data(st:nd));
                data_flash = v_data(st-gap_between_flashes:nd+gap_between_flashes);
                data_rep{flash_frame_num} = data_flash;
            end 

            % Add the rep data to "data"
            data(:, :, r) = data_rep(:, :);
        end 


        % Create "meanData" - mean across reps.
        [rows, cols, reps] = size(data);
        meanData = cell(rows, cols);
        
        for i = 1:rows
            for j = 1:cols
                ts = data(i,j,:);
                ts = ts(:);
                maxlen = max(cellfun(@length, ts));
                
                padded = NaN(maxlen, reps);
                for k = 1:reps
                    y = ts{k};
                    padded(1:length(y),k) = y;
                end
                
                meanData{i,j} = nanmean(padded, 2);
            end
        end

        if sp == 1
            data_slow = data;
            mean_slow = meanData;
        elseif sp ==2 
            data_fast = data;
            mean_fast = meanData;
        end 

    end 

end 



