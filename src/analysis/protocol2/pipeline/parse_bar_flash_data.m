function [data_slow, data_fast, mean_slow, mean_fast] = parse_bar_flash_data(f_data, v_data, prop_int)
% PARSE_BAR_FLASH_DATA  Extract bar flash responses from Protocol 2 data.
%
%   [DATA_SLOW, DATA_FAST, MEAN_SLOW, MEAN_FAST] = ...
%       PARSE_BAR_FLASH_DATA(F_DATA, V_DATA, PROP_INT)
%   parses the raw frame and voltage data to extract individual bar flash
%   responses at two speeds (slow = 80ms flash, fast = 14ms flash), across
%   3 repetitions, 8 orientations, and 11 bar positions per orientation.
%
%   INPUTS:
%     f_data   - 1xM double
%                Frame position data sampled at 10 kHz. Values are frame
%                indices; zero indicates grey screen / no stimulus.
%     v_data   - 1xM double
%                Voltage data sampled at 10 kHz (already scaled to mV).
%     prop_int - double (0 to 1)
%                Proportion of the inter-flash interval to include in each
%                extracted timeseries. E.g. 0.75 keeps 75% of the gap
%                before and after each flash. For slow flashes the full
%                gap is 10000 samples (1s); for fast, 5000 samples (0.5s).
%
%   OUTPUTS:
%     data_slow - 11x8x3 cell array
%                 Bar flash voltage timeseries for the slow (80ms) condition.
%                 Dimensions: positions x orientations x repetitions.
%                 Each cell contains a 1xT double voltage timeseries.
%     data_fast - 11x8x3 cell array
%                 Same format as data_slow for the fast (14ms) condition.
%     mean_slow - 11x8 cell array
%                 Mean timeseries across the 3 reps for the slow condition.
%                 Each cell contains a Tx1 double voltage vector (NaN-padded
%                 if reps differ in length).
%     mean_fast - 11x8 cell array
%                 Same format as mean_slow for the fast condition.
%
%   STIMULUS ORDER (within one repetition of the full protocol):
%     1.  10s grey screen
%     2.  4 pixel flashes
%     3.  6 pixel flashes
%     4.  3s grey screen
%     5.  Bar sweeps - 28 dps
%     6.  3s grey screen
%     7.  Bar sweeps - 56 dps
%     8.  3s grey screen
%     9.  Bar sweeps - 168 dps
%     10. 3s grey screen
%     11. Bar flashes - 80 ms (slow, ~28 dps equivalent)
%     12. 3s grey screen
%     13. Bar flashes - 14 ms (fast, ~168 dps equivalent)
%
%   INDEXING:
%     Repetition boundaries are found by detecting zero-value gaps >= 30000
%     samples (3s) in the frame data (idx_3). The slow bar flash block
%     falls between idx_3(5:6), idx_3(11:12), idx_3(17:18) for reps 1-3.
%     The fast block uses idx_3(6:7), idx_3(12:13), idx_3(18:19).
%     Individual flashes within each block are identified by transitions
%     in the frame data (frame > 0), and each flash is indexed into the
%     output cell array by its maximum frame number during the flash.
%
%   See also PROCESS_BAR_FLASHES_P2, PLOT_BAR_FLASH_DATA, PARSE_BAR_DATA

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

    % % Plot all of the 3s gaps
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
            gap_between_flashes = 10000*prop_int; % Gap between flashes is 1s (10000) but only clip half that for the timeseries.
        elseif sp == 2 % fast
            rep1_rng = idx_3(6:7);
            rep2_rng = idx_3(12:13);
            rep3_rng = idx_3(18:19);
            gap_between_flashes = 5000*prop_int;
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



