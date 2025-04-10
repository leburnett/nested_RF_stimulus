function [data_comb, cmap_id, var_across_reps, var_within_reps, diff_mean, max_data, min_data] = parse_flash_data(f_data, v_data, on_off, slow_fast, PROJECT_ROOT)

    diff_f_data = diff(f_data);
    median_v = median(v_data);
    v2_data = v_data - median_v;

    % Find where the flash stimuli end. 
    idx = find(diff_f_data == min(diff_f_data));

    % TEST - 'idx' values = end of the different groups of full flash stimuli. 
    % figure; 
    % plot(f_data);
    % hold on;
    % for iii = 1:numel(idx)
    %     plot([idx(iii), idx(iii)], [0 400], 'm');
    % end 
    
    % % 340 ms OFF - bkg - 160 ms FLASH. 500 ms = 0.5s. 
    % % 170 ms OFF - bkg - 80 ms FLASH - 0.25s
    if slow_fast == "slow"
        slow_flashes_dur = 5000; % 0.5s * sampling rate.
        flash_dur_ms_slow = 976700;
    elseif slow_fast == "fast"
        fast_flashes_dur = 2500;
        flash_dur_ms_fast = 976700/2;
    end 
    
    % Empty arrays to fill - 196 flashes/values. 
    data_comb = zeros(14, 14);
    cmap_id = zeros(14, 14);
    var_across_reps = zeros(14, 14);
    var_within_reps = zeros(14,14);
    diff_mean = zeros(14,14);
    max_data = zeros(14,14);
    min_data = zeros(14,14);
    % i_num = zeros(14,14);
    
    % Run through the responses to each flash.
    for i = 1:196

        if slow_fast == "slow"

            % Array to collect the data for each flash over the three repetitions.
            data_flash = ones(3, slow_flashes_dur); 
            data_flash_raw = ones(3, slow_flashes_dur); 
        
            for r = 1:3
        
                if r == 1 % rep 1
                    end_t = idx(1); 
                    edge_vals = end_t-flash_dur_ms_slow:slow_flashes_dur:end_t;
                elseif r == 2 % rep 2 
                    end_t = idx(3);
                    edge_vals = end_t-flash_dur_ms_slow:slow_flashes_dur:end_t;
                elseif r == 3 % rep3 
                    end_t = idx(5); 
                    edge_vals = end_t-flash_dur_ms_slow:slow_flashes_dur:end_t;
                end 
        
                if i < 196
                    d = f_data(edge_vals(i):edge_vals(i+1)-1); % frame data during flash. 
                    v = v2_data(edge_vals(i):edge_vals(i+1)-1); % median-subtracted voltage data
                    v_raw = v_data(edge_vals(i):edge_vals(i+1)-1);
                elseif i == 196
                    d = f_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
                    v = v2_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
                    v_raw = v_data(edge_vals(196):edge_vals(196)+slow_flashes_dur-1);
                end 
        
                data_flash(r, :) = v;
                data_flash_raw(r, :) = v_raw;
            end

        elseif slow_fast == "fast"

            % Array to collect the data for each flash over the three repetitions.
            data_flash = ones(3, fast_flashes_dur); 
            data_flash_raw = ones(3, fast_flashes_dur); 
        
            for r = 1:3
        
                if r == 1 % rep 1
                    end_t = idx(2); 
                    edge_vals = end_t-flash_dur_ms_fast:fast_flashes_dur:end_t;
                elseif r == 2 % rep 2 
                    end_t = idx(4);
                    edge_vals = end_t-flash_dur_ms_fast:fast_flashes_dur:end_t;
                elseif r == 3 % rep3 
                    end_t = idx(6); 
                    edge_vals = end_t-flash_dur_ms_fast:fast_flashes_dur:end_t;
                end 
        
                if i < 196
                    d = f_data(edge_vals(i):edge_vals(i+1)-1); % frame data during flash. 
                    v = v2_data(edge_vals(i):edge_vals(i+1)-1); % median-subtracted voltage data
                    v_raw = v_data(edge_vals(i):edge_vals(i+1)-1);
                elseif i == 196
                    d = f_data(edge_vals(196):edge_vals(196)+fast_flashes_dur-1);
                    v = v2_data(edge_vals(196):edge_vals(196)+fast_flashes_dur-1);
                    v_raw = v_data(edge_vals(196):edge_vals(196)+fast_flashes_dur-1);
                end 
        
                data_flash(r, :) = v;
                data_flash_raw(r, :) = v_raw;
            end

        end 
        
        % Check for the variance across reps:
        % var_X_rep = var(var(data_flash));
    
        % Standard deviation within each repetition (across time points)
        std_within_trial = std(data_flash_raw, 0, 2);  % std along rows (across time points)
        
        % Mean within each repetition (across time points)
        mean_within_trial = mean(data_flash_raw, 2);  % mean along rows (across time points)
        
        % Coefficient of Variation (CV) within each repetition
        cv_within_trial = std_within_trial ./ mean_within_trial;
        
        % Mean CV across all repetitions (mean within trial CV)
        mean_cv_within_trial = abs(mean(cv_within_trial));
    
        % Standard deviation across time points (across trials)
        std_across_trials = std(data_flash_raw, 0, 1);  % std along columns (across trials)
        
        % Mean across time points (across trials)
        mean_across_trials = mean(data_flash_raw, 1);  % mean along columns (across trials)
        
        % Coefficient of Variation (CV) across repetitions (across trials)
        cv_across_trials = std_across_trials ./ mean_across_trials;
    
        % max_cv_across_trials = prctile(cv_across_trials, 98);
        mean_cv_across_trials = abs(mean(cv_across_trials));
    
        mean_data_flash = mean(data_flash);
        n_vals = size(mean_data_flash, 2);
    
        % Max and min of the mean flash response:
        if slow_fast == "slow"
            max_val_flash = prctile(mean_data_flash(500:end), 98); % Max 
            % min_val_flash = prctile(mean_data_flash(1250:end), 2); % Min in the second half - ignore if min is early. 
            min_val_flash = prctile(mean_data_flash(2500:end), 2); % Min in the second half - ignore if min is early. 
        elseif slow_fast == "fast"
            max_val_flash = prctile(mean_data_flash(200:end), 98); % Max 
            % min_val_flash = prctile(mean_data_flash(1250:end), 2); % Min in the second half - ignore if min is early. 
            min_val_flash = prctile(mean_data_flash(1250:end), 2); % Min in the second half - ignore if min is early. 
        end 
    
        diff_resp = max_val_flash - min_val_flash; 
        diff_med = median_v - min_val_flash;
    
        if abs(max_val_flash)>=abs(min_val_flash) % larger excitatory peak
            
            if diff_resp>3
                    val = max_val_flash;
                    cm = 1;
            else 
                val = mean(mean_data_flash(n_vals*0.75:end));
                cm = 3;
            end 
        elseif abs(max_val_flash)<abs(min_val_flash) % larger inhibitory peak 
            if diff_resp>2.8
                val = min_val_flash;
                cm = 2;
            else 
                val = mean(mean_data_flash(n_vals*0.75:end));
                cm = 3;
            end 
        end 
    
        flash_frame_num = max(d)-1;
        % fnum(i) = flash_frame_num;
        % disp(flash_frame_num)
    
        if on_off == "on" % from 196
            rows = 14 - mod((flash_frame_num - 196), 14);   % Rows decrease from 14 to 1
            cols = floor((flash_frame_num - 196) / 14) + 1; % Columns increase normally
        elseif on_off == "off" % 1- 196
            rows = 14 - mod(flash_frame_num, 14);   % Rows decrease from 14 to 1
            cols = floor(flash_frame_num / 14) + 1; % Columns increase normally
        end
    
        data_comb(rows, cols) = val;
        cmap_id(rows, cols) = cm;
        var_across_reps(rows, cols) = mean_cv_across_trials;
        var_within_reps(rows, cols) = mean_cv_within_trial;
        diff_mean(rows, cols) = diff_resp;
        max_data(rows, cols) = max_val_flash;
        min_data(rows, cols) = min_val_flash;
        % i_num(rows, cols) = i;
    end 

    % Quality check figures:
    save_fig = 0;
    plot_quality_check_var_reps_prctile(var_across_reps, var_within_reps, diff_mean, max_data, min_data, save_fig, PROJECT_ROOT)
    
end 