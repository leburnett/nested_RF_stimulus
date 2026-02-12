function [data_slow, data_fast, mean_slow, mean_fast, debug_info] = parse_bar_flash_data(f_data, v_data, bar_flash_pattern)
% PARSE_BAR_FLASH_DATA  Extract bar flash stimulus responses from recording.
%
%   [DATA_SLOW, DATA_FAST, MEAN_SLOW, MEAN_FAST] = PARSE_BAR_FLASH_DATA(
%       F_DATA, V_DATA, BAR_FLASH_PATTERN) identifies bar flash stimulus
%   epochs and extracts voltage responses for all positions and orientations.
%
%   INPUTS:
%     f_data             - 1xN array of frame numbers at each time point (10kHz)
%     v_data             - 1xN array of voltage values at each time point
%     bar_flash_pattern  - pattern struct from the bar flash .mat file,
%                          containing .Pats (the pattern array) and .x_num
%                          (number of frames in pattern)
%
%   OUTPUTS:
%     data_slow  - 11x8x3 cell array: positions x orientations x reps
%                  (slow bar flashes, e.g., 80ms flash duration)
%     data_fast  - 11x8x3 cell array: same for fast bar flashes
%     mean_slow  - 11x8 cell array: mean across reps for slow
%     mean_fast  - 11x8 cell array: mean across reps for fast
%     debug_info - struct with epoch detection details for verification
%
%   BAR FLASH DETECTION:
%     Rather than using hardcoded gap indices, this function dynamically
%     detects bar flash epochs by looking for frame values in the range
%     of the bar flash pattern (typically 1-89 for 88 flash frames + 1 grey).
%     Bar flash epochs are distinguished from bar sweep epochs by their
%     characteristic frame number pattern: rapid alternation between 0
%     (grey interval) and flash frame numbers (2-89).
%
%   See also PROCESS_BAR_FLASHES_P2, GET_STIMULUS_METADATA

    num_reps = 3;
    num_speeds = 2;
    num_orient = 8;
    num_positions = 11;
    n_frames_pattern = num_orient * num_positions; % 88

    % Determine the max frame number for this pattern
    if nargin >= 3 && ~isempty(bar_flash_pattern)
        max_pattern_frame = bar_flash_pattern.x_num;
    else
        max_pattern_frame = n_frames_pattern + 1; % 89 (88 flashes + 1 grey)
    end

    %% Find all 3s+ grey gaps (frame == 0 for >= 30000 samples)
    zero_mask = f_data == 0;
    d = diff([0 zero_mask 0]);
    start_idx = find(d == 1);
    end_idx = find(d == -1) - 1;
    len = end_idx - start_idx + 1;

    % All gaps >= 3s
    gap_mask = len >= 30000;
    gap_starts = start_idx(gap_mask);
    gap_ends = end_idx(gap_mask);

    %% Identify bar flash epochs by examining the frame data between gaps
    % Bar flash epochs have these characteristics:
    % - Frame numbers alternate between 0 (interval) and small values (2-89)
    % - Many transitions (88 flashes per epoch)
    % - Max frame number <= max_pattern_frame (typically 89)
    %
    % Bar sweep epochs have:
    % - Sawtooth frame patterns sweeping through larger ranges
    % - Fewer transitions
    % - Max frame number up to 288

    n_gaps = numel(gap_starts);
    epoch_info = struct('start', {}, 'stop', {}, 'max_frame', {}, 'n_transitions', {}, 'is_bar_flash', {});

    for g = 1:(n_gaps - 1)
        epoch_start = gap_ends(g) + 1;
        epoch_stop = gap_starts(g + 1) - 1;

        if epoch_stop <= epoch_start
            continue;
        end

        epoch_frames = f_data(epoch_start:epoch_stop);
        max_frame = max(epoch_frames);

        % Count transitions from 0 to >0 (each flash is one transition)
        above_zero = epoch_frames > 0;
        transitions = sum(diff(above_zero) == 1);

        is_bar_flash = (max_frame <= max_pattern_frame) && (transitions >= 40);

        epoch_info(end+1) = struct('start', epoch_start, 'stop', epoch_stop, ...
            'max_frame', max_frame, 'n_transitions', transitions, 'is_bar_flash', is_bar_flash);
    end

    % Extract bar flash epochs
    bf_epochs = epoch_info([epoch_info.is_bar_flash]);
    n_bf_epochs = numel(bf_epochs);

    % We expect num_reps * num_speeds bar flash epochs
    expected_epochs = num_reps * num_speeds;
    if n_bf_epochs ~= expected_epochs
        warning('parse_bar_flash_data:epochCount', ...
            'Expected %d bar flash epochs, found %d. Check gap detection.', ...
            expected_epochs, n_bf_epochs);
    end

    %% Assign epochs to speeds and reps
    % Bar flash epochs should appear in order: slow_rep1, fast_rep1, slow_rep2, fast_rep2, slow_rep3, fast_rep3
    % Or they may be grouped differently. We'll determine speed from the
    % duration of individual flashes within each epoch.

    % Sort epochs by time
    [~, sort_idx] = sort([bf_epochs.start]);
    bf_epochs = bf_epochs(sort_idx);

    % Determine speed assignment by measuring flash duration in each epoch
    epoch_speeds = zeros(1, n_bf_epochs);
    for e = 1:n_bf_epochs
        ep_frames = f_data(bf_epochs(e).start:bf_epochs(e).stop);
        above_zero = ep_frames > 0;
        transitions_up = find(diff(above_zero) == 1);
        transitions_down = find(diff(above_zero) == -1);

        if numel(transitions_up) >= 2 && numel(transitions_down) >= 2
            % Measure typical flash duration (in samples)
            flash_durations = transitions_down(1:min(5, numel(transitions_down))) - ...
                              transitions_up(1:min(5, numel(transitions_up)));
            median_dur = median(flash_durations);
            % 80ms flash = 800 samples at 10kHz, 14ms flash = 140 samples
            if median_dur > 400
                epoch_speeds(e) = 1; % slow (80ms)
            else
                epoch_speeds(e) = 2; % fast (14ms)
            end
        end
    end

    %% Parse each epoch
    for sp = 1:num_speeds
        data = cell(num_positions, num_orient, num_reps);

        if sp == 1
            gap_between_flashes = 5000; % Half of 1s interval for slow
        else
            gap_between_flashes = 2500; % Half of 0.5s interval for fast
        end

        speed_epochs = bf_epochs(epoch_speeds == sp);
        n_speed_epochs = numel(speed_epochs);

        if n_speed_epochs ~= num_reps
            warning('parse_bar_flash_data:repCount', ...
                'Expected %d reps for speed %d, found %d.', ...
                num_reps, sp, n_speed_epochs);
        end

        for r = 1:min(n_speed_epochs, num_reps)
            data_rep = cell(num_positions, num_orient);

            rng_start = speed_epochs(r).start;
            rng_stop = speed_epochs(r).stop;

            % Find all flash onsets and offsets within this epoch
            ep_frames = f_data(rng_start:rng_stop);
            above_zero = ep_frames > 0;
            flash_ups = find(diff(above_zero) == 1) + rng_start - 1;
            flash_downs = find(diff(above_zero) == -1) + rng_start - 1;

            % Pair up onsets and offsets
            n_flashes = min(numel(flash_ups), numel(flash_downs));

            for f = 1:n_flashes
                st = flash_ups(f);
                nd = flash_downs(f);
                flash_frame_num = max(f_data(st:nd));

                % Ensure we don't go out of bounds
                clip_start = max(1, st - gap_between_flashes);
                clip_end = min(numel(v_data), nd + gap_between_flashes);
                data_flash = v_data(clip_start:clip_end);

                % Store using frame number as linear index into (positions x orientations)
                % Frame 1 = grey background (skip), frames 2-89 = 88 flash stimuli
                % Subtract 1 so frames 2-89 map to linear indices 1-88
                if flash_frame_num >= 2 && flash_frame_num <= n_frames_pattern + 1
                    data_rep{flash_frame_num - 1} = data_flash;
                end
            end

            data(:, :, r) = data_rep(:, :);
        end

        % Create mean across reps
        [rows, cols, reps] = size(data);
        meanData = cell(rows, cols);

        for i = 1:rows
            for j = 1:cols
                ts = data(i, j, :);
                ts = ts(:);

                % Skip if all empty
                non_empty = ~cellfun(@isempty, ts);
                if ~any(non_empty)
                    meanData{i, j} = [];
                    continue;
                end

                maxlen = max(cellfun(@length, ts(non_empty)));
                padded = NaN(maxlen, reps);
                for k = 1:reps
                    if ~isempty(ts{k})
                        y = ts{k};
                        padded(1:length(y), k) = y;
                    end
                end
                meanData{i, j} = mean(padded, 2, 'omitnan');
            end
        end

        if sp == 1
            data_slow = data;
            mean_slow = meanData;
        else
            data_fast = data;
            mean_fast = meanData;
        end
    end

    %% Build debug info for verification
    if nargout >= 5
        debug_info = struct();
        debug_info.all_epochs = epoch_info;
        debug_info.bar_flash_epochs = bf_epochs;
        debug_info.epoch_speeds = epoch_speeds;
        debug_info.gap_starts = gap_starts;
        debug_info.gap_ends = gap_ends;
        debug_info.n_gaps = n_gaps;
        debug_info.n_bar_flash_epochs = n_bf_epochs;
    end

end
