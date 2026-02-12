function [corrected, orient_labels_deg] = reindex_bar_flash_data(raw_data, stim_meta)
% REINDEX_BAR_FLASH_DATA  Correct off-by-one indexing in bar flash data.
%
%   [CORRECTED, ORIENT_LABELS_DEG] = REINDEX_BAR_FLASH_DATA(RAW_DATA, STIM_META)
%   takes the raw 11x8 cell array from parse_bar_flash_data and returns a
%   corrected version where columns truly correspond to orientations and
%   rows to positions.
%
%   INPUTS:
%     raw_data  - 11x8 cell array (e.g., mean_slow or mean_fast from
%                 parse_bar_flash_data). Due to frame indexing, the mapping
%                 from (row,col) to (position, orientation) is shifted.
%     stim_meta - Output of GET_STIMULUS_METADATA, containing the
%                 frame-to-orientation and frame-to-position mappings.
%
%   OUTPUTS:
%     corrected        - 11x8 cell array where column c = orientation c
%                        and row r = position r, properly aligned.
%     orient_labels_deg - 8x1 array of orientation angles (degrees) for
%                         each column in the corrected array.
%
%   THE PROBLEM:
%     In parse_bar_flash_data, data_rep{flash_frame_num} stores data using
%     the frame number as a linear index into an 11x8 cell. Since frame 1
%     is the grey background, flash frames start at frame 2. MATLAB's
%     column-major ordering means:
%       - Linear index 1 = (row 1, col 1) → empty (grey frame)
%       - Linear index 2 = (row 2, col 1) → position 1, orientation 1
%       - Linear index 12 = (row 1, col 2) → position 11, orientation 1
%       - Linear index 13 = (row 2, col 2) → position 1, orientation 2
%     This causes a systematic shift in the data layout.
%
%   See also GET_STIMULUS_METADATA, PARSE_BAR_FLASH_DATA

    num_positions = 11;
    num_orientations = 8;

    frame_to_orient = stim_meta.bar_flash.frame_to_orientation_idx;
    frame_to_pos = stim_meta.bar_flash.frame_to_position_idx;
    orient_labels_deg = stim_meta.bar_flash.orientation_deg;

    corrected = cell(num_positions, num_orientations);

    % Iterate over all cells in the raw data
    for col = 1:num_orientations
        for row = 1:num_positions
            linear_idx = (col - 1) * num_positions + row;

            if linear_idx == 1
                % Frame 1 = grey, this cell is empty
                continue;
            end

            % The frame number that was stored at this linear index
            frame_num = linear_idx;

            % Look up what this frame actually corresponds to
            fi = frame_num - 1; % Index into frame_to_orient (skips frame 1)

            if fi >= 1 && fi <= numel(frame_to_orient)
                actual_orient = frame_to_orient(fi);
                actual_pos = frame_to_pos(fi);

                if actual_orient >= 1 && actual_orient <= num_orientations && ...
                   actual_pos >= 1 && actual_pos <= num_positions
                    corrected{actual_pos, actual_orient} = raw_data{row, col};
                end
            end
        end
    end

    % Handle the last flash (position 11, orientation 8) which maps to
    % frame 89 = linear index 89. This is beyond the 88-element cell array
    % (11x8=88), so it was stored via the position function mapping frame 1
    % to frame 89 (wrapping). In practice, this flash may be missing from
    % the raw data. Check if (row 1, col 1) has data that should be
    % position 11 of orientation 8.
    % Note: frame 89 would wrap to linear index 89 which is out of bounds
    % for an 11x8 cell. The position function replaces frame 1 with 89,
    % meaning this flash is effectively lost. Flag this.
    if isempty(corrected{num_positions, num_orientations})
        warning('reindex_bar_flash_data:missingFlash', ...
            'Position %d of orientation %d is empty (frame 89 is out of bounds for 11x8 cell).', ...
            num_positions, num_orientations);
    end

end
