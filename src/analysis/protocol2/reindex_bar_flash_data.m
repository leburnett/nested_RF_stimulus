function [corrected, orient_labels_deg] = reindex_bar_flash_data(raw_data, stim_meta)
% REINDEX_BAR_FLASH_DATA  Reorder bar flash data to match true orientations.
%
%   [CORRECTED, ORIENT_LABELS_DEG] = REINDEX_BAR_FLASH_DATA(RAW_DATA, STIM_META)
%   takes the 11x8 cell array from parse_bar_flash_data and returns a
%   corrected version where columns correspond to empirically determined
%   orientations and rows to positions, using the metadata from
%   get_stimulus_metadata.
%
%   INPUTS:
%     raw_data  - 11x8 cell array (e.g., mean_slow or mean_fast from
%                 parse_bar_flash_data). After the frame-1 subtraction fix,
%                 linear index i corresponds to flash frame i+1 in the
%                 pattern, which is the i-th flash stimulus.
%     stim_meta - Output of GET_STIMULUS_METADATA, containing the
%                 frame-to-orientation and frame-to-position mappings.
%
%   OUTPUTS:
%     corrected        - 11x8 cell array where column c = orientation c
%                        and row r = position r, properly aligned using
%                        the empirical orientation labels.
%     orient_labels_deg - 8x1 array of orientation angles (degrees) for
%                         each column in the corrected array.
%
%   INDEXING:
%     parse_bar_flash_data stores data_rep{flash_frame_num - 1}, where
%     flash_frame_num ranges from 2 to 89. So:
%       - Linear index 1 = flash frame 2 = flash stimulus 1
%       - Linear index 88 = flash frame 89 = flash stimulus 88
%     The frame_to_orientation_idx and frame_to_position_idx arrays from
%     stim_meta map flash index (1-88) to the actual orientation and
%     position indices.
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

            % linear_idx corresponds to flash stimulus number (1-88)
            % which maps to flash frame (linear_idx + 1) in the pattern
            if linear_idx >= 1 && linear_idx <= numel(frame_to_orient)
                actual_orient = frame_to_orient(linear_idx);
                actual_pos = frame_to_pos(linear_idx);

                if actual_orient >= 1 && actual_orient <= num_orientations && ...
                   actual_pos >= 1 && actual_pos <= num_positions
                    corrected{actual_pos, actual_orient} = raw_data{row, col};
                end
            end
        end
    end

end
