function data_ordered = align_data_by_seq_angles(data)
% ALIGN_DATA_BY_SEQ_ANGLES  Reorder bar responses into sequential angular order.
%
%   DATA_ORDERED = ALIGN_DATA_BY_SEQ_ANGLES(DATA) rearranges timeseries
%   data from forward/backward presentation order into consecutive
%   angular order (0, 22.5, 45, ... degrees).
%
%   INPUT:
%     data - 32x4 cell array from PARSE_BAR_DATA containing:
%            Rows 1-16: Slow bar stimuli (forward/backward pairs)
%            Rows 17-32: Fast bar stimuli (forward/backward pairs)
%            Columns 1-3: Individual repetitions
%            Column 4: Mean across repetitions
%
%   OUTPUT:
%     data_ordered - 32x4 cell array with same structure but reordered
%                    so consecutive rows represent consecutive angles
%
%   REORDERING:
%     Original order: 0, 180, 22.5, 202.5, 45, 225, ... (forward/backward pairs)
%     New order: 0, 22.5, 45, 67.5, ... 180, 202.5, ... (sequential angles)
%
%     plot_order = [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16]
%     Maps odd indices (forward) first, then even indices (backward)
%
%   PURPOSE:
%     Sequential angular ordering facilitates visualization of the
%     directional tuning curve and comparison across cells.
%
%   See also PARSE_BAR_DATA, FIND_PD_AND_ORDER_IDX
    
    plot_order= [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16];

    data_ordered = cell(32, 4);

    n_dir_and_sp = height(data);
    n_speeds = 2;
    n_dir = n_dir_and_sp/n_speeds;

    for k = 1:n_dir_and_sp
   
        if k >n_dir % 56 dps
            data_idx = plot_order(k-n_dir);
            data_idx = data_idx+n_dir;
        else % 28 dps
            data_idx = plot_order(k);
        end 
        data_ordered(k, :) = data(data_idx, :);
    end 

end 




