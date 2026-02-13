function [lut_directions, lut_orientations, lut_patterns, lut_functions] = ...
    verify_lut_directions(Tbl, pattern_order, func_order, plot_order)
% VERIFY_LUT_DIRECTIONS  Check that LUT directions match plot_order assumptions.
%
%   [LUT_DIRECTIONS, LUT_ORIENTATIONS, LUT_PATTERNS, LUT_FUNCTIONS] = ...
%       VERIFY_LUT_DIRECTIONS(TBL, PATTERN_ORDER, FUNC_ORDER, PLOT_ORDER)
%   looks up each slow bar stimulus condition in the lookup table and prints
%   a comparison between the angular ordering assumed by plot_order and the
%   actual LUT directions. Warns if mismatches or duplicate directions are
%   found.
%
%   INPUTS:
%     Tbl           - Table from bar_lut.mat with columns: pattern, function,
%                     orientation, direction
%     pattern_order - 1xM array of pattern numbers from currentExp.mat
%     func_order    - 1xM array of function numbers from currentExp.mat
%     plot_order    - 1x16 index mapping from subplot position to data row
%                     (default convention: [1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16])
%
%   OUTPUTS:
%     lut_directions   - 16x1 array of motion directions in degrees, indexed
%                        by data row (i.e., lut_directions(i) is the direction
%                        for bar sweep data row i)
%     lut_orientations - 16x1 array of bar orientations in degrees
%     lut_patterns     - 16x1 array of pattern numbers
%     lut_functions    - 16x1 array of function numbers
%
%   VERIFICATION OUTPUT:
%     Prints a table to the command window comparing each data row's LUT
%     direction against the direction assumed by plot_order. Flags any
%     mismatches and warns about duplicate directions in the LUT.
%
%   See also ANALYZE_SINGLE_EXPERIMENT, FIND_PD_FROM_LUT

    n_dir = 16;

    % Extract slow bar conditions from the full ordering
    slow_func_mask = (func_order == 3 | func_order == 4);
    slow_pattern_nums = pattern_order(slow_func_mask);
    slow_func_nums = func_order(slow_func_mask);

    % Look up each data row in the LUT
    lut_directions   = zeros(n_dir, 1);
    lut_orientations = zeros(n_dir, 1);
    lut_patterns     = zeros(n_dir, 1);
    lut_functions    = zeros(n_dir, 1);

    for i = 1:n_dir
        p = slow_pattern_nums(i);
        f = slow_func_nums(i);
        row_mask = (Tbl.pattern == p) & (Tbl.function == f);
        lut_directions(i)   = Tbl.direction(row_mask);
        lut_orientations(i) = Tbl.orientation(row_mask);
        lut_patterns(i)     = p;
        lut_functions(i)    = f;
    end

    % Compare LUT directions to assumed angles from plot_order
    assumed_angles = linspace(0, 360 - 22.5, n_dir)';

    % fprintf('\n=== LUT-Based Direction Verification ===\n');
    % fprintf('%-8s %-8s %-8s %-12s %-12s %-16s %-6s\n', ...
    %     'DataRow', 'Pattern', 'Func', 'LUT_Dir', 'LUT_Orient', 'Assumed_Dir', 'Match');
    % fprintf('%s\n', repmat('-', 1, 72));

    % % all_match = true;
    % for subplot_idx = 1:n_dir
    %     data_row    = plot_order(subplot_idx);
    %     assumed_dir = assumed_angles(subplot_idx);
    %     actual_dir  = lut_directions(data_row);
    %     actual_orient = lut_orientations(data_row);
    %     match = abs(actual_dir - assumed_dir) < 0.1;
    %     if ~match
    %         all_match = false;
    %     end
    %     fprintf('%-8d %-8d %-8d %-12.1f %-12.1f %-16.1f %-6s\n', ...
    %         data_row, lut_patterns(data_row), lut_functions(data_row), ...
    %         actual_dir, actual_orient, assumed_dir, string(match));
    % end

    % if all_match
    %     fprintf('\nAll directions match plot_order assumptions.\n');
    % else
    %     fprintf('\nWARNING: Mismatches found! Using LUT directions for plotting.\n');
    % end

    % Check for duplicate directions
    unique_dirs = unique(lut_directions);
    if numel(unique_dirs) < n_dir
        fprintf('\nWARNING: LUT contains duplicate directions!\n');
        for d_idx = 1:numel(unique_dirs)
            rows_at_dir = find(lut_directions == unique_dirs(d_idx));
            if numel(rows_at_dir) > 1
                fprintf('  Direction %.1f° appears in data rows: %s\n', ...
                    unique_dirs(d_idx), mat2str(rows_at_dir'));
            end
        end
        missing = setdiff(assumed_angles, lut_directions);
        if ~isempty(missing)
            fprintf('  Missing directions: %s\n', mat2str(missing'));
        end
    end

end
