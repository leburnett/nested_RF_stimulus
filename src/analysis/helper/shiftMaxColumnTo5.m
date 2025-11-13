function [meanData_out, data_out, respMap, maxCol] = shiftMaxColumnTo5(meanData, data)
% shiftMaxColumnTo5
%   Inputs:
%     meanData : [11 x 8] cell array (each cell is a numeric vector)
%     data     : [11 x 8 x 3] cell array (same layout; 3 reps in dim 3)
%
%   Outputs:
%     meanData_out : meanData with columns circularly shifted so that the
%                    column containing the global maximum response is at col 5
%     data_out     : data shifted identically along dimension 2
%     respMap      : [11 x 8] double, per-cell response using the 98th percentile
%     maxCol       : original column index (1..8) where the global max occurred
%
%   Response per timeseries d (vector of length n):
%     resp = prctile(d( ceil(0.4*n) : ceil(0.8*n) ), 98)
%
%   Notes:
%     * If a cell is empty or too short to form a range (start >= stop),
%       falls back to prctile(d, 98) or NaN if empty.
%     * Columns are circularly shifted (rotated), preserving relative layout.

    % Basic checks
    assert(iscell(meanData) && isequal(size(meanData), [11, 8]), ...
        'meanData must be a [11 x 8] cell array.');
    assert(iscell(data) && isequal(size(data), [11, 8, 3]), ...
        'data must be a [11 x 8 x 3] cell array.');

    % ----- Build response map (11x8) -----
    respMap = nan(11, 8);
    for i = 1:11
        for j = 1:8
            d = meanData{i, j};
            if isempty(d)
                respMap(i, j) = NaN;
                continue;
            end
            if ~isvector(d) || ~isnumeric(d)
                error('Cell {%d,%d} in meanData must be a numeric vector.', i, j);
            end
            d = d(:);
            n = numel(d);

            s = ceil(0.4*n);
            e = ceil(0.8*n);
            % Clamp indices to valid range
            s = max(1, min(s, n));
            e = max(1, min(e, n));

            if s < e
                respMap(i, j) = prctile(d(s:e), 98);
            else
                % Fallback if window collapses (very short signals)
                respMap(i, j) = prctile(d, 98);
            end
        end
    end

    % ----- Find the column with the global maximum -----
    [~, linIdx] = max(respMap(:), [], 'omitnan');
    if isempty(linIdx)
        % All NaNs: nothing to shift
        meanData_out = meanData;
        data_out     = data;
        maxCol       = NaN;
        warning('All responses are NaN. No shift applied.');
        return;
    end
    [~, maxCol] = ind2sub(size(respMap), linIdx);  % (row, col) -> take col

    % ----- Circularly shift columns so maxCol moves to column 5 -----
    targetCol = 5;
    shift = targetCol - maxCol;   % positive = shift right, negative = left

    if shift ~= 0
        meanData_out = circshift(meanData, [0, shift]);     % shift columns (dim 2)
        data_out     = circshift(data,     [0, shift, 0]);  % shift columns (dim 2), keep reps
    else
        meanData_out = meanData;
        data_out     = data;
    end
end
