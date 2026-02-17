function [meanData_out, data_out, respMap, maxCol] = shiftMaxColumnTo5(meanData, data)
% SHIFTMAXCOLUMNTO5  Circularly shift bar flash data so preferred direction is at column 5.
%
%   [MEANDATA_OUT, DATA_OUT, RESPMAP, MAXCOL] = SHIFTMAXCOLUMNTO5(MEANDATA, DATA)
%   identifies the orientation column (1-8) containing the strongest bar
%   flash response, then circularly shifts all columns so that this
%   preferred direction lands at column 5. This aligns the data for
%   consistent cross-experiment comparison and plotting.
%
%   INPUTS:
%     meanData - 11x8 cell array
%                Mean bar flash timeseries (positions x orientations).
%                Each cell contains a numeric vector of voltage values (mV).
%     data     - 11x8x3 cell array
%                Individual rep timeseries (positions x orientations x reps).
%                Same layout as meanData with a 3rd dimension for reps.
%
%   OUTPUTS:
%     meanData_out - 11x8 cell array
%                    meanData with columns circularly shifted so the column
%                    with the global maximum response is at column 5.
%     data_out     - 11x8x3 cell array
%                    data shifted identically along dimension 2.
%     respMap      - 11x8 double
%                    Response magnitude map. Each entry is the 98th
%                    percentile of the timeseries in the window from 40%
%                    to 80% of the signal length:
%                      resp = prctile(d(ceil(0.4*n):ceil(0.8*n)), 98)
%     maxCol       - double (1-8)
%                    Original column index where the global max occurred
%                    (before shifting). NaN if all responses are NaN.
%
%   NOTES:
%     - If a cell is empty, its response is set to NaN.
%     - If the 40%-80% window collapses (very short signals), falls back
%       to prctile(d, 98) over the full signal.
%     - Columns wrap around: shifting preserves the relative circular order
%       of the 8 orientations.
%
%   See also PLOT_BAR_FLASH_DATA, PARSE_BAR_FLASH_DATA

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
