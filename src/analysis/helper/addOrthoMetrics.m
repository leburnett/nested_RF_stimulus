function T = addOrthoMetrics(T)
% ADDORTHOMETRICS  Add orthogonal response metrics and angle grouping to results table.
%
%   T = ADDORTHOMETRICS(T) computes derived metrics from bar sweep
%   direction selectivity data and appends them as new columns to the
%   input table.
%
%   INPUT:
%     T - table
%         Must contain the following columns (numeric, cell, or string):
%           .v_ortho1_slow    - response in orthogonal direction 1 (slow)
%           .v_ortho2_slow    - response in orthogonal direction 2 (slow)
%           .v_max_slow       - maximum (preferred direction) response (slow)
%           .v_null_slow      - null direction response (slow)
%           .resultant_angle  - preferred direction angle in radians
%
%   OUTPUT:
%     T - table (modified in-place with 4 new columns):
%         .v_ortho_slow       - double, mean of orthogonal responses:
%                               mean([v_ortho1_slow, v_ortho2_slow], 'omitnan')
%         .v_ortho_v_max_slow - double, orthogonal-to-max ratio:
%                               v_ortho_slow ./ v_max_slow
%         .v_null_v_max_slow  - double, null-to-max ratio:
%                               v_null_slow ./ v_max_slow
%         .angGroup           - double (1-4), angular quadrant grouping based
%                               on resultant_angle wrapped to [0, 2*pi):
%                                 1: [5.105, 2*pi) U [0, 0.393)   (~E/W axis)
%                                 2: [0.393, 1.963)                (~NE/SW)
%                                 3: [1.963, 3.534)                (~N/S axis)
%                                 4: [3.534, 5.105)                (~NW/SE)
%
%   EXAMPLE:
%     T = addOrthoMetrics(T);
%
%   See also COMPUTE_BAR_RESPONSE_METRICS, PLOTGROUPEDBOX, RUNGROUPEDSTATS

    % --- Coerce numeric inputs (supports numeric, cell, numeric strings) ---
    v1   = coerceToNumericVector(T.v_ortho1_slow, "v_ortho1_slow");
    v2   = coerceToNumericVector(T.v_ortho2_slow, "v_ortho2_slow");
    vmax = coerceToNumericVector(T.v_max_slow,    "v_max_slow");
    vnull = coerceToNumericVector(T.v_null_slow,   "v_null_slow");

    % --- Compute metrics ---
    v_ortho_slow       = mean([v1, v2], 2, 'omitnan');
    v_ortho_v_max_slow = v_ortho_slow ./ vmax;
    v_null_v_max_slow = vnull./ vmax;

    % --- Add metrics to table ---
    T.v_ortho_slow        = v_ortho_slow;
    T.v_ortho_v_max_slow  = v_ortho_v_max_slow;
    T.v_null_v_max_slow = v_null_v_max_slow;

    % --- Angle grouping ---
    theta = coerceToNumericVector(T.resultant_angle, "resultant_angle");
    theta = wrapTo2Pi(theta);  % [0, 2π)

    % Boundaries (radians)
    b1 = 0.3926991; %
    b2 = 1.9634954; % 
    b3 = 3.5342917; %
    b4 = 5.1050881; %

    angGroup = nan(size(theta));  % default NaN for missing/invalid

    % Group 1: [b4, 2π) U [0, b1)
    mask1 = (theta >= b4) | (theta < b1);
    % Group 2: [b1, b2)
    mask2 = (theta >= b1) & (theta < b2);
    % Group 3: [b2, b3)
    mask3 = (theta >= b2) & (theta < b3);
    % Group 4: [b3, b4)
    mask4 = (theta >= b3) & (theta < b4);

    angGroup(mask1) = 1;
    angGroup(mask2) = 2;
    angGroup(mask3) = 3;
    angGroup(mask4) = 4;

    T.angGroup = angGroup;
end

% ========= Helpers =========
function v = coerceToNumericVector(x, columnName)
    if isnumeric(x)
        v = double(x(:));
        return;
    end
    if iscell(x)
        v = nan(numel(x),1);
        for k = 1:numel(x)
            ck = x{k};
            if isempty(ck)
                v(k) = NaN;
            elseif isnumeric(ck) && isscalar(ck)
                v(k) = double(ck);
            elseif ischar(ck) || (isstring(ck) && isscalar(ck))
                v(k) = str2double(string(ck));
            else
                error('Unsupported cell content at row %d in column "%s".', k, columnName);
            end
        end
        return;
    end
    if isstring(x) || ischar(x) || iscellstr(x)
        v = str2double(string(x(:)));
        return;
    end
    error('Unsupported data type for column "%s".', columnName);
end

function ang = wrapTo2Pi(ang)
    ang = mod(ang, 2*pi);
end
