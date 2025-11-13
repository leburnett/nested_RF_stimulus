function T = addOrthoMetrics(T)
% addOrthoMetrics: Adds v_ortho_slow, v_ortho_v_max_slow, and angGroup.
%
%   v_ortho_slow       = mean([v_ortho1_slow, v_ortho2_slow], 2, 'omitnan')
%   v_ortho_v_max_slow = v_ortho_slow ./ v_max_slow
%
%   angGroup is assigned by resultant_angle (radians) into 4 groups:
%     Group 1: [5.8905, 2π) U [0, 1.1781)
%     Group 2: [1.1781, 2.7489)
%     Group 3: [2.7489, 4.3197)
%     Group 4: [4.3197, 5.8905)
%
% Edge handling (inclusive lower bounds):
%   theta == 5.8905 -> 1
%   theta == 1.1781 -> 2
%   theta == 2.7489 -> 3
%   theta == 4.3197 -> 4
%
% Usage:
%   T = addOrthoMetrics(T);

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
