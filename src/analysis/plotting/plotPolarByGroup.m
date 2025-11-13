function plotPolarByGroup(T, metric_type)
% plotPolarByGroup: Polar scatter of resultant_angle (rad) vs magnitude,
% colored by Strain × Type.
%
% Groups & colors:
%   Strain=="42F06_T4T5_control", Type=="on"  -> light gray
%   Strain=="42F06_T4T5_control", Type=="off" -> black
%   Strain=="42F06_T4T5_ttl",     Type=="on"  -> light pink
%   Strain=="42F06_T4T5_ttl",     Type=="off" -> dark magenta
%
% Input:
%   T : table with columns:
%       - Strain (string/char/cellstr/categorical)
%       - Type   (string/char/cellstr/categorical)
%       - resultant_angle (radians; numeric or cell with numeric/strings)
%       - magnitude       (numeric or cell with numeric/strings)
%
% Example:
%   plotPolarByGroup(T);

    arguments
        T table
        metric_type string
    end

    % ---------- Coerce columns ----------
    T.Strain = coerceToCategorical(T, 'Strain');
    T.Type   = coerceToCategorical(T, 'Type');

    theta = coerceToNumericVector(T.('resultant_angle'), "resultant_angle");
    rho   = coerceToNumericVector(T.(metric_type),       metric_type);

    % Normalize/validate
    theta = wrapTo2Pi(theta);                 % keep angles in [0, 2π)
    valid = isfinite(theta) & isfinite(rho) & (rho >= 0);
    T     = T(valid, :);
    theta = theta(valid);
    rho   = rho(valid);

    if isempty(theta)
        warning('No valid rows to plot (check angles/magnitudes).');
        return;
    end

    % ---------- Prepare groups ----------
    % normalize possible typo: 'tll' -> 'ttl'
    s = lower(string(T.Strain));
    s = strrep(s, 'tll', 'ttl');
    t = lower(string(T.Type));

    isCtrlOn  = (s == "42f06_t4t5_control") & (t == "on");
    isCtrlOff = (s == "42f06_t4t5_control") & (t == "off");
    isTtlOn   = (s == "42f06_t4t5_ttl")     & (t == "on");
    isTtlOff  = (s == "42f06_t4t5_ttl")     & (t == "off");

    % Colors
    colCtrlOn  = [0.80 0.80 0.80]; % light gray
    colCtrlOff = [0.00 0.00 0.00]; % black
    colTtlOn   = [1.00 0.75 0.85]; % light pink
    colTtlOff  = [0.60 0.00 0.60]; % dark magenta

    % Marker size
    msz = 100;

    % ---------- Plot ----------
    figure;
    pax = polaraxes; hold(pax, 'on');

    % Use polarplot with markers (no connecting line)
    h = gobjects(4,1);
    if any(isCtrlOn)
        h(1) = polarplot(pax, theta(isCtrlOn),  rho(isCtrlOn),  'o', ...
            'MarkerFaceColor', colCtrlOn,  'MarkerEdgeColor','k', ...
            'LineStyle','none', 'MarkerSize', sqrt(msz));
    end
    if any(isCtrlOff)
        h(2) = polarplot(pax, theta(isCtrlOff), rho(isCtrlOff), 'o', ...
            'MarkerFaceColor', colCtrlOff, 'MarkerEdgeColor','k', ...
            'LineStyle','none', 'MarkerSize', sqrt(msz));
    end
    if any(isTtlOn)
        h(3) = polarplot(pax, theta(isTtlOn),   rho(isTtlOn),   'o', ...
            'MarkerFaceColor', colTtlOn,   'MarkerEdgeColor','k', ...
            'LineStyle','none', 'MarkerSize', sqrt(msz));
    end
    if any(isTtlOff)
        h(4) = polarplot(pax, theta(isTtlOff),  rho(isTtlOff),  'o', ...
            'MarkerFaceColor', colTtlOff,  'MarkerEdgeColor','k', ...
            'LineStyle','none', 'MarkerSize', sqrt(msz));
    end

    % Cosmetics
    title(pax, strcat("Preferred Direction vs ", strrep(metric_type, "_", "-")));
    thetaticks(0:22.5:337.5);
    rtickangle(pax, 0);
    grid(pax, 'on');

    % Legend (only for groups that exist)
    legNames = {'control-on','control-off','ttl-on','ttl-off'};
    present  = arrayfun(@(hh) ~isempty(hh) && isvalid(hh), h);
    if any(present)
        legend(h(present), legNames(present), 'Location','bestoutside');
    end
end

% ===================== Helpers =====================

function C = coerceToCategorical(T, varName)
    x = T.(varName);
    if iscategorical(x)
        C = x; return;
    end
    if iscell(x) || ischar(x)
        x = string(x);
    end
    C = categorical(x);
end

function v = coerceToNumericVector(x, colName)
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
                error('Unsupported cell content at row %d in column "%s".', k, colName);
            end
        end
        return;
    end
    if isstring(x) || ischar(x) || iscellstr(x)
        v = str2double(string(x(:)));
        return;
    end
    error('Unsupported data type for column "%s".', colName);
end

function ang = wrapTo2Pi(ang)
    % Wrap angles (radians) to [0, 2π)
    ang = mod(ang, 2*pi);
end
