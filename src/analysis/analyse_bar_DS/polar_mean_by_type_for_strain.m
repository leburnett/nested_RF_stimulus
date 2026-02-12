function polar_mean_by_type_for_strain(T, strainSel)
% POLAR_MEAN_BY_TYPE_FOR_STRAIN  Single polar plot of mean responses (d_slow(:,2))
% for ONE selected Strain, overlaying two groups by Type: "on" and "off".
%
% Usage:
%   polar_mean_by_type_for_strain(T, 'control')
%   polar_mean_by_type_for_strain(T, 'ttl')
%
% Inputs:
%   T table with:
%       - T.Strain : string/char/categorical (contains 'control' or 'ttl')
%       - T.Type   : string/char/categorical ('on' or 'off')
%       - T.d_slow : cell, each a 16x2 numeric [theta(rad), response]
%   strainSel : char/string to select strain ('control' or 'ttl')
%
% Colors:
%   on  -> dark green   (line) with light green SEM shading
%   off -> black        (line) with light gray SEM shading

    % ---- validation ----
    validate_table_local(T);
    if nargin < 2 || isempty(strainSel)
        error('Provide strainSel as ''control'' or ''ttl''.');
    end

    % canonical theta from first valid entry
    theta = get_theta_from_table_local(T);
    if isempty(theta)
        error('Could not find a valid 16x2 d_slow to read angles from.');
    end

    % filter to selected strain (case-insensitive substring)
    s = string(T.Strain);
    maskStrain = contains(lower(s), lower(string(strainSel)));
    Tsub = T(maskStrain, :);
    if isempty(Tsub)
        error('No rows found for Strain containing "%s".', string(strainSel));
    end

    % stats by Type within this strain
    stats_on  = get_group_stats_by_type_local(Tsub, 'on');
    stats_off = get_group_stats_by_type_local(Tsub, 'off');

    % colors
    col_on_line  = [0.70 0.95 0.70];       % dark green
    col_on_fill  = [0.70 0.88 0.70]; % light green
    col_off_line = [0 0.35 0];          % black
    col_off_fill = [0.55 0.7 0.55]; % light gray

    % r-limit based on both groups
    rmax = estimate_rmax_local(stats_on, stats_off);

    % ---- single polar axes + overlay for shaded patches ----
    figure;
    axPolar = polaraxes;
    hold(axPolar, 'on');
    title(axPolar, sprintf('Strain: %s  (Type: on/off)', string(strainSel)));

    if ~isnan(rmax) && rmax > 0
        rlim(axPolar, [0 rmax]);
    end

    % overlay cartesian axes for patches
    axFill = axes('Position', axPolar.Position, 'Color', 'none', ...
                  'XColor','none','YColor','none', 'HitTest','off');
    axis(axFill, 'equal'); hold(axFill, 'on');
    match_overlay_limits(axPolar, axFill);

    % ---- plot ON (shade + mean line) ----
    % hOn = gobjects(1,1);
    if ~isempty(stats_on.mean)
        plot_with_shade_local(axPolar, axFill, theta, ...
                                    stats_on.mean, stats_on.sem, ...
                                    col_on_line, col_on_fill, 0.30, 2);
    else
        warning('No valid "on" entries found for "%s".', string(strainSel));
    end

    % % ---- plot OFF (shade + mean line) ----
    % hOff = gobjects(1,1);
    if ~isempty(stats_off.mean)
        plot_with_shade_local(axPolar, axFill, theta, ...
                                     stats_off.mean, stats_off.sem, ...
                                     col_off_line, col_off_fill, 0.30, 2);
    else
        warning('No valid "off" entries found for "%s".', string(strainSel));
    end

    % keep polar grid/ticks on top
    % uistack(axPolar, 'top');

    % legend for whichever groups exist
    hs = []; labs = {};
    if ~isempty(stats_on.mean)
        % hs(end+1)   = hOn;
        labs{end+1} = sprintf('on (n=%d)', stats_on.n);
    end
    if ~isempty(stats_off.mean)
        % hs(end+1)   = hOff;
        labs{end+1} = sprintf('off (n=%d)', stats_off.n);
    end
    if ~isempty(hs)
        legend(axPolar, hs, labs, 'Location', 'best');
    end
end

% ---------- helpers ----------
function validate_table_local(T)
    if ~istable(T)
        error('Input must be a table.');
    end
    needVars = {'Strain','Type','d_slow'};
    missing = setdiff(needVars, T.Properties.VariableNames);
    if ~isempty(missing)
        error('Table is missing required variables: %s', strjoin(missing, ', '));
    end
end

function theta = get_theta_from_table_local(T)
    theta = [];
    for i = 1:height(T)
        if ~isempty(T.d_slow{i})
            d = T.d_slow{i};
            if isnumeric(d) && isequal(size(d), [16 2])
                theta = d(:,1);
                return;
            end
        end
    end
end

function stats = get_group_stats_by_type_local(T, typeVal)
    TT = string(T.Type);
    maskType = strcmpi(strtrim(TT), string(typeVal));

    idx = find(maskType);
    vals = [];
    for k = idx'
        d = T.d_slow{k};
        if isnumeric(d) && isequal(size(d), [16 2])
            vals = [vals, d(:,2)];
        end
    end

    if isempty(vals)
        stats.mean = [];
        stats.sem  = [];
        stats.n    = 0;
        return;
    end

    stats.n    = size(vals, 2);
    stats.mean = mean(vals, 2, 'omitnan');
    sd         = std(vals, 0, 2, 'omitnan');
    stats.sem  = sd ./ sqrt(stats.n);
end

function rmax = estimate_rmax_local(sa, sb)
    allVals = [];
    if ~isempty(sa) && isfield(sa,'mean') && ~isempty(sa.mean)
        allVals = [allVals; sa.mean(:) + sa.sem(:)];
    end
    if ~isempty(sb) && isfield(sb,'mean') && ~isempty(sb.mean)
        allVals = [allVals; sb.mean(:) + sb.sem(:)];
    end
    if isempty(allVals)
        rmax = NaN;
    else
        rmax = max(allVals, [], 'omitnan');
        if ~isfinite(rmax) || rmax <= 0
            rmax = NaN;
        end
    end
end

function match_overlay_limits(axPolar, axFill)
    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);
end

function plot_with_shade_local(axPolar, axFill, theta, meanVals, bandVals, lineColor, fillColor, alphaVal, lw)
    th = theta(:);
    m  = meanVals(:);
    b  = bandVals(:);

    % non-negative radii
    upper = max(m + b, 0);
    lower = max(m - b, 0);

    % closed polygon for SEM patch on overlay axis
    th_poly = [th; th(1); flipud(th); th(1)];
    r_poly  = [upper; upper(1); flipud(lower); lower(1)];
    [x_poly, y_poly] = pol2cart(th_poly, r_poly);
    patch('XData', x_poly, 'YData', y_poly, ...
          'FaceColor', fillColor, 'FaceAlpha', alphaVal, ...
          'EdgeColor', 'none', 'Parent', axFill, 'HitTest','off');

    % mean line on polar axis (closed)
    th_line = [th; th(1)];
    r_line  = [m;  m(1)];
    polarplot(axPolar, th_line, r_line, '-', 'Color', lineColor, 'LineWidth', lw);
end
