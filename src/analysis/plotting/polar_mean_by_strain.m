function polar_mean_by_strain(T, speed)
% POLAR_MEAN_BY_STRAIN  Plot mean polar responses (d_slow(:,2)) grouped by Strain.
% Groups are defined by whether T.Strain contains:
%   - "control"  -> black (mean) with light gray SEM shading
%   - "ttl"      -> red (mean)   with light red SEM shading
%
% Input:
%   T table with:
%       - T.Strain : string/char/categorical (substrings "control" / "ttl")
%       - T.d_slow : cell, each a 16x2 numeric [theta(rad), response]

    validate_table(T);

    % Use a canonical theta from the first valid entry
    theta = get_theta_from_table(T);
    if isempty(theta)
        error('Could not find a valid 16x2 d_slow to read angles from.');
    end

    % Compute group stats (mean and SEM across cells) by Strain substrings
    stats_control = get_group_stats_by_strain(T, 'control', speed);
    stats_ttl     = get_group_stats_by_strain(T, 'ttl', speed);

    % Prepare polar axes (for grid + mean lines)
    figure;
    axPolar = polaraxes; 
    hold(axPolar, 'on');

    % Colors
    col_control_line = [0 0 0];          % black
    col_control_fill = [0.80 0.80 0.80]; % light gray
    col_ttl_line     = [1 0 0];          % red
    col_ttl_fill     = [1 0.70 0.70];    % light red

    % Determine a reasonable r-limit (based on all available data)
    rmax = estimate_rmax(stats_control, stats_ttl);
    if ~isnan(rmax) && rmax > 0
        % rlim(axPolar, [0 rmax]);
         rlim(axPolar, [0 32]);
    end

    % Create an overlay Cartesian axes for patches
    axFill = axes('Position', axPolar.Position, 'Color', 'none', ...
                  'XColor','none','YColor','none', 'HitTest','off');
    axis(axFill, 'equal');
    hold(axFill, 'on');

    % Match limits of overlay to polar r-limits
    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);

    % CONTROL shaded band + mean
    hCtrl = gobjects(1,1);
    if ~isempty(stats_control.mean)
        hCtrl = plot_with_shade(axPolar, axFill, theta, ...
                                stats_control.mean, stats_control.sem, ...
                                col_control_line, col_control_fill, 0.30, 2);
    else
        warning('No valid "control" strain entries found.');
    end

    % TTL shaded band + mean
    hTtl = gobjects(1,1);
    if ~isempty(stats_ttl.mean)
        hTtl = plot_with_shade(axPolar, axFill, theta, ...
                               stats_ttl.mean, stats_ttl.sem, ...
                               col_ttl_line, col_ttl_fill, 0.30, 2);
    else
        warning('No valid "ttl" strain entries found.');
    end

    % Keep polar axes on top so grid/ticks are visible
    uistack(axPolar, 'top');

    % Legend (lines only)
    labs = {};
    hs   = [];
    if ~isempty(stats_control.mean)
        hs(end+1)   = hCtrl;
        labs{end+1} = sprintf('control (n=%d)', stats_control.n);
    end
    if ~isempty(stats_ttl.mean)
        hs(end+1)   = hTtl;
        labs{end+1} = sprintf('ttl (n=%d)', stats_ttl.n);
    end
    if ~isempty(hs)
        legend(axPolar, hs, labs, 'Location', 'best');
    end
    % title(axPolar, 'Mean Polar Response by Strain');
end

% ---------- helpers ----------
function validate_table(T)
    if ~istable(T)
        error('Input must be a table.');
    end
    needVars = {'Strain','d_slow'};
    missing = setdiff(needVars, T.Properties.VariableNames);
    if ~isempty(missing)
        error('Table is missing required variables: %s', strjoin(missing, ', '));
    end
end

function theta = get_theta_from_table(T)
    theta = [];
    for i = 1:height(T)
        d = T.d_slow{i};
        if isnumeric(d) && isequal(size(d), [16 2])
            theta = d(:,1);
            return;
        end
    end
end

function stats = get_group_stats_by_strain(T, substr, speed)
    % Case-insensitive substring match on T.Strain (accepts string/char/categorical)
    S = T.Strain;
    s = string(S); % works for string, char, categorical
    mask = contains(lower(s), lower(string(substr)));

    idx = find(mask);
    vals = [];
    for k = idx'
        if speed == "slow"
            d = T.d_slow{k};
        elseif speed == "fast"
            d = T.d_fast{k};
        elseif speed == "vfast"
            d = T.d_vfast{k};
        end 
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

    % SEM band (change to "sd" if you want SD instead)
    stats.sem  = sd ./ sqrt(stats.n);
end

function rmax = estimate_rmax(sc, st)
    allVals = [];
    if ~isempty(sc) && isfield(sc,'mean') && ~isempty(sc.mean)
        allVals = [allVals; sc.mean(:) + sc.sem(:)]; 
    end
    if ~isempty(st) && isfield(st,'mean') && ~isempty(st.mean)
        allVals = [allVals; st.mean(:) + st.sem(:)]; 
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

function hLine = plot_with_shade(axPolar, axFill, theta, meanVals, bandVals, lineColor, fillColor, alphaVal, lw)
    % Close the curve for plotting as a loop
    th = theta(:);
    m  = meanVals(:);
    b  = bandVals(:);

    % Ensure non-negative radii
    upper = max(m + b, 0);
    lower = max(m - b, 0);

    % Build closed polygon in polar -> convert to cartesian for patch (on axFill)
    th_poly = [th; th(1); flipud(th); th(1)];
    r_poly  = [upper; upper(1); flipud(lower); lower(1)];
    [x_poly, y_poly] = pol2cart(th_poly, r_poly);

    patch('XData', x_poly, 'YData', y_poly, ...
          'FaceColor', fillColor, 'FaceAlpha', alphaVal, ...
          'EdgeColor', 'none', 'Parent', axFill, 'HitTest','off');

    % Draw the mean line on polar axes (closed for aesthetics)
    th_line = [th; th(1)];
    r_line  = [m;  m(1)];
    hLine   = polarplot(axPolar, th_line, r_line, '-', 'Color', lineColor, 'LineWidth', lw);
end
