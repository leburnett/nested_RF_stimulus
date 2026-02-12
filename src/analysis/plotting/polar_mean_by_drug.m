function polar_mean_by_drug(T, speed)
% POLAR_MEAN_BY_DRUG  Plot mean polar responses grouped by drug/mins_post.
%
% Groups:
%   - drug == 0                      -> black (mean) with light gray SEM shading
%   - drug == 1 & mins_post == 15    -> medium red (mean) with light red SEM shading
%   - drug == 1 & mins_post == 35    -> dark red (mean) with darker red SEM shading
%
% Input:
%   T table with:
%       - T.drug      : 0 or 1
%       - T.mins_post : 15 or 35 (for drug == 1)
%       - T.d_slow    : cell, each a 16x2 numeric [theta(rad), response]
%       - (optionally) T.d_fast, T.d_vfast if you use those speeds

    validate_table(T);

    % Use a canonical theta from the first valid entry
    theta = get_theta_from_table(T);
    if isempty(theta)
        error('Could not find a valid 16x2 d_slow to read angles from.');
    end

    % --- Define group masks ---
    sDrug = T.drug;          % numeric or logical
    sMins = T.mins_post;     % numeric

    mask_noDrug     = (sDrug == 0);
    mask_drug_15    = (sDrug == 1) & (sMins == 15);
    mask_drug_35    = (sDrug == 1) & (sMins == 35);

    % Compute group stats (mean and SEM across cells) for each mask
    stats_noDrug  = get_group_stats_by_mask(T, mask_noDrug,  speed);
    stats_15min   = get_group_stats_by_mask(T, mask_drug_15, speed);
    stats_35min   = get_group_stats_by_mask(T, mask_drug_35, speed);

    % Prepare polar axes (for grid + mean lines)
    figure;
    axPolar = polaraxes; 
    thetaticks([])
    hold(axPolar, 'on');

    % Colors
    col_noDrug_line = [0 0 0];          % black
    col_noDrug_fill = [0.80 0.80 0.80]; % light gray

    col_15_line     = [1.0 0.3 0.3];    % medium red
    col_15_fill     = [1.0 0.75 0.75];  % light red

    col_35_line     = [0.7 0   0];      % dark red
    col_35_fill     = [1.0 0.5 0.5];    % darker red-ish fill

    % Determine a reasonable r-limit (based on all available data)
    rmax = estimate_rmax(stats_noDrug, stats_15min, stats_35min);
    if ~isnan(rmax) && rmax > 0
        % rlim(axPolar, [0 rmax]);
        rlim(axPolar, [0 35]);  % keep your fixed limit if you prefer
    end

    % Create an overlay Cartesian axes for patches
    axFill = axes('Position', axPolar.Position, 'Color', 'none', ...
                  'XColor','none','YColor','none', 'HitTest','off');
    axis(axFill, 'equal');
    hold(axFill, 'on');

    % Match limits of overlay to polar r-limits
    rL = rlim(axPolar);
    set(axFill, 'XLim', [-rL(2) rL(2)], 'YLim', [-rL(2) rL(2)]);

    % ---- Plot groups ----
    % drug == 0
    hNoDrug = gobjects(1,1);
    if ~isempty(stats_noDrug.mean)
        hNoDrug = plot_with_shade(axPolar, axFill, theta, ...
                                  stats_noDrug.mean, stats_noDrug.sem, ...
                                  col_noDrug_line, col_noDrug_fill, 0.30, 2);
    else
        warning('No entries with drug == 0 found.');
    end

    % drug == 1 & mins_post == 15
    h15 = gobjects(1,1);
    if ~isempty(stats_15min.mean)
        h15 = plot_with_shade(axPolar, axFill, theta, ...
                              stats_15min.mean, stats_15min.sem, ...
                              col_15_line, col_15_fill, 0.30, 2);
    else
        warning('No entries with drug == 1 & mins_post == 15 found.');
    end

    % drug == 1 & mins_post == 35
    h35 = gobjects(1,1);
    if ~isempty(stats_35min.mean)
        h35 = plot_with_shade(axPolar, axFill, theta, ...
                              stats_35min.mean, stats_35min.sem, ...
                              col_35_line, col_35_fill, 0.30, 2);
    else
        warning('No entries with drug == 1 & mins_post == 35 found.');
    end

    % Keep polar axes on top so grid/ticks are visible
    % uistack(axPolar, 'top');

    % Legend (lines only)
    labs = {};
    hs   = [];
    if ~isempty(stats_noDrug.mean)
        hs(end+1)   = hNoDrug;
        labs{end+1} = sprintf('drug = 0 (n=%d)', stats_noDrug.n);
    end
    if ~isempty(stats_15min.mean)
        hs(end+1)   = h15;
        labs{end+1} = sprintf('drug = 1, %d min (n=%d)', 15, stats_15min.n);
    end
    if ~isempty(stats_35min.mean)
        hs(end+1)   = h35;
        labs{end+1} = sprintf('drug = 1, %d min (n=%d)', 35, stats_35min.n);
    end

    if ~isempty(hs)
        legend(axPolar, hs, labs, 'Location', 'best');
    end
    % title(axPolar, 'Mean Polar Response by Drug/Time');
end

% ---------- helpers ----------
function validate_table(T)
    if ~istable(T)
        error('Input must be a table.');
    end
    needVars = {'drug','mins_post','d_slow'};
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

function stats = get_group_stats_by_mask(T, mask, speed)
    % mask: logical index for rows belonging to this group

    idx = find(mask);
    vals = [];
    for k = idx'
        if speed == "slow"
            d = T.d_slow{k};
        elseif speed == "fast"
            d = T.d_fast{k};
        elseif speed == "vfast"
            d = T.d_vfast{k};
        else
            error('Unknown speed: %s (use "slow", "fast", or "vfast")', speed);
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

function rmax = estimate_rmax(sc, s15, s35)
    allVals = [];
    if ~isempty(sc) && isfield(sc,'mean') && ~isempty(sc.mean)
        allVals = [allVals; sc.mean(:) + sc.sem(:)]; 
    end
    if ~isempty(s15) && isfield(s15,'mean') && ~isempty(s15.mean)
        allVals = [allVals; s15.mean(:) + s15.sem(:)]; 
    end
    if ~isempty(s35) && isfield(s35,'mean') && ~isempty(s35.mean)
        allVals = [allVals; s35.mean(:) + s35.sem(:)]; 
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
    hLine   = polarplot(axPolar, th_line, r_line, '-', ...
                        'Color', lineColor, 'LineWidth', lw);
end
