function polar_cells_by_strain_subplots(T)
% POLAR_CELLS_BY_STRAIN_SUBPLOTS
% Make a [1,2] subplot figure:
%   (1) CONTROL: all cells (light gray) + mean (black)
%   (2) TURTLE : all cells (light red) + mean (red)
%
% Input table T:
%   - T.Strain : 'control'/'turtle' (string/char/categorical)
%   - T.d_slow : cell, each 16x2 numeric [theta(rad), response]

    validate_table_strain(T);

    % Canonical theta from the first valid entry
    theta = get_theta_from_table(T);
    if isempty(theta)
        error('Could not find a valid 16x2 d_slow to read angles from.');
    end

    % Gather per-strain values and means
    [vals_control, n_control] = get_group_values_by_strain(T, 'control');
    [vals_turtle,  n_turtle ] = get_group_values_by_strain(T, 'ttl');

    mean_control = mean(vals_control, 2, 'omitnan');
    mean_turtle  = mean(vals_turtle,  2, 'omitnan');

    % Colors
    col_ctrl_all = [0.80 0.80 0.80]; % light gray
    col_ctrl_mean= [0 0 0];          % black
    col_turt_all = [1.00 0.70 0.70]; % light red
    col_turt_mean= [1 0 0];          % red

    % Reasonable shared r-limit across both panels
    rmax = estimate_rmax_from_vals(vals_control, vals_turtle);
    if ~isfinite(rmax) || isnan(rmax) || rmax<=0
        rmax = [];
    end

    % --- Figure & [1,2] polar subplots (robust across MATLAB versions) ---
    fh = figure;
    axStub1 = subplot(1,2,1);
    axStub2 = subplot(1,2,2);
    ax1 = polaraxes('Position', axStub1.Position); delete(axStub1);
    ax2 = polaraxes('Position', axStub2.Position); delete(axStub2);

    % Panel 1: CONTROL
    hold(ax1, 'on');
    if ~isempty(vals_control)
        % individual cells
        for k = 1:size(vals_control,2)
            polarplot(ax1, [theta; theta(1)], [vals_control(:,k); vals_control(1,k)], ...
                '-', 'Color', col_ctrl_all, 'LineWidth', 0.75);
        end
        % mean
        polarplot(ax1, [theta; theta(1)], [mean_control; mean_control(1)], ...
            '-', 'Color', col_ctrl_mean, 'LineWidth', 2.0);
    else
        text(0.5,0.5,'No control data', 'Parent', axes('Position', ax1.Position, 'Color','none'), ...
            'HorizontalAlignment','center');
    end
    if ~isempty(rmax); rlim(ax1,[0 rmax]); end
    title(ax1, sprintf('control (n=%d)', n_control));

    % Panel 2: TURTLE
    hold(ax2, 'on');
    if ~isempty(vals_turtle)
        % individual cells
        for k = 1:size(vals_turtle,2)
            polarplot(ax2, [theta; theta(1)], [vals_turtle(:,k); vals_turtle(1,k)], ...
                '-', 'Color', col_turt_all, 'LineWidth', 0.75);
        end
        % mean
        polarplot(ax2, [theta; theta(1)], [mean_turtle; mean_turtle(1)], ...
            '-', 'Color', col_turt_mean, 'LineWidth', 2.0);
    else
        text(0.5,0.5,'No turtle data', 'Parent', axes('Position', ax2.Position, 'Color','none'), ...
            'HorizontalAlignment','center');
    end
    if ~isempty(rmax); rlim(ax2,[0 rmax]); end
    title(ax2, sprintf('turtle (n=%d)', n_turtle));

    % Optional: uniform ThetaZero/Dir so panels match
    ax1.ThetaZeroLocation = 'left'; ax1.ThetaDir = 'clockwise'; ax1.ThetaTickLabel = [];
    ax2.ThetaZeroLocation = 'left'; ax2.ThetaDir = 'clockwise'; ax2.ThetaTickLabel = [];
    ax1.RLim = [0 40]; ax2.RLim = [0 40];

    % Overall title
    sgtitle(fh, 'Polar Responses by Strain');
    fh.Position = [258  537  1158  472];
end

% ---------- helpers ----------
function validate_table_strain(T)
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

function [vals, n] = get_group_values_by_strain(T, target)
    % Return a 16 x N matrix of responses for the given strain
    S = T.Strain;
    if iscellstr(S) || isstring(S)
        mask = strcmpi(string(S), target);
    elseif iscategorical(S)
        mask = strcmpi(string(S), target);
    else
        error('Unsupported type for T.Strain.');
    end

    idx  = find(mask);
    vals = [];
    for k = idx'
        d = T.d_slow{k};
        if isnumeric(d) && isequal(size(d), [16 2])
            vals = [vals, d(:,2)]; %#ok<AGROW>
        end
    end
    n = size(vals,2);
end

function rmax = estimate_rmax_from_vals(vc, vt)
    allv = [];
    if ~isempty(vc), allv = [allv; vc(:)]; end
    if ~isempty(vt), allv = [allv; vt(:)]; end
    if isempty(allv)
        rmax = NaN;
        return;
    end
    rmax = max(allv, [], 'omitnan');
    if ~isfinite(rmax) || rmax<=0
        rmax = NaN;
    end
end
