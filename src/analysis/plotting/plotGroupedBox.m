function plotGroupedBox(T, columnName)
    % plotGroupedBox: Create boxplots grouped by Strain and Type.
    % T: input table
    % columnName: string specifying the column to plot (numeric data)
    %
    % Expected columns: 'Strain' and 'Type'

    % ---------- Coerce Strain & Type to categorical safely ----------
    T.Strain = coerceToCategorical(T, 'Strain');
    T.Type   = coerceToCategorical(T, 'Type');

    % ---------- Coerce the Y column to numeric vector ----------
    data = coerceToNumericVector(T.(columnName), columnName);

    % Remove rows with missing Strain/Type or completely missing/NaN data
    good = ~ismissing(T.Strain) & ~ismissing(T.Type) & isfiniteOrNaN(data);
    T = T(good,:);
    data = data(good);

    % ---------- Build groups ----------
    groups = findgroups(T.Type, T.Strain);
    [uniqueGroups, ~, ~] = unique(groups);
    groupLabels = cell(numel(uniqueGroups),1);

    data = T.(columnName);

    % Check if data is a cell array and if so, convert to double.
    if iscell(data)
        data = cell2mat(data);
    end 

    % Create boxchart
    figure;
    hold on;

    % Define color scheme lookup table
    colorMap = containers.Map;
    colorMap('42F06_T4T5_control_on')  = [0.8 0.8 0.8]; % light grey
    colorMap('42F06_T4T5_control_off') = [0.4 0.4 0.4]; % dark grey
    colorMap('42F06_T4T5_ttl_on')      = [0.6 0.9 0.6]; % light green
    colorMap('42F06_T4T5_ttl_off')     = [0.2 0.5 0.2]; % dark green
    

    % Loop through each group and plot separately (so we can assign colors)
    for i = 1:numel(uniqueGroups)
        mask = groups == uniqueGroups(i);
        thisStrain = string(T.Strain(find(mask, 1)));
        thisType   = string(T.Type(find(mask, 1)));
        key = sprintf('%s_%s', thisStrain, thisType);
        color = colorMap(key);

        % Plot box for this group
        b = boxchart(repmat(i,sum(mask),1), data(mask), 'BoxFaceColor', color, ...
                     'MarkerStyle', 'none'); %#ok<NASGU>

        % Plot individual data points
        scatter(repmat(i,sum(mask),1), data(mask), 40, 'k', 'o', ...
            'MarkerFaceColor', 'none', 'LineWidth', 1.2);

        % Label
        groupLabels{i} = sprintf('%s-%s', strrep(thisStrain, "_", "-"), strrep(thisType, "_", "-"));
    end

    % Aesthetics
    set(gca, 'XTick', 1:numel(uniqueGroups), 'XTickLabel', groupLabels);
    xlabel('Group (Strain-Type)');
    colStr = strrep(columnName, '_', '-');
    ylabel(colStr);
    title(sprintf('%s', colStr));
    box off;
    hold off;
    ax = gca;
    ax.TickDir = "out";
    ax.LineWidth = 1.2;
    ax.TickLength = [0.015 0.015];
    ax.FontSize = 14;
    f = gcf;
    f.Position = [620   392   450   575];

    % ===================== Helpers =====================

    function C = coerceToCategorical(T, varName)
        % Coerce a table variable to categorical, tolerating cellstr/char/string/categorical.
        x = T.(varName);
        if iscell(x)
            x = string(x);                 % cellstr or mixed -> string
        elseif ischar(x)
            x = string(x);                 % char row -> string
        elseif isstring(x)
            % ok
        elseif iscategorical(x)
            C = x; return;
        else
            % Any other type (e.g., logical/numeric) -> string to be safe
            x = string(x);
        end
        C = categorical(x);
    end
    
    function v = coerceToNumericVector(x, columnName)
        % Convert input column into a numeric column vector (double).
        % Accepts numeric arrays, cell arrays containing numerics/strings/empties,
        % string arrays, char arrays, and cellstr. Empty -> NaN.
    
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
                    tmp = str2double(string(ck));
                    if isnan(tmp) && ~(ischar(ck) || isstring(ck))
                        error('Cannot convert cell(%d) in column "%s" to numeric.', k, columnName);
                    end
                    v(k) = tmp;
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

    function tf = isfiniteOrNaN(v)
        % True for finite numbers or NaN (keeps NaNs so boxchart can ignore them gracefully)
        tf = isnumeric(v) & (isfinite(v) | isnan(v));
    end

end
