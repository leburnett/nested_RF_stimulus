function results = runGroupedStats(T, columnName)
% runGroupedStats: Tests assumptions and compares 4 groups (Strain x Type)
% Chooses one-way ANOVA or Kruskal–Wallis based on normality + equal variances.

    arguments
        T table
        columnName {mustBeTextScalar}
    end
    columnName = string(columnName);

    % ---------- Coerce group columns ----------
    T.Strain = coerceToCategorical(T, 'Strain');
    T.Type   = coerceToCategorical(T, 'Type');

    % ---------- Coerce Y column to numeric ----------
    x = coerceToNumericVector(T.(columnName), columnName);

    % Keep only rows with defined groups and numeric (NaN allowed)
    good = ~ismissing(T.Strain) & ~ismissing(T.Type) & isnumeric(x);
    T = T(good,:);
    x = x(good);

    % ---------- Build groups ----------
    G = findgroups(T.Type, T.Strain);

    % Robust group names (avoid relying on 2nd output of findgroups)
    nGroups = max(G);
    groupNames = strings(nGroups,1);
    for i = 1:nGroups
        idx = find(G==i, 1, 'first');
        groupNames(i) = string(T.Strain(idx)) + "-" + string(T.Type(idx));
    end

    % Split by group
    dataCell = cell(nGroups,1);
    Ns = zeros(nGroups,1);
    for i = 1:nGroups
        xi = x(G == i);
        xi = xi(isfinite(xi));   % remove NaNs for stats
        dataCell{i} = xi;
        Ns(i) = numel(xi);
    end

    results.groupsTable = table(groupNames, Ns, 'VariableNames', {'Group','n'});
    results.groupNames  = groupNames;

    % ---------- Normality per group (Lilliefors if n>=5) ----------
    hVec = nan(nGroups,1);
    pVec = nan(nGroups,1);
    call = strings(nGroups,1);
    hasLillie = exist('lillietest','file') == 2;

    for i = 1:nGroups
        xi = dataCell{i};
        if numel(xi) >= 5 && hasLillie
            try
                [h,p] = lillietest(xi);
            catch
                h = NaN; p = NaN;
            end
            hVec(i) = h;  pVec(i) = p;
            if ~isnan(h) && h==0
                call(i) = "normal (fail to reject)";
            elseif ~isnan(h) && h==1
                call(i) = "non-normal (reject)";
            else
                call(i) = "unknown";
            end
        else
            hVec(i) = NaN; pVec(i) = NaN;
            if ~hasLillie
                call(i) = "unknown (no lillietest)";
            elseif numel(xi) < 5
                call(i) = "unknown (n<5)";
            else
                call(i) = "unknown";
            end
        end
    end

    results.normality = table(groupNames, Ns, hVec, pVec, call, ...
        'VariableNames', {'Group','n','h','p','call'});

    % ---------- Homogeneity of variances (Levene) ----------
    if sum(Ns > 0) >= 2 && all(~cellfun(@isempty, dataCell))
        try
            finiteMask = isfinite(x);
            pVar = vartestn(x(finiteMask), G(finiteMask), ...
                'TestType','LeveneAbsolute', 'Display','off');
        catch
            pVar = NaN;
        end
    else
        pVar = NaN;
    end
    results.leveneP = pVar;

    % ---------- Choose omnibus test ----------
    allNok   = all(Ns >= 5);
    allNorm  = all(hVec == 0);   % unknown counts against "all normal"
    equalVar = ~isnan(pVar) && pVar >= 0.05;
    useANOVA = allNok && allNorm && equalVar;

    if sum(Ns > 0) < 2
        warning('Fewer than two groups with data. No omnibus test performed.');
        results.chosenTest   = 'None';
        results.omnibusP     = NaN;
        results.posthocTable = table();
        return;
    end

    % ---------- Run chosen test + post-hoc ----------
    if useANOVA
        [p, anovaTbl, stats] = anova1(x, G, 'off');
        c = multcompare(stats, 'Display','off');  % Tukey-Kramer
        results.chosenTest = 'ANOVA';
        results.omnibusP   = p;
        results.anovaTbl   = anovaTbl;
    else
        [p, kwTbl, stats] = kruskalwallis(x, G, 'off');
        c = multcompare(stats, 'Display','off');  % mean-rank comps (Dunn–Šidák control)
        results.chosenTest = 'KruskalWallis';
        results.omnibusP   = p;
        results.kwTbl      = kwTbl;
    end

    % Pairwise table
    pairI = results.groupNames(c(:,1));
    pairJ = results.groupNames(c(:,2));
    posthocTbl = table(pairI, pairJ, c(:,3), c(:,4), c(:,5), c(:,6), ...
        'VariableNames', {'GroupA','GroupB','LowerCI','Estimate','UpperCI','pValue'});
    results.posthocTable = posthocTbl;

    % ---------- Console summary ----------
    fprintf('\n=== Grouped comparison for %s ===\n', columnName);
    disp(results.groupsTable);
    fprintf('Normality by group (Lilliefors):\n');
    disp(results.normality(:,{'Group','n','p','call'}));
    fprintf('Levene''s p (equal variances): %s\n', num2str(results.leveneP,3));
    fprintf('Chosen omnibus test: %s | p = %s\n', results.chosenTest, num2str(results.omnibusP,3));
    fprintf('Pairwise post-hoc comparisons (adjusted):\n');
    disp(results.posthocTable);
end

% ===================== Helpers (same coercion as plotting fn) =====================

function C = coerceToCategorical(T, varName)
    x = T.(varName);
    if iscell(x), x = string(x);
    elseif ischar(x), x = string(x);
    elseif iscategorical(x)
        C = x; return;
    end
    C = categorical(x);
end

function v = coerceToNumericVector(x, columnName)
    if isnumeric(x), v = double(x(:)); return; end
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
