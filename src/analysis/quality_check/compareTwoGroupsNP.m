function S = compareTwoGroupsNP(x, y, varargin)
%COMPARETWOGROUPSNP Nonparametric comparison and boxplot with colored boxes and means.
%   S = compareTwoGroupsNP(x, y) runs a Wilcoxon rank-sum test (Mann–Whitney)
%   comparing vectors x and y (independent samples), and produces a figure with
%   two colored boxplots (red & blue) and open-circle datapoints.
%
%   Adds group mean values as text on the plot.
%
%   Optional name-value pairs:
%     'paired'   (false)
%     'alpha'    (0.05)
%     'labels'   ({'Group 1','Group 2'})
%     'title'    ('')
%     'showStats'(true)

% -------- Parse & validate inputs
p = inputParser;
p.addRequired('x', @(v) isnumeric(v) && isvector(v));
p.addRequired('y', @(v) isnumeric(v) && isvector(v));
p.addParameter('paired', false, @(b) islogical(b) && isscalar(b));
p.addParameter('alpha', 0.05, @(a) isnumeric(a) && isscalar(a) && a>0 && a<1);
p.addParameter('labels', {'Group 1','Group 2'}, @(c) iscellstr(c) && numel(c)==2);
p.addParameter('title', '', @(s) isstring(s) || ischar(s));
p.addParameter('ylabel', '', @(s) isstring(s) || ischar(s));
p.addParameter('xlabel', '', @(s) isstring(s) || ischar(s));
p.addParameter('showStats', true, @(b) islogical(b) && isscalar(b));
p.addParameter('colors', [0.8 0 0.7; 0 0 0.7], @(c) isnumeric(c) && isequal(size(c), [2,3]));
p.parse(x, y, varargin{:});
opt = p.Results;

x = x(:); y = y(:);
x = x(~isnan(x)); y = y(~isnan(y));
n1 = numel(x); n2 = numel(y);

if opt.paired && n1 ~= n2
    error('For paired design, x and y must have same number of elements.');
end

% -------- Nonparametric test
if opt.paired
    [pval,~,st] = signrank(x, y, 'alpha', opt.alpha);
    testName = 'signrank (Wilcoxon signed-rank)';
    if isfield(st, 'zval')
        eff = st.zval / sqrt(n1);
    else
        eff = NaN;
    end
else
    [pval,~,st] = ranksum(x, y, 'alpha', opt.alpha);
    testName = 'ranksum (Mann–Whitney U)';
    U = st.ranksum - n1*(n1+1)/2;
    A = U / (n1*n2);
    eff = 2*A - 1;
end
h = pval < opt.alpha;

% -------- Plot setup
figure('Color','w');
data = [x; y];
grp  = [ones(n1,1); 2*ones(n2,1)];
colors = opt.colors;

useBoxchart = ~isempty(which('boxchart'));
if useBoxchart
    b1 = boxchart(ones(n1,1), x, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0.15, 'MarkerStyle','none');
    hold on
    b2 = boxchart(2*ones(n2,1), y, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0.15, 'MarkerStyle','none');
    set(gca, 'XTick', [1 2], 'XTickLabel', opt.labels);
else
    boxplot(data, grp, 'Labels', opt.labels, 'Colors', 'k', 'Symbol', '', 'MarkerStyle','none');
    hB = findobj(gca,'Tag','Box');
    for j = 1:2
        patch(get(hB(j),'XData'), get(hB(j),'YData'), colors(3-j,:), 'FaceAlpha', 0.3);
    end
end

% Overlay datapoints (open circles)
rng('default');
jitter = 0.12;
scatter(1 + jitter*randn(n1,1), x, 45, 'o', 'MarkerEdgeColor', colors(1,:), 'LineWidth', 1);
scatter(2 + jitter*randn(n2,1), y, 45, 'o', 'MarkerEdgeColor', colors(2,:), 'LineWidth', 1);

if ~isempty(opt.ylabel)
    ylb = opt.ylabel;
else
    ylb = "";
end
ylabel(ylb);

if ~isempty(opt.title)
    ttl = opt.title;
else
    ttl = sprintf('Nonparametric comparison: %s', testName);
end
title(ttl);
grid on; box on;

% Add mean text
m1 = mean(x); m2 = mean(y);
yl = ylim;
text(1, yl(2) - 0.05*range(yl), sprintf('Mean = %.3f', m1), 'Color', colors(1,:), ...
     'HorizontalAlignment','center', 'FontWeight','bold');
text(2, yl(2) - 0.05*range(yl), sprintf('Mean = %.3f', m2), 'Color', colors(2,:), ...
     'HorizontalAlignment','center', 'FontWeight','bold');

% Add p-value annotation
if opt.showStats
    ystar = yl(2) + 0.05*range(yl);
    plot([1 2], [yl(2) yl(2)], 'k-', 'LineWidth', 1);
    text(1.5, ystar, sprintf('p = %.4g', pval), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
    ylim([yl(1), yl(2) + 0.12*range(yl)]);
end

ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.2;
box off
grid off
f = gcf;
f.Position = [364   357   308   507];

% -------- Summary struct
S = struct();
S.test = testName;
S.p = pval;
S.h = h;
S.alpha = opt.alpha;
S.n = [n1 n2];
S.median = [median(x) median(y)];
S.mean = [m1 m2];
S.effect = eff;
S.details = st;
S.labels = opt.labels;

end
