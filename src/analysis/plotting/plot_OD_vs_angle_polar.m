function plot_OD_vs_angle_polar()
%
% For each .mat file in folderPath, this:
%   - Loads struct bar_results
%   - Finds PD = max(bar_results.max_v{1,1}) over 16 elements
%   - Finds OD = mean of values 4 steps away on either side (orthogonal directions), circularly
%   - Computes r = (PD - OD) / (PD + OD)
%   - Uses theta = bar_results.resultant_angle (wrapped to [0, 2*pi))
%   - Plots a single point on polar axes for each file
%
% Example:
%   plot_bar_results_polar('C:\data\my_mats')

% If no input, use current directory
% if nargin < 1 || ~isfolder(folderPath)
%     error('Provide a valid folder path containing .mat files.');
% end
folderPath = cd;

files = dir(fullfile(folderPath, '*.mat'));
if isempty(files)
    warning('No .mat files found in %s', folderPath);
    return;
end

% Figure & axes
figure('Name','bar\_results Polar Plot','Color','w');
pax = polaraxes;
hold(pax, 'on');
title(pax, 'Polar plot of (PD-OD)/(PD+OD) per file');
% rlim(pax, [0 1]);  % initial radial limit; adjusted later below if desired
legendEntries = {};

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);

    % Load only the needed variable to be safe/fast
    S = load(fpath, 'bar_results');
    if ~isfield(S, 'bar_results')
        warning('Skipping %s: missing bar_results.', files(k).name);
        continue;
    end
    br = S.bar_results;

    % Basic validation
    if ~isfield(br, 'resultant_angle') || ~isfield(br, 'max_v')
        warning('Skipping %s: bar_results missing required fields.', files(k).name);
        continue;
    end

    % Extract the 16x1 vector
    v = br.max_v;
    if iscell(v)
        try
            v = v{1,1};
        catch
            warning('Skipping %s: max_v cell not at {1,1}.', files(k).name);
            continue;
        end
    end
    v = v(:);  % ensure column
    if numel(v) < 16
        warning('Skipping %s: expected >=16 elements in max_v{1,1}.', files(k).name);
        continue;
    end
    v = v(1:16);  % take the first 16 if longer

    % Find PD and its index
    [PD, idxPD] = max(v);

    % Orthogonal indices are 90 degrees away: ±4 steps (circular)
    idxOD1 = mod(idxPD + 4 - 1, 16) + 1;   % +4 with wrap
    idxOD2 = mod(idxPD - 4 - 1, 16) + 1;   % -4 with wrap (equiv. +12)
    OD = mean([v(idxOD1), v(idxOD2)]);

    % Compute ratio, guard against divide-by-zero
    denom = PD + OD;
    if denom == 0
        r = NaN;
    else
        r = (PD - OD) / denom;
    end

    % Angle handling: wrap to [0, 2*pi)
    theta = br.resultant_angle;
    theta = mod(theta, 2*pi);

    if ~isnan(r)
        % Plot a single point
        polarplot(pax, theta, r, 'o', 'MarkerSize', 12, ...
                  'DisplayName', files(k).name, 'LineWidth', 1.5, 'MarkerEdgeColor', 'r');
        legendEntries{end+1} = files(k).name; 
    else
        warning('Skipping plot for %s: r is NaN (PD+OD==0).', files(k).name);
    end
end

if ~isempty(legendEntries)
    legend(pax, 'Interpreter', 'none', 'Location', 'bestoutside');
else
    warning('No valid data points were plotted.');
end

% Your preferred display range and figure size
rlim([0 70])
f = gcf;
f.Position = [192 379 854 644];

end
