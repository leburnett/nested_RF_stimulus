function plot_PD_vs_angle_polar()
% 
% For each .mat file in folderPath, this:
%   - Loads struct bar_results
%   - Finds PD = max(bar_results.max_v{1,1}) over 16 elements
%   - Finds ND = element 8 positions away, circularly (wrap-around)
%   - Computes r = (PD - ND) / (PD + ND)
%   - Uses theta = bar_results.resultant_angle (auto-detects degrees vs radians)
%   - Plots a single point on polar axes for each file
%
% Notes:
% - If theta > 2*pi, it will be treated as degrees and converted to radians.
% - Files missing required fields are skipped with a warning.
%
% Example:
%   plot_bar_results_polar('C:\data\my_mats')

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
title(pax, 'Polar plot of (PD-ND)/(PD+ND) per file');
rlim(pax, [0 1]);  % the ratio is naturally in [0,1]; NaNs will be skipped
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

    % ND is 8 positions away, circularly (wrap-around)
    idxND = mod(idxPD + 8 - 1, 16) + 1;  % 1-based indexing
    ND = v(idxND);

    % Compute ratio, guard against divide-by-zero
    denom = PD + ND;
    if denom == 0
        r = NaN;
    else
        r = (PD - ND) / denom;
    end

    % Angle handling: treat as degrees if clearly > 2*pi
    theta = br.resultant_angle;
    theta = mod(theta, 2*pi);

    if ~isnan(r)
        % Plot a single point and a radial guide line (optional)
        polarplot(pax, theta, r, 'o', 'MarkerSize', 12, ...
                  'DisplayName', files(k).name, 'LineWidth', 1.5, 'MarkerEdgeColor', 'k');
        % Uncomment the next line if you'd also like a radial line from origin to the point:
        % polarplot(pax, [theta theta], [0 r], '-', 'HandleVisibility', 'off');
        legendEntries{end+1} = files(k).name;
    else
        warning('Skipping plot for %s: r is NaN (PD+ND==0).', files(k).name);
    end
end

if ~isempty(legendEntries)
    legend(pax, 'Interpreter', 'none', 'Location', 'bestoutside');
else
    warning('No valid data points were plotted.');
end

rlim([0 0.4])
f = gcf;
f.Position = [192   379   854   644];

end
