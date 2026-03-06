function fig_path = find_matching_p1_gridplot(p1_base, p2_date, p2_time)
% FIND_MATCHING_P1_GRIDPLOT  Find the Protocol 1 GridPlots.fig for a P2 cell.
%
%   FIG_PATH = FIND_MATCHING_P1_GRIDPLOT(P1_BASE, P2_DATE, P2_TIME)
%   matches a Protocol 2 experiment (by date and time) to the correct
%   Protocol 1 GridPlots.fig file using closest-preceding-time heuristic.
%
%   INPUTS:
%     p1_base  - Base path to Protocol 1 data directory
%     p2_date  - P2 date string in 'YYYY_MM_DD' format
%     p2_time  - P2 time string in 'HH_MM' format
%
%   OUTPUT:
%     fig_path - Full path to GridPlots.fig, or empty string if not found
%
%   MATCHING LOGIC:
%     1. Convert P2 date (YYYY_MM_DD) to P1 format (MM_DD_YYYY)
%     2. List all fly folders for that date in P1
%     3. Parse time from fly folder names (StrainCross-HH_MM_CS)
%     4. Select the folder with the closest preceding time to P2 time

fig_path = '';

% Convert P2 date format (YYYY_MM_DD) to P1 format (MM_DD_YYYY)
parts = strsplit(p2_date, '_');
if numel(parts) ~= 3
    warning('Invalid date format: %s', p2_date);
    return;
end
p1_date_folder = sprintf('%s_%s_%s', parts{2}, parts{3}, parts{1});

p1_date_path = fullfile(p1_base, p1_date_folder);
if ~isfolder(p1_date_path)
    return;
end

% List subdirectories (fly folders)
d = dir(p1_date_path);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));

% Filter to fly folders (contain a dash separator with time info)
fly_folders = {};
fly_times_min = [];

for i = 1:numel(d)
    name = d(i).name;
    % Pattern: StrainCross-HH_MM_CS (e.g., 42F06XJFRC28_T4T5-11_40_91)
    dash_idx = strfind(name, '-');
    if isempty(dash_idx)
        continue;
    end
    time_part = name(dash_idx(end)+1:end); % e.g., '11_40_91'
    time_tokens = strsplit(time_part, '_');
    if numel(time_tokens) < 2
        continue;
    end
    hh = str2double(time_tokens{1});
    mm = str2double(time_tokens{2});
    if isnan(hh) || isnan(mm)
        continue;
    end

    % Check that GridPlots.fig exists in this folder
    gp = fullfile(p1_date_path, name, 'GridPlots.fig');
    if ~isfile(gp)
        continue;
    end

    fly_folders{end+1} = name; %#ok<AGROW>
    fly_times_min(end+1) = hh * 60 + mm; %#ok<AGROW>
end

if isempty(fly_folders)
    return;
end

% Parse P2 time
p2_parts = strsplit(p2_time, '_');
p2_hh = str2double(p2_parts{1});
p2_mm = str2double(p2_parts{2});
p2_time_min = p2_hh * 60 + p2_mm;

% Find closest preceding time (P1 time <= P2 time)
time_diffs = p2_time_min - fly_times_min;
valid = time_diffs >= 0;

if ~any(valid)
    % No preceding P1 recording; pick the closest overall
    [~, best_idx] = min(abs(time_diffs));
else
    % Pick the closest preceding one
    time_diffs(~valid) = Inf;
    [~, best_idx] = min(time_diffs);
end

fig_path = fullfile(p1_date_path, fly_folders{best_idx}, 'GridPlots.fig');

end
