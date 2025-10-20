% Scan date-structured folders, collect metadata, and build a summary table.
% Pattern: "DD_MM_YYYY" (e.g., "07_09_2025")

% ---- Find date folders in the current directory ----
listing = dir('.');
isDir = [listing.isdir];
names = {listing(isDir).name};
names = names(~ismember(names, {'.', '..'}));

datePattern = '^\d{2}_\d{2}_\d{4}$';
dateFolders = names(~cellfun('isempty', regexp(names, datePattern, 'once')));

records = struct( ...
    'date_folder', string.empty(0,1), ...
    'subdir', string.empty(0,1), ...
    'experiment_name', string.empty(0,1), ...
    'fly_name', string.empty(0,1), ...
    'fly_genotype', string.empty(0,1), ...
    'fly_age', string.empty(0,1), ...
    'fly_sex', string.empty(0,1), ...
    'comments', string.empty(0,1), ...
    'timestamp_raw', string.empty(0,1), ...
    'timestamp_date', string.empty(0,1), ...
    'timestamp_time', string.empty(0,1), ...
    'g4_files_found', zeros(0,1), ...
    'mean_vm', NaN(0,1), ...
    'std_vm', NaN(0,1), ...
    'min_vm', NaN(0,1), ...
    'max_vm', NaN(0,1) );

for i = 1:numel(dateFolders)
    df = dateFolders{i};
    subs = dir(fullfile(df, '*'));
    subs = subs([subs.isdir]);                            % only directories
    subnames = setdiff({subs.name}, {'.', '..'});         % skip dot dirs

    for j = 1:numel(subnames)
        sd = subnames{j};
        metaPath = fullfile(df, sd, 'metadata.mat');

        exp_name = ""; fly_name = ""; fly_geno = ""; fly_age = ""; fly_sex = "";
        comments = ""; ts_raw = ""; ts_date = ""; ts_time = "";

        % ---- Read metadata (if present) ----
        if exist(metaPath, 'file') == 2
            try
                S = load(metaPath);                       % load .mat
            catch ME
                warning('Failed to load %s: %s', metaPath, ME.message);
                S = struct(); % continue with empty metadata
            end

            if isfield(S, 'metadata')
                md = S.metadata;
            else
                md = S;
            end

            exp_name = getFieldAsString(md, 'experiment_name');
            fly_name = getFieldAsString(md, 'fly_name');
            fly_geno = getFieldAsString(md, 'fly_genotype');
            fly_age  = getFieldAsString(md, 'fly_age');
            fly_sex  = getFieldAsString(md, 'fly_sex');
            comments = getFieldAsString(md, 'comments');

            ts_raw = getFieldAsString(md, 'timestamp');
            [ts_date, ts_time] = splitTimestamp(ts_raw);
        end

        % ---- Find and process G4 files in this experiment folder ----
        g4List = dir(fullfile(df, sd, 'G4_TDMS*'));
        v_all = [];   % will hold concatenated voltage across all matches
        for k = 1:numel(g4List)
            g4Path = fullfile(g4List(k).folder, g4List(k).name);
            try
                Sg4 = load(g4Path);  % load into struct
                v = extractVoltageVector(Sg4); % returns [] if not found/invalid
                if ~isempty(v)
                    v_all = [v_all, v]; %#ok<AGROW>
                end
            catch ME
                warning('Failed to load/parse %s: %s', g4Path, ME.message);
            end
        end

        % Compute stats over the entire recording(s) in this subdir
        if ~isempty(v_all)
            mean_vm = mean(v_all, 'omitnan');
            std_vm  = std(v_all, 0, 'omitnan');

            % NaN-safe min/max (works on older MATLABs too)
            vn = v_all(~isnan(v_all));
            if isempty(vn)
                min_vm = NaN; max_vm = NaN;
            else
                min_vm = min(vn);
                max_vm = max(vn);
            end
        else
            mean_vm = NaN;
            std_vm  = NaN;
            min_vm  = NaN;
            max_vm  = NaN;
        end

        % Append a record
        records(end+1) = struct( ... %#ok<SAGROW>
            'date_folder', string(df), ...
            'subdir', string(sd), ...
            'experiment_name', exp_name, ...
            'fly_name', fly_name, ...
            'fly_genotype', fly_geno, ...
            'fly_age', fly_age, ...
            'fly_sex', fly_sex, ...
            'comments', comments, ...
            'timestamp_raw', ts_raw, ...
            'timestamp_date', ts_date, ...
            'timestamp_time', ts_time, ...
            'g4_files_found', numel(g4List), ...
            'mean_vm', mean_vm, ...
            'std_vm', std_vm, ...
            'min_vm', min_vm, ...
            'max_vm', max_vm );
    end
end

% ---- Convert to table ----
if isempty(records)
    warning('No records found beneath date-structured folders.');
    summaryTable = table;
else
    summaryTable = struct2table(records, 'AsArray', true);
end

% Display the table in the Command Window
disp(summaryTable);

% Optionally save to a CSV for convenience:
% if ~isempty(summaryTable)
%     writetable(summaryTable, 'metadata_summary.csv');
%     fprintf('Saved summary to metadata_summary.csv\n');
% end

% ---- Helper functions (local) ----
function out = getFieldAsString(s, fname)
    if isstruct(s) && isfield(s, fname)
        val = s.(fname);
        if isstring(val)
            out = val;
        elseif ischar(val)
            out = string(val);
        elseif isnumeric(val) && isscalar(val)
            out = string(val);
        else
            try
                out = string(evalc('disp(val)')); % captures readable form
                out = strtrim(erase(out, newline));
            catch
                out = "";
            end
        end
    else
        out = "";
    end
end

function [d, t] = splitTimestamp(ts)
    d = ""; t = "";
    if strlength(ts) == 0
        return;
    end
    tok = regexp(char(ts), '^(\d{2}-\d{2}-\d{4})(\d{2}_\d{2}_\d{2})$', 'tokens', 'once');
    if ~isempty(tok)
        d = string(tok{1});
        t = string(tok{2});
        return;
    end
    if strlength(ts) >= 10
        d = extractBetween(ts, 1, 10);
        t = extractAfter(ts, 10);
        if strlength(t) == 0
            t = "";
        end
    end
end

function v = extractVoltageVector(Sg4)
% Try to extract voltage vector from a loaded G4 file struct.
% Expected: Sg4.Log.ADC.Volts with size [channels x samples]; use channel 2.
% Multiply by 10 (as in your example).
    v = [];
    try
        if isfield(Sg4,'Log') && isfield(Sg4.Log,'ADC') && isfield(Sg4.Log.ADC,'Volts')
            volts = Sg4.Log.ADC.Volts;
            if isnumeric(volts)
                if ndims(volts) == 2 && size(volts,1) >= 2
                    v = volts(2,:)*10;
                elseif isvector(volts)
                    v = volts(:).'; % row vector
                    v = v*10;
                elseif ndims(volts) == 2 && size(volts,1) == 1
                    v = volts(1,:)*10;
                end
            end
        end
    catch
        v = [];
    end
end
