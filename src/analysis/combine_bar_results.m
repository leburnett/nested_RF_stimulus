function T = combine_bar_results()
% Create a summary table from all .mat files in the current directory.
% Compatible with two data layouts:
%  (A) Old format:
%      - bar_results: struct with fields
%          .Date, .Time, .Strain, .Type, .median_voltage, .resultant_angle
%          .slow: struct with
%             .magnitude, .angle_rad, .fwhm, .cv, .thetahat, .kappa, ...
%             .sym_ratio, .vector_sum, .DSI_vector, .DSI_pdnd
%      - d_slow: [16 x 2] numeric array at top level
%  (B) New format:
%      - bar_results: TABLE with columns
%          Date, Time, Strain, Type, median_voltage
%          and slow metrics in columns with "_slow" suffix, e.g.:
%              magnitude_slow, angle_rad_slow, fwhm_slow, cv_slow,
%              thetahat_slow, kappa_slow, sym_ratio_slow,
%              vector_sum_slow, DSI_vector_slow, DSI_pdnd_slow
%      - d: [16 x 2] numeric array at top level
%
% In all cases, we also compute
%   v_max    = d(5,2)
%   v_null   = d(13,2)
%   v_ortho1 = d(1,2)
%   v_ortho2 = d(9,2)

files = dir('*.mat');
n = numel(files);

% Preallocate table columns
Date            = strings(n,1);
Time            = strings(n,1);
Strain          = strings(n,1);
Type            = strings(n,1);
median_voltage  = nan(n,1);
resultant_angle = nan(n,1);

magnitude   = nan(n,1);
angle_rad   = nan(n,1);
fwhm        = nan(n,1);
cv          = nan(n,1);
thetahat    = nan(n,1);
kappa       = nan(n,1);
sym_ratio   = nan(n,1);

vector_sum  = cell(n,1);   % may be vector/array
DSI_vector  = cell(n,1);   % may be vector/array or scalar (store as cell safely)
DSI_pdnd    = nan(n,1);

d_slow      = cell(n,1);   % store [16 x 2] per row (new format uses "d", but we keep column name 'd_slow' for continuity)
d_fast      = cell(n,1);   % store [16 x 2] per row (new format uses "d", but we keep column name 'd_slow' for continuity)
d_vfast     = cell(n,1);   % store [16 x 2] per row (new format uses "d", but we keep column name 'd_slow' for continuity)

% Derived from d (or d_slow)
v_max_slow    = nan(n,1);
v_null_slow   = nan(n,1);
v_ortho1_slow = nan(n,1);
v_ortho2_slow = nan(n,1);
v_max_fast    = nan(n,1);
v_null_fast   = nan(n,1);
v_ortho1_fast = nan(n,1);
v_ortho2_fast = nan(n,1);
v_max_vfast    = nan(n,1);
v_null_vfast   = nan(n,1);
v_ortho1_vfast = nan(n,1);
v_ortho2_vfast = nan(n,1);


for k = 1:n
    fname = files(k).name;
    S = load(fname);

    if ~isfield(S, 'bar_results')
        warning('File "%s" has no "bar_results".', fname);
    else
        br = S.bar_results;

        % --- Common fields (Date/Time/Strain/Type/median_voltage/resultant_angle) ---
        Date(k)            = asString(getFromBar(br, 'Date',            false));
        Time(k)            = asString(getFromBar(br, 'Time',            false));
        Strain(k)          = asString(getFromBar(br, 'Strain',          false));
        Type(k)            = asString(getFromBar(br, 'Type',            false));
        median_voltage(k)  = asScalar(getFromBar(br, 'median_voltage',  false));
        resultant_angle(k) = asScalar(getFromBar(br, 'resultant_angle', false));

        % --- "slow" metrics ---
        magnitude(k)   = asScalar(getFromBar(br, 'magnitude',  true));
        angle_rad(k)   = asScalar(getFromBar(br, 'angle_rad',  true));
        fwhm(k)        = asScalar(getFromBar(br, 'fwhm',       true));
        cv(k)          = asScalar(getFromBar(br, 'cv',         true));
        thetahat(k)    = asScalar(getFromBar(br, 'thetahat',   true));
        kappa(k)       = asScalar(getFromBar(br, 'kappa',      true));
        sym_ratio(k)   = asScalar(getFromBar(br, 'sym_ratio',  true));

        % Array-capable fields (store in cells)
        vector_sum{k}  = asArray(getFromBar(br, 'vector_sum',  true));
        DSI_vector{k}  = asArray(getFromBar(br, 'DSI_vector',  true));

        % This one is typically scalar; keep numeric column
        DSI_pdnd(k)    = asScalar(getFromBar(br, 'DSI_pdnd',   true));
    end

    % --- d_slow ---
    D_slow = [];
    if isfield(S, 'd') && isnumeric(S.d)
        D_slow = S.d;
    elseif isfield(S, 'd_slow') && isnumeric(S.d_slow)
        D_slow = S.d_slow;
    end

    if ~isempty(D_slow)
        d_slow{k}  = D_slow;
        v_max_slow(k)    = getDSValue(D_slow, 5,  2);
        v_null_slow(k)   = getDSValue(D_slow, 13, 2);
        v_ortho1_slow(k) = getDSValue(D_slow, 1,  2);
        v_ortho2_slow(k) = getDSValue(D_slow, 9,  2);
    else
        d_slow{k} = [];
    end

    % --- d_fast ---
    D_fast = [];
    if isfield(S, 'd_fast') && isnumeric(S.d_fast)
        D_fast = S.d_fast;
    end

    if ~isempty(D_fast)
        d_fast{k}  = D_fast;
        v_max_fast(k)    = getDSValue(D_fast, 5,  2);
        v_null_fast(k)   = getDSValue(D_fast, 13, 2);
        v_ortho1_fast(k) = getDSValue(D_fast, 1,  2);
        v_ortho2_fast(k) = getDSValue(D_fast, 9,  2);
    else
        d_fast{k} = [];
    end

    % --- d_vfast ---
    D_vfast = [];
    if isfield(S, 'd_vfast') && isnumeric(S.d_vfast)
        D_vfast = S.d_vfast;
    end

    if ~isempty(D_vfast)
        d_vfast{k}  = D_vfast;
        v_max_vfast(k)    = getDSValue(D_vfast, 5,  2);
        v_null_vfast(k)   = getDSValue(D_vfast, 13, 2);
        v_ortho1_vfast(k) = getDSValue(D_vfast, 1,  2);
        v_ortho2_vfast(k) = getDSValue(D_vfast, 9,  2);
    else
        d_vfast{k} = [];
    end

end

T = table( ...
    Date, Time, Strain, Type, ...
    median_voltage, resultant_angle, ...
    magnitude, angle_rad, fwhm, cv, thetahat, kappa, sym_ratio, ...
    vector_sum, DSI_vector, DSI_pdnd, ...
    v_max_slow, v_null_slow, v_ortho1_slow, v_ortho2_slow, ...
    d_slow, ...
    v_max_fast, v_null_fast, v_ortho1_fast, v_ortho2_fast, ...
    d_fast, ...
    v_max_vfast, v_null_vfast, v_ortho1_vfast, v_ortho2_vfast, ...
    d_vfast, ...
    'VariableNames', { ...
        'Date','Time','Strain','Type', ...
        'median_voltage','resultant_angle', ...
        'magnitude','angle_rad','fwhm','cv','thetahat','kappa','sym_ratio', ...
        'vector_sum','DSI_vector','DSI_pdnd', ...
        'v_max_slow','v_null_slow','v_ortho1_slow','v_ortho2_slow', ...
        'd_slow' ...
        'v_max_fast', 'v_null_fast', 'v_ortho1_fast', 'v_ortho2_fast', ...
        'd_fast', ...
        'v_max_vfast', 'v_null_vfast', 'v_ortho1_vfast', 'v_ortho2_vfast', ...
        'd_vfast', ...
    } ...
);

% Optional: preview
disp(T);

% Optional: save
% save('summary_table.mat','T');

%% ----------------- Helper functions -----------------
function val = getFromBar(br, baseField, isSlow)
% Pull value from bar_results regardless of layout:
%   - struct:       br.(baseField) or br.slow.(baseField)
%   - table layout: br.(baseField) or br.([baseField '_slow'])
% Handles single-row tables by unwrapping their single entry.

    val = [];  % default if not found

    if istable(br)
        if isSlow
            col = baseField + "_slow";
        else
            col = baseField;
        end
        if ismember(col, br.Properties.VariableNames)
            v = br.(col);
            % Unwrap if it's a single-row table variable
            if istable(v) || istimetable(v)
                % unlikely, but keep for completeness
                val = v{1,:};
            else
                if iscell(v)
                    if ~isempty(v)
                        val = v{1};
                    end
                elseif isstring(v) || ischar(v) || isnumeric(v) || islogical(v)
                    if isscalar(v)
                        val = v;
                    elseif size(v,1) == 1
                        val = v(1,:);  % row vector
                    else
                        % multi-row variable; try first entry
                        val = v(1,:);
                    end
                else
                    % unknown type; leave as []
                    val = [];
                end
            end
        end

    elseif isstruct(br)
        if ~isSlow
            if isfield(br, baseField)
                val = br.(baseField);
            end
        else
            if isfield(br, 'slow') && isstruct(br.slow) && isfield(br.slow, baseField)
                val = br.slow.(baseField);
            end
        end

    else
        % Unsupported br type
        val = [];
    end
end

function out = asScalar(x)
    if isnumeric(x) && isscalar(x)
        out = x;
    elseif islogical(x) && isscalar(x)
        out = double(x);
    elseif isstring(x) || ischar(x)
        num = str2double(string(x));
        if isnan(num)
            out = NaN;
        else
            out = num;
        end
    else
        % If a 1x1 cell, unwrap and retry
        if iscell(x) && numel(x)==1
            out = asScalar(x{1});
        else
            out = NaN;
        end
    end
end

function out = asArray(x)
    % Return numeric arrays as-is; unwrap 1x1 cells; otherwise [].
    if isnumeric(x)
        out = x;
    elseif iscell(x) && numel(x)==1 && isnumeric(x{1})
        out = x{1};
    else
        out = [];
    end
end

function s = asString(x)
    if isstring(x) || ischar(x)
        s = string(x);
    elseif isdatetime(x)
        s = string(x);
    elseif isnumeric(x) && isscalar(x)
        try
            s = string(datetime(x, 'ConvertFrom', 'datenum'));
        catch
            s = string(x);
        end
    elseif iscell(x) && ~isempty(x)
        s = asString(x{1});
    else
        s = "";
    end
end

function v = getDSValue(d, r, c)
    if isnumeric(d) && ndims(d)==2 && size(d,1) >= r && size(d,2) >= c
        v = d(r,c);
    else
        v = NaN;
    end
end

end
