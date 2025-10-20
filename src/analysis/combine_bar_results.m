% function T = combine_bar_results()
% % Create a summary table from all .mat files in the current directory.
% % Each .mat file is expected to contain:
% %   - bar_results (a struct)
% %       .Date, .Time, .Strain, .Type, .median_voltage, .resultant_angle
% %       .slow (a struct)
% %            .magnitude, .angle_rad, .fwhm, .cv, .thetahat, .kappa, ...
% %            .sym_ratio, .vector_sum, .DSI_vector, .DSI_pdnd
% %   - d_slow (a [16 x 2] array)  <-- stored at top level of the .mat file
% %
% % This version also adds four columns derived from d_slow(:,2):
% %   v_max   = d_slow(5,2)
% %   v_null  = d_slow(13,2)
% %   v_ortho1= d_slow(1,2)
% %   v_ortho2= d_slow(9,2)
% 
% files = dir('*.mat');
% n = numel(files);
% 
% % Preallocate columns
% Date            = strings(n,1);
% Time            = strings(n,1);
% Strain          = strings(n,1);
% Type            = strings(n,1);
% median_voltage  = nan(n,1);
% resultant_angle = nan(n,1);
% 
% magnitude   = nan(n,1);
% angle_rad   = nan(n,1);
% fwhm        = nan(n,1);
% cv          = nan(n,1);
% thetahat    = nan(n,1);
% kappa       = nan(n,1);
% sym_ratio   = nan(n,1);
% vector_sum  = cell(n,1);   % arrays -> store in cells
% DSI_vector  = nan(n,1);    % kept scalar to match your current code
% DSI_pdnd    = nan(n,1);
% 
% d_slow      = cell(n,1);   % each row holds a [16 x 2] numeric array
% 
% % New columns derived from d_slow(:,2)
% v_max    = nan(n,1);  % d_slow(5,2)
% v_null   = nan(n,1);  % d_slow(13,2)
% v_ortho1 = nan(n,1);  % d_slow(1,2)
% v_ortho2 = nan(n,1);  % d_slow(9,2)
% 
% for k = 1:n
%     fname = files(k).name;
%     S = load(fname);  % load everything
% 
%     if ~isfield(S, 'bar_results')
%         warning('File "%s" has no "bar_results" struct. Skipping fields.', fname);
%         % still attempt to capture d_slow if present
%     else
%         br = S.bar_results;
% 
%         % Top-level fields in bar_results
%         Date(k)            = asString(getFieldOr([], br, 'Date'));
%         Time(k)            = asString(getFieldOr([], br, 'Time'));
%         Strain(k)          = asString(getFieldOr([], br, 'Strain'));
%         Type(k)            = asString(getFieldOr([], br, 'Type'));
%         median_voltage(k)  = asScalar(getFieldOr(NaN, br, 'median_voltage'));
%         resultant_angle(k) = asScalar(getFieldOr(NaN, br, 'resultant_angle'));
% 
%         % Nested "slow" fields
%         if isfield(br, 'slow') && isstruct(br.slow)
%             sl = br.slow;
%             magnitude(k)  = asScalar(getFieldOr(NaN, sl, 'magnitude'));
%             angle_rad(k)  = asScalar(getFieldOr(NaN, sl, 'angle_rad'));
%             fwhm(k)       = asScalar(getFieldOr(NaN, sl, 'fwhm'));
%             cv(k)         = asScalar(getFieldOr(NaN, sl, 'cv'));
%             thetahat(k)   = asScalar(getFieldOr(NaN, sl, 'thetahat'));
%             kappa(k)      = asScalar(getFieldOr(NaN, sl, 'kappa'));
%             sym_ratio(k)  = asScalar(getFieldOr(NaN, sl, 'sym_ratio'));
%             vector_sum{k} = asArray(getFieldOr([],  sl, 'vector_sum'));
%             DSI_vector(k) = asScalar(getFieldOr(NaN,  sl, 'DSI_vector'));
%             DSI_pdnd(k)   = asScalar(getFieldOr(NaN, sl, 'DSI_pdnd'));
%         end
%     end
% 
%     % d_slow is top-level in the .mat file
%     if isfield(S, 'd_slow') && isnumeric(S.d_slow)
%         d_slow{k} = S.d_slow;
% 
%         % Safely extract requested values from column 2
%         v_max(k)    = getDSValue(S.d_slow, 5,  2);
%         v_null(k)   = getDSValue(S.d_slow, 13, 2);
%         v_ortho1(k) = getDSValue(S.d_slow, 1,  2);
%         v_ortho2(k) = getDSValue(S.d_slow, 9,  2);
%     else
%         d_slow{k} = [];
%         % leave v_* as NaN
%     end
% end
% 
% T = table( ...
%     Date, Time, Strain, Type, ...
%     median_voltage, resultant_angle, ...
%     magnitude, angle_rad, fwhm, cv, thetahat, kappa, sym_ratio, ...
%     vector_sum, DSI_vector, DSI_pdnd, ...
%     v_max, v_null, v_ortho1, v_ortho2, ...  % <- new columns
%     d_slow, ...
%     'VariableNames', { ...
%         'Date','Time','Strain','Type', ...
%         'median_voltage','resultant_angle', ...
%         'magnitude','angle_rad','fwhm','cv','thetahat','kappa','sym_ratio', ...
%         'vector_sum','DSI_vector','DSI_pdnd', ...
%         'v_max','v_null','v_ortho1','v_ortho2', ...
%         'd_slow' ...
%     } ...
% );
% 
% % Optional: preview
% disp(T);
% 
% % Optional: save as MAT
% % save('summary_table.mat', 'T');
% 
% %% --------- Helper functions (can live at end of script) ----------------
% function val = getFieldOr(defaultVal, S, fieldName)
%     if isstruct(S) && isfield(S, fieldName)
%         val = S.(fieldName);
%     else
%         val = defaultVal;
%     end
% end
% 
% function out = asScalar(x)
%     if isnumeric(x) && isscalar(x)
%         out = x;
%     elseif islogical(x) && isscalar(x)
%         out = double(x);
%     elseif isstring(x) || ischar(x)
%         num = str2double(string(x));
%         if isnan(num)
%             out = NaN;
%         else
%             out = num;
%         end
%     else
%         out = NaN;
%     end
% end
% 
% function out = asArray(x)
%     if isnumeric(x)
%         out = x;
%     else
%         out = [];
%     end
% end
% 
% function s = asString(x)
%     if isstring(x) || ischar(x)
%         s = string(x);
%     elseif isdatetime(x)
%         s = string(x);
%     elseif isnumeric(x) && isscalar(x)
%         try
%             s = string(datetime(x, 'ConvertFrom', 'datenum'));
%         catch
%             s = string(x);
%         end
%     else
%         s = "";
%     end
% end
% 
% function v = getDSValue(d, r, c)
%     % Safely fetch d(r,c) if it exists; otherwise return NaN.
%     if isnumeric(d) && ndims(d) == 2 && size(d,1) >= r && size(d,2) >= c
%         v = d(r,c);
%     else
%         v = NaN;
%     end
% end
% 
% end


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

% Derived from d (or d_slow)
v_max    = nan(n,1);
v_null   = nan(n,1);
v_ortho1 = nan(n,1);
v_ortho2 = nan(n,1);

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

    % --- d array (new format) or d_slow (old format) ---
    D = [];
    if isfield(S, 'd') && isnumeric(S.d)
        D = S.d;
    elseif isfield(S, 'd_slow') && isnumeric(S.d_slow)
        D = S.d_slow;
    end

    if ~isempty(D)
        d_slow{k}  = D;
        v_max(k)    = getDSValue(D, 5,  2);
        v_null(k)   = getDSValue(D, 13, 2);
        v_ortho1(k) = getDSValue(D, 1,  2);
        v_ortho2(k) = getDSValue(D, 9,  2);
    else
        d_slow{k} = [];
    end
end

T = table( ...
    Date, Time, Strain, Type, ...
    median_voltage, resultant_angle, ...
    magnitude, angle_rad, fwhm, cv, thetahat, kappa, sym_ratio, ...
    vector_sum, DSI_vector, DSI_pdnd, ...
    v_max, v_null, v_ortho1, v_ortho2, ...
    d_slow, ...
    'VariableNames', { ...
        'Date','Time','Strain','Type', ...
        'median_voltage','resultant_angle', ...
        'magnitude','angle_rad','fwhm','cv','thetahat','kappa','sym_ratio', ...
        'vector_sum','DSI_vector','DSI_pdnd', ...
        'v_max','v_null','v_ortho1','v_ortho2', ...
        'd_slow' ...
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
