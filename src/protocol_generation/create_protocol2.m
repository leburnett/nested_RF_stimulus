function [pattern_order, func_order, trial_dur] = create_protocol2(exp_folder)

%% check for correct folders
pattern_folder = fullfile(exp_folder, 'Patterns');
function_folder = fullfile(exp_folder, 'Functions');

assert(exist(pattern_folder,'dir')==7,'did not detect pattern folder')
assert(exist(function_folder,'dir')==7,'did not detect position function folder')

%% read patterns
matinfo = dir(fullfile(pattern_folder, "*.mat"));
patinfo = dir(fullfile(pattern_folder, '*.pat'));
% num_files is the number of patterns.
num_files = length({patinfo.name});

flash_patt = [1, 2];
bar_patt = flash_patt(end)+1:num_files;
n_bar_patt = numel(bar_patt);
n_dir = 2;
n_speeds = 2;
% 2 reps = forward and reverse direction
% 3 reps = n_speeds
pattern_order = [flash_patt, repmat(repelem(bar_patt, n_dir), [1, n_speeds])];
n_patts = numel(pattern_order);

for p = 1:n_patts
    f = pattern_order(p);
    currentExp.pattern.pattNames{p} = matinfo(f).name;
    currentExp.pattern.patternList{p} = patinfo(f).name;
    patternIN = load(fullfile(pattern_folder, matinfo(f).name));
    currentExp.pattern.x_num(p) = patternIN.pattern.x_num; % Matlab now has its own built in class "pattern" which conflicts here.
    currentExp.pattern.y_num(p) = patternIN.pattern.y_num;
    currentExp.pattern.gs_val(p) = patternIN.pattern.gs_val;
    currentExp.pattern.arena_pitch(p) = round(rad2deg(patternIN.pattern.param.arena_pitch));
end
currentExp.pattern.num_patterns = num_files;


%% read position functions
matinfo = dir(fullfile(function_folder, "*.mat"));
pfninfo = dir(fullfile(function_folder, "*.pfn"));
num_files = length({pfninfo.name});
total_dur = 0;
trial_dur = nan([1, n_patts]);

fns_flash = [1,2];
fns_28dps = [3,4];
fns_56dps = [5,6];

func_order = [fns_flash, repmat(fns_28dps, [1,n_bar_patt]), repmat(fns_56dps, [1,n_bar_patt])];

for p = 1:n_patts
    f = func_order(p);
    currentExp.function.functionName{p} = matinfo(f).name;
    currentExp.function.functionList{p} = pfninfo(f).name;
    load(fullfile(function_folder, matinfo(f).name))
    currentExp.function.functionSize(p) = pfnparam.size;
    currentExp.function.trial_dur(p) = pfnparam.dur;
    total_dur = total_dur + pfnparam.dur;
    trial_dur(p) = pfnparam.dur;
end
currentExp.totalDuration = total_dur;
currentExp.function.numFunc = num_files;

%% save currentExp structure
save(fullfile(exp_folder, 'currentExp.mat'),'currentExp', 'pattern_order', 'func_order', 'trial_dur');

end