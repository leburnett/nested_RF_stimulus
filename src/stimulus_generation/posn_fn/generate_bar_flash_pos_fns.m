function generate_bar_flash_pos_fns(save_dir, n_flank)

% 80ms flashes + 920ms interval - - - match 28 dps moving bars. 

% For each bar pattern, generate a random order 

n_orientations = 8;
n_frames_orient = (n_flank*2 + 1);
n_frames_pattern = n_frames_orient * n_orientations;

flash_dur = 0.08;
flash_str = "0-08";
int_dur = 0.92; 
int_str = "0-92";

% Make 3 different position functions, one for each rep.
for r = 1:n_reps

    % Random order in which to present flashes.
    rand_order = randperm(88);

    value_array = zeros(1, 2 * numel(rand_order)); 
    value_array(1:2:end) = 1; 
    value_array(2:2:end) = rand_order; 

    % user-defined function parameters
    pfnparam.type = 'pfn'; 
    pfnparam.frames = n_frames_pattern; %number of frames in pattern
    pfnparam.gs_val = 4; %brightness bits in pattern
    pfnparam.section = repmat({'static'}, 1, n_frames_pattern*2); 
    pfnparam.dur = repmat([int_dur flash_dur], 1, n_frames_pattern); %section duration (in s)
    pfnparam.val = value_array; %function value for static sections
    
    %% generate function
    func = G4_Function_Generator(pfnparam);
    
    %% save function
    pfnparam.ID = get_function_ID('pfn',save_dir);
    filename = strcat("bar_flashes_", flash_str, 's_', int_str, 's_', string(n_frames_pattern), 'frames-Rep', string(r));
    save_function_G4(func, pfnparam, save_dir, filename);

end 


end 

