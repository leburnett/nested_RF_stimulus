function generate_bar_pos_fns(save_dir, px_crop)
% For bar patterns - want position function to move from frames 1 - 68.
% Create functions for the forward and 'flip' direction. 
% 3 different speeds - 14 dps, 28 dps and 56 dps.

% July 2025 - added 1s of static first frame before each
% sawtooth.

    for dps = [28, 56]
    
        if px_crop == 45 

            if dps == 14
                dur = 6.08;
                dur_str = '6-08';
                freq = 0.16471;
            elseif dps == 28
                dur = 3.04;
                dur_str = '3-04';
                freq = 0.32941;
            elseif dps == 56 
                dur = 1.52;
                dur_str = '1-52';
                freq = 0.65882;
            end
            
            fr_low = 1;
            fr_high = 68;

        elseif px_crop == 30

            if dps == 28
                dur = 2.328;
                dur_str = '2-328';
                freq = 0.43077;
            elseif dps == 56 
                dur = 1.164;
                dur_str = '1-164';
                freq = 0.86154;
            end

            fr_low = 11;
            fr_high = 62;

        end 
    
        for flip = [0, 1]
            if flip == 1
                flip_str = 'FLIP_';
            else 
                flip_str = '';
            end
    
            % user-defined function parameters
            pfnparam.type = 'pfn'; %number of frames in pattern
            pfnparam.frames = 288; %number of frames in pattern
            pfnparam.gs_val = 4; %brightness bits in pattern
            pfnparam.section = {'static' 'sawtooth'}; %static, sawtooth, traingle, sine, cosine, or square
            pfnparam.dur = [1 dur]; %section duration (in s)
            pfnparam.val = [1 1]; %function value for static sections
            pfnparam.high = [62 fr_high]; %high end of function range {for non-static sections}
            pfnparam.low = [1 fr_low ]; %low end of function range {for non-static sections}
            pfnparam.freq = [1 freq ]; %frequency of section {for non-static sections}
            pfnparam.size_speed_ratio = [40 40]; %size/speed ratio {for looms}
            pfnparam.flip = [0 flip]; %flip the range of values of function {for non-static sections}
            
            %% generate function
            func = G4_Function_Generator(pfnparam);
            
            %% save function
            pfnparam.ID = get_function_ID('pfn',save_dir);
            filename = strcat(string(dps), 'dps_', flip_str, dur_str, 's_', string(fr_low), '-', string(fr_high), '_288fr_1-25step_3s_static');
            save_function_G4(func, pfnparam, save_dir, filename);
        end
    end 

end 