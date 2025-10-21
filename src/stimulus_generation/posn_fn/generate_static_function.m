function generate_static_function(func_save_dir, dur)
% Generate a static position function for presenting the grey screen at the
% beginning of the protocol.
    
    pfnparam.type = 'pfn'; %number of frames in pattern
    pfnparam.frames = 24; %number of frames in pattern
    pfnparam.gs_val = 4; %brightness bits in pattern
    pfnparam.section = { 'static' }; %static, sawtooth, traingle, sine, cosine, or square
    pfnparam.dur = [ dur ]; %section duration (in s)
    pfnparam.val = [ 1 ]; %function value for static sections
    pfnparam.high = [ 24 ]; %high end of function range {for non-static sections}
    pfnparam.low = [ 1 ]; %low end of function range {for non-static sections}
    pfnparam.freq = [ 3 ]; %frequency of section {for non-static sections}
    pfnparam.size_speed_ratio = [ 40 ]; %size/speed ratio {for looms}
    pfnparam.flip = [ 0 ]; %flip the range of values of function {for non-static sections}
    
    %% generate function
    func = G4_Function_Generator(pfnparam);
    
    pfnparam.ID = get_function_ID('pfn',func_save_dir);
    filename = strcat(string(dur), 's_static_frame1_24frames');
    save_function_G4(func, pfnparam, func_save_dir, filename);

end 