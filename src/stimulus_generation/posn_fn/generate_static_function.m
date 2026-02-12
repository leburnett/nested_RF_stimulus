function generate_static_function(func_save_dir, dur)

% GENERATE_STATIC_FUNCTION  Create 10s static background function for Protocol 2.
%
%   GENERATE_STATIC_FUNCTION(FUNC_SAVE_DIR) generates a position function
%   that displays frame 1 (gray background) for 10 seconds. Used at the
%   start of Protocol 2 to provide a baseline period.
%
%   INPUT:
%     func_save_dir - Directory to save the position function files
%
%   PURPOSE:
%     Called by CREATE_PROTOCOL2 to add a static background display at
%     the beginning of the experiment. This allows the neuron's baseline
%     activity to stabilize before stimulus presentation.
%
%   OUTPUT FILES:
%     - 10s_static_frame1_24frames.pfn - Binary function file
%     - 10s_static_frame1_24frames_G4.mat - MATLAB structure
%
%   PARAMETERS:
%     Duration: 10 seconds
%     Frame value: 1 (background frame)
%     Grayscale: 4-bit (gs_val = 4)
%
%   See also CREATE_PROTOCOL2, GENERATE_STATIC_POS_FN
    
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