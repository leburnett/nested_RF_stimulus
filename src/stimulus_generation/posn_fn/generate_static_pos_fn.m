function generate_static_pos_fn(save_dir, dur_t)
% GENERATE_STATIC_POS_FN  Create position function for static background display.
%
%   GENERATE_STATIC_POS_FN(SAVE_DIR, DUR_T) generates a position function
%   that displays frame 1 (background) for a specified duration. Used to
%   present a uniform gray screen before starting the flash sequence.
%
%   INPUTS:
%     save_dir - Directory path to save the position function files
%     dur_t    - Duration in seconds to display the static background
%
%   PURPOSE:
%     Provides a baseline period at the start of each protocol where only
%     the background is shown. This allows the neural response to stabilize
%     before stimulus presentation begins.
%
%   OUTPUT FILES:
%     - <dur_t>s_static_frame1_24frames.pfn - Binary function file
%     - <dur_t>s_static_frame1_24frames_G4.mat - MATLAB structure
%
%   PARAMETERS:
%     Frame value: 1 (background frame)
%     Grayscale: 4-bit (gs_val = 4)
%     Section type: 'static' (constant frame display)
%
%   See also GENERATE_STATIC_FUNCTION, GENERATE_FLASH_FUNCTION,
%            GENERATE_PROTOCOL1_STIMULI

    %% user-defined function parameters
    pfnparam.type = 'pfn'; %number of frames in pattern
    pfnparam.frames = 24; %number of frames in pattern
    pfnparam.gs_val = 4; %brightness bits in pattern
    pfnparam.section = { 'static' }; %static, sawtooth, traingle, sine, cosine, or square
    pfnparam.dur = [ dur_t ]; %section duration (in s)
    pfnparam.val = [ 1 ]; %function value for static sections
    pfnparam.high = [ 24 ]; %high end of function range {for non-static sections}
    pfnparam.low = [ 1 ]; %low end of function range {for non-static sections}
    pfnparam.freq = [ 3 ]; %frequency of section {for non-static sections}
    pfnparam.size_speed_ratio = [ 40 ]; %size/speed ratio {for looms}
    pfnparam.flip = [ 0 ]; %flip the range of values of function {for non-static sections}
    
    %% generate function
    func = G4_Function_Generator(pfnparam);
    
    %% save function
    pfnparam.ID = get_function_ID('pfn',save_dir);
    filename = strcat(string(dur_t), 's_static_frame1_24frames');
    save_function_G4(func, pfnparam, save_dir, filename);
end 