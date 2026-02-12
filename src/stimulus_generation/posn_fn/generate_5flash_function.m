function total_dur = generate_5flash_function(peak_frame, n_flashes, flash_dur, bkg_dur, function_folder)
% GENERATE_5FLASH_FUNCTION  Create position function for test flash sequence.
%
%   TOTAL_DUR = GENERATE_5FLASH_FUNCTION(PEAK_FRAME, N_FLASHES, FLASH_DUR, ...
%       BKG_DUR, FUNCTION_FOLDER)
%   generates a simple position function that repeatedly presents a single
%   flash at the specified frame location. Used for quick testing of RF
%   location before running the full protocol.
%
%   INPUTS:
%     peak_frame      - Frame number from Protocol 1 pattern to display
%     n_flashes       - Number of flash repetitions to present
%     flash_dur       - Duration of each flash in seconds
%     bkg_dur         - Duration of interval between flashes in seconds
%     function_folder - Directory to save the position function
%
%   OUTPUT:
%     total_dur - Total duration of the function in seconds
%
%   PURPOSE:
%     Allows quick verification that the identified RF location is correct
%     by presenting repeated flashes at a single position. Useful for
%     confirming neural responses before committing to the full Protocol 2.
%
%   FUNCTION STRUCTURE:
%     [bkg_dur] - [flash_dur] - [bkg_dur] - [flash_dur] - ... (repeated n_flashes times)
%
%   NOTES:
%     - Uses SAVE_FUNCTION_G4_TEST which overwrites existing function
%     - Function ID is always set to 1
%     - Called by PRESENT_FLASHES for quick RF verification
%
%   See also PRESENT_FLASHES, SAVE_FUNCTION_G4_TEST, GENERATE_FLASH_FUNCTION

%% Make func

% 1 - Generate array of ones for the 'interval' when the background frame
% is shown. Background frame in the pattern should always be frame 1.

% 0.5s = 500ms = [1 x 250] array of numbers. 
n_frames_bkg = (bkg_dur * 1000) / 2; % convert from s to ms. % one entry every 2ms.
bkg_array = ones(1, n_frames_bkg);

% 2 - Generate func array for when an individual flash stimulus is being 
% shown. This is just an array of ones, of the size of however many 2ms
% steps are within the desired flash duration. For example, if you want to
% present a flash for 500ms, then this would be a [1 x 1000] array of 1s.
n_frames_flash = (flash_dur * 1000)/2; 
flash_array = ones(1, n_frames_flash);

% Start with the background frame
func = bkg_array;

for fl = 1:n_flashes
    array_to_add = flash_array*peak_frame;
    func = cat(2, func, array_to_add, bkg_array);
end 

%% Save the position function.
param.func = func;
function_type = 'pfn';
ID = 1; % Always set to 1 - overwrite every time.
param.ID = ID;

param.frames = numel(func); % - - -- CHECK THIS. 
param.gs_val = 4;
param.type = function_type;

% add duration of the stimulus
total_dur = (numel(func)*2)/1000; % func consists of 1 element every 2ms. n_elements*2 = time in ms. Divide by 100 to get time in seconds. %(n_flashes*flash_dur + (n_flashes)*interval_dur); 
param.dur = total_dur;

fname = "Func";
% Uses a different version of "save_function_G4" that will overwrite the
% function in the folder.
save_function_G4_test(func, param, function_folder, fname);

end