function total_dur = generate_5flash_function(peak_frame, n_flashes, flash_dur, bkg_dur, function_folder)

% peak_frame - frame number from pattern displayed in protocol 1 for which
% you would like to display the flash. 
% n_flashes - the number of flashes to be presented.
% flash_dur - duration of flash in seconds.
% bkg_dur - Duration of interval between flashes in seconds. 
% function_folder - where to save the position function.

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