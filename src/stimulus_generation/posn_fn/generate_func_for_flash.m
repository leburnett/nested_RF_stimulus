
function func = generate_func_for_flash(bkg_frame, bkg_dur, flash_seq, flash_dur)

% Generate the 'func' position function to be read by the controller. 
% 'bkg_frame' - frame of pattern which contains the background frame.
% Should normally be the first frame of the pattern.

% 'bkg_dur' - duration in seconds for how long to present this background
% frame between flash stimuli

% 'flash_seq' - a [1 x n_flashes] array of the order in which to display
% the flashes. The index number refers to which frame of the pattern to read in. 

% 'flash_dur' - duratio in seconds for which to display each flash.

% The system when display 4-bit patterns runs at 500Hz (updates every 2ms).

% ______________________________________________________________

% Generate part of func array for background.
% 0.5s = 500ms = [1 x 250] array of numbers. 
n_frames_bkg = (bkg_dur * 1000) / 2; % convert from s to ms. % one entry every 2ms.
bkg_array = ones(1, n_frames_bkg)*bkg_frame;

% Generate func array for flashes. 
n_frames_flash = (flash_dur * 1000)/2; 
flash_array = ones(1, n_frames_flash);

% Number of flashes 
n_flashes = numel(flash_seq);

% Generate 'func'
% func_sz = (n_frames_bkg*n_flashes+1)+(n_frames_flash+n_flashes);

% Start with the background frame
func = bkg_array;

% Then run through all of the ON and OFF flashes
for fl = 1:n_flashes*2

    if fl<=n_flashes
        frame_idx = flash_seq(fl);
    elseif fl >n_flashes
        frame_idx = flash_seq(fl-n_flashes)+n_flashes; % Should be  up to 64! 
    end 

    array_to_add = flash_array*frame_idx;

    func = cat(2, func, array_to_add, bkg_array);
end 

























