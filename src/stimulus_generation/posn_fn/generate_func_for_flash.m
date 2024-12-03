function func = generate_func_for_flash(bkg_frame, bkg_dur, flash_seq, flash_dur, on_off)
% Generate the 'func' position function array to be read by the controller.
% This array consists of integer values which refer to the pattern id to be
% shown at each 2ms interval because the system runs at 500Hz (updates 
% every 2ms) when showing 4-bit patterns.

% Inputs
% ______

% 'bkg_frame' - frame of pattern which contains the background frame.
% Should normally be the first frame of the pattern.

% 'bkg_dur' - duration in seconds for how long to present this background
% frame between flash stimuli

% 'flash_seq' - a [1 x n_flashes] array of the order in which to display
% the flashes. The index number refers to which frame of the pattern to read in. 

% 'flash_dur' - duration in seconds for which to display each flash.

% 'on_off' - whether to show just ON flashes, just OFF flashes or both.
% Possible options are 'on', 'off' or 'both'.
% ______________________________________________________________

% 1 - Generate array of ones for the 'interval' when the background frame
% is shown. Background frame in the pattern should always be frame 1.

% 0.5s = 500ms = [1 x 250] array of numbers. 
n_frames_bkg = (bkg_dur * 1000) / 2; % convert from s to ms. % one entry every 2ms.
bkg_array = ones(1, n_frames_bkg)*bkg_frame;

% 2 - Generate func array for when an individual flash stimulus is being 
% shown. This is just an array of ones, of the size of however many 2ms
% steps are within the desired flash duration. For example, if you want to
% present a flash for 500ms, then this would be a [1 x 1000] array of 1s.
n_frames_flash = (flash_dur * 1000)/2; 
flash_array = ones(1, n_frames_flash);

% 3 - Generate the overall 'func' array.

% Number of flashes - for only ON or OFF flashes. The total number of 
% flashes will be 2 * n_flashes if 'on_off' = 'both'. 
n_flashes = numel(flash_seq);

% Start with the background frame
func = bkg_array;

% Then run through all of the ON and OFF flashes. Add 'bkg_array' between
% flashes which will be the interval stimulus shown for 'bkg_dur' seconds.
if on_off == "off"

    for fl = 1:n_flashes
        frame_idx = flash_seq(fl)+1;
        array_to_add = flash_array*frame_idx;
        func = cat(2, func, array_to_add, bkg_array);
    end 

elseif on_off == "on"

    for fl = 1:n_flashes
        frame_idx = flash_seq(fl-n_flashes)+n_flashes+1; 
        array_to_add = flash_array*frame_idx;
        func = cat(2, func, array_to_add, bkg_array);
    end 

elseif on_off == "both"

    % %  Displaying all OFF then all ON flashes in order:
    % for fl = 1:n_flashes*2
    % 
    %     if fl<=n_flashes
    %         frame_idx = flash_seq(fl)+1;
    %     elseif fl >n_flashes
    %         frame_idx = flash_seq(fl-n_flashes)+n_flashes+1;
    %     end 
    % 
    %     array_to_add = flash_array*frame_idx;
    % 
    %     func = cat(2, func, array_to_add, bkg_array);
    % end 

    % % Displaying alternate ON and OFF flashes.

    % Alternate between on and off flashes in sequence
    a = 1:2:n_flashes;                      % Odd indices
    b = n_flashes+2:2:n_flashes*2;          % Even indices (second sequence)
    c = 2:2:n_flashes;                      % Even indices
    d = n_flashes+1:2:n_flashes*2;          % Odd indices (second sequence)
    % Interleave and combine sequences
    v1 = reshape([a; b], 1, []);            % Interleaving for the first set
    v2 = reshape([d; c], 1, []);            % Interleaving for the second set
    v = [v1, v2];                           % Combine interleaved arrays
    % Generate full flash sequence for quadrants
    flash_seq2 = [flash_seq, flash_seq + n_flashes];  % Combine original and offset sequences
    % Loop through and generate the full sequence
    for fl = 1:n_flashes*2
        idx = v(fl);
        frame_idx = flash_seq2(idx) + 1;          
        array_to_add = flash_array * frame_idx; 
        func = cat(2, func, array_to_add, bkg_array);
    end

end 



























