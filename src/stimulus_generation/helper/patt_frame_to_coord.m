function [x, y, on_off] = patt_frame_to_coord(peak_frame, arena_side)
% PATT_FRAME_TO_COORD  Convert Protocol 1 peak frame to screen coordinates.
%
%   [X, Y, ON_OFF] = PATT_FRAME_TO_COORD(PEAK_FRAME, ARENA_SIDE)
%   determines the screen position [x,y] and contrast polarity from a
%   Protocol 1 frame number. This is the primary coordinate conversion
%   function used when setting up Protocol 2.
%
%   INPUTS:
%     peak_frame  - Frame number (1-256) from Protocol 1 that elicited
%                   the strongest neural response
%     arena_side  - 'L' or 'R' indicating which arena hemisphere was used
%                   'L' = left hemisphere recording (fly's left visual field)
%                   'R' = right hemisphere recording (fly's right visual field)
%
%   OUTPUTS:
%     x      - Horizontal pixel coordinate (column) of flash center
%     y      - Vertical pixel coordinate (row) of flash center
%     on_off - 'on' if flash was bright, 'off' if flash was dark
%
%   ALGORITHM:
%     1. Loads the 6px Protocol 1 pattern file for the specified arena side
%     2. Extracts the specified frame from the pattern
%     3. Finds pixels different from background to locate the flash
%     4. Calculates centroid of flash location
%     5. Determines ON/OFF based on whether flash pixels are brighter
%        or darker than background
%
%   PATTERN FILES:
%     Left hemisphere:  protocols/LHS/protocol1_.../Patterns/0002_*.mat
%     Right hemisphere: protocols/RHS2/protocol1_.../Patterns/0002_*.mat
%
%   EXAMPLE:
%     % Convert frame 45 from left hemisphere recording
%     [x, y, on_off] = patt_frame_to_coord(45, 'L');
%     % Returns: x=85, y=24, on_off='on'
%
%   See also GENERATE_PROTOCOL2, FRAME_TO_COORD, FIND_INV_PEAK_FRAME
%_______________________________________________________________________

%% Load the patterns used depending on which arena half the pattern was presented on: 
assert(ismember(arena_side, {'L', 'R'}), 'screen_hemi must be either "L" or "R"')

% Assuming using left at the moment:
if arena_side == "L"
    pattern_path = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\LHS\protocol1_10kHz_4reps_12px_6px_LHS_2sbkg_200msfl_50msint_12-13-24_14-33-03\Patterns';
elseif arena_side == "R"
    % pattern_path = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\RHS\protocol1_10kHz_4reps_12px_6px_RHS_2sbkg_200msfl_50msint_04-08-25_08-08-42\Patterns'; 
    pattern_path = "C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\RHS2\protocol1_10kHz_4reps_12px_6px_RHS2_2sbkg_200msfl_50msint_81_180_05-05-25_16-18-66\Patterns";
end 

cd(pattern_path)
% Load the pattern with the smaller 6px square flashes:
pat2 = dir('0002_*');
pattern2 = load(pat2.name, 'pattern');
allf2 = pattern2.pattern.Pats;

%% Find the centre coordinate of the flash that was presented during the peak frame.

% If the peak frame found is not NaN.
if ~isnan(peak_frame)

    % Use the pattern for condition 2 with the smaller flashes. 
    f = allf2(:, :, peak_frame);

    bkg_color = mode(mode(f)); % The background colour will be the most common pixel value in the frame. 

    [a, b] = find(f~=bkg_color); 
    max_col = max(max(f));
    if max_col>bkg_color % contains pixels higher than bkg - ON 
        on_off = 'on';
    else 
        on_off = 'off';
    end 
    y = int16(median(a));
    x = int16(median(b));
end

disp(['Final coordinate to centre stimuli on: [', num2str(x), ',', num2str(y), '] and ', on_off, ' flashes.'])

end 