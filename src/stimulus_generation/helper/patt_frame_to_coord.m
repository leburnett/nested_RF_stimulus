function [x, y, on_off] = patt_frame_to_coord(peak_frame, arena_side)

% Find the [x,y] coordinate of the screen from the identifying which frames 
% of protocol 1 the cell responded to best.

% 'peak_frame' - [n_condition, n_rep, [peak_frame, peak_voltage]]

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