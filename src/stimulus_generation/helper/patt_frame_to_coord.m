function [x, y, on_off] = patt_frame_to_coord(peak_frame, bkg_color)
% Find the [x,y] coordinate of the screen from the identifying which frames 
% of protocol 1 the cell responded to best.

% 'peak_frame' - [n_condition, n_rep, [peak_frame, peak_voltage]]
% bkg_color - int - pixel intensityvalue used as the background color in the
% patterns. 
%_______________________________________________________________________

%% Load the patterns used: 

% Assuming using left at the moment:
pattern_path = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\LHS\protocol1_10kHz_4reps_12px_6px_LHS_2sbkg_200msfl_50msint_12-13-24_14-33-03\Patterns';
cd(pattern_path)

pat2 = dir('0002_*');
pattern2 = load(pat2.name, 'pattern');
allf2 = pattern2.pattern.Pats;

%% Find the centre coordinate of the flash that was presented during the peak frame.
peakf = peak_frame;

% If the peak frame found is not NaN.
if ~isnan(peakf)

    % Use the pattern for condition 2 with the smaller flashes. 
    f = allf2(:, :, peakf);

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