% GENERATE_FAKE_PEAK_FRAMES  Create synthetic peak frame data for testing.
%
%   This script generates random peak_frames data for testing Protocol 2
%   generation without requiring actual Protocol 1 recordings.
%
%   OUTPUT:
%     peak_frames - 4D array [cond_dim, rep_dim, 1, 2] containing:
%                   peak_frames(:,:,:,1) = frame numbers
%                   peak_frames(:,:,:,2) = peak voltage values
%
%   DATA STRUCTURE:
%     cond_dim = 2: Two conditions (12px and 6px flashes)
%     rep_dim = 4:  Four repetitions per condition
%     Condition 1: Frame numbers randomly drawn from 2-65
%     Condition 2: Frame numbers randomly drawn from 2-257
%     Peak voltages: Random values 1-90 (arbitrary units)
%
%   USAGE:
%     Run this script before testing FRAME_TO_COORD or other functions
%     that require peak_frames input. The generated data simulates
%     noisy RF measurements across multiple repetitions.
%
%   See also FRAME_TO_COORD, PATT_FRAME_TO_COORD

% Dimensions
cond_dim = 2;  % 1 or 2
rep_dim = 4;   % 1, 2, or 3
% Generate fake data
peak_frames = zeros(cond_dim, rep_dim, 1, 1); % Preallocate for 2 peak values (frame & voltage)
for cond = 1:cond_dim
    for rep = 1:rep_dim
            if cond == 1
                peak_frame = randi([2, 65]);      % Random peak frame
            elseif cond == 2
                peak_frame = randi([2, 257]);      % Random peak frame
            end 
            peak_voltage = randi([1, 90]);    % Random peak voltage
            % Store in 4D array
            peak_frames(cond, rep, 1) = peak_frame;
            peak_frames(cond, rep, 2) = peak_voltage;
    end
end
