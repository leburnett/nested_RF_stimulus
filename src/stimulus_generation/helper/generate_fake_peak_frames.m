% generate fake peak_data

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
