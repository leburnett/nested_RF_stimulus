% generate fake peak_data

% Dimensions
cond_dim = 2;  % 1 or 2
rep_dim = 4;   % 1, 2, or 3
% Generate fake data
peak_frames = zeros(cond_dim, rep_dim, 1, 1); % Preallocate for 2 peak values (frame & voltage)
for cond = 1:cond_dim
    for rep = 1:rep_dim
            peak_frame = randi([1, 192]);      % Random peak frame
            peak_voltage = randi([40, 60]);    % Random peak voltage
            % Store in 4D array
            peak_frames(cond, rep, 1) = peak_frame;
            peak_frames(cond, rep, 2) = peak_voltage;
    end
end
