function [x, y, on_off] = frame_to_coord(peak_frames, bkg_color, threshold_distance, hemi)
% FRAME_TO_COORD  Convert multiple peak frames to consensus RF coordinates.
%
%   [X, Y, ON_OFF] = FRAME_TO_COORD(PEAK_FRAMES, BKG_COLOR, THRESHOLD_DISTANCE, HEMI)
%   analyzes peak frames across multiple conditions and repetitions to
%   determine a robust estimate of the receptive field center position.
%   Includes outlier detection and removal.
%
%   INPUTS:
%     peak_frames        - 3D array [n_condition, n_rep, 2] containing:
%                          peak_frames(:,:,1) = frame numbers
%                          peak_frames(:,:,2) = peak voltage values
%     bkg_color          - Background pixel intensity value (0-15)
%     threshold_distance - [dist1, dist2] maximum allowed distance in pixels
%                          between flash centroids across repetitions
%                          dist1 = threshold for condition 1 (12px flashes)
%                          dist2 = threshold for condition 2 (6px flashes)
%     hemi               - 'left' or 'right' recording hemisphere
%
%   OUTPUTS:
%     x      - Consensus horizontal pixel coordinate (column)
%     y      - Consensus vertical pixel coordinate (row)
%     on_off - 'on' or 'off' based on majority of responses
%
%   ALGORITHM:
%     1. Loads Protocol 1 patterns for both flash sizes
%     2. Extracts flash centroids for each rep/condition
%     3. Calculates pairwise distances between centroids
%     4. Identifies and removes outliers beyond threshold
%     5. Computes centroid of remaining valid points
%     6. Averages centroids from both conditions
%     7. Determines ON/OFF preference from majority vote
%
%   DIAGNOSTICS:
%     Prints variance in peak voltage, outlier information, and
%     final centroid coordinates to console.
%
%   NOTE:
%     This function is more sophisticated than PATT_FRAME_TO_COORD and
%     handles multiple repetitions with outlier removal. Use when
%     Protocol 1 data shows variable responses across repetitions.
%
%   See also PATT_FRAME_TO_COORD, GENERATE_FAKE_PEAK_FRAMES
%_______________________________________________________________________

%% Load the patterns used: 
if hemi == "left"
    pattern_path = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\LHS\protocol1_4reps_12px_6px_LHS_2sbkg_200msfl_50msint_12-03-24_15-11-40\Patterns';
elseif hemi == "right"
    pattern_path = 'C:\matlabroot\G4_Protocols\nested_RF_stimulus\protocols\RHS\protocol1_4reps_12px_6px_RHS_2sbkg_200msfl_50msint_12-03-24_15-25-60\Patterns';
end 
cd(pattern_path)

pat1 = dir('0001_*');
pattern1 = load(pat1.name, 'pattern');
allf1 = pattern1.pattern.Pats;

pat2 = dir('0002_*');
pattern2 = load(pat2.name, 'pattern');
allf2 = pattern2.pattern.Pats;

%% Print the variances in the peak voltage across the reps per condition.

% Variance in peak voltage between reps for condition 1.
varV1 = var(peak_frames(1, :, 2));
disp(strcat("Variance in voltage across reps for condition 1: ", string(varV1)))
% Variance in peak voltage between reps for condition 2.
varV2 = var(peak_frames(2, :, 2));
disp(strcat("Variance in voltage across reps for condition 2: ", string(varV2)))

%% Find the centre coordinate of the flash that was presented during the peak frame.

n_reps = numel(peak_frames(1, :, 1));

% Go through each rep - extract the frame that was presented - find the
% centre of the flash.
coords = zeros(2, n_reps, 2);
on_off_array = NaN(2, n_reps); % array to store whether the flash that 
% elicited the maximal response was an ON or an OFF flash. 

for c = 1:2 % per condition
    for r = 1:n_reps % per rep
        peakf = peak_frames(c, r, 1);

        % If the peak frame found is not NaN.
        if ~isnan(peakf)
            % Find the coord of the centre of the flash.
            if c == 1
                f = allf1(:, :, peakf);
            elseif c == 2
                f = allf2(:, :, peakf);
            end 
    
            [a, b] = find(f~=bkg_color);
            max_col = max(f);
            if max_col>bkg_color % contains pixels higher than bkg - ON 
                on_off_array(c, r) = 1;
            else 
                on_off_array(c, r) = 0;
            end 
            coords(c,r,1) = ceil(median(a));
            coords(c,r,2) = ceil(median(b));
        end

    end 
end 

% Determine whether the cell is responding to on / off stimuli.
mean_on_off = mean(mean(on_off_array));
disp(["Mean on-off score:", mean_on_off])
if mean_on_off > 0.5
    on_off = "on";
else 
    on_off = "off";
end

%% Find the euclidean distance between the centre of the flashes that were
% presented during the 'peak frames'. 

% Determine whether any of the central points are further than
% 'threshold_distance' away from the other centroids. 
% If there are 'outlier' points, remove them and find the middle point of
% the remaining points to get [x,y]. 

% Should I do this for condition 1, then condition 2, or altogether????????????
% ________________________________________________________________________

for i = 1:2
    th_dist = threshold_distance(i);
    if i == 1
        crds = reshape(coords(1, :, :), [4,2]);
    elseif i == 2
        crds = reshape(coords(2, :, :), [4,2]);
    end 

    % Calculate pairwise distances
    distances = pdist(crds, 'euclidean');      % Compute pairwise distances
    dist_matrix = squareform(distances);         % Convert to symmetric matrix form
    % Additional: plot as figure: 
    % figure; imagesc(dist_matrix);
    
    % Find if any distance exceeds the threshold
    [max_dist, row_idx] = max(dist_matrix, [], 2);   % Max distance from each point
    outlier_indices = find(max_dist > th_dist);  % Identify points exceeding the threshold

    inlier_indices = max(dist_matrix, [], 2) <= th_dist; % Indices of points within the threshold
    remaining_coords = crds(inlier_indices, :);
    centroid = floor(mean(remaining_coords, 1));
    if i == 1
        centroid1 = centroid;
    elseif i == 2 
        centroid2 = centroid;
    end 
    
    % Check if there are any outliers
    disp(['For condition ', num2str(i), ':'])
    if ~isempty(outlier_indices)
        disp(['Points further than ', num2str(th_dist), ' pixels away:']);
        disp(['Indices of points: ', num2str(outlier_indices')]);
        disp('Coordinates of outlier points:');
        disp(crds(outlier_indices, :));
    else
        disp(['All points are within ', num2str(th_dist), ' pixels of each other .']);
    end

    disp(['Central position of flashes presented during peak frames - condition ', num2str(i), ': [', num2str(centroid), ']'])
end 

dist_centroids = pdist2(centroid1, centroid2);
disp(['Distance between centroids from condition 1 & 2: ', num2str(dist_centroids), ' pixels.'])

% Take the mean of the centroids found for condition 1 and condition 2.
x = int16(mean([centroid1(1), centroid2(1)]));
y = int16(mean([centroid1(2), centroid2(2)]));
disp(['Final coordinate to centre stimuli on: [', num2str(x), ',', num2str(y), ']'])

end 