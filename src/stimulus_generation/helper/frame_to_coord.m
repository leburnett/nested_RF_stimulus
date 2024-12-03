function [x, y, on_off] = frame_to_coord(peak_frames, bkg_color, threshold_distance, hemi)
% Find the [x,y] coordinate of the screen from the identifying which frames 
% of protocol 1 the cell responded to best.

% 'peak_frames' - [n_condition, n_rep, [peak_frame, peak_voltage]]
% bkg_color - int - pixel intensityvalue used as the background color in the
% patterns. 
% 'threshold_distance' - tuple - the maximum desired distance (in pixels) between
% flash centroids across all of the peak frames found across the
% repeitions. Tuple because the first number is for condition 1, 12 pixel
% flashes, and the second is for condition 2 with pixel flashes.
% 'hemi' - hemisphere that is being recorded from. 
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
x = mean(centroid1(1), centroid2(1));
y = mean(centroid1(2), centroid2(2));
disp(['Final coordinate to centre stimuli on: [', num2str(x), ',', num2str(y), ']'])

end 