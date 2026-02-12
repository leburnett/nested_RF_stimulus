function [d, ord, magnitude, angle_rad, fwhm, cv, thetahat, kappa] = find_PD_and_order_idx(max_v_polar, median_voltage)
% FIND_PD_AND_ORDER_IDX  Compute preferred direction and tuning metrics.
%
%   [D, ORD, MAGNITUDE, ANGLE_RAD, FWHM, CV, THETAHAT, KAPPA] = ...
%       FIND_PD_AND_ORDER_IDX(MAX_V_POLAR, MEDIAN_VOLTAGE)
%   calculates the preferred direction using vector sum and computes
%   multiple metrics characterizing the directional tuning curve.
%
%   INPUTS:
%     max_v_polar    - 17x1 array of max voltage responses for 16 directions
%                      (last element duplicates first for circular plotting)
%     median_voltage - Baseline voltage for subtraction
%
%   OUTPUTS:
%     d         - 16x2 array: [aligned_angles, median-subtracted_responses]
%                 Data rotated so preferred direction aligns to pi/2
%     ord       - 16x1 array of indices for reordering data to align PD
%     magnitude - Normalized vector sum magnitude (0-1, higher=more selective)
%     angle_rad - Preferred direction in radians (0 to 2*pi)
%     fwhm      - Full-width half-maximum of tuning curve (degrees)
%     cv        - Circular variance (0=sharp, 1=broad tuning)
%     thetahat  - Von Mises mean direction parameter
%     kappa     - Von Mises concentration parameter (higher=sharper tuning)
%
%   VECTOR SUM METHOD:
%     Each direction contributes a vector with magnitude=response and
%     angle=direction. The sum of all vectors gives the preferred direction
%     (angle) and selectivity strength (normalized magnitude).
%
%   ALIGNMENT:
%     The output 'd' is rotated so the preferred direction appears at pi/2.
%     This allows comparison across cells with different preferred directions.
%     'ord' provides the index mapping for this rotation.
%
%   DEPENDENCIES:
%     Uses circ_vmpar() from the Circular Statistics Toolbox
%
%   See also COMPUTE_FWHM, COMPUTE_CIRCULAR_VAR, COMPUTE_BAR_RESPONSE_METRICS
% ________________________________________________________________________________________

    % Angles in radians to be used.
    angls = linspace(0, 2*pi, 17)';
    angls = angls(1:end-1); 

    % Median-subtracted peak voltages - remove last value to not double
    % weight.
    responses = max_v_polar(1:end-1); %-median_voltage;  % peak neural responses
    
    %% 1 - Find vector sum of responses to find PD. 

    % Compute vector components
    x_component = responses.*cos(angls);
    y_component = responses.*sin(angls);

    % Compute vector sum
    vector_sum_x = sum(x_component)';
    vector_sum_y = sum(y_component)';

    % Compute magnitude and direction
    magnitude = sqrt(vector_sum_x.^2 + vector_sum_y.^2)/sum(responses); % Normalise between 0 and 1. 

    % This is the angle in which the vector sum points.
    angle_rad = atan2(vector_sum_y, vector_sum_x);

    % Different methods of finding out how wide the response is:
    fwhm = compute_FWHM(angls, responses);
    cv = compute_circular_var(angls, responses);
    [thetahat, kappa] = circ_vmpar(angls, responses);

    % Find the closest direction in 'directions' to the vector sum angle. 
    if angle_rad < 0 
        angle_rad = angle_rad+2*pi;
    end

    diff_angl = abs(angls - angle_rad);
    loc_min = find(diff_angl == min(diff_angl));

    % Find order of directions to plot - if the 'peak' direction is in
    % subplot (position) 8. 
    ord = nan(16, 1);
    for j = 1:16
        id = (j - loc_min)+8;
        if id<1
            id = id+16;
        elseif id>16
            id = id-16;
        end 
        ord(j, 1) = id;
    end 

    target_direction = pi/2;  % Target direction (π/2 radians)
    current_peak_direction = angls(loc_min);
    rotation_offset = target_direction - current_peak_direction;
    
    aligned_directions = mod(angls + rotation_offset, 2*pi);
    % Combine the aligned angles and the responses
    d = [aligned_directions, responses];
    % Sort by angle so that arrays can be combined across cells. 
    d = sortrows(d, 1);

end 