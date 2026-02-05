function fwhm = compute_FWHM(angls, responses)
% COMPUTE_FWHM  Calculate full-width at half-maximum of tuning curve.
%
%   FWHM = COMPUTE_FWHM(ANGLS, RESPONSES) computes the angular width
%   of the directional tuning curve at half its maximum value.
%
%   INPUTS:
%     angls     - 16x1 array of stimulus directions in radians
%     responses - 16x1 array of neural response values
%
%   OUTPUT:
%     fwhm - Full-width at half-maximum in degrees
%            Smaller values indicate sharper directional tuning
%
%   ALGORITHM:
%     1. Find maximum response and calculate half-maximum threshold
%     2. Identify all angles where response exceeds half-maximum
%     3. FWHM = difference between first and last crossing angles
%     4. Adjust for circular angle wrap-around if needed
%
%   INTERPRETATION:
%     FWHM of 90 degrees indicates the neuron responds strongly to
%     a 90-degree range of directions (sharp tuning).
%     FWHM of 180 degrees indicates broad directional tuning.
%
%   See also FIND_PD_AND_ORDER_IDX, COMPUTE_CIRCULAR_VAR

    % Find maximum response and half-maximum
    R_max = max(responses);
    half_max = R_max / 2;

    % Find indices where responses cross the half-maximum threshold
    crossIdx = find(responses >= half_max);

    % Get the corresponding angles at those indices
    ang1 = angls(crossIdx(1)); % First crossing point
    ang2 = angls(crossIdx(end)); % Last crossing point

    % Compute FWHM as the difference
    fwhm = rad2deg(ang2 - ang1);

    % Ensure correct handling of circular angles (0 to 360 degrees)
    if fwhm < 0
        fwhm = fwhm + 360; % Adjust for circular range
    end

end