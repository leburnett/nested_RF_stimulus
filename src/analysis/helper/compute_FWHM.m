function fwhm = compute_FWHM(angls, responses)
    % Ensure angles are sorted in ascending order
    % [angls, sortIdx] = sort(angls);
    % responses = responses(sortIdx);

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