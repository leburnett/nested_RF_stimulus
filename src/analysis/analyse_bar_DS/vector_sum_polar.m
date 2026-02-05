function [resultant_magnitude, resultant_angle] = vector_sum_polar(rho, theta)
% VECTOR_SUM_POLAR  Compute vector sum of polar coordinates.
%
%   [RESULTANT_MAGNITUDE, RESULTANT_ANGLE] = VECTOR_SUM_POLAR(RHO, THETA)
%   calculates the vector sum of polar coordinate data, typically used
%   to determine preferred direction from directional tuning responses.
%
%   INPUTS:
%     rho   - 1xN or Nx1 array of magnitudes (response values)
%     theta - 1xN or Nx1 array of angles in radians
%             Must be same size as rho
%
%   OUTPUTS:
%     resultant_magnitude - Magnitude of the vector sum
%     resultant_angle     - Direction of vector sum in radians [0, 2*pi]
%
%   ALGORITHM:
%     1. Converts polar coordinates to Cartesian (x = rho*cos(theta))
%     2. Sums x and y components separately
%     3. Converts back to polar (angle and magnitude)
%
%   CIRCULAR DATA HANDLING:
%     If the first and last theta values are identical (completing a
%     circle), the last value is excluded to avoid double-weighting
%     that direction.
%
%   USAGE:
%     The resultant angle indicates the preferred direction.
%     The resultant magnitude indicates the strength of directional
%     preference (larger = more selective).
%
%   See also FIND_PD_AND_ORDER_IDX, PLOT_POLAR_WITH_ARROW

    % Ensure input vectors are the same size. 
    if size(rho) ~= size(theta)
        error(strcat("Error: rho and theta must be the same size." + ...
            " rho = [", string(size(rho, 1)),",", string(size(rho, 2)), "]" + ...
            " and theta = [", string(size(theta, 1)), ",", string(size(theta, 2)),"]" ...
            ));
    end

    % Convert polar coordinates to Cartesian coordinates
    % [x, y] = pol2cart(theta, rho);

    % Exclude the last value if it's the same as the first, don't want to
    % double weight that direction.
    if theta(1) == theta(end)
        x = sum(rho(1:end-1) .* cos(theta(1:end-1)));
        y = sum(rho(1:end-1) .* sin(theta(1:end-1)));
    else
        x = sum(rho .* cos(theta));
        y = sum(rho .* sin(theta));
    end 

    % % Compute the vector sum in Cartesian coordinates
    % sum_x = sum(x);
    % sum_y = sum(y);
    % 
    % % Convert back to polar coordinates
    % [resultant_angle, resultant_magnitude] = cart2pol(sum_x, sum_y); % cart2pol returns angle first

    resultant_angle = atan2(y, x);
    resultant_angle = mod(resultant_angle, 2*pi); % ensure it's in [0, 2pi]
    resultant_magnitude = sqrt(x^2 + y^2);

end
