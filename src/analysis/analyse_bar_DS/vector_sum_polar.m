function [resultant_magnitude, resultant_angle] = vector_sum_polar(rho, theta)
    
    % Ensure input vectors are the same size. 
    if size(rho) ~= size(theta)
        error(strcat("Error: rho and theta must be the same size." + ...
            " rho = [", string(size(rho, 1)),",", string(size(rho, 2)), "]" + ...
            " and theta = [", string(size(theta, 1)), ",", string(size(theta, 2)),"]" ...
            ));
    end

    % Convert polar coordinates to Cartesian coordinates
    % [x, y] = pol2cart(theta, rho);
    x = sum(rho .* cos(theta));
    y = sum(rho .* sin(theta));

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
