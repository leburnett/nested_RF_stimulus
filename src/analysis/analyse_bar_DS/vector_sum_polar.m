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
