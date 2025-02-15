function circular_variance = compute_circular_var(angls, responses)
% If CV = 0, all responses are in one direction (sharp tuning).
% If CV = 1, responses are uniformly spread (broad tuning).

vector_sum_x = sum(responses .* cos(angls));
vector_sum_y = sum(responses .* sin(angls));
magnitude = sqrt(vector_sum_x^2 + vector_sum_y^2);
circular_variance = 1 - (magnitude / sum(responses));

end 