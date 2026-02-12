function circular_variance = compute_circular_var(angls, responses)
% COMPUTE_CIRCULAR_VAR  Calculate circular variance of directional tuning.
%
%   CIRCULAR_VARIANCE = COMPUTE_CIRCULAR_VAR(ANGLS, RESPONSES)
%   computes the circular variance metric indicating how broadly
%   the neural response is distributed across directions.
%
%   INPUTS:
%     angls     - 16x1 array of stimulus directions in radians
%     responses - 16x1 array of neural response values
%
%   OUTPUT:
%     circular_variance - Value between 0 and 1
%                         0 = all responses in one direction (sharp tuning)
%                         1 = responses uniformly distributed (no selectivity)
%
%   FORMULA:
%     CV = 1 - (|vector_sum| / sum(responses))
%     where vector_sum = sum(response * exp(i*angle))
%
%   INTERPRETATION:
%     Circular variance is the complement of the direction selectivity
%     index (DSI). Low CV indicates strong directional preference.
%     High CV indicates the neuron responds similarly to all directions.
%
%   See also FIND_PD_AND_ORDER_IDX, COMPUTE_FWHM, COMPUTE_BAR_RESPONSE_METRICS

vector_sum_x = sum(responses .* cos(angls));
vector_sum_y = sum(responses .* sin(angls));
magnitude = sqrt(vector_sum_x^2 + vector_sum_y^2);
circular_variance = 1 - (magnitude / sum(responses));

end 