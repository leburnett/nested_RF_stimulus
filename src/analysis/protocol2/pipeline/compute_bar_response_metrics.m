function [sym_ratio, DSI, DSI_pdnd, vector_sum] = compute_bar_response_metrics(d)
% COMPUTE_BAR_RESPONSE_METRICS  Calculate direction selectivity metrics.
%
%   [SYM_RATIO, DSI, DSI_PDND, VECTOR_SUM] = COMPUTE_BAR_RESPONSE_METRICS(D)
%   computes direction selectivity indices and symmetry metrics from
%   responses to 16 bar directions.
%
%   INPUT:
%     d - 16x2 matrix where:
%         d(:,1) = bar angles in radians (0 to 2*pi, in 22.5 degree steps)
%         d(:,2) = max voltage response for each angle
%         Data should be aligned so row 5 = 90 degrees (1.5708 rad)
%
%   OUTPUTS:
%     sym_ratio  - Symmetry ratio (0-1), higher = more symmetric tuning
%                  Compares responses on opposite sides of preferred direction
%     DSI        - Direction Selectivity Index via vector sum method
%                  |vector_sum| / sum(all responses)
%                  Range 0-1, higher = more selective
%     DSI_pdnd   - Direction Selectivity Index via PD/ND method
%                  (PD - ND) / (PD + ND)
%                  PD = preferred direction (row 5), ND = null (row 13)
%     vector_sum - Complex vector sum of all direction responses
%                  Magnitude indicates selectivity, angle indicates PD
%
%   SYMMETRY CALCULATION:
%     Compares responses at paired angles equidistant from PD:
%       Pairs: (4,6), (3,7), (2,8), (1,9), (16,10), (15,11), (14,12)
%     sym_ratio = 1 - (sum of absolute differences / sum of all responses)
%
%   VECTOR SUM:
%     vector_sum = sum(response * exp(i * angle))
%     The resulting complex number's angle gives the preferred direction
%     and magnitude indicates tuning strength
%
%   See also FIND_PD_AND_ORDER_IDX, PROCESS_BARS_P2

    %% Symmetry index:
    % 'd' - column 1 = angle in radians.
    % column 2 = max_v for that angle. Adjusted so that the peak is 90 deg.
    % (1.5708 rad). 
    % Row 5 = the row corresponding to 90 deg. 
    % 4-6
    % 3-7
    % 2-8
    % 1-9
    % 16-10
    % 15-11
    % 14-12
    
    vals1 = d([4,3,2,1,16,15,14], 2); % Increases in voltage from median per angle. 
    vals2 = d([6,7,8,9,10,11,12], 2);
    diff_vals = abs(vals1-vals2); % Find the diff voltage between the peak responses from the 2 halves.
    
    % Symmetry ratio:
    sym_val = sum(diff_vals)/sum(d(:, 2)); % Sum these differences and normalise by the sum of all of the responses.
    sym_ratio = 1 - sym_val; % Closer to 1 = more symmetric.
    
    % Compute the vector sum
    vector_sum = sum(d(:,2) .* exp(1i * d(:,1)));
    % Compute the DSI
    DSI = abs(vector_sum) / sum(d(:, 2));
    % DSI: - PD - ND
    DSI_pdnd = (d(5,2)-d(13,2))/(d(5,2)+d(13,2));

end 