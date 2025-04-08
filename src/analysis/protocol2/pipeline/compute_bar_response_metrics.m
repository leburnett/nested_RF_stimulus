function [sym_ratio, DSI, DSI_pdnd, vector_sum] = compute_bar_response_metrics(d)

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