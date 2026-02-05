function optParams = fitGaussian(xData, yData, zData, idx)
% FITGAUSSIAN  Fit a 2D Gaussian to selected data points.
%
%   OPTPARAMS = FITGAUSSIAN(XDATA, YDATA, ZDATA, IDX) fits a non-rotated
%   2D Gaussian function to a subset of response data specified by idx.
%
%   INPUTS:
%     xData - 1D array of x coordinates (column positions)
%     yData - 1D array of y coordinates (row positions)
%     zData - 1D array of response values
%     idx   - Logical or numeric indices specifying which points to fit
%
%   OUTPUT:
%     optParams - 6-element array of optimized Gaussian parameters:
%                 [A, x0, y0, sigma_x, sigma_y, B]
%                 A       = amplitude
%                 (x0,y0) = center coordinates
%                 sigma_x = x standard deviation
%                 sigma_y = y standard deviation
%                 B       = baseline offset
%
%   MODEL:
%     z = A * exp(-((x-x0)^2/(2*sigma_x^2) + (y-y0)^2/(2*sigma_y^2))) + B
%
%   PREPROCESSING:
%     Response data is log-transformed: sign(z) * log(1 + |z|)
%     This reduces the influence of outliers and improves fit stability.
%
%   NOTE:
%     This is a simplified version without rotation. For rotated Gaussian
%     fitting, see GAUSSIAN_RF_ESTIMATE.
%
%   See also GAUSSIAN_RF_ESTIMATE, LSQCURVEFIT

    zData = sign(zData) .* log(1 + abs(zData)); % Log-transform responses

    % Initial parameter guess
    A_init = max(zData(idx)) - min(zData(idx));
    x0_init = mean(xData(idx));
    y0_init = mean(yData(idx));
    sigma_x_init = std(xData(idx));
    sigma_y_init = std(yData(idx));
    B_init = min(zData(idx)); % mean(zData(idx)); 
    initParams = [A_init, x0_init, y0_init, sigma_x_init, sigma_y_init, B_init];

    % Bounds (optional)
    lb = [0, min(xData(idx)), min(yData(idx)), 0, 0, min(zData(idx))]; % Lower bounds
    ub = [Inf, max(xData(idx)), max(yData(idx)), Inf, Inf, max(zData(idx))]; % Upper bounds

    % Define the 2D Gaussian function
    gauss2D = @(params, xy) params(1) * exp(-((xy(:,1) - params(2)).^2 / (2*params(4)^2) + ...
                                         (xy(:,2) - params(3)).^2 / (2*params(5)^2))) + params(6);

    % Optimize using least squares
    optParams = lsqcurvefit(gauss2D, initParams, [xData(idx), yData(idx)], zData(idx), lb, ub);
end

% 
% function optParams = fitGaussian(xData, yData, zData, idx)
%     % Select data points in the inhibitory region
%     xSel = xData(idx);
%     ySel = yData(idx);
%     zSel = zData(idx);
% 
%     % Improved initial parameter guess
%     A_init = abs(min(zSel)) - mean(zSel);  % Ensure A is positive
%     x0_init = sum(xSel .* abs(zSel)) / sum(abs(zSel));  % Weighted mean
%     y0_init = sum(ySel .* abs(zSel)) / sum(abs(zSel));  % Weighted mean
%     sigma_x_init = sqrt(sum((xSel - x0_init).^2 .* abs(zSel)) / sum(abs(zSel))); 
%     sigma_y_init = sqrt(sum((ySel - y0_init).^2 .* abs(zSel)) / sum(abs(zSel))); 
%     B_init = mean(zSel);
% 
%     initParams = [A_init, x0_init, y0_init, sigma_x_init, sigma_y_init, B_init];
% 
%     % Bounds (adjust for stability)
%     lb = [0, min(xSel), min(ySel), 0, 0, min(zSel)];
%     ub = [Inf, max(xSel), max(ySel), Inf, Inf, max(zSel)];
% 
%     gauss2D = @(params, xy) params(1) * exp(-((xy(:,1) - params(2)).^2 / (2*params(4)^2) + ...
%                                          (xy(:,2) - params(3)).^2 / (2*params(5)^2))) + params(6);
% 
%     % Optimize using least squares
%     optParams = lsqcurvefit(gauss2D, initParams, [xSel, ySel], zSel, lb, ub);
% end
