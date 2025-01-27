% Function to fit Gaussian to given indices
function optParams = fitGaussian(xData, yData, zData, idx)
    % Initial parameter guess
    A_init = max(zData(idx)) - min(zData(idx));
    x0_init = mean(xData(idx));
    y0_init = mean(yData(idx));
    sigma_x_init = std(xData(idx));
    sigma_y_init = std(yData(idx));
    B_init = min(zData(idx));
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
