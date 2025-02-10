function [optEx, R_squared, optInh, R_squaredi, f1, f2] = gaussian_RF_estimate(response, min_data)
%% Fit and plot a rotated 2D Gaussian around excitatory and inhibitory lobes

% Define the x and y coordinate grid
[xGrid, yGrid] = meshgrid(1:size(response,2), 1:size(response,1));
xData = xGrid(:);
yData = yGrid(:);
zData = response(:); 
zData = sign(zData) .* log(1 + abs(zData));

% Define the rotated 2D Gaussian function
gauss2D = @(params, xy) params(1) * exp(-((( (xy(:,1) - params(2)) * cos(params(6)) + (xy(:,2) - params(3)) * sin(params(6)) ).^2 / (2*params(4)^2)) + ...
                                         (( (- (xy(:,1) - params(2)) * sin(params(6)) + (xy(:,2) - params(3)) * cos(params(6)) ).^2 / (2*params(5)^2))))) + params(7);

% Initial parameter guess: [A, x0, y0, sigma_x, sigma_y, theta, B]
A_init = max(zData) - min(zData);
x0_init = mean(xData);
y0_init = mean(yData);
sigma_x_init = std(xData);
sigma_y_init = std(yData);
theta_init = 0; % Assume no initial rotation
B_init = min(zData);
initParams = [A_init, x0_init, y0_init, sigma_x_init, sigma_y_init, theta_init, B_init];

% Bounds for parameters
lb = [0, min(xData), min(yData), 0, 0, -pi, min(zData)]; % Lower bounds
ub = [Inf, max(xData), max(yData), Inf, Inf, pi, max(zData)]; % Upper bounds

% Optimize using least squares for excitatory lobe
optEx = lsqcurvefit(gauss2D, initParams, [xData, yData], zData, lb, ub);

% Compute fitted values
zFit = gauss2D(optEx, [xData, yData]);

% Compute R-squared
SS_res = sum((zData - zFit).^2);
SS_tot = sum((zData - mean(zData)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Reshape fitted values into grid form
% zFitGrid = reshape(zFit, size(response));

% Plot original and fitted response
% figure;
% subplot(1,3,1);
% imagesc(response);
% med_val = median(response(:));
% max_val = prctile(response(:), 98);
% clim([med_val-max_val med_val+max_val])
% title('Original Response Data');
% axis image; colorbar; colormap redblue;

%% Fit inhibitory region

% Prepare inhibitory response data
zData2 = (min_data(:)) * -1; 

% Fit Gaussian to inhibitory region
optInh = lsqcurvefit(gauss2D, initParams, [xData, yData], zData2, lb, ub);

% Compute fitted values for inhibition
zFiti = gauss2D(optInh, [xData, yData]);

% Compute R-squared
SS_resi = sum((zData2 - zFiti).^2);
SS_toti = sum((zData2 - mean(zData2)).^2);
R_squaredi = 1 - (SS_resi / SS_toti);

% Create meshgrid for contour plotting
[Xq, Yq] = meshgrid(1:size(response,2), 1:size(response,1));
excFitValues = gauss2D(optEx, [Xq(:), Yq(:)]);
inhFitValues = gauss2D(optInh, [Xq(:), Yq(:)]);

% Reshape into grid form
excFitGrid = reshape(excFitValues, size(response));
inhFitGrid = reshape(inhFitValues, size(response));

subplot(1,3,2);
imagesc(excFitGrid);
title('Fitted Excitatory Gaussian');
axis image; colorbar; colormap redblue;
med_val = median(excFitGrid(:));
max_val = prctile(excFitGrid(:), 98);
clim([med_val-max_val med_val+max_val])

subplot(1,3,3);
imagesc(inhFitGrid*-1);
title('Fitted Inhibitory Gaussian');
axis image; colorbar; colormap redblue;
med_val = median(inhFitGrid(:));
max_val = prctile(inhFitGrid(:), 98);
clim([med_val-max_val med_val+max_val])
f1 = gcf; 
f1.Position = [58, 795, 646, 203];

%% Compute and plot 2σ contours with rotation
theta_vals = linspace(0, 2*pi, 100); % Define angles for contour
R_exc = [cos(optEx(6)), -sin(optEx(6)); sin(optEx(6)), cos(optEx(6))]; % Rotation matrix for excitatory lobe
R_inh = [cos(optInh(6)), -sin(optInh(6)); sin(optInh(6)), cos(optInh(6))]; % Rotation matrix for inhibitory lobe

% Compute rotated contour points
exc_ellipse = R_exc * [1.5 * optEx(4) * cos(theta_vals); 1.5 * optEx(5) * sin(theta_vals)];
exc_x2sig = optEx(2) + exc_ellipse(1, :);
exc_y2sig = optEx(3) + exc_ellipse(2, :);

inh_ellipse = R_inh * [1.5 * optInh(4) * cos(theta_vals); 1.5 * optInh(5) * sin(theta_vals)];
inh_x2sig = optInh(2) + inh_ellipse(1, :);
inh_y2sig = optInh(3) + inh_ellipse(2, :);

% Plot results
figure;
imagesc(response); hold on;
colormap redblue; colorbar;
title('Receptive Field with 1.5σ Contours');

med_val = median(response(:));
max_val = prctile(response(:), 98);
clim([med_val-max_val med_val+max_val])

% Plot excitatory contour
plot(exc_x2sig, exc_y2sig, 'r', 'LineWidth', 2);

% Plot inhibitory contour
plot(inh_x2sig, inh_y2sig, 'k', 'LineWidth', 2);
axis square;

f2 = gcf;
f2.Position = [712, 576, 560, 420];

% Display the results from both the excitatory and inhibitory lobes
disp('Excitatory lobe:')
fprintf('Optimized Parameters: A = %.2f, x0 = %.2f, y0 = %.2f, sigma_x = %.2f, sigma_y = %.2f, theta = %.2f, B = %.2f\n', optEx);
fprintf('R-squared: %.4f\n', R_squared);

disp('Inhibitory lobe:')
fprintf('Optimized Parameters: A = %.2f, x0 = %.2f, y0 = %.2f, sigma_x = %.2f, sigma_y = %.2f, theta = %.2f, B = %.2f\n', optInh);
fprintf('R-squared: %.4f\n', R_squaredi);

end 
