function [optEx, R_squared, optInh, R_squaredi, f1, f2] = gaussian_RF_estimate(response, min_data)
%% Plot gaussian fit of excitatory and inhibitory lobes. 

% Response data:
% response = data_comb;

% Define the x and y coordinate grid
[xGrid, yGrid] = meshgrid(1:size(response,2), 1:size(response,1));
xData = xGrid(:);
yData = yGrid(:);
zData = response(:); % Flatten response matrix into a vector
zData = sign(zData) .* log(1 + abs(zData));

% Define the 2D Gaussian function
gauss2D = @(params, xy) params(1) * exp(-((xy(:,1) - params(2)).^2 / (2*params(4)^2) + ...
                                         (xy(:,2) - params(3)).^2 / (2*params(5)^2))) + params(6);

% Initial parameter guess: [A, x0, y0, sigma_x, sigma_y, B]
A_init = max(zData) - min(zData);
x0_init = mean(xData);
y0_init = mean(yData);
sigma_x_init = std(xData);
sigma_y_init =  std(yData);
B_init = min(zData);
initParams = [A_init, x0_init, y0_init, sigma_x_init, sigma_y_init, B_init];

% Bounds (optional, can be adjusted)
lb = [0, min(xData), min(yData), 0, 0, min(zData)]; % Lower bounds
ub = [Inf, max(xData), max(yData), Inf, Inf, max(zData)]; % Upper bounds

% Optimize using least squares
optEx = lsqcurvefit(gauss2D, initParams, [xData, yData], zData, lb, ub);

% Compute fitted values
zFit = gauss2D(optEx, [xData, yData]);

% Compute R-squared
SS_res = sum((zData - zFit).^2);
SS_tot = sum((zData - mean(zData)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Display results
% fprintf('Optimized Parameters: A = %.2f, x0 = %.2f, y0 = %.2f, sigma_x = %.2f, sigma_y = %.2f, B = %.2f\n', optParams);
% fprintf('R-squared: %.4f\n', R_squared);

% Reshape fitted values into grid form
zFitGrid = reshape(zFit, size(response));

% Plot original and fitted response
figure;
subplot(1,3,1);
imagesc(response);
med_val = median(response(:));
max_val = prctile(response(:), 98);
clim([med_val-max_val med_val+max_val])
title('Original Response Data');
axis image; colorbar; colormap redblue;


%% Gaussian fit - plot 2SD contour around data. 

% Define the x and y coordinate grid
[xGrid, yGrid] = meshgrid(1:size(response,2), 1:size(response,1));
xData = xGrid(:);
yData = yGrid(:);
zData = response(:); %a(:); %response(:); % Flatten response matrix into a vector
zData2 = (min_data(:))*-1; %b(:)*-1; %min_data(:)*-1; 

% Fit separate Gaussians to excitatory and inhibitory regions
optExc = fitGaussian(xData, yData, zData, 1:numel(xData));
optInh = fitGaussian(xData, yData, zData2, 1:numel(xData));

% Compute fitted values for inhibition
zFiti = gauss2D(optInh, [xData, yData]);

% Compute R-squared
SS_resi = sum((zData2 - zFiti).^2);
SS_toti = sum((zData2 - mean(zData2)).^2);
R_squaredi = 1 - (SS_resi / SS_toti);

% Define the 2D Gaussian function
gauss2D = @(params, xy) params(1) * exp(-((xy(:,1) - params(2)).^2 / (2*params(4)^2) + ...
                                         (xy(:,2) - params(3)).^2 / (2*params(5)^2))) + params(6);

% Create meshgrid for contour plotting
[Xq, Yq] = meshgrid(1:size(response,2), 1:size(response,1));
excFitValues = gauss2D(optExc, [Xq(:), Yq(:)]);
inhFitValues = gauss2D(optInh, [Xq(:), Yq(:)]);

% % Reshape into grid form
excFitGrid = reshape(excFitValues, size(response));
inhFitGrid = reshape(inhFitValues, size(response));

subplot(1,3,2);
imagesc(excFitGrid);
title('Fitted 2D Gaussian');
axis image; colorbar; colormap redblue;
med_val = median(excFitGrid(:));
max_val = prctile(excFitGrid(:), 98);
clim([med_val-max_val med_val+max_val])

subplot(1,3,3);
imagesc(inhFitGrid*-1);
title('Fitted 2D Gaussian');
axis image; colorbar; colormap redblue;
med_val = median(inhFitGrid(:));
max_val = prctile(inhFitGrid(:), 98);
clim([med_val-max_val med_val+max_val])
f1 = gcf; f1.Position = [58   795   646   203];

% Compute 2-sigma contours
theta = linspace(0, 2*pi, 100); % Circle for contour
exc_x2sig = optExc(2) + 1.5*optExc(4) * cos(theta);
exc_y2sig = optExc(3) + 1.5*optExc(5) * sin(theta);
inh_x2sig = optInh(2) + 1.5*optInh(4) * cos(theta);
inh_y2sig = optInh(3) + 1.5*optInh(5) * sin(theta);

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
axis square

f2 = gcf;
f2.Position =[712   576   560   420];

% Display the results from both the excitatory and inhibitory lobes:
disp('Excitatory lobe:')
fprintf('Optimized Parameters: A = %.2f, x0 = %.2f, y0 = %.2f, sigma_x = %.2f, sigma_y = %.2f, B = %.2f\n', optEx);
fprintf('R-squared: %.4f\n', R_squared);

disp('Inhibitory lobe:')
fprintf('Optimized Parameters: A = %.2f, x0 = %.2f, y0 = %.2f, sigma_x = %.2f, sigma_y = %.2f, B = %.2f\n', optInh);
fprintf('R-squared: %.4f\n', R_squaredi);

end 
