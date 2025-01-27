
%% Plot gaussian fit of excitatory and inhibitory lobes. 

% Response data:
response = data_comb2;

% Define the x and y coordinate grid
[xGrid, yGrid] = meshgrid(1:size(response,2), 1:size(response,1));
xData = xGrid(:);
yData = yGrid(:);
zData = response(:); % Flatten response matrix into a vector

% Define the 2D Gaussian function
gauss2D = @(params, xy) params(1) * exp(-((xy(:,1) - params(2)).^2 / (2*params(4)^2) + ...
                                         (xy(:,2) - params(3)).^2 / (2*params(5)^2))) + params(6);

% Initial parameter guess: [A, x0, y0, sigma_x, sigma_y, B]
A_init = max(zData) - min(zData);
x0_init = mean(xData);
y0_init = mean(yData);
sigma_x_init = std(xData);
sigma_y_init = std(yData);
B_init = min(zData);
initParams = [A_init, x0_init, y0_init, sigma_x_init, sigma_y_init, B_init];

% Bounds (optional, can be adjusted)
lb = [0, min(xData), min(yData), 0, 0, min(zData)]; % Lower bounds
ub = [Inf, max(xData), max(yData), Inf, Inf, max(zData)]; % Upper bounds

% Optimize using least squares
optParams = lsqcurvefit(gauss2D, initParams, [xData, yData], zData, lb, ub);

% Compute fitted values
zFit = gauss2D(optParams, [xData, yData]);

% Compute R-squared
SS_res = sum((zData - zFit).^2);
SS_tot = sum((zData - mean(zData)).^2);
R_squared = 1 - (SS_res / SS_tot);

% Display results
fprintf('Optimized Parameters: A = %.2f, x0 = %.2f, y0 = %.2f, sigma_x = %.2f, sigma_y = %.2f, B = %.2f\n', optParams);
fprintf('R-squared: %.4f\n', R_squared);

% Reshape fitted values into grid form
zFitGrid = reshape(zFit, size(response));

% Plot original and fitted response
figure;
subplot(1,2,1);
imagesc(response);
title('Original Response Data');
axis image; colorbar; colormap redblue;

subplot(1,2,2);
imagesc(zFitGrid);
title('Fitted 2D Gaussian');
axis image; colorbar; colormap redblue;



%% Gaussian fit - plot 2SD contour around data. 

% Define the x and y coordinate grid
[xGrid, yGrid] = meshgrid(1:size(response,2), 1:size(response,1));
xData = xGrid(:);
yData = yGrid(:);
zData = response(:); % Flatten response matrix into a vector

% Define a threshold to separate excitatory and inhibitory responses
threshold = 2 ; %median(zData);  
excIdx = zData > threshold; % Excitatory indices
inhIdx = zData <= threshold; % Inhibitory indices

% Fit separate Gaussians to excitatory and inhibitory regions
optExc = fitGaussian(xData, yData, zData, excIdx);
optInh = fitGaussian(xData, yData, zData, inhIdx);

% Define the 2D Gaussian function
gauss2D = @(params, xy) params(1) * exp(-((xy(:,1) - params(2)).^2 / (2*params(4)^2) + ...
                                         (xy(:,2) - params(3)).^2 / (2*params(5)^2))) + params(6);

% Create meshgrid for contour plotting
[Xq, Yq] = meshgrid(1:size(response,2), 1:size(response,1));
excFitValues = gauss2D(optExc, [Xq(:), Yq(:)]);
inhFitValues = gauss2D(optInh, [Xq(:), Yq(:)]);

% Reshape into grid form
excFitGrid = reshape(excFitValues, size(response));
inhFitGrid = reshape(inhFitValues, size(response));

% Compute 2-sigma contours
theta = linspace(0, 2*pi, 100); % Circle for contour
exc_x2sig = optExc(2) + 2*optExc(4) * cos(theta);
exc_y2sig = optExc(3) + 2*optExc(5) * sin(theta);
inh_x2sig = optInh(2) + 2*optInh(4) * cos(theta);
inh_y2sig = optInh(3) + 2*optInh(5) * sin(theta);

% Plot results
figure;
imagesc(response); hold on;
colormap redblue; colorbar;
title('Receptive Field with 2σ Contours');

% Plot excitatory contour
plot(exc_x2sig, exc_y2sig, 'r', 'LineWidth', 2);

% Plot inhibitory contour
plot(inh_x2sig, inh_y2sig, 'k', 'LineWidth', 2);


%% 
figure; 
subplot(2,1,1)
imagesc(max_data); colorbar; title('max - 98th prctile')
subplot(2,1,2)
histogram(max_data)
f=gcf;
f.Position = [620   501   275   466];

figure; 
subplot(2,1,1)
imagesc(min_data); colorbar; title('min - 98th prctile')
subplot(2,1,2)
histogram(min_data)
f=gcf;
f.Position = [620   501   275   466];