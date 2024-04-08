% Example: Replace these with your actual data
xData = randn(16, 8); % 16 points for 8 circles
yData = randn(16, 8); % 16 points for 8 circles

% Combine x and y data into a single matrix for optimization
allData = [xData(:); yData(:)];

% Number of circles
numCircles = size(xData, 2);

% Initial guesses for centers and radii
initialGuesses = zeros(1, 3 * numCircles); % [centerX1, centerY1, radius1, centerX2, centerY2, radius2, ...]

% Optimize using lsqnonlin
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
estimatedParameters = lsqnonlin(@(params) circleResidual(params, allData, numCircles), initialGuesses, [], [], options);

% Extract estimated centers and radii
estimatedCenters = reshape(estimatedParameters(1:2*numCircles), [], 2);
estimatedRadii = estimatedParameters(2*numCircles+1:end);

% Plot the data and fitted circles
figure;
for i = 1:numCircles
    
    scatter(xData(:, i), yData(:, i), 'filled');
    hold on;
    
    % Plot fitted circle
    thetaFit = linspace(0, 2*pi, 100);
    xFit = estimatedCenters(i, 1) + estimatedRadii(i) * cos(thetaFit);
    yFit = estimatedCenters(i, 2) + estimatedRadii(i) * sin(thetaFit);
    plot(xFit, yFit, 'LineWidth', 2);
end

axis equal;
title('Fitted Circles');

% Residual function for lsqnonlin
function residual = circleResidual(params, data, numCircles)
    residual = [];
    for i = 1:numCircles
        centerX = params(2*i - 1);
        centerY = params(2*i);
        radius = params(numCircles*2 + i);
        
        xData = data(1:end/2);
        yData = data(end/2 + 1:end);
        
        residual = [residual; sqrt((xData - centerX).^2 + (yData - centerY).^2) - radius];
    end
end
