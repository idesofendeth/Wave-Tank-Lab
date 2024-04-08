clear all
close all

%%
scale=3.004*10^-4; %scaling factor
image=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
%load('C:\Users\Ides\OneDrive\Dokument\MATLAB\ProjectCourse\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image; %desktop
YOff=0.0105;
XOff=0;
Rmini=250;
thetaStep=45/2;
centerIm=size(image,1)/2;
X=(-centerIm+1:centerIm)*scale;
Y=[centerIm:-1:-centerIm+1]*scale;
Y=Y+YOff; %offset found by experimentation, this is a decent center of image
X=X-XOff;
V=image;
%%
[XPeakVec,YPeakVec,Centers,Radius] = CrestFinderV3(image,scale,X,Y,Rmini,thetaStep);
%%

numPoints = size(XPeakVec,1);
numCircles = size(XPeakVec,2);

xData = XPeakVec;
yData = YPeakVec;

% for i = 1:numCircles
%     theta = linspace(0, 2*pi, numPoints);
%     noise = 0.1 * randn(1, numPoints);
%     radius = randi([5, 10]); % Random radius for each circle
%     xData(:, i) = radius * cos(theta) + noise;
%     yData(:, i) = radius * sin(theta) + noise;
% end
RadiiVec=[];
CenterVec=[];
figure;
fig=imagesc(X,Y,V./max(V)) %,255,'linecolor','none')

for jj=1:numCircles
    % Combine x and y data into a single matrix for optimization
    allData = [xData(:,jj) yData(:,jj)];

    % Initial guesses for centers and radii
    initialGuesses = zeros(numCircles,3 ); % [centerX1, centerY1, radius1, centerX2, centerY2, radius2, ...]

    % Optimize using lsqnonlin
    options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
    estimatedParameters = lsqnonlin(@(params) circleResidual(params, allData, numPoints, numCircles), initialGuesses, [], [], options);

    % Extract estimated centers and radii
    estimatedCentersX = reshape(estimatedParameters(1:2*numCircles), [], 2);
    estimatedRadii = estimatedParameters(2*numCircles+1:end);
    meanCenter=mean(estimatedCenters);
    meanRadii=mean(estimatedRadii);
    RadiiVec=[RadiiVec meanRadii];
    CenterVec=[CenterVec meanCenter];

    hold on
    scatter(xData(:, i), yData(:, i), 'filled');
    plot(estimatedCenters(i, 1), estimatedCenters(i, 2), 'o', 'MarkerSize', 10, 'LineWidth', 2);
    thetaFit = linspace(0, 2*pi, 100);
    xFit = estimatedCenters(i, 1) + estimatedRadii(i) * cos(thetaFit);
    yFit = estimatedCenters(i, 2) + estimatedRadii(i) * sin(thetaFit);
    plot(xFit, yFit, 'LineWidth', 2);
    pbaspect([1 1 1])
end
hold off
%%
% Plot the data and fitted circles
figure;
imagesc(X,Y,V./max(V)) %,255,'linecolor','none')
hold on
% Residual function for lsqnonlin
function residual = circleResidual(params, data, numPoints, numCircles)
residual = [];
for i = 1:numCircles
    centerX = params(2*i - 1);
    centerY = params(2*i);
    radius = params(numCircles*2 + i);

    xData = data(1:numPoints,1);
    yData = data(1:numPoints,2);

    residual = [residual; sqrt((xData - centerX).^2 + (yData - centerY).^2) - radius];
end
end
