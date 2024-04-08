clear all
close all
%% testing ideas many parts of this starting code are taking from Albins main code file.
I = load( sprintf('Images_Marble_%dcm.mat',6) ).I;
%I = load( sprintf('Images_MarbleADJUSTED_%dcm.mat',10) ).I;
%I = load( sprintf('Images_Droplet_%dcm.mat',6) ).I;


% Create an averaged image-------------------------------------------------
% NOTES: MAYBE MAKE A IF STATEMENT HERE SUCH THAT ONE CHOOSES AN IMAGE
% SERIES NUMBER (eg 6, 10 ect) and the avg and thres are chosen
% automatically

% furthermore, the following work quite well and give good results when
% matching the correct average and thresholding for the given dataset.
% Which is awesome!
%Iavg = AverageImageFunc(I(1:15)); %works for marble 10cm, this works in this code aswell
Iavg = AverageImageFunc(I(:));%works for marble depth 1cm, 6cm. Also works for 1cm droplet!
%Iavg = AverageImageFunc(I(185:200));%works for marble depth 1cm, 6cm
%Iavg=AverageImageFunc(I(45:end));
%Threshold level for binarizing--------------------------------------------
%thres = 18; %for droplet with depth 6cm
%thres = 23; %for marble with depth 1cm
thres = 25; %for marble with depth 6cm


%Remove background, adjust contrast, threshold, edge-detection-------------
I2 = cell(1,height(I));
for i = 1:height(I)
    %Remove background from images-----------------------------------------
    I2{i} = I{i}-Iavg;%imsubtract(I{i},Iavg) ;

    %Adjust image contrast-------------------------------------------------
    I3{i} = imadjust(I2{i});

    %Filter image----------------------------------------------------------
   % I3{i} = imdiffusefilt(I3{i});

    % % %Threshold image-------------------------------------------------------
    % I3{i} = I3{i} > thres;
    % 
    % %Clean up image--------------------------------------------------------
    % I3{i} = bwareaopen(I3{i},50) ;
    % I3{i} = imclearborder(I3{i});
    %edgesmoothing
    % windowSize = 51;
    % kernel = ones(windowSize) / windowSize ^ 2;
    % blurryImage = conv2(single(I3{i}), kernel, 'same');
    % I3{i} = blurryImage > thres; % Rethreshold

end
%%
scale=3.004*10^-4; %scaling factor
image=I2{125}%load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
image=im2double(image);%load('C:\Users\Ides\OneDrive\Dokument\MATLAB\ProjectCourse\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image; %desktop
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
    [XPeakVec,YPeakVec,Centers,Radius] = CrestFinderV3(image,scale,X,Y,Rmini,thetaStep);


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

for jj=1:numCircles-3
    % Combine x and y data into a single matrix for optimization
    allData = [xData(:,jj) yData(:,jj)];

    % Initial guesses for centers and radii
    initialGuesses = zeros(numCircles,3 ); % [centerX1, centerY1, radius1, centerX2, centerY2, radius2, ...]

    % Optimize using lsqnonlin
    options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
    estimatedParameters = lsqnonlin(@(params) circleResidual(params, allData, numPoints, numCircles), initialGuesses, [], [], options);

    % Extract estimated centers and radii
    estimatedCenters = reshape(estimatedParameters(1:2*numCircles), [], 2);
    estimatedRadii = estimatedParameters(2*numCircles+1:end);
    meanCenter=mean(estimatedCenters);
    meanRadii=mean(estimatedRadii);
    RadiiVec=[RadiiVec meanRadii];
    CenterVec=[CenterVec meanCenter];

    hold on
    scatter(xData(:, jj), yData(:, jj), 'filled');
    plot(estimatedCenters(jj, 1), estimatedCenters(jj, 2), 'o', 'MarkerSize', 10, 'LineWidth', 2);
    thetaFit = linspace(0, 2*pi, 100);
%     xFit = estimatedCenters(jj, 1) + estimatedRadii(jj) * cos(thetaFit);
%     yFit = estimatedCenters(jj, 2) + estimatedRadii(jj) * sin(thetaFit);
     xFit = CenterVec(1, 1) + estimatedRadii(jj) * cos(thetaFit);
    yFit = CenterVec(1, 2) + estimatedRadii(jj) * sin(thetaFit);
    plot(xFit, yFit, 'LineWidth', 2);
    pbaspect([1 1 1])
end
hold off

bestCenters=mean(estimatedCenters);
%%
StartIm=125;
EndIm=150;
intervall=abs(EndIm-StartIm)+1;

CellRadius = cell(intervall,2);
cellindex=1;
for imdex=StartIm:EndIm
    image=I2{imdex};
    image=im2double(image);
    [XPeakVec,YPeakVec,Centers,Radius] = CrestFinderV3(image,scale,X,Y,Rmini,thetaStep);
V=image;
CellRadius{cellindex,1}=Radius;
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
figure(150),clf;
imagesc(X,Y,V./max(V)) %,255,'linecolor','none')


for jj=1:numCircles-3
    % Combine x and y data into a single matrix for optimization
    allData = [xData(:,jj) yData(:,jj)];

    % Initial guesses for centers and radii
    initialGuesses = zeros(numCircles,3 ); % [centerX1, centerY1, radius1, centerX2, centerY2, radius2, ...]

    % Optimize using lsqnonlin
    options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);
    estimatedParameters = lsqnonlin(@(params) circleResidual(params, allData, numPoints, numCircles), initialGuesses, [], [], options);

    % Extract estimated centers and radii
    estimatedCenters = reshape(estimatedParameters(1:2*numCircles), [], 2);
    estimatedRadii = estimatedParameters(2*numCircles+1:end);
    meanCenter=mean(estimatedCenters);
    meanRadii=mean(estimatedRadii);
    RadiiVec=[RadiiVec meanRadii];
    CenterVec=[CenterVec meanCenter];

    hold on
    scatter(xData(:, jj), yData(:, jj), 'filled');
    plot(estimatedCenters(jj, 1), estimatedCenters(jj, 2), 'o', 'MarkerSize', 10, 'LineWidth', 2);
    thetaFit = linspace(0, 2*pi, 100);
%     xFit = estimatedCenters(jj, 1) + estimatedRadii(jj) * cos(thetaFit);
%     yFit = estimatedCenters(jj, 2) + estimatedRadii(jj) * sin(thetaFit);
     xFit = bestCenters(1, 1) + estimatedRadii(jj) * cos(thetaFit);
    yFit = bestCenters(1, 2) + estimatedRadii(jj) * sin(thetaFit);
    plot(xFit, yFit, 'LineWidth', 2);
    pbaspect([1 1 1])
end
hold off
CellRadius{cellindex,2}=RadiiVec;
cellindex=cellindex+1;
end
%%
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
