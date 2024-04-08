%% loop time find the circle peaks
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

%%
[XPeakVec,YPeakVec,Centers,Radius] = CrestFinderV3(image,scale,X,Y,Rmini,thetaStep);


%%

% Number of circles and datasets
numCircles = size(XPeakVec,2);
numDatasets = size(XPeakVec,1);
allX = XPeakVec;
allY = YPeakVec;
% Initial guesses for the centers and radii
initialGuesses = zeros(numCircles, 3); % [centerX, centerY, radius]

% Optimize for each circle and dataset
estimatedParameters = zeros(numCircles, 3);
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-6);

for i = 1:numCircles
    %for j = 1:numDatasets
        % Optimize using lsqnonlin
        currentDataX = allX(:,i);
        currentDataY = allY(:,i);
        initialGuess = initialGuesses(i, :);
        estimatedParameters(i, :) = lsqnonlin(@(params) circleResidual(params, currentDataX, currentDataY), initialGuess,[],[],options);
    %end
end

% Extract estimated radii
estimatedRadii = estimatedParameters(:, 3);

% Plot the data and fitted circles
figure;
for i = 1:numDatasets
    for j = 1:numCircles
        scatter(allX(i, j), allY(i, j), 'filled');
        hold on;
    end
end

for i = 1:numCircles
    plot(estimatedParameters(i, 1), estimatedParameters(i, 2), 'o', 'MarkerSize', 10, 'LineWidth', 2);
    thetaFit = linspace(0, 2*pi, 100);
    xFit = estimatedParameters(i, 1) + estimatedRadii(i) * cos(thetaFit);
    yFit = estimatedParameters(i, 2) + estimatedRadii(i) * sin(thetaFit);
    plot(xFit, yFit, 'LineWidth', 2);
end

axis equal;
title('Fitted Circles');

function residual = circleResidual(center, x, y)
    % Residual function for least squares fitting of a circle
    radius = sqrt((x - center(1)).^2 + (y - center(2)).^2);
    residual = radius - mean(radius);
end
