%% Crest Finder script for project course
clc
clear all
close all

%% notes
% albin stored each video as a series of images
% in a .mat file. Each mat file is a struct, which contains cells
% each position in the cell relates to a frame in the video
% IMPORTANT NOTE: THE CM for droplet and marble refer to the DEPTH! Albin
% had a constant marble size of 1.6cm and unknown droplet size
% 3cm depth droplet seems to have issues, will check later

%% constants
scale=3.004*10^-4; %scaling factor
dt=1/300; %framerate time delta
centers=[446.497245615476 486.383816964470]; % For MARBLE FOUND FROM PREVIOUS EXPERIMENTATION using findCircles function
%centers=[406.0999  501.9339]; %for DROPLET 1cm
%centers=[435.9794  485.6637]; % for DROPLET 6cm
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

%% curve fit test
image=I2{125};
image=im2double(image);

x1=image(467,686:end)';
imax=max(x1);
x1=x1/imax;
x2=x1*imax;
L=length(image(467,686:end));
xgrid=0:1:L-1;
xgrid=xgrid'*scale;
xlower=0.005;
xupper=0.01;
[fitresult, gof] = createPoly2FitV2(xgrid,x1,peakPos,deltax)

coeffvals= coeffvalues(fitresult);
k=coeffvals(2);
lambda=(2*pi)/k; %meters

%ok cool, we now have a lambda for the first delta r. Lets make a loop
res=length(unique(x1))

resavg=length(unique(Iavg))
%% loop time
xlower=0.005;
xupper=0.0075;
deltar=abs(xlower-xupper);
lambdavec=[];
rgrid=xlower:deltar:length(xgrid)*scale;
tol=10^-2;
lamdex=[];
for i=0:length(rgrid)
    
    [fitresult, gof] = createFitSine(xgrid,x1,xlower+i*deltar,xupper+i*deltar);
    g=struct2cell(gof);
    if g{1}<tol
    coeffvals= coeffvalues(fitresult);
    k=coeffvals(2);
    lambda=(2*pi)/k; %meters
    lambdavec=[lambdavec; lambda];
    lamdex=[lamdex; i*deltar];
    end



    %end loop
end

%ok this works sort of, needs some conditions probably to ignore some
%points.
%% plot
figure;
plot(lamdex*100,lambdavec*100) %(1:end-1)
ylabel('Lambda (cm)')
xlabel('r (cm)')

%% findpeaks test

[peaks, locs]=findpeaks(x1,xgrid,'MinPeakHeight',0.25);
xlower=0.0075;
xupper=0.01;
deltar=abs(xlower-xupper);
%ok this works

figure;
plot(locs*scale,peaks,'*')
xlim([0 160*scale])
ylim([0 1])

%% test of function poly2fitV2
scalelocs=locs;%*scale;
peakPos=scalelocs(5);
deltax=3*scale; % is there a way to systematically identify x number of points? instead of a window, yes integer*scale.
% scale is the delta x for each pixel. We get discrete points at each pixel
% see below
[fitresult, gof] = createPoly2FitV2(xgrid,x1,peakPos,deltax)

%% loop poly2fitV2
deltax=3*scale;
fitPeakLoc=[];
for j=1:length(locs)
    peakPos=scalelocs(j);
    [fitresult, gof] = createPoly2FitV2(xgrid,x1,peakPos,deltax)
    coeffvals= coeffvalues(fitresult);
    b=coeffvals(2);
    fitPeakLoc=[fitPeakLoc b];

end

lambFit=diff(fitPeakLoc);


%% plot test poly2fitV2

% this looks promising

figure;
plot(fitPeakLoc(2:end)*100,lambFit*100,'*') %(1:end-1)
ylabel('Lambda (cm)')
xlabel('r (cm)')

%% interp2 test DEREK THIS IS FOR YOU
scale=3.004*10^-4; %scaling factor
imagedata=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
image=imagedata;
X=(1:892)*scale;
Y=X;
V=image;
Rmax=(892/2)*scale;
Rmin=(220)*scale;
R=linspace(Rmin,Rmax,1000);
theta=0;
Xq=R*cosd(theta)+(892/2)*scale;
Yq=R*sind(theta)+(892/2)*scale;
intervall=((892/2)+220:892)';
%Xq=((892/2)+220:892)';
%Yq=ones(length(Xq),1)*462;
Vq = interp2(X,Y,V,Xq,Yq);
% Vqmax=max(Vq);
% VqNorm=Vq/imax;
% [Vqpeaks Vqlocs]=findpeaks(VqNorm,Xq,'MinPeakHeight',0.25);

figure;
hold on
plot(R,Vq)
%plot(Vqlocs-((892/2-1)+220),Vqpeaks,'*')
xpos=(892/2)+220:892;
xpos1=(xpos-xpos(1)+220)*scale;
plot(xpos1,image(467,intervall))
legend('polyfit','raw')
hold off









%% image line rotate test chatgpt

% Define image size
imageSize = 892; % Adjust the size as needed

% Create a blank image matrix
imageMatrix = image;

% Define center of the image
center = imageSize / 2;

% Define the length of the line
lineLength = imageSize / 2;

% Define angles for rotation (in degrees)
angles = 0:1:360;

% Initialize arrays to store x and y values
xValues = zeros(size(angles));
yValues = zeros(size(angles));

% Create a figure to display the rotating line
figure;

for i = 1:length(angles)
    % Convert angle to radians
    angleRad = deg2rad(angles(i));
    
    % Calculate x and y coordinates of the line
    x = center + lineLength * cos(angleRad);
    y = center + lineLength * sin(angleRad);
    
    % Store x and y values in arrays
    xValues(i) = x;
    yValues(i) = y;
    
    % Plot the image matrix
    imagesc(imageMatrix);
    colormap(gray);
    hold on;
    
    % Plot the rotating line
    plot([center, x], [center, y], 'r', 'LineWidth', 2);
    
    % Set axis limits
    axis([1, imageSize, 1, imageSize]);
    
    % Set aspect ratio to be equal
    axis equal;
    
    % Pause to visualize the rotation
    pause(0.01);
    
    % Clear the plot for the next iteration
    clf;
end

% Display the x and y values on the line as it rotates
disp('Angle   X Value   Y Value');
disp([angles' xValues' yValues']);
