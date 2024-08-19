function [XPeakVec,YPeakVec,Centers,Radius,Radiuskmeans,RadiusDerek] = CrestFinderV5(image,scale,X,Y,Rmini,thetaStep)
%CrestFinderV4 COPY: 
% INPUT:
% image = image data in greyscale
% scale = scaling factor from pixels to meters
% center = center of image IMPORTANT, we are assuming a square image
% XOff = X offset for the real or approximate center of the wave circles
% YOff = Y Offset for the real or approximate center of the wave circles
% Rmini = Minimum radius to consider, use this to cut away center of image
% dtheta = angle step to step around the image with
% OUTPUT:
% XPeakVec = X location of reported peaks of the crests of the waves
% YPeakVec = Y location of reported peaks of the cressts of the waves
% Centers = reported centers of the wave formations calculated from the
% crest data
% Radius = reported radius of the wave formations calculated from the crest
% data

% What needs to be done: Break this function up into a series of
% subfunctions such that troubleshooting can be done easily
%% initialization
%scale=3.004*10^-4; %scaling factor
%imagedata=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
%load('C:\Users\Ides\OneDrive\Dokument\MATLAB\ProjectCourse\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image; %desktop
%YOff=-0.0105;
%XOff=0;
%initialize center and image data
centerIm=size(image,1)/2;
%image=imagedata;

%x and y grid in image


V=image;

%define a radius from the center, cut anything below 220 pixels (can be
%changed
Rmax=(centerIm)*scale; % edge of image is exactly equal to the center position after one shifts the center to zero
Rmin=(Rmini)*scale; % minimum radius to be examined

% define the interpolation line
R=linspace(Rmin,Rmax,1000);

%angle of line
theta=0;

%angle step
dtheta=thetaStep;

% x and y peak location vectors, initialize with zeros as one does not know
% a priori how many peaks the program will extract, this is a very basic
%solution
XPeakVec=ones(360/dtheta,10)*(-100);%[]; multiplied by -100 as it is a value that will not exist in the solution
%thus, it can be removed later to get rid of non physical points. This is
%done because the vector/matrix sizes must be the same in the caluclations
%that follow for this code. Perhaps it can be done in a smarter way.
YPeakVec=ones(360/dtheta,10)*(-100);%[];

%index for inserting peak data into vector
i=1;
%space threshold for filtering peaks too close to one another
spaceThresh=0.0025;
% initialize variable for filtered peak locations and values
VqPeaksFixed=[];
VqlocsFixed=[];
VqwidthsFixed=[];
VqpromsFixed=[];
%feel free to change dtheta as one sees fit

%contourf(V)
%% looping through the image in a circle, using dtheta as a angle step. 
for theta=0:dtheta:360-dtheta %270-dtheta*2
    % x and y positions to be interpolated as function of R and theta
    Xq=R*cosd(theta);%+(892/2)*scale;
    Yq=R*sind(theta);%(892/2)*scale;

    % interpolated values from the data.
    Vq = interp2(X,Y,V,Xq,Yq);
   

    %Normalize
    VqMax=max(Vq);
    VqNorm=Vq/VqMax;
   
    VqNormMean=movmean(Vq/VqMax,10); %moving mean 15 is default
    %troubleshooting figure
%     figure(97),clf
%     plot(R,VqNorm,'-x')
%     hold on
%     plot(R,VqNormMean)
    %Find Peaks threshold usually .15, maybe need to play with this value
    [Vqpeaks, Vqlocs, Vqwidths, Vqproms ]=findpeaks(VqNormMean,R,'MinPeakHeight',0.085,'annotate','extents','MinPeakWidth',0.0001);
    %findpeaks(VqNormMean,R,'MinPeakHeight',0.085,'annotate','extents')
    [XPeakLoc,YPeakLoc] = VqSort(Vqpeaks, Vqlocs, Vqwidths, Vqproms,VqNorm,theta,spaceThresh,R);

    XPeakVec(i,1:length(XPeakLoc))=XPeakLoc;
    YPeakVec(i,1:length(YPeakLoc))=YPeakLoc;
    % index
    i=i+1;

%end of for loop
end

%% SOMEWHERE HERE THERE NEEDS TO BE DATA SORTING DONE IE KMEANS or OTHER methods.


%% fitting circles to data
[Centers,Radius,Radiuskmeans, RadiusDerek] = circleFitterV2(XPeakVec,YPeakVec,scale,centerIm);

%% visualization
figure(100),clf
hold on
imagesc(X,Y,V./max(V)) %,255,'linecolor','none')
colorbar
colormap('bone')
grid on
plot(XPeakVec,YPeakVec,'*','MarkerSize',10)
plot(0,0,'+','MarkerSize',10)
pbaspect([1 1 1])
xlim([-centerIm*scale centerIm*scale])
ylim([-centerIm*scale centerIm*scale])
circcenters=zeros(length(Radius),2);
viscircles(circcenters,Radius,'LineWidth',1)
title('Previous Implementation')


figure(1000),clf
hold on
imagesc(X,Y,V./max(V)) %,255,'linecolor','none')
colorbar
colormap('bone')
grid on
plot(XPeakVec,YPeakVec,'*','MarkerSize',10,'LineWidth', 2)
plot(0,0,'+','MarkerSize',10,'LineWidth', 2)
pbaspect([1 1 1])
xlim([-centerIm*scale centerIm*scale])
ylim([-centerIm*scale centerIm*scale])
circcenters=zeros(length(Radiuskmeans),2);
viscircles(circcenters,Radiuskmeans,'LineWidth',1)
title('Kmeans Solution')


figure(1001),clf
hold on
imagesc(X,Y,V./max(V)) %,255,'linecolor','none')
colorbar
colormap('bone')
grid on
plot(XPeakVec,YPeakVec,'*','MarkerSize',10,'LineWidth', 2)
plot(0,0,'+','MarkerSize',10,'LineWidth', 2)
pbaspect([1 1 1])
xlim([-centerIm*scale centerIm*scale])
ylim([-centerIm*scale centerIm*scale])
circcenters=zeros(length(RadiusDerek),2);
viscircles(circcenters,RadiusDerek,'LineWidth',1)
title('Derek Threshold Solution')


%end of function
end 