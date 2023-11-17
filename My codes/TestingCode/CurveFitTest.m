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
%IMAGE CENTER IS 446, 892/2=446! MATH IS HARD SOMETIMES
scale=3.004*10^-4; %scaling factor
imagedata=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
%imagedata=load('C:\Users\Ides\OneDrive\Dokument\MATLAB\ProjectCourse\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
centerIm=892/2;

image=imagedata;
X=(1:892)*scale;
Y=X;
V=image;
Rmax=(892/2)*scale;
Rmin=(250)*scale;
R=linspace(Rmin,Rmax,1000);
theta=-90;
Xq=R*cosd(theta)+(892/2)*scale;
Yq=R*sind(theta)+(892/2)*scale;
intervall=((892/2)+220:892)';
%Xq=((892/2)+220:892)';
%Yq=ones(length(Xq),1)*462;
Vq = interp2(X,Y,V,Xq,Yq);
% Vqmax=max(Vq);
% VqNorm=Vq/imax;


figure;
hold on
plot(R,Vq,'--','LineWidth',2)
%plot(Vqlocs-((892/2-1)+220),Vqpeaks,'*')
xpos=(892/2)+220:892;
xpos1=(xpos-xpos(1)+220)*scale;
plot(xpos1,image(446,intervall))
legend('polyfit','raw')
hold off
%% find those peaks of the interp data
%normalize
VqMax=max(Vq);
VqNorm=Vq/VqMax;
[Vqpeaks Vqlocs]=findpeaks(VqNorm,R,'MinPeakHeight',0.25);

figure;
hold on
plot(R,VqNorm,'--','LineWidth',2)
plot(Vqlocs,Vqpeaks,'*')
hold off
% very good it works
figure;
plot(Xq,Yq)
xlim([1*scale 892*scale])
ylim([1*scale 892*scale])

%% Turn peak location into x,y points

XPeakLoc=Vqlocs*cosd(theta);%+(892/2)*scale;
YPeakLoc=Vqlocs*sind(theta); %+(892/2)*scale;

figure;
plot(XPeakLoc,YPeakLoc,'*')
xlim([-(892/2)*scale (892/2)*scale])
ylim([-(892/2)*scale (892/2)*scale])

% ok this seems to work, lets make a loop

%% loop time find the circle peaks

scale=3.004*10^-4; %scaling factor
imagedata=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
%load('C:\Users\Ides\OneDrive\Dokument\MATLAB\ProjectCourse\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image; %desktop

%initialize center and image data
centerIm=892/2;
image=imagedata;

%x and y grid in image
X=(-centerIm+1:centerIm)*scale;
Y=[centerIm:-1:-centerIm+1]*scale;
Y=Y+0.0105;
V=image;

%define a radius from the center, cut anything below 220 pixels (can be
%changed
Rmax=(892/2)*scale;
Rmin=(250)*scale;

% define the interpolation line
R=linspace(Rmin,Rmax,1000);

%angle of line
theta=0;

%angle step
dtheta=90;

% x and y peak location vectors, initialize with zeros as one does not know
% a priori how many peaks the program will extract, this is a very basic
%solution
XPeakVec=zeros(360/dtheta,5);%[];
YPeakVec=zeros(360/dtheta,5);%[];

%index for inserting peak data into vector
i=1;

%feel free to change dtheta as one sees fit
for theta=0:dtheta:360
    % x and y positions to be interpolated as function of R and theta
    Xq=R*cosd(theta);%+(892/2)*scale;
    Yq=R*sind(theta);%(892/2)*scale;
   
   % interpolated values
    Vq = interp2(X,Y,V,Xq,Yq);

    %Normalize
    VqMax=max(Vq);
    VqNorm=Vq/VqMax;
    %Find Peaks
    [Vqpeaks Vqlocs]=findpeaks(VqNorm,R,'MinPeakHeight',0.15);
    % Extract x and y positions from the peak radial locations
    XPeakLoc=Vqlocs*cosd(theta);%+(892/2)*scale;
    YPeakLoc=Vqlocs*sind(theta); %+(892/2)*scale;
    % add values to x and y peak vectors
    XPeakVec(i,1:length(XPeakLoc))=XPeakLoc;
    YPeakVec(i,1:length(YPeakLoc))=YPeakLoc;
    % index
    i=i+1;

    % troubleshooting polyfit
%     figure;
%     hold on
%     plot(Vqlocs,Vqpeaks,'*')
%     plot(R,VqNorm,'LineWidth',2)
%     hold off
end
%% plot results
figure;
plot(XPeakVec,YPeakVec,'*')
xlim([-(892/2)*scale (892/2)*scale])
ylim([-(892/2)*scale (892/2)*scale])
%%
figure(100),clf
hold on
contourf(X,Y,V,255,'linecolor','none')
grid on
plot(XPeakVec,YPeakVec,'*')
plot(0,0,'+','MarkerSize',10)
