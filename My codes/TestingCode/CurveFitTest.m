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
[fitresult, gof] = createFitSine(xgrid,x1,xlower,xupper);

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

[peaks, locs]=findpeaks(x1,'MinPeakHeight',0.25);
xlower=0.0075;
xupper=0.01;
deltar=abs(xlower-xupper);
%ok this works



