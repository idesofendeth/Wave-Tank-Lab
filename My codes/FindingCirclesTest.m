clc
clear all
close all

%% start
%scaling factor
f=3.004*10^-4;
%%
A=load('ImageDataAlbin\Images_Droplet_1cm.mat');

I = load( sprintf('Images_Marble_%dcm.mat',6) ).I; 
%I = load( sprintf('Images_MarbleADJUSTED_%dcm.mat',NR) ).I; 
%I = load( sprintf('Images_Droplet_%dcm.mat',NR) ).I; 
%%

% Create an averaged image-------------------------------------------------
%Iavg = AverageImageFunc(I(1:15)); %works for marble 10cm,
Iavg = AverageImageFunc(I(:));%works for marble depth 1cm, 6cm
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
    I2{i} = imsubtract(I{i},Iavg) ;

    % %Adjust image contrast-------------------------------------------------
     I3{i} = imadjust(I2{i});
    % 
    % %Filter image----------------------------------------------------------
     I3{i} = imdiffusefilt(I3{i});
    % 
    % %Threshold image-------------------------------------------------------
      I3{i} = I3{i} > thres;
    % % 
    % % %Clean up image--------------------------------------------------------
    % I3{i} = bwareaopen(I3{i},50) ;
    I3{i} = imclearborder(I3{i});
    %edgesmoothing
    % windowSize = 10;
    % kernel = ones(windowSize) / windowSize ^ 2;
    % blurryImage = conv2(single(I3{i}), kernel, 'same');
    % I3{i} = blurryImage >0.5; % Rethreshold
end
%I3=I2;
%% stuff
image1=cell2mat(I3(150));


figure;
image1=I3{100};
imagehalf=image1(:,446:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);
figure;
imshow(I{100},[])

% figure;
% imshow(mirror,[])
%%
radiusRange=[50 250]
[centers1, radii1, metric1] = imfindcircles(image1,radiusRange,'ObjectPolarity','bright','EdgeThreshold',0.4)%'Sensitivity',0.95,'Method','twostage')%"Method","twostage")
figure;
imshow(I{100},[])

figure;
imshow(image1,[])
centersStrong1 = centers1(1:end,:); 
radiiStrong1 = radii1(1:end);
metricStrong1 = metric1(1:end);
viscircles(centersStrong1, radiiStrong1,'EdgeColor','b');

%% next frame
figure;
image2=I3{110};
imagehalf=image2(:,446:end);
imageflip=flipdim(imagehalf,2);
mirror2=horzcat(imageflip,imagehalf);
figure;
imshow(I{110},[])
figure;
imshow(image2,[])
% figure;
% imshow(mirror,[])


[centers2, radii2, metric2] = imfindcircles(image2,radiusRange,'ObjectPolarity','bright','EdgeThreshold',0.4)%'Sensitivity',0.95,'Method','twostage')%,'Sensitivity',0.95,'EdgeThreshold',0.5,"Method","twostage")

centersStrong2 = centers2(1:end,:); 
radiiStrong2 = radii2(1:end);
metricStrong2 = metric2(1:end);
viscircles(centersStrong2, radiiStrong2,'EdgeColor','b');
%% loop frame 1
Cvec1=[];
Rvec1=[];
Mvec1=[];
for k=1:25
radiusRange=[50+k*10 200+k*10];
[centers1, radii1, metric1] = imfindcircles(image1,radiusRange,'ObjectPolarity','bright','EdgeThreshold',0.4,'Sensitivity',0.87);%,'Method','twostage')%"Method","twostage")

Cvec1=[Cvec1; centers1];
Rvec1=[Rvec1; radii1];
Mvec1=[Mvec1; metric1];
end

%% loop frame 2
Cvec2=[];
Rvec2=[];
Mvec2=[];
for k=1:25
radiusRange=[50+k*10 200+k*10];
[centers2, radii2, metric2] = imfindcircles(image2,radiusRange,'ObjectPolarity','bright','EdgeThreshold',0.4,'Sensitivity',0.87);%,'Method','twostage')%"Method","twostage")

Cvec2=[Cvec2; centers2];
Rvec2=[Rvec2; radii2];
Mvec2=[Mvec2; metric2];
end


%%
figure;
imshow(image1,[])


viscircles(Cvec1, Rvec1,'EdgeColor','b');
 


figure;
imshow(image2,[])

viscircles(Cvec2, Rvec2,'EdgeColor','b');



%%

r1avg=mean(Rvec1);
r2avg=mean(Rvec2);
dr=abs(radii2-radii1)*f
dravg=abs(r1avg-r2avg)*f
dt=1/300;
Dt=10*dt;

v=dr/Dt

vavg=dravg/Dt %ms

g=9.81;
lambtest=(284.78-268.85)*f
k=2*pi/lambtest
c=sqrt(g/k)*10^2 %cm/s
lambcm=lambtest*10^2

lambtest2=(315-300)*f

k2=2*pi/lambtest2
c2=sqrt(g/k2)*10^2 %cm/s
lambcm2=lambtest2*10^2