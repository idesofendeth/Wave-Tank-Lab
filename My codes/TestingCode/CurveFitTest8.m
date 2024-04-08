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
imagedata=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
%load('C:\Users\Ides\OneDrive\Dokument\MATLAB\ProjectCourse\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image; %desktop
YOff=-0.0105;
XOff=0;
%initialize center and image data
centerIm=size(imagedata,1)/2;
image=I2{125};%imagedata;
image=im2double(image);

%x and y grid in image
X=(-centerIm+1:centerIm)*scale;
Y=[centerIm:-1:-centerIm+1]*scale;
Y=Y-YOff; %offset found by experimentation, this is a decent center of image
X=X-XOff;
V=image;

%define a radius from the center, cut anything below 220 pixels (can be
%changed
Rmax=(centerIm)*scale;
Rmin=(250)*scale;

% define the interpolation line
R=linspace(Rmin,Rmax,1000);

%angle of line
theta=0;

%angle step
dtheta=45/2;

% x and y peak location vectors, initialize with zeros as one does not know
% a priori how many peaks the program will extract, this is a very basic
%solution
XPeakVec=zeros(360/dtheta,5);%[];
YPeakVec=zeros(360/dtheta,5);%[];

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
for theta=0:dtheta:360-dtheta
    % x and y positions to be interpolated as function of R and theta
    Xq=R*cosd(theta);%+(892/2)*scale;
    Yq=R*sind(theta);%(892/2)*scale;

    % interpolated values
    Vq = interp2(X,Y,V,Xq,Yq);

    %Normalize
    VqMax=max(Vq);
    VqNorm=Vq/VqMax;
    %Find Peaks
    [Vqpeaks, Vqlocs, Vqwidths, Vqproms ]=findpeaks(VqNorm,R,'MinPeakHeight',0.15);
    %while index
    k=2;
    while k<=length(Vqlocs)
        %diff of location positions
        Vqlocdiff=Vqlocs(k)-Vqlocs(k-1);
        % windows of location positions and peak positions one is comparing
        VqLocWindow=[Vqlocs(k-1) Vqlocs(k)];
        VqPeakWindow=[Vqpeaks(k-1) Vqpeaks(k)];
        VqwidthsWindow=[Vqwidths(k-1) Vqwidths(k)];
        VqpromsWindow=[Vqproms(k-1) Vqproms(k)];
        % filtering out locations of peaks that are too close to eachother
        % taking the largest peak to save, throwing away the other
        % there is some bug here, some of the points arent disappearing..
        if abs(Vqlocdiff)<spaceThresh
            %finding the value and index of the maximum between the
            %investigated points
            [VqPeakFilter I]=max(VqPeakWindow);
            %edge cases for start and end of vector
            if k==2
                VqlocsFixed=[VqLocWindow(I) Vqlocs(k+1:end)];
                VqPeaksFixed=[VqPeakFilter Vqpeaks(k+1:end)];
                VqwidthsFixed=[VqwidthsWindow(I) Vqwidths(k+1:end)];
                VqpromsFixed=[VqpromsWindow(I) Vqproms(k+1:end)];
            elseif k==length(Vqlocs)
                VqlocsFixed=[Vqlocs(1:k-2) VqLocWindow(I)];
                VqPeaksFixed=[Vqpeaks(1:k-2) VqPeakFilter];
                VqwidthsFixed=[Vqwidths(1:k-2) VqwidthsWindow(I)];
                VqpromsFixed=[Vqproms(1:k-2) VqpromsWindow(I)];
            else %otherwise do this
                VqlocsFixed=[Vqlocs(1:k-2) VqLocWindow(I) Vqlocs(k+1:end)];
                VqPeaksFixed=[Vqpeaks(1:k-2) VqPeakFilter Vqpeaks(k+1:end)];
                VqwidthsFixed=[Vqwidths(1:k-2) VqwidthsWindow(I) Vqwidths(k+1:end)];
                VqpromsFixed=[Vqproms(1:k-2) VqpromsWindow(I) Vqproms(k+1:end)];

            end
        else
            VqlocsFixed=Vqlocs;
            VqPeaksFixed=Vqpeaks;
            VqwidthsFixed=Vqwidths;
            VqpromsFixed=Vqproms;


        end
        Vqlocs=VqlocsFixed;
        Vqpeaks=VqPeaksFixed;
        Vqwidths=VqwidthsFixed;
        Vqproms=VqpromsFixed

        k=k+1;
    end
    PolyPeakLocs=[];
    for ii=1:length(VqlocsFixed)
        [fitresult, gof] = createPoly2FitV3(R, VqNorm,VqlocsFixed(ii),VqwidthsFixed(ii)/2)

        PolyPeakLocs=[PolyPeakLocs fitresult.b];

    end
    VqlocsFixed=PolyPeakLocs;

    % Extract x and y positions from the peak radial locations
    XPeakLoc=VqlocsFixed*cosd(theta);
    YPeakLoc=VqlocsFixed*sind(theta);
    % add values to x and y peak vectors
    XPeakVec(i,1:length(XPeakLoc))=XPeakLoc;
    YPeakVec(i,1:length(YPeakLoc))=YPeakLoc;
    % index
    i=i+1;

    % troubleshooting polyfit plot
%     figure(98),clf
%     hold on
%     plot(VqlocsFixed,VqPeaksFixed,'*')
%     plot(R,VqNorm,'LineWidth',2)
%     hold off
end
%% plot results
figure(99),clf
plot(XPeakVec,YPeakVec,'*')
xlim([-(892/2)*scale (892/2)*scale])
ylim([-(892/2)*scale (892/2)*scale])
pbaspect([1 1 1])
%%
figure(100),clf
hold on
contourf(X,Y,V,255,'linecolor','none')
grid on
plot(XPeakVec,YPeakVec,'*','MarkerSize',10)
plot(0,0,'+','MarkerSize',10)
pbaspect([1 1 1])

%% circle testing circlefit 1991
xCentVec=[];
yCentVec=[];
radiusVec=[];
for l=1:size(XPeakVec,2)
    [xCenter, yCenter, radiusFit, a] = circlefit(XPeakVec(:,l), YPeakVec(:,l))
    xCentVec=[xCentVec; xCenter];
    yCentVec=[yCentVec; yCenter];
    if l<size(XPeakVec,2)
        radiusVec=[radiusVec; radiusFit];
    end
end
meanYCent=mean(rmoutliers(yCentVec));
meanXCent=mean(rmoutliers(xCentVec));


%% attempt to center based on data


scale=3.004*10^-4; %scaling factor
imagedata=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
%load('C:\Users\Ides\OneDrive\Dokument\MATLAB\ProjectCourse\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image; %desktop

%initialize center and image data
%centerIm=892/2;
image=imagedata;

%x and y grid in image
%X=(-centerIm+1:centerIm)*scale;
%Y=[centerIm:-1:-centerIm+1]*scale;

Y=Y-meanYCent; %offset found by experimentation, this is a decent center of image
X=X-meanXCent;
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
dtheta=45/2;

% x and y peak location vectors, initialize with zeros as one does not know
% a priori how many peaks the program will extract, this is a very basic
%solution
XPeakVec=zeros(360/dtheta,5);%[];
YPeakVec=zeros(360/dtheta,5);%[];

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
for theta=0:dtheta:360-dtheta
    % x and y positions to be interpolated as function of R and theta
    Xq=R*cosd(theta);%+(892/2)*scale;
    Yq=R*sind(theta);%(892/2)*scale;

    % interpolated values
    Vq = interp2(X,Y,V,Xq,Yq);

    %Normalize
    VqMax=max(Vq);
    VqNorm=Vq/VqMax;
    %Find Peaks
    [Vqpeaks, Vqlocs, Vqwidths, Vqproms ]=findpeaks(VqNorm,R,'MinPeakHeight',0.15);
    %while index
    k=2;
    while k<=length(Vqlocs)
        %diff of location positions
        Vqlocdiff=Vqlocs(k)-Vqlocs(k-1);
        % windows of location positions and peak positions one is comparing
        VqLocWindow=[Vqlocs(k-1) Vqlocs(k)];
        VqPeakWindow=[Vqpeaks(k-1) Vqpeaks(k)];
        VqwidthsWindow=[Vqwidths(k-1) Vqwidths(k)];
        VqpromsWindow=[Vqproms(k-1) Vqproms(k)];
        % filtering out locations of peaks that are too close to eachother
        % taking the largest peak to save, throwing away the other
        % there is some bug here, some of the points arent disappearing..
        if abs(Vqlocdiff)<spaceThresh
            %finding the value and index of the maximum between the
            %investigated points
            [VqPeakFilter I]=max(VqPeakWindow);
            %edge cases for start and end of vector
            if k==2
                VqlocsFixed=[VqLocWindow(I) Vqlocs(k+1:end)];
                VqPeaksFixed=[VqPeakFilter Vqpeaks(k+1:end)];
                VqwidthsFixed=[VqwidthsWindow(I) Vqwidths(k+1:end)];
                VqpromsFixed=[VqpromsWindow(I) Vqproms(k+1:end)];
            elseif k==length(Vqlocs)
                VqlocsFixed=[Vqlocs(1:k-2) VqLocWindow(I)];
                VqPeaksFixed=[Vqpeaks(1:k-2) VqPeakFilter];
                VqwidthsFixed=[Vqwidths(1:k-2) VqwidthsWindow(I)];
                VqpromsFixed=[Vqproms(1:k-2) VqpromsWindow(I)];
            else %otherwise do this
                VqlocsFixed=[Vqlocs(1:k-2) VqLocWindow(I) Vqlocs(k+1:end)];
                VqPeaksFixed=[Vqpeaks(1:k-2) VqPeakFilter Vqpeaks(k+1:end)];
                VqwidthsFixed=[Vqwidths(1:k-2) VqwidthsWindow(I) Vqwidths(k+1:end)];
                VqpromsFixed=[Vqproms(1:k-2) VqpromsWindow(I) Vqproms(k+1:end)];

            end
        else
            VqlocsFixed=Vqlocs;
            VqPeaksFixed=Vqpeaks;
            VqwidthsFixed=Vqwidths;
            VqpromsFixed=Vqproms;


        end
        Vqlocs=VqlocsFixed;
        Vqpeaks=VqPeaksFixed;
        Vqwidths=VqwidthsFixed;
        Vqproms=VqpromsFixed

        k=k+1;
    end
    PolyPeakLocs=[];
    for ii=1:length(VqlocsFixed)
        [fitresult, gof] = createPoly2FitV3(R, VqNorm,VqlocsFixed(ii),VqwidthsFixed(ii)/2)

        PolyPeakLocs=[PolyPeakLocs fitresult.b];

    end
    VqlocsFixed=PolyPeakLocs;

    % Extract x and y positions from the peak radial locations
    XPeakLoc=VqlocsFixed*cosd(theta);
    YPeakLoc=VqlocsFixed*sind(theta);
    % add values to x and y peak vectors
    XPeakVec(i,1:length(XPeakLoc))=XPeakLoc;
    YPeakVec(i,1:length(YPeakLoc))=YPeakLoc;
    % index
    i=i+1;

    % troubleshooting polyfit plot
%     figure(98),clf
%     hold on
%     plot(VqlocsFixed,VqPeaksFixed,'*')
%     plot(R,VqNorm,'LineWidth',2)
%     hold off
end
%%
xCentVec=[];
yCentVec=[];
radiusVec=[];
for l=1:size(XPeakVec,2)
    zcheck=numel(find(~XPeakVec(:,l)));
    if zcheck<=2
    [xCenter, yCenter, radiusFit, a] = circlefit(XPeakVec(:,l), YPeakVec(:,l))
    xCentVec=[xCentVec; xCenter];
    yCentVec=[yCentVec; yCenter];
    radiusVec=[radiusVec; radiusFit];
    end
end
meanYCent=mean(yCentVec);
meanXCent=mean(xCentVec);
%%
figure(1000),clf
hold on
contourf(X,Y,V,255,'linecolor','none')
grid on
plot(XPeakVec,YPeakVec,'*','MarkerSize',10)
plot(0,0,'+','MarkerSize',10)
pbaspect([1 1 1])
circcenters=zeros(length(radiusVec),2);
% circcenters(:,1)=circcenters(:,1)-meanXCent;
% circcenters(:,2)=circcenters(:,2)-meanYCent;
viscircles(circcenters,radiusVec)

%% prominence test
figure(101),clf
findpeaks(VqNorm,R,'MinPeakHeight',0.15,'annotate','extents')
[pks,locs,widths,proms]=findpeaks(VqNorm,R,'MinPeakHeight',0.15,'annotate','extents')


%% fit using width from prominence using poly2fitv2

[fitresult, gof] = createPoly2FitV2(R, VqNorm,VqlocsFixed(1),widths(1)/2)


%% loop fit test
PolyPeakLocs=[];
for ii=1:length(VqlocsFixed)
    [fitresult, gof] = createPoly2FitV2(R, VqNorm,VqlocsFixed(ii),VqwidthsFixed(ii)/2)

    PolyPeakLocs=[PolyPeakLocs fitresult.b];

end


