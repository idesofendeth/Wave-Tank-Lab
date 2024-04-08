function [XPeakVec,YPeakVec,Centers,Radius] = CrestFinderV4(image,scale,X,Y,Rmini,thetaStep,InitialSwitch)
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
%% looping through the image in a circle, using dtheta as a angle step. This should probably be its own function.
% perhaps remove the loop from the function? could be cleaner.

% input: 

% output: 
for theta=0:dtheta:360-dtheta %270-dtheta*2
    % x and y positions to be interpolated as function of R and theta
    Xq=R*cosd(theta);%+(892/2)*scale;
    Yq=R*sind(theta);%(892/2)*scale;

    % interpolated values
    Vq = interp2(X,Y,V,Xq,Yq);

    %Normalize
    VqMax=max(Vq);
    VqNorm=Vq/VqMax;
    VqNormMean=movmean(Vq/VqMax,15); %moving mean
    %Find Peaks threshold usually .15, maybe need to play with this value
    [Vqpeaks, Vqlocs, Vqwidths, Vqproms ]=findpeaks(VqNormMean,R,'MinPeakHeight',0.085,'annotate','extents')
    %while index
    k=2;

    %% this should be its own function? VqSort? 
    % input: line 85, output from findpeaks
    % output: fixed version of locations
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
        Vqproms=VqpromsFixed;

        k=k+1;
    
    %end of while loop
    end
    PolyPeakLocs=[];
    for ii=1:length(VqlocsFixed)
        [fitresult, gof] = createPoly2FitV3(R, VqNorm,VqlocsFixed(ii),VqwidthsFixed(ii)/2);

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
    figure(98),clf
    hold on
    a=length(R(1:end))
    plot(VqlocsFixed,VqPeaksFixed,'*')
    plot(R,VqNorm,'LineWidth',2)
    %plot(R,VqNorm,'g--','LineWidth',2)
    hold off
    %pause

%end of for loop
end
%output goes here


%% circle testing circlefit 1991 this should also be its own function: circlefitter?
% input: XPeakVec and YPeakVec from previous
% function. These are then cleaned up and
% non-physical points are removed
% output: x and y center positions with radii for each reported circle
xCentVec=[];
yCentVec=[];
radiusVec=[];
for l=1:size(XPeakVec,2)
    A=[XPeakVec(:,l) YPeakVec(:,l)];
    indA=find(A(:,1)>=-50); % values that are not found are set to -100, this index finds anything larger, could be set to -99, doesnt matter
    B=A(indA,:);
    radiusVectemp=sqrt(B(:,1).^2+B(:,2).^2);
    
    outliers=radiusVectemp>=1.02*mean(radiusVectemp) | radiusVectemp<=0.98*mean(radiusVectemp); % setting outliers to be +-2 procent of data
    B=B(~outliers,:); %remove outliers
    figure(101),clf
    hold on
    grid on
    %plot(B(:,1),B(:,2),'*','MarkerSize',10)
    %[xCenter, yCenter, radiusFit, a] = circlefit(XPeakVec(:,l), YPeakVec(:,l))
    [xCenter, yCenter, radiusFit, a] = circlefit(B(:,1), B(:,2));

    xCentVec=[xCentVec; xCenter];
    yCentVec=[yCentVec; yCenter];
    if l<size(XPeakVec,2)
        if radiusFit<=centerIm*scale
        radiusVec=[radiusVec; radiusFit];
        end
    end
    clear A B radiusVectemp
    %viscircles([0 0],radiusFit,'LineWidth',1)
    %pause
end

%output goes here

meanYCent=mean(rmoutliers(yCentVec));
meanXCent=mean(rmoutliers(xCentVec));

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
circcenters=zeros(length(radiusVec),2);
viscircles(circcenters,radiusVec,'LineWidth',1)
Radius=radiusVec;
Centers=[meanXCent meanYCent];


%end of function
end 