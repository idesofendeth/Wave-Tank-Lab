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
thres = 23; %for marble with depth 1cm
%thres = 25; %for marble with depth 6cm


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
image=I2{127};
image=im2double(image);
%image=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
YOff=0.0105;
XOff=0;
Rmini=250;
thetaStep=45/2;
centerIm=size(image,1)/2;
X=(-centerIm+1:centerIm)*scale;
Y=[centerIm:-1:-centerIm+1]*scale;
Y=Y+YOff; %offset found by experimentation, this is a decent center of image
X=X-XOff;

XCenterDelta=XOff;
YCenterDelta=YOff;




%% initial center finding RUN ONCE
[XPeakVec,YPeakVec,Centers,Radius] = CrestFinderV5(image,scale,X,Y,Rmini,thetaStep);

%% real center RUN ONCE
XOff=Centers(1,1);
YOff=Centers(1,2);
XCenterDelta=XCenterDelta-XOff;
YCenterDelta=YCenterDelta-YOff;
Y=Y-YOff; %offset found by experimentation, this is a decent center of image
X=X-XOff;
[XPeakVec2,YPeakVec2,Centers2,Radius2] = CrestFinderV5(image,scale,X,Y,Rmini,thetaStep);


%% real center
XOff=Centers2(1,1);
YOff=Centers2(1,2);
XCenterDelta=XCenterDelta-XOff;
YCenterDelta=YCenterDelta-YOff;
Y=Y-YOff; %offset found by experimentation, this is a decent center of image
X=X-XOff;
[XPeakVec2,YPeakVec2,Centers2,Radius2] = CrestFinderV5(image,scale,X,Y,Rmini,thetaStep);

%% loop as much as you want center offsets xoff=0.0026 y=0.0129 marble 6cm

scale=3.004*10^-4; %scaling factor
image=I2{125};
image=im2double(image);
%image=load('C:\Users\ideso\OneDrive\Dokument\MATLAB\Project course\ProjectCourse\My codes\TestingCode\imageForDerek.mat').image;
YOff=0.0129;
XOff=0.0026;
Rmini=250;
thetaStep=45/5;
centerIm=size(image,1)/2; 
X=(-centerIm+1:centerIm)*scale;
Y=[centerIm:-1:-centerIm+1]*scale;
Y=Y+YOff; %offset found by experimentation, this is a decent center of image
X=X+XOff;
%%
%image=I2{125};
%image=im2double(image);
[XPeakVec2,YPeakVec2,Centers2,Radius2] = CrestFinderV4(image,scale,X,Y,Rmini,thetaStep);
%% edgecase testing

%Vq = interp2(X,Y,V,Xq,Yq);
%figure(95),clf
%contourf(X,Y,image,256,'LineColor','none')

echo off;
thetaStep=45/5;
for imdex=137
    image=I2{imdex};
    image=im2double(image);
    [XPeakVec,YPeakVec,Centers,Radius,Radiuskmeans,RadiusDerek] = CrestFinderV5(image,scale,X,Y,Rmini,thetaStep)
    %contourf(X,Y,image,256,'LineColor','none')
    V=image;
    figure(95), clf
    imagesc(X,-Y,V./max(V)) %,255,'linecolor','none')
    pbaspect([1 1 1])

colorbar
colormap('bone')
    Radius 
    Radiuskmeans
    RadiusDerek
    %pause
 
end



%% lambda testing


 %% Dispersion relation plot
 istart=100;
 iend=210;
    fig=figure(150);
    g = 9.82;
    %lambda = linspace(0,0.15);
    lambda = linspace(0,0.35,10000);
    rho = 997;
    sigma = 0.07275;
    H = linspace(1,1,length(lambda));

    mainaxes=axes(fig);
    hold(mainaxes,"on") 
    %Dispersion relation
    c = sqrt( ( g* lambda /(2* pi) + 2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
    disp('The limiting cases')
    %pure capillary wave
    c_cap = sqrt( (  2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
    %Pure gravity wave
    c_grav = sqrt( ( g* lambda /(2* pi) ) .*tanh( 2*pi* H./lambda ) );
   
    % plot(lambda,c,lambda,c_cap,'-.',lambda,c_grav,'-.','LineWidth',1)

    plot(lambda*100,c*100,lambda*100,c_cap*100,'-.',lambda*100,c_grav*100,'-.','LineWidth',1.5)
for i=istart:iend-1
    V=image;
    %Vq = interp2(X,Y,V,Xq,Yq);
    %Rmini=150;
    Rmini=200;
    echo off;
    thetaStep=45/5;
     imdex=i
        image=I2{imdex};
        image=im2double(image);
        [XPeakVec,YPeakVec,Centers,Radius1,Radiuskmeans1,RadiusDerek1] = CrestFinderV5(image,scale,X,Y,Rmini,thetaStep);

        %Radius
        %pause

    
    V=image;
    figure(95), clf
    imagesc(X,flip(Y),V./max(V)) %,255,'linecolor','none')
    colorbar
    colormap('bone')
    title('image #'+string(imdex))

     imdex=i+1
        image=I2{imdex};
        image=im2double(image);
        [XPeakVec,YPeakVec,Centers,Radius2,Radiuskmeans2,RadiusDerek2] = CrestFinderV5(image,scale,X,Y,Rmini,thetaStep);

        %Radius
        %pause

    
    V=image;

%     figure(96), clf
%     imagesc(X,flip(Y),V./max(V)) %,255,'linecolor','none')
%     colorbar
%     colormap('bone')
    %% phase speeds
    dt=1/300;
    lambvec1=diff(RadiusDerek1)*100;
    lambvec2=diff(RadiusDerek2)*100;
    %lambavg=(lambvec1+lambvec2)/2;
    cderek=phaseSpeedCalc(RadiusDerek1,RadiusDerek2,dt);
    cphasederek=cderek*100;
    cphasederek=cphasederek(cphasederek>0);
    ckmeans=phaseSpeedCalc(Radiuskmeans1,Radiuskmeans2,dt);
    cphasekmeans=ckmeans*100;
    cphasekmeans=cphasekmeans(cphasekmeans>0);

    hold on
    % making vectors equal for plotting for kmeans solution
    cphasekmeans=rmoutliers(cphasekmeans);
    desired_length=min(numel(cphasekmeans),numel(lambvec1));
    lambvec1=lambvec1(1:desired_length);
    cphasekmeans=cphasekmeans(1:desired_length);
    %plot(mainaxes,lambvec1,cphasekmeans,'x','LineWidth',1.5)
    % making vectors equal for plotting for dereks solution
     cphasederek=rmoutliers(cphasederek);
    desired_length=min(numel(cphasederek),numel(lambvec1));
    lambvec1=lambvec1(1:desired_length);
    cphasederek=cphasederek(1:desired_length);
    plot(mainaxes,lambvec1,cphasederek,'.','LineWidth',1.5,'MarkerSize',10)
    

   
    

   
    
    %plot(lamdavec,phasevec(:,2:end),'*') %lamda will always have 1 less column, thus one must cut some data from phase vec.
    %this might be losing some data, but it works for the most part
    



    ylabel(mainaxes,'$c$ [cm/s]','Interpreter','latex')
    xlabel(mainaxes,'$\lambda$ [cm]','Interpreter','latex')

    legend(mainaxes,'$c$','$c_{capillary}$','$c_{gravity}$','Data','Interpreter','latex')
    title(mainaxes,'Dispersion relation: Theory vs captured data')
    ylim(mainaxes,[0 100])
    xlim(mainaxes,[0 6])
    % disp('Shallow water approx, kH<<1, leads to tanh(kH)~kH')
    % disp('Deep water approx, kH>>1, leads to tanh(kH)~1')
    % %Deep water surface waves ( kH>>1 )
    % c_cap_deep = sqrt( (2* pi* sigma)./ (rho*lambda)); %Pure capillary wave  (dispersive)
    % c_grav_deep = sqrt( (g* lambda)./ (2*pi) ); %Pure gravity wave (dispersive)
    % %Shallow water surface waves ( kH<<1 )
    % c_cap_shallow = 2* pi./lambda .*sqrt( sigma.*H/ rho); %Pure capillary wave (dispersive)
    % c_grav_shallow = sqrt(g* H); % Pure gravity wave (non dispersive)
    lambda_min = 2*pi*sqrt(sigma/(rho*g))*10^2;
    lambda_cap = sqrt(sigma/(rho*g))*10^2;
    %2*pi*sqrt(lambda_)


    set(mainaxes,'fontsize',15)

end
hold off