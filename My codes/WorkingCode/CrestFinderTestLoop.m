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
%centers=[446.497245615476 486.383816964470]; % For MARBLE FOUND FROM PREVIOUS EXPERIMENTATION using findCircles function
%centers=[406.0999  501.9339]; %for DROPLET 1cm
centers=[435.9794  485.6637]; % for DROPLET 6cm
%% testing ideas many parts of this starting code are taking from Albins main code file.
%I = load( sprintf('Images_Marble_%dcm.mat',1) ).I;
%I = load( sprintf('Images_MarbleADJUSTED_%dcm.mat',10) ).I;
I = load( sprintf('Images_Droplet_%dcm.mat',6) ).I;


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
thres = 18; %for droplet with depth 6cm
%thres = 23; %for marble with depth 1cm
%thres = 25; %for marble with depth 6cm


%Remove background, adjust contrast, threshold, edge-detection-------------
I2 = cell(1,height(I));
for i = 1:height(I)
    %Remove background from images-----------------------------------------
    I2{i} = imsubtract(I{i},Iavg) ;

    %Adjust image contrast-------------------------------------------------
    I3{i} = imadjust(I2{i});

    %Filter image----------------------------------------------------------
    I3{i} = imdiffusefilt(I3{i});

    % %Threshold image-------------------------------------------------------
    I3{i} = I3{i} > thres;

    %Clean up image--------------------------------------------------------
    I3{i} = bwareaopen(I3{i},50) ;
    I3{i} = imclearborder(I3{i});
    %edgesmoothing
    % windowSize = 51;
    % kernel = ones(windowSize) / windowSize ^ 2;
    % blurryImage = conv2(single(I3{i}), kernel, 'same');
    % I3{i} = blurryImage > thres; % Rethreshold

end

%% finding centers experimental stuff
%c=imfindcircles(I3{67},[150 500]);
%visCircles(centers,250)

%% loop
% indexes for start and end frames, change these as one sees fit

% for 10cm marble start 30 end 210
% for 6cm marble start 80 end 285
% for 1cm marble start 60 end 235
% for 6cm droplet start 20 end 200
% for 1cm droplet start 67 end 250
%halfdex is the vertical position of the line to mirror the image across
%the position of the circles is different for both marble and droplet

%halfdex=447; %for marble
%halfdex=402; %for droplet 1cm
halfdex=436; %for droplet 6cm
imstart=20;
imend=200;
% initializing variables
size=[(abs(imstart-imend)) length(I)];
phasevec=[];
lamdavec=[];
phasevec2=zeros(size);
lamdavec2=zeros(size);
j=1;
m=1;
figure;
%begin loop
for k=imstart:imend

    %frame 1
    image=I3;
    imdex1=k;
    image1=cell2mat(image(imdex1));
    imagehalf=image1(:,halfdex:end);
    imageflip=flipdim(imagehalf,2);
    mirror=horzcat(imageflip,imagehalf);
  
    minradius=150;

    [Crest1,lamda1, innerdistvec1,outerdistvec1] = CrestFinder(mirror,centers,minradius);
    drawnow


    % frame 2
    imdex1=k+1;
    image1=cell2mat(image(imdex1));
    imagehalf=image1(:,halfdex:end);
    imageflip=flipdim(imagehalf,2);
    mirror=horzcat(imageflip,imagehalf);
    
    minradius=150;

    [Crest2,lamda2, innerdistvec2,outerdistvec2] = CrestFinder(mirror,centers,minradius);



    %Checking that the number of crests found are equal in frame 1 and 2
    %for comparison. This is to prevent errors in the vector size when
    %calculating phase speed.
    %This is a simple solution, which could use some looking into as I
    %believe I may be throwing away some data which could be used for the
    %sake of making the code run
    
    if length(Crest2)==length(Crest1)
        phase=abs(Crest2-Crest1)/dt;

        phaseScaled=phase*scale;

        phaseCMs=phaseScaled*10^2;

        lamda1Scaled=lamda1*scale*10^2;
        lamda2Scaled=lamda2*scale*10^2;
        phasevec(j,1:length(phaseCMs))=phaseCMs;
        lamdavec(j,1:length(lamda2Scaled))=lamda2Scaled;
        j=j+1;

    else
        ldex=min([length(Crest1) length(Crest2)]);
        phase=abs(Crest2(1:ldex)-Crest1(1:ldex))/dt;

        phaseScaled=phase*scale;

        phaseCMs=phaseScaled*10^2;
        phasevec2(m,1:length(phaseCMs))=phaseCMs;

        lamda1Scaled=lamda1*scale*10^2;
        lamda2Scaled=lamda2*scale*10^2;
        lamdavec2(m,1:length(lamda2Scaled))=lamda2Scaled;
    end

end

%% Dispersion relation plot
g = 9.82;
lambda = linspace(0,0.15);
lambda = linspace(0,0.35,10000);
rho = 997;
sigma = 0.07275;
H = linspace(1,1,length(lambda));


%Dispersion relation
c = sqrt( ( g* lambda /(2* pi) + 2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
disp('The limiting cases')
%pure capillary wave
c_cap = sqrt( (  2*pi*sigma./ (rho *lambda) ) .*tanh( 2*pi* H./lambda ) );
%Pure gravity wave
c_grav = sqrt( ( g* lambda /(2* pi) ) .*tanh( 2*pi* H./lambda ) );
%Plots

figure(1);
% plot(lambda,c,lambda,c_cap,'-.',lambda,c_grav,'-.','LineWidth',1)

plot(lambda*100,c*100,lambda*100,c_cap*100,'-.',lambda*100,c_grav*100,'-.','LineWidth',1.5)
hold on
plot(lamdavec,phasevec(:,2:end),'*')
hold off



ylabel('$c$ [cm/s]','Interpreter','latex')
xlabel('$\lambda$ [cm]','Interpreter','latex')

legend('$c$','$c_{capillary}$','$c_{gravity}$','data','Interpreter','latex')
title('Dispersion relation: Theory vs captured data')
ylim([0 100])
xlim([0 10])
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


set(gca,'fontsize',15)