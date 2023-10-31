%% scratch code for project course
clc
clear all
close all

%% notes
% albin stored each video as a series of images 
% in a .mat file. Each mat file is a struct, which contains cells
% each position in the cell relates to a frame in the video
%constants
scale=3.004*10^-4; %scaling factor
dt=1/300; %framerate time delta
centers=[446.497245615476 486.383816964470]; %FOUND FROM PREVIOUS EXPERIMENTATION
%% testing ideas
I = load( sprintf('Images_Marble_%dcm.mat',6) ).I; 
%I = load( sprintf('Images_MarbleADJUSTED_%dcm.mat',NR) ).I; 
%I = load( sprintf('Images_Droplet_%dcm.mat',NR) ).I; 


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
%% loop

imstart=60;
imend=190;
size=[(abs(imstart-imend)) length(I)];
phasevec=[];
lamdavec=[];
phasevec2=zeros(size);
lamdavec2=zeros(size);
j=1;
m=1;
figure;
for k=imstart:imend

%frame 1
image=I3;
imdex1=k;
image1=cell2mat(image(imdex1));
imagehalf=image1(:,447:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);
% figure;
% imshow(I{imdex1},[])
% figure;
 %imshow(mirror,[])
%radiusRange=[50 300]
% finding the center of the image
%[centers, radii1, metric1] = imfindcircles(mirror,radiusRange,'ObjectPolarity','bright','EdgeThreshold',0.4)
%viscircles(centers, radii1,'EdgeColor','b');
minradius=150;

[Crest1,lamda1, innerdistvec1,outerdistvec1] = CrestFinder(mirror,centers,minradius);
drawnow


% frame 2
imdex1=k+1;
image1=cell2mat(image(imdex1));
imagehalf=image1(:,447:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);
% figure;
% imshow(I{imdex1},[])
% figure;
% imshow(mirror,[])
%radiusRange=[50 300]
% finding the center of the image
%[centers, radii1, metric1] = imfindcircles(mirror,radiusRange,'ObjectPolarity','bright','EdgeThreshold',0.4)
%viscircles(centers, radii1,'EdgeColor','b');
minradius=150;

[Crest2,lamda2, innerdistvec2,outerdistvec2] = CrestFinder(mirror,centers,minradius);



%phase speed = abs(rt2-rt1)/dt
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
%% albin plot
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

legend('$c$','$c_{capillary}$','$c_{gravity}$','diff1','diff2','Interpreter','latex')

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