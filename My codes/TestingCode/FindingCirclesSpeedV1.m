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
%%
% A=load('ImageDataAlbin\Images_Droplet_6cm.mat');
% 
 image=I3; %load cell from the mat file
% image1=cell2mat(image(150)); %choose a frame, convert to matrix
% %figure(1) %display
% %imshow(image1,[]);
[nx ny]=size(I3{1});




%% various stuff
yCentroid=486;
xCentroid=nx/2;
% range taken from albins code. before 90 there is too much noise
imstart=90;
imend=140;
avginner1vec=[];
avginner2vec=[];
avgouter1vec=[];
avgouter2vec=[];
phasevec=[];
diffvec1=[];
diffvec2=[];
%% loop
for j=imstart:imend
% frame 1
imdex1=j;
image1=cell2mat(image(imdex1));
imagehalf=image1(:,446:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);


[B,L,N,A] = bwboundaries(mirror);
innerdistvec1=[];
outerdistvec1=[];
minrad=150;
hold on
for k = 1:length(B)
   boundary = B{k};
   if(k>N)
        plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
        innerdistances1 = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if innerdistances1>=minrad
            innerdistvec1=[innerdistvec1 innerdistances1];
        end
   else
        plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
        outerdistances1 = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if outerdistances1>=minrad
            outerdistve1c=[outerdistvec1 outerdistances1];
            
        end
   end
    
end
sortouterdist1=sort(outerdistvec1);
outerdiff1=(diff(sortouterdist1))*scale*10^2;
%outerdiff=outerdiff>10^-3

sortinnerdist1=sort(innerdistvec1);
innerdiff1=(diff(sortinnerdist1))*scale*10^2;
%innerdiff=innerdiff>10^-3
avginner1=mean(innerdiff1);
avgouter1=mean(outerdiff1);

avginner1vec=[avginner1vec avginner1];
avgouter1vec=[avgouter1vec avgouter1];
% frame 2
imdex2=imdex1+1;
image1=cell2mat(image(imdex2));
imagehalf=image1(:,446:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);


[B,L,N,A] = bwboundaries(mirror);
innerdistvec2=[];
outerdistvec2=[];
minrad=150;
hold on
for k = 1:length(B)
   boundary = B{k};
   if(k>N)
        plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
        innerdistances2 = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if innerdistances2>=minrad
            innerdistvec2=[innerdistvec2 innerdistances2];
        end
   else
        plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
        outerdistances2 = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if outerdistances2>=minrad
            outerdistvec2=[outerdistvec2 outerdistances2];
            
        end
   end
    
end
sortouterdist2=sort(outerdistvec2);
outerdiff2=(diff(sortouterdist2))*scale*10^2;
%outerdiff=outerdiff>10^-3

sortinnerdist2=sort(innerdistvec2);
innerdiff2=(diff(sortinnerdist2))*scale*10^2;
%innerdiff=innerdiff>10^-3
avginner2=mean(innerdiff2);
avginner2vec=[avginner2vec avginner2];
avgouter2=mean(outerdiff2);
avgouter2vec=[avgouter2vec avgouter2];
%
%phase=abs((avginner2-avginner1))/2*dt;


ldex=min(length(innerdistvec1),length(innerdistvec2));
for ii=1:ldex-1
    phase=abs(innerdistvec2(ii)-innerdistvec1(ii));
    phasevec=[phasevec phase];
    
end
indiff2=abs(diff(innerdistvec2(1:ldex)));
indiff1=abs(diff(innerdistvec1(1:ldex)));
diffvec1=[diffvec1 indiff1];
diffvec2=[diffvec2 indiff2];

end
%%
sortp=sort(phasevec)/dt*scale*100;
revp=flip(sortp);
sortdiff=sort(diffvec1);%sort(mean([diffvec1' diffvec2'],2));
% figure;
% plot(sortdiff,revp)
% xlim([0 20])

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
hold on
plot(lambda*100,c*100,lambda*100,c_cap*100,'-.',lambda*100,c_grav*100,'-.','LineWidth',1.5)
plot(sortdiff,revp,'-*')
hold off



ylabel('$c$ [cm/s]','Interpreter','latex')
xlabel('$\lambda$ [cm]','Interpreter','latex')

legend('$c$','$c_{capillary}$','$c_{gravity}$','ExpData','Interpreter','latex')
%     ylim([0 100])
     xlim([0 6])
        ylim([0 100])
    %xlim([0 6])
% disp('Shallow water approx, kH<<1, leads to tanh(kH)~kH')
% disp('Deep water approx, kH>>1, leads to tanh(kH)~1')
% %Deep water surface waves ( kH>>1 )
% c_cap_deep = sqrt( (2* pi* sigma)./ (rho*lambda)); %Pure capillary wave  (dispersive)
% c_grav_deep = sqrt( (g* lambda)./ (2*pi) ); %Pure gravity wave (dispersive)
% %Shallow water surface waves ( kH<<1 )
% c_cap_shallow = 2* pi./lambda .*sqrt( sigma.*H/ rho); %Pure capillary wave (dispersive)
% c_grav_shallow = sqrt(g* H); % Pure gravity wave (non dispersive)
lambda_min = 2*pi*sqrt(sigma/(rho*g))*10^2
lambda_cap = sqrt(sigma/(rho*g))*10^2
%2*pi*sqrt(lambda_)


set(gca,'fontsize',15) 