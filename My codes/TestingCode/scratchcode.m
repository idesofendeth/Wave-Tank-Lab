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

%using a loop to reshape the data from the image files into a
%data matrix to be able to use SVD(POD)
datamatrix=[];
for ii=90:140%length(image)
    image1=cell2mat(image(ii));
    image1=reshape(image1,[],1);
    datamatrix=[datamatrix image1];



end

datamatrix=double(datamatrix);

dt=1/300;
[x t0]=size(datamatrix)
t=0:dt:dt*(t0-1); %time vector, FPS=300. dt=0.0033

[U,S,V]=svd(datamatrix,'econ');
%% spatial modes? Might be thinking of this wrong.

test=reshape(U(:,1),[nx nx]);
figure(1)
subplot(2,1,1)
%image1=cell2mat(image(1));
h=surf(test);
set(h,'LineStyle','none');
subplot(2,1,2)
%imshow(image1,[]);
contourf(test)

test=reshape(U(:,2),[nx nx]);
figure(2)
subplot(2,1,1)
%image1=cell2mat(image(50));
h=surf(test);
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(test)
%imshow(image1,[]);

test=reshape(U(:,3),[nx nx]);
figure(3)
subplot(2,1,1)
%image1=cell2mat(image(100));
h=surf(test);
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(test)
%imshow(image1,[]);

test=reshape(U(:,4),[nx nx]);
figure(4)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(test);
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(test)
%imshow(image1,[]);

test=reshape(U(:,5),[nx nx]);
figure(5)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(test);
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(test)
%imshow(image1,[]);
test=reshape(U(:,6),[nx nx]);
figure(6)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(test);
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(test)


figure(7)
hold on
plot(t,V(:,1))
plot(t,V(:,2))
plot(t,V(:,3))
plot(t,V(:,4))
plot(t,V(:,5))
plot(t,V(:,6))
hold off
title('Temporal modes POD')
xlabel('t')
ylabel('magnitude?')
legend('1','2','3','4','5','6')

figure(8)
image1=cell2mat(image(50));
imshow(image1,[])

figure(9)
image1=cell2mat(image(60));
imshow(image1,[])

figure(10)
image1=cell2mat(image(75));
imshow(image1,[])

figure(11)
image1=cell2mat(image(100));
imshow(image1,[])

figure(12)
image1=cell2mat(image(150));
imshow(image1,[])

%% lets try DMD

r=50;
X1=datamatrix(:,1:end-1);
X2=datamatrix(:,2:end);
% Reduce rank
% SVD rank
[Ud,Sd,Vd]=svd(X1,'econ');
Ur=Ud(:,1:r);
Sr=Sd(1:r,1:r);
Vr=Vd(:,1:r);

% Build Atilde and DMD modes
Atilde=Ur'*X2*Vr/Sr;

%% spectrum
% alternate scaling of DMD modes
Ahat = (Sr^(-1/2))*Atilde*(Sr^(1/2));
[What,D]=eig(Ahat); %solve eigenvalue problem
W_r = Sr^(1/2)*What;
Phi = X2*Vr/Sr*W_r;

P = diag(Phi'*Phi);


%% PLOT RITZ spectrum

figure(1000)
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(D)),imag(diag(D)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

%%

% DMD spectra
dt=t(2)-t(1);
lambda=diag(D);
omega=log(lambda)/dt;

% %Reconstruct the solution with DMD
% X1=datamatrix(:,1); % Initial condition
% b=Phi\X1; % This is an inverse regression.
% %But we can solve and optimization problem to obtain b: Optimized DMD
% time_dynamics=zeros(r,length(t));
% for iter=1:length(t)
%     time_dynamics(:,iter)=(b.*exp(omega*t(iter)));
% end
% X_dmd=Phi*time_dynamics;
% 
% figure(5)
% surfl(real(X_dmd'));
% shading interp;


%% Order DMD modes and frequencies as function of P
dummy1=[omega P Phi'];
dummy2=sortrows(dummy1,-2);
omega=dummy2(:,1);
P=dummy2(:,2);
Phi=dummy2(:,3:end)';

%% RRMS error
% dif=X-real(X_dmd);
% RRMSE=norm(dif,2)/norm(X,2)

%%
% DMD spectrum
figure(120)
f = abs(imag(omega(2:end))); % PLOT Freq vs. amplitude
stem(f, P(2:end), 'k');
xlabel('Freq')
ylabel('Amplitude')

axis square;

% nx=892;
% ny=892;
% for i=1:2:11
%     surf(reshape(real(Phi(:,i)),nx,ny),nx,ny);
%     plotCylinder(reshape(imag(Phi(:,i)),nx,ny),nx,ny);
% end


test=reshape(Phi(:,1),[nx nx]);
figure(13)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,2),[nx nx]);
figure(14)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,3),[nx nx]);
figure(15)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,4),[nx nx]);
figure(16)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,5),[nx nx]);
figure(17)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,6),[nx nx]);
figure(18)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,10),[nx nx]);
figure(19)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,15),[nx nx]);
figure(20)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,20),[nx nx]);
figure(21)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,25),[nx nx]);
figure(22)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))

test=reshape(Phi(:,30),[nx nx]);
figure(23)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(test))
%% analysis test
%mode3
peaks = imregionalmax(real(reshape(Phi(:,3),[nx nx])));

testpeaks=reshape(Phi(:,3),[nx nx]);
figure;
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(testpeaks));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(testpeaks))

[TF,Prom] = islocalmax(real(testpeaks));

figure;
contourf(Prom)


%mode4
peaks = imregionalmax(real(reshape(Phi(:,4),[nx nx])));

testpeaks=reshape(Phi(:,4),[nx nx]);
figure;
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(testpeaks));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(testpeaks))

[TF,Prom] = islocalmax(real(testpeaks));

figure;
contourf(Prom)


%mode5
peaks = imregionalmax(real(reshape(Phi(:,5),[nx nx])));

testpeaks=reshape(Phi(:,5),[nx nx]);
figure;
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(testpeaks));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(testpeaks))

[TF,Prom] = islocalmax(real(testpeaks));

figure;
contourf(Prom)


%mode6
peaks = imregionalmax(real(reshape(Phi(:,6),[nx nx])));

testpeaks=reshape(Phi(:,6),[nx nx]);
figure;
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(testpeaks));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(testpeaks))

[TF,Prom] = islocalmax(real(testpeaks));

figure;
contourf(Prom)

%mode7
peaks = imregionalmax(real(reshape(Phi(:,7),[nx nx])));

testpeaks=reshape(Phi(:,7),[nx nx]);
figure;
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(testpeaks));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(testpeaks))

[TF,Prom] = islocalmax(real(testpeaks));

figure;
contourf(Prom)

%mode15
peaks = imregionalmax(real(reshape(Phi(:,15),[nx nx])));

testpeaks=reshape(Phi(:,15),[nx nx]);
figure;
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(real(testpeaks));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(real(testpeaks))

[TF,Prom] = islocalmax(real(testpeaks));

figure;
contourf(Prom)


%%
yCentroid=486;
xCentroid=nx/2;
boundaries = bwboundaries(Prom);
x = boundaries{1}(:, 2);
y = boundaries{1}(:, 1);
distances = sqrt((x-xCentroid).^2 + (y-yCentroid).^2);

%% frame 1
imdex1=100;
image1=cell2mat(image(imdex1));
imagehalf=image1(:,446:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);
figure;
imshow(I{imdex1},[])
figure;
imshow(mirror,[])

[B,L,N,A] = bwboundaries(mirror);
innerdistvec=[];
outerdistvec=[];
minrad=150;
hold on
for k = 1:length(B)
   boundary = B{k};
   if(k>N)
        plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
        innerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if innerdistances>=minrad
            innerdistvec=[innerdistvec innerdistances];
        end
   else
        plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
        outerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if outerdistances>=minrad
            outerdistvec=[outerdistvec outerdistances];
            
        end
   end
    
end
sortouterdist=sort(outerdistvec)
outerdiff=(diff(sortouterdist))*scale*10^2
%outerdiff=outerdiff>10^-3

sortinnerdist=sort(innerdistvec)
innerdiff=(diff(sortinnerdist))*scale*10^2
%innerdiff=innerdiff>10^-3
avginner1=mean(innerdiff)
avgouter1=mean(outerdiff)



%% frame 2
imdex2=imdex1+1;
image1=cell2mat(image(imdex2));
imagehalf=image1(:,446:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);
figure;
imshow(I{imdex2},[])
figure;
imshow(mirror,[])

[B,L,N,A] = bwboundaries(mirror);
innerdistvec=[];
outerdistvec=[];
minrad=150;
hold on
for k = 1:length(B)
   boundary = B{k};
   if(k>N)
        plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
        innerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if innerdistances>=minrad
            innerdistvec=[innerdistvec innerdistances];
        end
   else
        plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
        outerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
        if outerdistances>=minrad
            outerdistvec=[outerdistvec outerdistances];
            
        end
   end
    
end
sortouterdist=sort(outerdistvec)
outerdiff=(diff(sortouterdist))*scale*10^2
%outerdiff=outerdiff>10^-3

sortinnerdist=sort(innerdistvec)
innerdiff=(diff(sortinnerdist))*scale*10^2
%innerdiff=innerdiff>10^-3
avginner2=mean(innerdiff)
avgouter2=mean(outerdiff)
%%
phase=(avginner2-avginner1)/dt