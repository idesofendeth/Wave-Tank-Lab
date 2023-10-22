%% scratch code for project course
clc
clear all
close all

%% notes
% albin stored each video as a series of images 
% in a .mat file. Each mat file is a struct, which contains cells
% each position in the cell relates to a frame in the video
%
%% testing ideas
A=load('ImageDataAlbin\Images_Droplet_1cm.mat');

image=A.I; %load cell from the mat file
image1=cell2mat(image(150)); %choose a frame, convert to matrix
%figure(1) %display
%imshow(image1,[]);
[nx ny]=size(cell2mat(image(1)));

%using a loop to reshape the data from the image files into a
%data matrix to be able to use SVD(POD)
datamatrix=[];
for ii=1:length(image)
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

r=100;
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

test=reshape(Phi(:,6),[nx nx]);
figure(19)
subplot(2,1,1)
%image1=cell2mat(image(150));
h=surf(imag(test));
set(h,'LineStyle','none');
subplot(2,1,2)
contourf(imag(test))
