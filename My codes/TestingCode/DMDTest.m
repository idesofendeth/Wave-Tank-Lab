%% scratch code for project course
clc
clear all
close all

%% notes
% albin stored each video as a series of images 
% in a .mat file. Each mat file is a struct, which contains cells
% each position in the cell relates to a frame in the video
%% constants
scale=3.004*10^-4; %scaling factor
scale=3.004*10^-4; %scaling factor
dt=1/300; %framerate time delta
centers=[446.497245615476 486.383816964470]; %FOUND FROM PREVIOUS EXPERIMENTATION using findCircles function

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
    % 
    % %Threshold image-------------------------------------------------------
    I3{i} = I3{i} > thres;

    %Clean up image--------------------------------------------------------
    I3{i} = bwareaopen(I3{i},50) ;
    I3{i} = imclearborder(I3{i});
    % %edgesmoothing
    % % windowSize = 51;
    % % kernel = ones(windowSize) / windowSize ^ 2;
    % % blurryImage = conv2(single(I3{i}), kernel, 'same');
    % % I3{i} = blurryImage > thres; % Rethreshold

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

%[U,S,V]=svd(datamatrix,'econ');

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

%%
% test=reshape(Phi(:,1),[nx nx]);
% figure(13)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,2),[nx nx]);
% figure(14)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,3),[nx nx]);
% figure(15)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,4),[nx nx]);
% figure(16)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,5),[nx nx]);
% figure(17)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,6),[nx nx]);
% figure(18)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,10),[nx nx]);
% figure(19)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,15),[nx nx]);
% figure(20)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,20),[nx nx]);
% figure(21)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,25),[nx nx]);
% figure(22)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
% 
% test=reshape(Phi(:,30),[nx nx]);
% figure(23)
% subplot(2,1,1)
% %image1=cell2mat(image(150));
% h=surf(real(test));
% set(h,'LineStyle','none');
% subplot(2,1,2)
% contourf(real(test))
%% analysis test
%mode1
peaks = imregionalmax(real(reshape(Phi(:,1),[nx nx])));

testpeaks=reshape(Phi(:,1),[nx nx]);
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
%mode2
peaks = imregionalmax(real(reshape(Phi(:,2),[nx nx])));

testpeaks=reshape(Phi(:,2),[nx nx]);
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

%mode9
peaks = imregionalmax(real(reshape(Phi(:,9),[nx nx])));

testpeaks=reshape(Phi(:,9),[nx nx]);
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

%% test idea: Find distance from center of rings to each point in PROM
% This should, in theory give a radius or radii that can be compared giving
% the wavelengths captured in the DMD algorithm 
% see documentation for islocalmax
close all

% test for mode 3
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

% chatgpt idea

% Find the non-zero values and their positions
[row, col] = find(Prom>0.04);

% Calculate distances from each non-zero value to the center
distances = sqrt((row - centers(1)).^2 + (col - centers(2)).^2);
minradius=150;
% Define a tolerance for considering points with similar distances
tolerance = 25;  % You can adjust the tolerance as needed

% Calculate the average distance of all non-zero points
avg_distance = mean(distances);
difftances=sort(abs(distances));
unique_distances = unique(distances);
% Initialize cell arrays to store distances in different vectors
distances_by_avg_distance = cell(1, length(unique(distances)));

for i = 1:length(unique_distances)
    current_distance = unique_distances(i);
    points_within_tolerance = find(abs(distances - current_distance) >= tolerance);
    distances_by_avg_distance{i} = points_within_tolerance;
end


%% binarization idea

% Define the radius of the rings to be ignored
ignore_radius = 100;  % You can adjust the radius value as needed
matrix=Prom;
% Create a binary mask of non-zero values
binary_mask = matrix > 0;

% Use the bwlabel function to label connected components in the binary mask
labeled_matrix = bwlabel(binary_mask);

% Exclude rings with the specified radius and calculate remaining ring radii
unique_labels = unique(labeled_matrix);
ignore_labels = [];
ring_radii = [];

for i = 1:length(unique_labels)
    label = unique_labels(i);
    if label == 0
        continue;  % Skip background label
    end
    
    [r, c] = find(labeled_matrix == label);
    
    % Calculate the centroid of the ring
    centroid = [mean(r), mean(c)];
    
    % Calculate the distance from the centroid to the center
    distance_to_center = sqrt((centroid(1) - (size(matrix, 1) / 2))^2 + (centroid(2) - (size(matrix, 2) / 2))^2);
    
    if distance_to_center <= ignore_radius
        % Mark this label to be ignored
        ignore_labels = [ignore_labels, label];
    else
        % Store the radius of the remaining ring
        ring_radii = [ring_radii, distance_to_center];
    end
end

% Display the labeled matrix
disp('Labeled Matrix (after excluding specified radius):');
disp(labeled_matrix);

% Count the number of remaining rings
remaining_labels = unique(labeled_matrix);
num_remaining_rings = length(remaining_labels) - 1;  % Subtract 1 to exclude 0 label

% Display the number of remaining rings
%disp('Number of Remaining Rings: ', num2str(num_remaining_rings);

% Display the radii of the remaining rings
%disp('Radii of Remaining Rings:');
%disp(ring_radii);

%%

% Define the centroid coordinates
given_centroid = [446, 486];

% Define the search radius for the rings
search_radius = 2;  % You can adjust the radius value as needed
matrix=Prom>0.04;
% Create a binary mask of non-zero values
binary_mask = matrix > 0;

% Use the bwlabel function to label connected components in the binary mask
labeled_matrix = bwlabel(binary_mask);

% Find the label of the ring containing the given centroid
centroid_label = labeled_matrix(given_centroid(1), given_centroid(2));

% Identify rings within the search radius of the given centroid
rings_within_radius = unique(labeled_matrix(abs(labeled_matrix - centroid_label) <= search_radius));

% Exclude the background label (0)
rings_within_radius = rings_within_radius(rings_within_radius ~= 0);
% Calculate the radii of the identified rings
ring_radii = zeros(1, length(rings_within_radius));

for i = 1:length(rings_within_radius)
    current_ring_label = rings_within_radius(i);
    [row, col] = find(labeled_matrix == current_ring_label);
    
    % Calculate the radius based on the maximum distance from the centroid
    distances = sqrt((row - given_centroid(1)).^2 + (col - given_centroid(2)).^2);
    ring_radii(i) = max(distances);
end

% Display the labeled matrix
disp('Labeled Matrix:');
disp(labeled_matrix);

% Display the labels of the rings within the search radius
disp('Labels of Rings within Search Radius:');
disp(rings_within_radius);

% Display the radii of the identified rings
disp('Radii of Identified Rings:');
disp(ring_radii);
%%
% Create a binary mask of non-zero values
binary_mask = matrix > 0;

% Use the bwlabel function to label connected components in the binary mask
labeled_matrix = bwlabel(binary_mask);

% Calculate the radii of the identified rings
unique_labels = unique(labeled_matrix);
ring_radii = [];

for i = 1:length(unique_labels)
    current_ring_label = unique_labels(i);
    
    if current_ring_label == 0
        continue;  % Skip background label
    end
    
    [row, col] = find(labeled_matrix == current_ring_label);
    
    % Calculate the radius based on the maximum distance from the centroid
    distances = sqrt((row - given_centroid(1)).^2 + (col - given_centroid(2)).^2);
    ring_radius = max(distances);
    
    % Check if the radius is within the desired range
    if ring_radius >= 100 && ring_radius <= 600
        ring_radii = [ring_radii, ring_radius];
    end
end

% Plot the matrix with rings
figure;
imshow(matrix, 'InitialMagnification', 'fit');
hold on;

% Plot circles around the rings
for i = 1:length(ring_radii)
    ring_radius = ring_radii(i);
    viscircles(given_centroid, ring_radius, 'EdgeColor', 'r');
end

title('Matrix with Identified Rings');
hold off;

%%

given_centroid = [446, 486];

% Define the minimum and maximum ring radii
min_radius = 100;
max_radius = 600;

% Create a binary mask of non-zero values
binary_mask = matrix > 0;

% Use the regionprops function to find and analyze connected components (rings)
stats = regionprops(binary_mask, 'Centroid', 'MajorAxisLength', 'MinorAxisLength');

% Filter dominant rings based on their properties
dominant_rings = [];
for i = 1:length(stats)
    centroid = stats(i).Centroid;
    major_axis_length = stats(i).MajorAxisLength;
    minor_axis_length = stats(i).MinorAxisLength;
    ring_radius = (major_axis_length + minor_axis_length) / 4;
    
    % Check if the ring meets the criteria
    if ring_radius >= min_radius && ring_radius <= max_radius
        dominant_rings = [dominant_rings; centroid, ring_radius];
    end
end

% Plot the matrix with dominant rings
figure;
imshow(matrix, 'InitialMagnification', 'fit');
hold on;

% Plot circles around the dominant rings
for i = 1:size(dominant_rings, 1)
    centroid = dominant_rings(i, 1:2);
    ring_radius = dominant_rings(i, 3);
    viscircles(centroid, ring_radius, 'EdgeColor', 'r');
end

title('Matrix with Dominant Rings');
hold off;

%%
binary_mask = Prom > 0;
inputMatrix=binary_mask;
% Define the center and minimum radius
center = [446, 486];
minRadius = 100;

% Create a meshgrid to compute distances from the center
[x, y] = meshgrid(1:size(inputMatrix, 2), 1:size(inputMatrix, 1));
distances = sqrt((x - center(1)).^2 + (y - center(2)).^2);

% Initialize an empty matrix for the interpolated result
interpolatedMatrix = zeros(size(inputMatrix));

% Find the ring regions
ringRegion = (distances > minRadius);
ringMatrix = inputMatrix .* ringRegion;

% Interpolate the ring regions
for r = minRadius:max(distances(:))
    ringMask = (distances == r) & ringRegion;
    
    % Check if any non-zero elements exist in the ring
    if any(ringMask(:))
        % Interpolate the ring using a 2D convolution with a kernel
        interpolatedRing = conv2(double(ringMask), ones(3), 'same') > 0;
        interpolatedMatrix(interpolatedRing) = interp2(x(ringMask), y(ringMask), inputMatrix(ringMask), x(interpolatedRing), y(interpolatedRing), 'linear');
    end
end

% Display or save the interpolated matrix
imshow(interpolatedMatrix, [])