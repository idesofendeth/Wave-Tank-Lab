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
%yCentroid=487.580461180086;
%xCentroid=436.827061208240;
%centers=[xCentroid yCentroid];

%% frame 1
imdex1=100;
image1=cell2mat(image(imdex1));
imagehalf=image1(:,447:end);
imageflip=flipdim(imagehalf,2);
mirror=horzcat(imageflip,imagehalf);
figure;
imshow(I{imdex1},[])
figure;
imshow(mirror,[])
radiusRange=[50 300]
% finding the center of the image
[centers, radii1, metric1] = imfindcircles(mirror,radiusRange,'ObjectPolarity','bright','EdgeThreshold',0.4)
viscircles(centers, radii1,'EdgeColor','b');
yCentroid=centers(1,2);
xCentroid=centers(1,1);
% [B,L,N,A] = bwboundaries(mirror);
% innerdistvec=[];
% outerdistvec=[];
% minrad=150;
% hold on
% for k = 1:length(B)
%    boundary = B{k};
%    if(k>N)
%         plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
%         innerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
%         if innerdistances>=minrad
%             innerdistvec=[innerdistvec innerdistances];
%         end
%    else
%         plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
%         outerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
%         if outerdistances>=minrad
%             outerdistvec=[outerdistvec outerdistances];
% 
%         end
%    end
% 
% end
% sortouterdist1=sort(outerdistvec)
% outerdiff1=(diff(sortouterdist1))*scale*10^2
% %outerdiff=outerdiff>10^-3
% 
% sortinnerdist1=sort(innerdistvec)
% %Inner Diff gives lambdas
% innerdiff1=(diff(sortinnerdist1))*scale*10^2
% %innerdiff=innerdiff>10^-3
% avginner1=mean(innerdiff1)
% avgouter1=mean(outerdiff1)
% 
% 
% 
% %% frame 2
% imdex2=imdex1+1;
% image1=cell2mat(image(imdex2));
% imagehalf=image1(:,446:end);
% imageflip=flipdim(imagehalf,2);
% mirror=horzcat(imageflip,imagehalf);
% figure;
% imshow(I{imdex2},[])
% figure;
% imshow(mirror,[])
% 
% [B,L,N,A] = bwboundaries(mirror);
% innerdistvec=[];
% outerdistvec=[];
% minrad=150;
% hold on
% for k = 1:length(B)
%    boundary = B{k};
%    if(k>N)
%         plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
%         innerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
%         if innerdistances>=minrad
%             innerdistvec=[innerdistvec innerdistances];
%         end
%    else
%         plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
%         outerdistances = mean(sqrt((boundary(:,2)-xCentroid).^2 + (boundary(:,1)-yCentroid).^2));
%         if outerdistances>=minrad
%             outerdistvec=[outerdistvec outerdistances];
% 
%         end
%    end
% 
% end
% sortouterdist2=sort(outerdistvec)
% outerdiff2=(diff(sortouterdist2))*scale*10^2
% %outerdiff=outerdiff>10^-3
% 
% sortinnerdist2=sort(innerdistvec)
% innerdiff2=(diff(sortinnerdist2))*scale*10^2
% %innerdiff=innerdiff>10^-3
% avginner2=mean(innerdiff2)
% avgouter2=mean(outerdiff2)
% %%
% phase=(avginner2-avginner1)/dt




%% parent children test
minrad=150;
innerdistvec=[];
outerdistvec=[];
BW=mirror;
[B,L,N,A] = bwboundaries(BW); 
figure; imshow(BW); hold on; 
% Loop through object boundaries, this will "hopefully" find only the boundaries of proper rings  
for k = 1:N 
    % Boundary k is the parent of a hole if the k-th column 
    % of the adjacency matrix A contains a non-zero element 
    if (nnz(A(:,k)) > 0) 
        boundaryP = B{k}; 
        plot(boundaryP(:,2),... 
            boundaryP(:,1),'r','LineWidth',2); 
         outerdistances = mean(sqrt((boundaryP(:,2)-xCentroid).^2 + (boundaryP(:,1)-yCentroid).^2));
        if outerdistances>=minrad
            outerdistvec=[outerdistvec outerdistances];
            
        end
        % Loop through the children of boundary k 
        for l = find(A(:,k))' 
            boundaryC = B{l}; 
            plot(boundaryC(:,2),... 
                boundaryC(:,1),'g','LineWidth',2);
            innerdistances = mean(sqrt((boundaryC(:,2)-xCentroid).^2 + (boundaryC(:,1)-yCentroid).^2));
        if innerdistances>=minrad
            innerdistvec=[innerdistvec innerdistances];
        end
        end 
    end
    



end
hold off
%%
Crest=abs(outerdistvec+innerdistvec)/2
lamda=abs(diff(Crest))*scale*100
% figure;
% imshow(image1,[])
% hold on
viscircles(centers,Crest(1),'Color','b');
viscircles(centers,Crest(2),'Color','b');
viscircles(centers,Crest(3),'Color','b');