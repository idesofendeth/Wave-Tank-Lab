function [Crest,lamda,innerdistvec,outerdistvec] = CrestFinderv2(image,centers,minradius)
%CRESTFINDER: V1, This function goes through a binarized image and finds
%the inner and outer borders of the wave crest. With these two positions
%found a average ring is formed using all position of the inner and outer
%borders. 
% 
% As input the function takes in an image, the defined centers of
%the rings in the image (found through experimentation with findCircles
%function), and a minimum radius to be considered. Anything below the
%minimum radius will not be calculated.
% 
% The output is the crest radius in vector "crest" for each ring
%found. 
% Lamda is the wave lengths found within the image between the crests 
% inner and outer dist vecs are for troubleshooting and/or more analysis if
% needed. These may not be necessary but are kept in the current version of
% the code.

% defining constants, initializing variables
yCentroid=centers(1,2);
xCentroid=centers(1,1);
minrad=minradius;
innerdistvec=[];
outerdistvec=[];
%Using bwboundaries to calculate the boundaries in the image which are then
%stored in a data structure. 
% B gives coordinates of boundary verticies in a
%cell array
%L gives labels to contiguous regions (not used here, I might not need this
%variable)
% N is the number of objects found
% A is the parent child dependencies between boundaries and holes. In this
% case I am leveraging the algorithms ability to find "holes", which in
% reality are the inner boundaries found in green. 
% The issue with the previous algorithm using findCircles was that it could
% not account for concentric circles, and had issues with finding circles
% that were not perfectly circular.
BW=image;
[B,L,N,A] = bwboundaries(BW); 
imshow(BW); hold on; 
% Loop through object boundaries, this will find only the boundaries of proper rings
% The boundaries seem to be defined as pixel positions, thus one can
% leverage this to find a radial distance from a given centroid
for k = 1:N 
    % Boundary k is the parent of a hole if the k-th column 
    % of the adjacency matrix A contains a non-zero element. Boundary k is
    % the outer boundary of the circle/crest
    if length(B{k})>500
        boundaryP = B{k}; 
        plot(boundaryP(:,2),... 
            boundaryP(:,1),'r','LineWidth',2); 
         outerdistances = mean(sqrt((boundaryP(:,2)-xCentroid).^2 + (boundaryP(:,1)-yCentroid).^2)); %this gives a mean outer radius
        if outerdistances>=minrad
            outerdistvec=[outerdistvec outerdistances];
            
        end
        % Loop through the children of boundary k, the children of a given
        % boundary of k is the inner boundary of that circle
        % for l = find(A(:,k))' 
        %     boundaryC = B{l}; 
        %     plot(boundaryC(:,2),... 
        %         boundaryC(:,1),'g','LineWidth',2);
        %     innerdistances = mean(sqrt((boundaryC(:,2)-xCentroid).^2 + (boundaryC(:,1)-yCentroid).^2)); %This gives the mean inner radius
        % if innerdistances>=minrad
        %     innerdistvec=[innerdistvec innerdistances];
        % end
        % end 
    end
    



end
hold off
%% Here a check is made to see that there are equal number of inner and outer boundaries found. Otherwise one cannot find a proper crest.
% This however does miss some data, eg when there is enough noise in the
% data that causes a ring to not be properly identified, even though one
% might be able to see a ring with ones eye. Eg, there is a break in the
% ring or the ring is noisy and the algorithm does not detect inner and
% outer edges

%outerdistvec=sort(outerdistvec);
% if length(outerdistvec)==length(innerdistvec)
%     Crest=abs(outerdistvec+innerdistvec)/2; %averaging between the outer and inner radius to find the "true" crest position, as it should be in the middle of
%     %the outer and inner boudary
    Crest=outerdistvec;
     lamda=abs(diff(Crest)); %finding the reported wavelengths
     lamda=sort(lamda);
% else
    
   % lamda=[];
% end
% 
% end

