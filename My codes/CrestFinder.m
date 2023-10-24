function [Crest,lamda,innerdistvec,outerdistvec] = CrestFinder(image,centers,minradius)
%CRESTFINDER Summary of this function goes here
%   Detailed explanation goes here
yCentroid=centers(1,2);
xCentroid=centers(1,1);
minrad=minradius;
innerdistvec=[];
outerdistvec=[];
BW=image;
[B,L,N,A] = bwboundaries(BW); 
imshow(BW); hold on; 
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
if length(outerdistvec)==length(innerdistvec)
    Crest=abs(outerdistvec+innerdistvec)/2
    lamda=abs(diff(Crest))
else
    Crest=[];
    lamda=[];
end

end

