function [Centers,Radius,Radiuskmeans,radiusSortDerek] = circleFitterV2(XPeakVec,YPeakVec,scale,centerIm)
%% Least squares circle fitting and data processing
% input: XPeakVec and YPeakVec from previous
% function. These are then cleaned up and
% non-physical points are removed
% output: Centers  and Radius: x and y center positions with radii for each
% reported circle, Kmeans and Derek Threshold Solution are two sorting
% methods that sort the data for each circle into their own vectors.
% Previous implementation could not discern between datapoints for
% different circles.
xCentVec=[];
yCentVec=[];
radiusVec=[];
radMatrix=ones(size(XPeakVec,1),size(XPeakVec,2))*-100;
for l=1:size(XPeakVec,2)
    A=[XPeakVec(:,l) YPeakVec(:,l)];
    indA=find(A(:,1)>=-50); % values that are not found are set to -100, this index finds anything larger, could be set to -99, doesnt matter
    B=A(indA,:);
    radiusVectemp=sqrt(B(:,1).^2+B(:,2).^2);
    %[cluster1,cluster2]=clusterRadii(radiusVectemp);
    outliers=radiusVectemp>=1.02*mean(radiusVectemp) | radiusVectemp<=0.98*mean(radiusVectemp); % setting outliers to be +-2 procent of data
    B=B(~outliers,:); %remove outliers
%     figure(101),clf
%     hold on
%     grid on
    %plot(B(:,1),B(:,2),'*','MarkerSize',10)
    %[xCenter, yCenter, radiusFit, a] = circlefit(XPeakVec(:,l), YPeakVec(:,l))
    [xCenter, yCenter, radiusFit, a] = circlefit(B(:,1), B(:,2));

    xCentVec=[xCentVec; xCenter];
    yCentVec=[yCentVec; yCenter];
    if l<size(XPeakVec,2)
        if radiusFit<=centerIm*scale
            radiusVec=[radiusVec; radiusFit];
        end
    end
    radMatrix(1:size(radiusVectemp,1),l)=radiusVectemp;
    clear A B radiusVectemp
    %viscircles([0 0],radiusFit,'LineWidth',1)
    %pause
end

%output goes here
radMatrixSort=radMatrix~=-100;
[radMatrixIndx,radMatrixIndy]=find(radMatrixSort);
dimMax=max(radMatrixIndy);
radVec=reshape(radMatrix,[numel(radMatrix),1]);
radVecred=radVec(radVec~=-100,:);
radVecred=radVecred(radVecred~=0,:);
[sorted_radii, unnormalized_centroids] = sort_radii_with_kmeans(radVecred, dimMax);
%isnan removes NaN entries from bad data
meanYCent=mean(rmoutliers(yCentVec(~isnan(yCentVec))));
meanXCent=mean(rmoutliers(xCentVec(~isnan(xCentVec))));
% kmeans stuff
for ii=1:dimMax
    
    Radiuskmeans(ii,1)=mean(cell2mat(sorted_radii(ii)));
end
Radiuskmeans=sort(Radiuskmeans);

% dereks idea
radVecred=sort(radVecred);
radiusIndex=[1];
radiusSortDerek=[];
for jj=1:length(radVecred)-1
    if diff([radVecred(jj) radVecred(jj+1)])>=0.02*radVecred(jj)
        radiusIndex=[radiusIndex,jj];
        radiusSortDerek=[radiusSortDerek; mean(radVecred(radiusIndex(end-1):radiusIndex(end)))];

    end


end


%output

Radius=radiusVec;
Centers=[meanXCent meanYCent];
end