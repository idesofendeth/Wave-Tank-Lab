%% Calculate the wavecrests radii
%Choose the all x-values for vertical coordinates y0+-10 for all matrices in I3
Xmid = cellfun(@(x) x(y0-10:y0+10,:), I3,'UniformOutput',false);

%Cleaning up the binary image again, removes components smaller than 20 pixels
Xmid = cellfun(@(x) bwareaopen(x,20), Xmid,'UniformOutput',false);

%Find the centroid of the waves,. i.e. the x-position of the crests
props = cellfun(@(x) regionprops(x, 'Centroid'), Xmid,'UniformOutput',false); %properties given in structs
centroids = cellfun(@(x) cat(1,x.Centroid), props,'UniformOutput',false); %convert structs to double
centroids(cellfun(@isempty,centroids)) = {nan}; %This is needed, if there are any empty cells then r_crests does´nt work

x_crests = cellfun(@(x) x(:,1), centroids,'UniformOutput',false); %The x-position of the crests
x_crests_lhs = cellfun(@(x) x(x<x0-r_min), x_crests, 'UniformOutput',false); %Make sure to not keep too small radii´s. And sort the radii in ascending order
x_crests_rhs = cellfun(@(x) x(x>x0+r_min), x_crests, 'UniformOutput',false); %Make sure to not keep too small radii´s. And sort the radii in ascending order

r_crests = cellfun(@(x) abs(x(:,1)-x0), centroids,'UniformOutput',false); %The radial distance from origin to crest
r_crests = cellfun(@(x) x(x>r_min), r_crests, 'UniformOutput',false); %Make sure to not keep too small radii´s. And sort the radii in ascending order

%Split into rhs and lhs:
r_crests_lhs = cellfun(@(x) abs(x-x0), x_crests_lhs, 'UniformOutput',false); %Make sure to not keep too small radii´s. And sort the radii in ascending order
r_crests_rhs = cellfun(@(x) abs(x-x0), x_crests_rhs, 'UniformOutput',false); %Make sure to not keep too small radii´s. And sort the radii in ascending order

%Assemble rhs and lhs into a Nx2 cell array
R = [r_crests_lhs(t1:t2)',r_crests_rhs(t1:t2)'];


