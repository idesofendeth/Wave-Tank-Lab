function [cluster1_radii,cluster2_radii] = clusterRadii(radii,k)
N=length(radii);

% Normalize the data


% Perform k-means clustering
%k = 2; % Number of clusters
% Normalize the data
normalized_radii = normalize(radii);

% Perform k-means clustering
k = 2; % Number of clusters
[idx, centroids_normalized] = kmeans(normalized_radii, k);

% Un-normalize centroids
centroids = centroids_normalized * std(radii) + mean(radii);

% Un-normalize radii
denormalized_radii = normalized_radii * std(radii) + mean(radii);

% Sort denormalized radii into two distinct columns based on cluster assignment
cluster1_radii = denormalized_radii(idx == 1);
cluster2_radii = denormalized_radii(idx == 2);

% Plot radii and centroids together
figure;
scatter(denormalized_radii, zeros(N,1), 50, idx, 'filled'); % Radii colored by cluster
hold on;
scatter(centroids, zeros(k,1), 100, 'r', 'filled'); % Centroids
xlabel('Radii');
ylabel('Centroids');
title('Radii and Centroids');
legend('Radii', 'Centroids');
grid on;

end
