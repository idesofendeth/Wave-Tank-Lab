function [sorted_radii, unnormalized_centroids] = sort_radii_with_kmeans(radii, k)
    % Normalize the radii
    normalized_radii = (radii - min(radii)) / (max(radii) - min(radii));

    % Perform k-means clustering
    [idx, centroids] = kmeans(normalized_radii, k,'Replicates',5);

    % Unnormalize the centroids
    unnormalized_centroids = centroids * (max(radii) - min(radii)) + min(radii);

    % Sort radii into k distinct columns
    sorted_radii = cell(1, k);
    for i = 1:k
        cluster_indices = find(idx == i);
        sorted_radii{i} = radii(cluster_indices);
    end

    % Sort radii within each cluster
    for i = 1:k
        sorted_radii{i} = sort(sorted_radii{i});
    end

    % Plot unnormalized centroids and unnormalized data together
% figure;
% hold on;
% scatter(radii, zeros(size(radii)), 'bo', 'DisplayName', 'Data');
% scatter(unnormalized_centroids, zeros(size(unnormalized_centroids)), 'rx', 'DisplayName', 'Centroids');
% legend;
% xlabel('Radii');
% title('Unnormalized Centroids and Unnormalized Data');
% hold off;
end
