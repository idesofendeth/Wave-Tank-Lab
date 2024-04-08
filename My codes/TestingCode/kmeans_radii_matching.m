function [cluster_radii, unmatched_radii] = kmeans_radii_matching(radii1, radii2)
    % Determine the maximum number of clusters
    max_clusters = 2;
    
    % Concatenate radii into a single array and keep track of their original indices
    radii_combined = [radii1(:); radii2(:)];
    idx_radii1 = 1:length(radii1);
    idx_radii2 = (1:length(radii2)) + length(radii1);
    
    % Perform k-means clustering with a maximum of two clusters
    [centroids, ~] = kmeans(radii_combined, max_clusters);
    
    % Use pdist2 to calculate pairwise distances between data points and centroids
    distances = pdist2(radii_combined, centroids);
    
    % Find the nearest centroid for each data point
    [~, cluster_indices] = min(distances, [], 2);
    
    % Initialize cell arrays to store radii for each cluster and unmatched radii
    cluster_radii = cell(1, max_clusters);
    unmatched_radii = cell(1, 2);
    
    % Process results
    for idx = 1:length(radii_combined)
        radius = radii_combined(idx);
        cluster_idx = cluster_indices(idx);
        if ismember(idx, idx_radii1)
            radii_type = 1;
        else
            radii_type = 2;
        end
        if cluster_idx <= max_clusters
            if isempty(cluster_radii{cluster_idx})
                cluster_radii{cluster_idx} = [];
            end
            cluster_radii{cluster_idx}(end+1) = radius;
        else
            if isempty(unmatched_radii{radii_type})
                unmatched_radii{radii_type} = [];
            end
            unmatched_radii{radii_type}(end+1) = radius;
        end
    end
    
    % Remove radii from clusters that have more than two matching radii
    for i = 1:max_clusters
        if ~isempty(cluster_radii{i}) && length(cluster_radii{i}) > 2
            unmatched_radii{radii_type}(end+1:end+length(cluster_radii{i})-2) = cluster_radii{i}(3:end);
            cluster_radii{i} = cluster_radii{i}(1:2);
        end
    end
end

