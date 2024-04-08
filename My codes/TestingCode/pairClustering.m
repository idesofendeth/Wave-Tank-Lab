function clusterPairs = pairClustering(vector1, vector2)
    % Determine the number of clusters
    numClusters = min(length(vector1), length(vector2));

    % Normalize the vectors
    vector1 = normalize(vector1);
    vector2 = normalize(vector2);

    % Combine the vectors into a single matrix
    data = [vector1(:), vector2(:)];

    % Check if the number of clusters is greater than the number of data points
    if numClusters > size(data, 1)
        error('Number of clusters cannot be greater than the number of data points.');
    end

    % Perform k-means clustering
    [~, centroids] = kmeans(data, numClusters);

    % Assign each data point to its nearest centroid
    [~, clusterIndices] = pdist2(centroids, data, 'euclidean', 'Smallest', 1);

    % Separate the pairs and unmatched values into separate vectors
    clusterPairs = cell(numClusters, 1);
    unmatchedValues = [];
    for i = 1:numClusters
        clusterPairs{i} = data(clusterIndices == i, :);
        if isempty(clusterPairs{i})
            unmatchedValues = [unmatchedValues; centroids(i, :)];
        end
    end

    % Plot the scatter plot with centroids
    scatter(data(:, 1), data(:, 2), 'filled');
    hold on;
    scatter(centroids(:, 1), centroids(:, 2), 'r', 'filled');
    hold off;
    legend('Data Points', 'Centroids');
    xlabel('Feature 1');
    ylabel('Feature 2');
    title('Clustering Results');
end
