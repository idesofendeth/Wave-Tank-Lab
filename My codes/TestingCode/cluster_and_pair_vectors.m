function [paired_values, unpaired_values] = pair_between_vectors(vector1, vector2)
    % Combine vectors
    combined_vector = [vector1(:); vector2(:)];

    % Normalize the combined vector
    normalized_vector = (combined_vector - mean(combined_vector)) / std(combined_vector);

    % Determine the maximum number of clusters
    max_clusters = min(length(vector1), length(vector2));

    % Perform clustering
    idx = kmeans(normalized_vector, max_clusters);

    % Separate values into clusters
    clusters = cell(max_clusters, 1);
    for i = 1:max_clusters
        clusters{i} = combined_vector(idx == i);
    end

    % Pair values between clusters
    paired_values = [];
    for i = 1:max_clusters
        cluster1 = clusters{i};
        if ~isempty(cluster1)
            % Find the closest value in other clusters
            for j = i+1:max_clusters
                cluster2 = clusters{j};
                if ~isempty(cluster2)
                    % Pair values between clusters
                    [~, min_index] = min(abs(cluster2 - cluster1(1)));
                    paired_values(end+1, :) = [cluster1(1), cluster2(min_index)];
                    % Remove paired value from cluster2
                    cluster2(min_index) = [];
                    clusters{j} = cluster2;
                    break;
                end
            end
        end
    end

    % Values that do not have a match
    unpaired_values = [];
    for i = 1:max_clusters
        unpaired_values = [unpaired_values; clusters{i}];
    end
end
