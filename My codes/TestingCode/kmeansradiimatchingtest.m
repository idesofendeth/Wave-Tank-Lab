% Example usage
radii1 = RadiusDerek1;
radii2 = RadiusDerek2;

%% gpt test one
[cluster_radii, unmatched_radii] = kmeans_radii_matching(radii1, radii2);

% Plot clusters and radii
colors = {'r', 'b', 'g', 'm', 'c'};
figure;
hold on;
for i = 1:length(cluster_radii)
    scatter(cluster_radii{i}, ones(size(cluster_radii{i})) * i, 50, colors{i}, 'filled');
end
for i = 1:2
    scatter(unmatched_radii{i}, ones(size(unmatched_radii{i})) * (length(cluster_radii) + i), 50, colors{length(cluster_radii) + i}, 'filled');
end
hold off;
xlabel('Radius');
ylabel('Cluster');
title('Clustered Radii');
legend('Cluster 1', 'Cluster 2', 'Unmatched Radii from Radii 1', 'Unmatched Radii from Radii 2');
%% gpt test two

[paired_values, unpaired_values] = cluster_and_pair_vectors(radii1, radii2)
