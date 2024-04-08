% Sample radii matrix
radii_matrix = [1.2, 3.4, -100, 2.1, 5.6;
                1.5, -100, 4.2, 6.3, 3.9;
                1.5, -100, 4.2, 6.3, 3.9];

% Number of clusters
k = size(radii_matrix,1)-1;

% Call the function
radii_sorted = sort_radii_with_kmeans(radii_matrix, k);

% Display the sorted radii
disp('Sorted radii:');
disp(radii_sorted);
