vec1 = RadiusDerek1;
vec2 = RadiusDerek2;
[matchedPairs, unmatchedValues] = clusterAndPair(vec1, vec2);
%% matlabgpt

% Generate two vectors with arbitrary values and lengths
vector1 = randi([1, 10], [1, 8]);
vector2 = randi([1, 10], [1, 6]);

% Call the pairClustering function
clusterPairs = pairClustering(vector1, vector2);
