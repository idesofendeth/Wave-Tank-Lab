% Generate random vectors
vector1 = randn(1, 10);  % Arbitrary vector 1
vector2 = randn(1, 15);  % Arbitrary vector 2

% Call the function to pair values between vectors
[paired_values, unpaired_values] = pair_between_vectors(vector1, vector2);

% Display results
disp('Paired Values:');
disp(paired_values);
disp('Unpaired Values:');
disp(unpaired_values);
