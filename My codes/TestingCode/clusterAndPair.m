function [matchedPairs, unmatchedValues] = pairValues(vec1, vec2)
    % Normalize the vectors
    minVec1 = min(vec1);
    maxVec1 = max(vec1);
    minVec2 = min(vec2);
    maxVec2 = max(vec2);
    vec1 = (vec1 - minVec1) / (maxVec1 - minVec1);
    vec2 = (vec2 - minVec2) / (maxVec2 - minVec2);

    % Calculate the absolute differences between all pairs of values
    diffMatrix = abs(vec1 - vec2');

    % Initialize the output variables
    matchedPairs = cell(min(length(vec1), length(vec2)), 1);
    unmatchedValues = [];

    % Loop until all values are either matched or marked as unmatched
    while ~isempty(diffMatrix)
        % Find the pair with the smallest difference
        [minDiff, linearIndices] = min(diffMatrix(:));
        [i, j] = ind2sub(size(diffMatrix), linearIndices);

        % If the smallest difference is less than a threshold, match the pair
        if minDiff < 0.01
            matchedPairs{end+1} = [vec1(i), vec2(j)];
            diffMatrix(i, :) = [];
            diffMatrix(:, j) = [];
        else
            % Otherwise, mark the value with the smaller number of remaining pairs as unmatched
            if size(diffMatrix, 1) <= size(diffMatrix, 2)
                unmatchedValues = [unmatchedValues; vec1(i)];
                diffMatrix(i, :) = [];
            else
                unmatchedValues = [unmatchedValues; vec2(j)];
                diffMatrix(:, j) = [];
            end
        end
    end

    % Unnormalize the output values
    for i = 1:length(matchedPairs)
        matchedPairs{i} = matchedPairs{i} * (maxVec1 - minVec1) + minVec1;
    end
    unmatchedValues = unmatchedValues * (maxVec1 - minVec1) + minVec1;
end

