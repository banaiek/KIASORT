function [optimalEps, k] = estimate_dbscan_par(data, k)
% ESTIMATE_DBSCAN_PAR - Estimate optimal DBSCAN epsilon using k-distance graph
%
% This function uses the "elbow method" on the k-distance graph to find
% the optimal epsilon parameter for DBSCAN clustering.
%
% INPUTS:
%   data - Data matrix [N x D]
%   k    - Number of neighbors (optional, default: D+1)
%
% OUTPUTS:
%   optimalEps - Estimated optimal epsilon value
%   k          - Number of neighbors used
%
% ALGORITHM:
%   1. Compute k-th nearest neighbor distance for each point
%   2. Sort distances in ascending order
%   3. Find the "knee" point using perpendicular distance method
%   4. Return the distance at the knee as optimal epsilon

[N, D] = size(data);

% Default k based on dimensionality
if nargin < 2 || isempty(k)
    k = min(D + 1, max(3, round(log2(N))));
end

% Ensure k is valid
k = max(2, min(k, N - 1));

%% ========================================================================
%  COMPUTE K-NEAREST NEIGHBOR DISTANCES
%  ========================================================================
try
    % Use Minkowski distance (L1 norm) for robustness
    [~, distances] = knnsearch(data, data, 'K', k + 1, 'Distance', 'minkowski', 'P', 1);
catch
    % Fallback to Euclidean if Minkowski fails
    [~, distances] = knnsearch(data, data, 'K', k + 1, 'Distance', 'euclidean');
end

% k-th nearest neighbor distance (exclude self at position 1)
kDistances = distances(:, end);

%% ========================================================================
%  SORT AND PREPARE FOR KNEE DETECTION
%  ========================================================================
sortedKDistances = sort(kDistances);
n = length(sortedKDistances);

% Remove any infinite or NaN values
validIdx = isfinite(sortedKDistances);
sortedKDistances = sortedKDistances(validIdx);
n = length(sortedKDistances);

if n < 3
    optimalEps = median(kDistances(isfinite(kDistances)));
    return;
end

%% ========================================================================
%  KNEE DETECTION USING PERPENDICULAR DISTANCE METHOD
%  ========================================================================
x = (1:n)';
y = sortedKDistances;

% Normalize to [0,1] for better numerical stability
x_norm = (x - 1) / (n - 1);
y_norm = (y - y(1)) / (y(end) - y(1) + eps);

% Line from first to last point
lineVec = [1, y_norm(end) - y_norm(1)];
lineVecNorm = norm(lineVec);

% Compute perpendicular distance from each point to the line
distancesToLine = abs((x_norm - 0) * lineVec(2) - (y_norm - y_norm(1)) * lineVec(1)) / lineVecNorm;

% Find the knee point (maximum perpendicular distance)
[~, kneeIdx] = max(distancesToLine);

%% ========================================================================
%  VALIDATE AND REFINE KNEE DETECTION
%  ========================================================================
% Apply smoothing to avoid noise-induced false knees
windowSize = max(3, round(n * 0.05));
smoothedDistances = movmean(distancesToLine, windowSize);
[~, smoothedKneeIdx] = max(smoothedDistances);

% Use the more conservative (larger epsilon) knee if they differ significantly
if abs(kneeIdx - smoothedKneeIdx) > n * 0.1
    kneeIdx = max(kneeIdx, smoothedKneeIdx);
end

% Ensure knee is not at extreme ends
kneeIdx = max(round(n * 0.05), min(kneeIdx, round(n * 0.95)));

optimalEps = sortedKDistances(kneeIdx);

%% ========================================================================
%  SANITY CHECKS
%  ========================================================================
% Ensure epsilon is not too small (would result in many noise points)
minEps = prctile(kDistances, 5);
optimalEps = max(optimalEps, minEps);

% Ensure epsilon is not too large (would merge everything)
maxEps = prctile(kDistances, 95);
optimalEps = min(optimalEps, maxEps);

end