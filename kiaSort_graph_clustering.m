function labels = kiaSort_graph_clustering(X, varargin)
% 
%        labels = advancedGraphClustering(X, 'existingLabels', dbscan_labels, ...)
%
% Inputs:
%   X - Data points
%
% Optional Parameters:
%   'existingLabels' labels from dbscan
%   'k' - number of nearest neighbors
% Output:
%   labels - Cluster labels (1 to K for clusters, -1 for noise points)


p = inputParser;
addRequired(p, 'X', @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'existingLabels', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'k', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'minClusterSize', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'qualityThreshold', 0.8, @(x) isnumeric(x) && isscalar(x));
parse(p, X, varargin{:});

[N, D] = size(X);

if isempty(p.Results.k)
    k = max(3, min(20, round(log2(N) + 1)));
else
    k = p.Results.k;
end

if isempty(p.Results.minClusterSize)
    minClusterSize = max(5, round(0.005*N));
else
    minClusterSize = p.Results.minClusterSize;
end

existingLabels = p.Results.existingLabels;
qualityThreshold = p.Results.qualityThreshold;

%  clustering quality 
if ~isempty(existingLabels) && length(existingLabels) == N
    quality = evaluateClusteringQuality(X, existingLabels, k);
    
    if quality >= qualityThreshold
        labels = existingLabels;
        return;
    end
end

% graph-based clustering

[W, localDensity] = buildOptimizedGraph(X, k);


noisePoints = detectAdvancedNoise(X, localDensity, k);

validPoints = ~noisePoints;
if sum(validPoints) < minClusterSize
    labels = -ones(N, 1);
    return;
end

X_clean = X(validPoints, :);
W_clean = W(validPoints, validPoints);

% optimal number of clusters using stability analysis
maxClusters = min(round(sqrt(sum(validPoints))), 20);
optimalK = findOptimalClustersAdvanced(W_clean, X_clean, maxClusters, minClusterSize);

% spectral clustering
cleanLabels = spectralClusteringRobust(W_clean, optimalK);

%  post-processing 
cleanLabels = mergeSmallClusters(X_clean, cleanLabels, minClusterSize);
cleanLabels = mergeOverSplitClusters(X_clean, cleanLabels);

% final labels
labels = -ones(N, 1);
labels(validPoints) = cleanLabels;

% fewer clusters use initial labels
if (~isempty(existingLabels) && max(labels) < max(existingLabels(existingLabels > 0))) || (~isempty(existingLabels) && (sum(labels==-1)/length(labels))>.2)

    labels = existingLabels;
end
end

function quality = evaluateClusteringQuality(X, labels, k)
% clustering quality evaluation 

if isempty(labels) || all(labels <= 0)
    quality = 0;
    return;
end

validPoints = labels > 0;
if sum(validPoints) < 10
    quality = 0;
    return;
end

X_valid = X(validPoints, :);
labels_valid = labels(validPoints);
uniqueLabels = unique(labels_valid);
numClusters = length(uniqueLabels);

if numClusters < 2
    quality = 0.5;
    return;
end

%  silhouette 
sampleSize = min(300, size(X_valid, 1));
if size(X_valid, 1) > sampleSize
    idx = randperm(size(X_valid, 1), sampleSize);
    X_sample = X_valid(idx, :);
    labels_sample = labels_valid(idx);
else
    X_sample = X_valid;
    labels_sample = labels_valid;
end

try
    silScore = mean(silhouette(X_sample, labels_sample));
    silScore = max(0, silScore);
catch
    silScore = 0;
end

% density-based
densityScore = computeDensityScore(X_valid, labels_valid, k);

quality = 0.7 * silScore + 0.3 * densityScore;
quality = max(0, min(1, quality));
end

function densityScore = computeDensityScore(X, labels, k)
%  density-based quality 

try
    [~, D] = knnsearch(X, X, 'K', k+1);
    kDist = D(:, end);
    
    uniqueLabels = unique(labels);
    intraScore = 0;
    interScore = 0;
    
    for i = 1:length(uniqueLabels)
        clusterMask = labels == uniqueLabels(i);
        if sum(clusterMask) < 2, continue; end
        
        clusterKDist = kDist(clusterMask);
        intraScore = intraScore + 1 / (mean(clusterKDist) + eps);
        
        if i < length(uniqueLabels)
            clusterPoints = X(clusterMask, :);
            otherPoints = X(~clusterMask, :);
            
            nSample = min(15, size(clusterPoints, 1));
            if size(clusterPoints, 1) > nSample
                sampleIdx = randperm(size(clusterPoints, 1), nSample);
                samplePoints = clusterPoints(sampleIdx, :);
            else
                samplePoints = clusterPoints;
            end
            
            if ~isempty(otherPoints) && ~isempty(samplePoints)
                minDists = min(pdist2(samplePoints, otherPoints), [], 2);
                interScore = interScore + mean(minDists);
            end
        end
    end
    
    densityScore = (intraScore / length(uniqueLabels)) / (interScore / length(uniqueLabels) + 1);
    densityScore = min(1, densityScore / 8);
    
catch
    densityScore = 0;
end
end

function [W, localDensity] = buildOptimizedGraph(X, k)
% k-NN graph 

[N, ~] = size(X);

[knnIdx, knnDist] = knnsearch(X, X, 'K', k+1);
knnIdx = knnIdx(:, 2:end); 
knnDist = knnDist(:, 2:end);

% local density using k-distance
localDensity = 1 ./ (mean(knnDist, 2) + eps);
localDensity = localDensity / max(localDensity);

I = repmat((1:N)', 1, k);
J = knnIdx;
weights = zeros(N, k);

for i = 1:N
    neighbors = knnIdx(i, :);
    dists = knnDist(i, :);
    
    sigma = mean(dists);
    distWeight = exp(-dists.^2 / (2 * sigma^2));
    
    densityWeight = sqrt(localDensity(i) * localDensity(neighbors)');
    
    weights(i, :) = distWeight .* densityWeight;
end

W = sparse([I(:); J(:)], [J(:); I(:)], [weights(:); weights(:)], N, N);
W = max(W, W'); 
end

function noisePoints = detectAdvancedNoise(X, localDensity, k)
% noise detection

[N, ~] = size(X);

% Local density threshold
densityThreshold = quantile(localDensity, 0.08);
lowDensity = localDensity < densityThreshold;

% isolation-based detection
[~, knnDist] = knnsearch(X, X, 'K', k+1);
kDist = knnDist(:, end); % k-th nearest neighbor distance
isolationScore = kDist / median(kDist);
isolationThreshold = quantile(isolationScore, 0.92);
isolated = isolationScore > isolationThreshold;

% outlier 
lof = computeOptimizedLOF(X, k);
lofThreshold = quantile(lof, 0.95);
outliers = lof > lofThreshold;


noisePoints = lowDensity & (isolated | outliers);

extremeIsolation = isolationScore > quantile(isolationScore, 0.98);
noisePoints = noisePoints | extremeIsolation;
end

function lof = computeOptimizedLOF(X, k)

[N, ~] = size(X);
[knnIdx, knnDist] = knnsearch(X, X, 'K', k+1);
knnIdx = knnIdx(:, 2:end); 
kDist = knnDist(:, end); 

reachDist = max(knnDist(:, 2:end), kDist(knnIdx));

lrd = zeros(N, 1);
for i = 1:N
    lrd(i) = 1 / (mean(reachDist(i, :)) + eps);
end

lof = zeros(N, 1);
for i = 1:N
    neighbors = knnIdx(i, :);
    lof(i) = mean(lrd(neighbors)) / (lrd(i) + eps);
end
end

function optimalK = findOptimalClustersAdvanced(W, X_clean, maxClusters, minClusterSize)
%  stability analysis

N = size(W, 1);
if N < 2 * minClusterSize
    optimalK = 1;
    return;
end

%  normalized Laplacian
d = sum(W, 2);
d(d == 0) = eps; 
D_inv_sqrt = spdiags(1./sqrt(d), 0, N, N);
L_norm = speye(N) - D_inv_sqrt * W * D_inv_sqrt;

L_norm = L_norm + 1e-10 * speye(N);

try
    sigma = 1e-8; 
    [V, eigenvals] = eigs(L_norm, min(maxClusters + 3, N-1), sigma);
    eigenvals = real(diag(eigenvals));
    [eigenvals, idx] = sort(eigenvals);
    V = real(V(:, idx));
catch
    
    try
        [V, D_eig] = eig(full(L_norm));
        eigenvals = real(diag(D_eig));
        [eigenvals, idx] = sort(eigenvals);
        V = real(V(:, idx));
    catch
        optimalK = max(1, min(3, maxClusters));
        return;
    end
end

if length(eigenvals) < 3
    optimalK = 2;
    return;
end

eigengaps = diff(eigenvals(1:min(maxClusters, length(eigenvals)-1)));
[~, gapOptimalK] = max(eigengaps);
gapOptimalK = min(gapOptimalK + 1, maxClusters);

silScores = zeros(maxClusters, 1);
for testK = 2:min(maxClusters, 8) 
    if testK > size(V, 1), break; end
    try
        tempLabels = robustSpectralClustering(W, testK, V);
        if max(tempLabels) == testK && min(tempLabels) == 1
            
            X_embed = V(:, 1:min(testK, size(V, 2)));
            silScore = mean(silhouette(X_embed, tempLabels));
            silScores(testK) = max(0, silScore);
        end
    catch
        silScores(testK) = 0;
    end
end

[maxSilScore, silOptimalK] = max(silScores);


if maxSilScore > 0.4 %  silhouette threshold
    if silOptimalK < gapOptimalK
        optimalK = silOptimalK; % fewer clusters if silhouette is good
    else
        if silScores(min(gapOptimalK, length(silScores))) > 0.2
            optimalK = gapOptimalK;
        else
            optimalK = silOptimalK;
        end
    end
else
    optimalK = max(1, min(gapOptimalK, 4)); 
end

optimalK = max(1, min(optimalK, maxClusters));
end

function labels = spectralClusteringRobust(W, k)
%  spectral clustering

N = size(W, 1);
if k >= N || k < 1
    labels = ones(N, 1);
    return;
end

% normalized Laplacian
d = sum(W, 2);
d(d == 0) = eps;
D_inv_sqrt = spdiags(1./sqrt(d), 0, N, N);
L_norm = speye(N) - D_inv_sqrt * W * D_inv_sqrt;
L_norm = L_norm + 1e-10 * speye(N);

try
    sigma = 1e-8;
    [V, ~] = eigs(L_norm, k, sigma);
    V = real(V);
catch
    try
        [V, D_eig] = eig(full(L_norm));
        eigenvals = real(diag(D_eig));
        [~, idx] = sort(eigenvals);
        V = real(V(:, idx(1:k)));
    catch
        labels = ones(N, 1);
        return;
    end
end

V_norm = V ./ (sqrt(sum(V.^2, 2)) + eps);

labels = robustKmeans(V_norm, k);
end

function labels = robustSpectralClustering(W, k, V)
% silhouette validation

if k > size(V, 2)
    labels = ones(size(W, 1), 1);
    return;
end

V_k = V(:, 1:k);
V_norm = V_k ./ (sqrt(sum(V_k.^2, 2)) + eps);
labels = robustKmeans(V_norm, k);
end

function labels = robustKmeans(X, k)

bestLabels = [];
bestInertia = inf;

for trial = 1:3
    try
        [tempLabels, ~, sumd] = kmeans(X, k, 'MaxIter', 100, 'Replicates', 1, ...
                                      'Start', 'plus');
        inertia = sum(sumd);
        if inertia < bestInertia
            bestInertia = inertia;
            bestLabels = tempLabels;
        end
    catch
        
        [~, tempLabels] = max(abs(X), [], 2);
        if ~isempty(tempLabels)
            bestLabels = tempLabels;
        end
    end
end

if isempty(bestLabels)
    bestLabels = ones(size(X, 1), 1);
end

labels = bestLabels;
end

function labels = mergeSmallClusters(X, labels, minClusterSize)
% merge  small clusters

uniqueLabels = unique(labels);
clusterSizes = histcounts(labels, [uniqueLabels; max(uniqueLabels)+1]);

smallClusters = uniqueLabels(clusterSizes < minClusterSize);
largeClusters = uniqueLabels(clusterSizes >= minClusterSize);

if isempty(smallClusters) || isempty(largeClusters)
    return;
end

for i = 1:length(smallClusters)
    smallClusterIdx = find(labels == smallClusters(i));
    smallClusterCentroid = mean(X(smallClusterIdx, :), 1);
    
    minDist = inf;
    targetCluster = largeClusters(1);
    
    for j = 1:length(largeClusters)
        largeClusterIdx = find(labels == largeClusters(j));
        largeClusterCentroid = mean(X(largeClusterIdx, :), 1);
        
        dist = norm(smallClusterCentroid - largeClusterCentroid);
        if dist < minDist
            minDist = dist;
            targetCluster = largeClusters(j);
        end
    end
    
    % Merge
    labels(smallClusterIdx) = targetCluster;
end

uniqueLabels = unique(labels);
newLabels = labels;
for i = 1:length(uniqueLabels)
    newLabels(labels == uniqueLabels(i)) = i;
end
labels = newLabels;
end

function labels = mergeOverSplitClusters(X, labels)

uniqueLabels = unique(labels);
numClusters = length(uniqueLabels);

if numClusters <= 2
    return;
end

centroids = zeros(numClusters, size(X, 2));
clusterSizes = zeros(numClusters, 1);
intraClusterDist = zeros(numClusters, 1);
clusterDensities = zeros(numClusters, 1);

for i = 1:numClusters
    clusterPoints = X(labels == uniqueLabels(i), :);
    centroids(i, :) = mean(clusterPoints, 1);
    clusterSizes(i) = size(clusterPoints, 1);
    
    if size(clusterPoints, 1) > 1
        intraClusterDist(i) = mean(pdist(clusterPoints));
        clusterSpread = max(pdist(clusterPoints));
        clusterDensities(i) = size(clusterPoints, 1) / (clusterSpread + eps);
    else
        intraClusterDist(i) = 0;
        clusterDensities(i) = 1;
    end
end

dataSpread = range(X); 
globalScale = mean(dataSpread);
avgIntraDist = mean(intraClusterDist(intraClusterDist > 0));

centroidDist = pdist2(centroids, centroids);
separationScores = zeros(numClusters, numClusters);

% separation analysis
for i = 1:numClusters-1
    for j = i+1:numClusters
        cluster1_points = X(labels == uniqueLabels(i), :);
        cluster2_points = X(labels == uniqueLabels(j), :);
        
        separation = computeClusterSeparation(cluster1_points, cluster2_points, ...
                                            centroids(i,:), centroids(j,:), globalScale);
        separationScores(i, j) = separation;
        separationScores(j, i) = separation;
    end
end

baseMergeThreshold = avgIntraDist * 0.5; 
adaptiveFactor = min(2.0, globalScale / avgIntraDist); 
adaptiveThreshold = baseMergeThreshold * adaptiveFactor;

mergeCandidates = [];

for i = 1:numClusters-1
    for j = i+1:numClusters
        cluster1_idx = i;
        cluster2_idx = j;
        
        closeEnough = centroidDist(i, j) < adaptiveThreshold;
        
        separationOK = separationScores(i, j) < 0.3; % Threshold for oversplitting
        
        sizeRatio = min(clusterSizes(i), clusterSizes(j)) / max(clusterSizes(i), clusterSizes(j));
        sizeCompatible = sizeRatio > 0.3; % At least 30% size ratio
        
        densityRatio = min(clusterDensities(i), clusterDensities(j)) / ...
                      max(clusterDensities(i), clusterDensities(j));
        densityCompatible = densityRatio > 0.4; % At least 40% density ratio
        
        medianSize = median(clusterSizes);
        sizeBasedMerge = clusterSizes(i) <= medianSize || clusterSizes(j) <= medianSize;
        
        shouldMerge = closeEnough && separationOK && sizeCompatible && ...
                     densityCompatible && sizeBasedMerge;
        
        if shouldMerge
            mergeCandidates = [mergeCandidates; i, j];
        end
    end
end

%  merging
for k = 1:size(mergeCandidates, 1)
    cluster1 = uniqueLabels(mergeCandidates(k, 1));
    cluster2 = uniqueLabels(mergeCandidates(k, 2));
    
    if cluster1 ~= cluster2
        labels(labels == cluster2) = cluster1;
        uniqueLabels(uniqueLabels == cluster2) = cluster1;
    end
end

uniqueLabels = unique(labels);
newLabels = labels;
for i = 1:length(uniqueLabels)
    newLabels(labels == uniqueLabels(i)) = i;
end
labels = newLabels;
end

function separation = computeClusterSeparation(cluster1, cluster2, centroid1, centroid2, globalScale)
%  cluster separation score

centroidDist = norm(centroid1 - centroid2);
cluster1_spread = 0;
cluster2_spread = 0;

if size(cluster1, 1) > 1
    cluster1_spread = mean(pdist(cluster1));
end
if size(cluster2, 1) > 1  
    cluster2_spread = mean(pdist(cluster2));
end

avgSpread = (cluster1_spread + cluster2_spread) / 2;
centroidSeparation = centroidDist / (avgSpread + eps);

nSample = min(15, size(cluster1, 1));
if size(cluster1, 1) > nSample
    sample1_idx = randperm(size(cluster1, 1), nSample);
    sample1 = cluster1(sample1_idx, :);
else
    sample1 = cluster1;
end

nSample = min(15, size(cluster2, 1));
if size(cluster2, 1) > nSample
    sample2_idx = randperm(size(cluster2, 1), nSample);
    sample2 = cluster2(sample2_idx, :);
else
    sample2 = cluster2;
end

if ~isempty(sample1) && ~isempty(sample2)
    interClusterDist = pdist2(sample1, sample2);
    minInterDist = min(interClusterDist(:));
    relativeSeparation = minInterDist / (avgSpread + eps);
else
    relativeSeparation = 1;
end

overlapScore = computeClusterOverlap(cluster1, cluster2);

separation = 0.4 * centroidSeparation + 0.4 * relativeSeparation + 0.2 * (1 - overlapScore);
separation = max(0, separation); 
end

function overlap = computeClusterOverlap(cluster1, cluster2)

if size(cluster1, 1) < 2 || size(cluster2, 1) < 2
    overlap = 0;
    return;
end

min1 = min(cluster1, [], 1);
max1 = max(cluster1, [], 1);
min2 = min(cluster2, [], 1);
max2 = max(cluster2, [], 1);

dim_overlaps = zeros(1, size(cluster1, 2));
for d = 1:size(cluster1, 2)
    range1_start = min1(d);
    range1_end = max1(d);
    range2_start = min2(d);
    range2_end = max2(d);
    
    overlap_start = max(range1_start, range2_start);
    overlap_end = min(range1_end, range2_end);
    overlap_length = max(0, overlap_end - overlap_start);
    
    range1_length = range1_end - range1_start;
    range2_length = range2_end - range2_start;
    min_range_length = min(range1_length, range2_length);
    
    if min_range_length > eps
        dim_overlaps(d) = overlap_length / min_range_length;
    else
        dim_overlaps(d) = 0;
    end
end

overlap = mean(dim_overlaps);
end