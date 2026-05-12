function finalLabels = noise_exclude_dbscan(dataAll, ampVals, thr, epsilon, numPt, useSpecialLabels)
% NOISE_EXCLUDE_DBSCAN - DBSCAN clustering with threshold-based noise identification
%
% This function performs DBSCAN clustering and identifies clusters that are
% likely noise artifacts (spikes near the detection threshold).
%
% INPUTS:
%   dataAll         - Feature matrix [N x D]
%   ampVals         - Amplitude values for each point [N x 1]
%   thr             - Detection threshold value
%   epsilon         - DBSCAN epsilon parameter
%   numPt           - DBSCAN minimum points parameter
%   useSpecialLabels - If true, mark noise clusters with special labels (-2, -3)
%
% OUTPUT:
%   finalLabels - Cluster labels (-1 = noise, -2 = near negative threshold,
%                 -3 = near positive threshold, positive integers = clusters)

if nargin < 6
    useSpecialLabels = true;
end

%% ========================================================================
%  INITIAL DBSCAN CLUSTERING
%  ========================================================================
labels = dbscan(dataAll, epsilon, numPt, 'Distance', 'minkowski', 'P', 1);

%% ========================================================================
%  IDENTIFY NEAR-THRESHOLD CLUSTERS
%  ========================================================================
uniqueLabels = unique(labels(labels > 0));
nearThresholdClusters = [];
clusterAmps = cell(length(uniqueLabels), 1);
clusterPolarities = cell(length(uniqueLabels), 1);

for i = 1:length(uniqueLabels)
    clusterMask = labels == uniqueLabels(i);
    clusterAmps{i} = ampVals(clusterMask);
    
    meanAmp = mean(clusterAmps{i});
    stdAmp = std(clusterAmps{i});
    ampRange = range(clusterAmps{i});
    
    % Distance to negative and positive thresholds (relative to threshold)
    distToNegThr = abs(meanAmp + thr) / thr;
    distToPosThr = abs(meanAmp - thr) / thr;
    minDistToThr = min(distToNegThr, distToPosThr);
    
    % Criteria for near-threshold cluster:
    % 1. Mean amplitude is close to threshold (within 30%)
    % 2. Amplitude range is small (< 40% of threshold)
    isNearThreshold = minDistToThr < 0.3 && ampRange < 0.4 * thr;
    
    if isNearThreshold
        % Additional check: most points should be near threshold
        nearNegThr = sum(abs(clusterAmps{i} + thr) < 0.2 * thr);
        nearPosThr = sum(abs(clusterAmps{i} - thr) < 0.2 * thr);
        fractionNearThr = (nearNegThr + nearPosThr) / length(clusterAmps{i});
        
        if fractionNearThr > 0.7
            nearThresholdClusters(end+1) = uniqueLabels(i); %#ok<AGROW>
            
            % Determine polarity
            if meanAmp < 0
                clusterPolarities{length(nearThresholdClusters)} = 'neg';
            else
                clusterPolarities{length(nearThresholdClusters)} = 'pos';
            end
        end
    end
end

%% ========================================================================
%  ASSIGN NOISE POINTS TO NEAR-THRESHOLD CLUSTERS
%  ========================================================================
finalLabels = labels;

if ~isempty(nearThresholdClusters)
    % Compute cluster prototypes (centroids)
    nThresholdClusters = length(nearThresholdClusters);
    prototypes = zeros(nThresholdClusters, size(dataAll, 2));
    protoTypes = cell(nThresholdClusters, 1);
    
    for i = 1:nThresholdClusters
        clusterMask = labels == nearThresholdClusters(i);
        prototypes(i, :) = mean(dataAll(clusterMask, :), 1);
        
        meanAmp = mean(ampVals(clusterMask));
        if meanAmp < 0
            protoTypes{i} = 'neg';
        else
            protoTypes{i} = 'pos';
        end
    end
    
    % Find unclustered points near threshold
    noisePoints = find(labels == -1);
    
    for i = 1:length(noisePoints)
        pt = dataAll(noisePoints(i), :);
        amp = ampVals(noisePoints(i));
        
        % Check if amplitude is in near-threshold range
        nearNegThr = abs(amp + thr) < 0.3 * thr;
        nearPosThr = abs(amp - thr) < 0.3 * thr;
        
        if nearNegThr || nearPosThr
            % Find nearest prototype with matching polarity
            distances = sum(abs(prototypes - pt), 2);  % L1 distance
            
            % Prefer prototypes with matching polarity
            for j = 1:nThresholdClusters
                if (nearNegThr && strcmp(protoTypes{j}, 'pos')) || ...
                   (nearPosThr && strcmp(protoTypes{j}, 'neg'))
                    distances(j) = distances(j) * 2;  % Penalize polarity mismatch
                end
            end
            
            [minDist, nearest] = min(distances);
            
            % Assign if close enough
            if minDist < epsilon * 1.5
                finalLabels(noisePoints(i)) = nearThresholdClusters(nearest);
            end
        end
    end
    
    % Apply special labels for near-threshold clusters
    if useSpecialLabels
        for i = 1:nThresholdClusters
            clusterMask = finalLabels == nearThresholdClusters(i);
            if strcmp(protoTypes{i}, 'neg')
                finalLabels(clusterMask) = -2;  % Near negative threshold
            else
                finalLabels(clusterMask) = -3;  % Near positive threshold
            end
        end
    end
end

%% ========================================================================
%  EXTEND MATCHING: MERGE SIMILAR CLUSTERS ACROSS POLARITIES
%  ========================================================================
finalLabels = extendMatching(dataAll, finalLabels, ampVals, epsilon);

end

%% ========================================================================
%  HELPER FUNCTION: EXTEND MATCHING
%  ========================================================================

function labels = extendMatching(dataAll, labels, ampVals, epsilon)
% Find and merge similar clusters that span the negative/positive divide
%
% This handles cases where a single unit's spikes are split into two clusters
% based on amplitude polarity.

uniqueLabels = unique(labels(labels > 0));
nClusters = length(uniqueLabels);

if nClusters < 2
    return;
end

% Compute cluster statistics
centroids = zeros(nClusters, size(dataAll, 2));
ampMeans = zeros(nClusters, 1);
clusterSizes = zeros(nClusters, 1);

for i = 1:nClusters
    mask = labels == uniqueLabels(i);
    centroids(i, :) = mean(dataAll(mask, :), 1);
    ampMeans(i) = mean(ampVals(mask));
    clusterSizes(i) = sum(mask);
end

% Track which clusters have been merged
merged = false(nClusters, 1);

for i = 1:nClusters
    if merged(i)
        continue;
    end
    
    % Find potential merge candidates (opposite polarity)
    if ampMeans(i) < 0
        candidates = find(ampMeans > 0 & ~merged);
    else
        candidates = find(ampMeans < 0 & ~merged);
    end
    
    if isempty(candidates)
        continue;
    end
    
    % Compute distances to candidates
    distances = sum(abs(centroids(candidates, :) - centroids(i, :)), 2);
    
    % Find best match
    [minDist, idx] = min(distances);
    
    % Merge threshold: consider feature space dimensionality
    mergeThreshold = epsilon * size(dataAll, 2) * 0.5;
    
    if minDist < mergeThreshold
        % Merge clusters - keep the one with more samples
        candidateIdx = candidates(idx);
        
        if clusterSizes(i) >= clusterSizes(candidateIdx)
            labels(labels == uniqueLabels(candidateIdx)) = uniqueLabels(i);
        else
            labels(labels == uniqueLabels(i)) = uniqueLabels(candidateIdx);
        end
        
        merged(candidateIdx) = true;
    end
end
end