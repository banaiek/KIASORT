function finalLabels = noise_exclude_dbscan(dataAll, ampVals, thr, epsilon, numPt, useSpecialLabels)

if nargin < 6
        useSpecialLabels = true;
    end
    
    labels = dbscan(dataAll, epsilon, numPt, 'Distance', 'minkowski', 'P', 1);
    
    uniqueLabels = unique(labels(labels > 0));
    nearThresholdClusters = [];
    clusterAmps = cell(length(uniqueLabels), 1);
    
    for i = 1:length(uniqueLabels)
        clusterMask = labels == uniqueLabels(i);
        clusterAmps{i} = ampVals(clusterMask);
        
        meanAmp = mean(clusterAmps{i});
        stdAmp = std(clusterAmps{i});
        

        relDist = min(abs(meanAmp + thr), abs(meanAmp - thr)) / thr;
        ampRange = range(clusterAmps{i});
        
        if relDist < 0.3 && ampRange < 0.4 * thr
            % Additional check: percentage of points near threshold
            nearNegThr = sum(abs(clusterAmps{i} + thr) < 0.2 * thr);
            nearPosThr = sum(abs(clusterAmps{i} - thr) < 0.2 * thr);
            
            if (nearNegThr + nearPosThr) / length(clusterAmps{i}) > 0.7
                nearThresholdClusters(end+1) = uniqueLabels(i);
            end
        end
    end
    
    finalLabels = labels;
    
    if ~isempty(nearThresholdClusters)
        % Build KD-tree for efficient nearest neighbor search
        kdtree = KDTreeSearcher(dataAll, 'Distance', 'cityblock');
        
        % Compute cluster prototypes
        prototypes = zeros(length(nearThresholdClusters), size(dataAll, 2));
        protoTypes = cell(length(nearThresholdClusters), 1);
        
        for i = 1:length(nearThresholdClusters)
            clusterMask = labels == nearThresholdClusters(i);
            prototypes(i, :) = mean(dataAll(clusterMask, :), 1);
            
            meanAmp = mean(ampVals(clusterMask));
            if meanAmp < 0
                protoTypes{i} = 'neg';
            else
                protoTypes{i} = 'pos';
            end
        end
        
        % Find unclustered points that match near-threshold characteristics
        noisePoints = find(labels == -1);
        
        for i = 1:length(noisePoints)
            pt = dataAll(noisePoints(i), :);
            amp = ampVals(noisePoints(i));
            
            % Check if amplitude is in near-threshold range
            if (abs(amp + thr) < 0.3 * thr) || (abs(amp - thr) < 0.3 * thr)
                % Find nearest prototype
                distances = sum(abs(prototypes - pt), 2);
                [minDist, nearest] = min(distances);
                
                % Assign if close enough
                if minDist < epsilon * 1.5
                    finalLabels(noisePoints(i)) = nearThresholdClusters(nearest);
                end
            end
        end
        
        if useSpecialLabels
            for i = 1:length(nearThresholdClusters)
                clusterMask = finalLabels == nearThresholdClusters(i);
                if strcmp(protoTypes{i}, 'neg')
                    finalLabels(clusterMask) = -2;
                else
                    finalLabels(clusterMask) = -3;
                end
            end
        end
    end
    
    finalLabels = extendMatching(dataAll, finalLabels, ampVals, epsilon);
end

function labels = extendMatching(dataAll, labels, ampVals, epsilon)
    % Find and merge similar clusters across negative/positive divide
    uniqueLabels = unique(labels(labels > 0));
    
    % Compute all cluster centroids
    centroids = zeros(length(uniqueLabels), size(dataAll, 2));
    ampMeans = zeros(length(uniqueLabels), 1);
    
    for i = 1:length(uniqueLabels)
        mask = labels == uniqueLabels(i);
        centroids(i, :) = mean(dataAll(mask, :), 1);
        ampMeans(i) = mean(ampVals(mask));
    end
    
    % Find matching pairs
    merged = false(length(uniqueLabels), 1);
    
    for i = 1:length(uniqueLabels)
        if merged(i), continue; end
        
        % Only match if one is negative and one is positive
        if ampMeans(i) < 0
            candidates = find(ampMeans > 0 & ~merged);
        else
            candidates = find(ampMeans < 0 & ~merged);
        end
        
        if ~isempty(candidates)
            distances = sum(abs(centroids(candidates, :) - centroids(i, :)), 2);
            [minDist, idx] = min(distances);
            
            if minDist < epsilon * size(dataAll, 2) * 0.5
                % Merge clusters
                labels(labels == uniqueLabels(candidates(idx))) = uniqueLabels(i);
                merged(candidates(idx)) = true;
            end
        end
    end
end