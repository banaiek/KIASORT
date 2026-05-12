function merge_matrix = kiaSort_merge_features(labels, pca, cfg)
% KIASORT_MERGE_FEATURES - Compute merge matrix based on PCA feature overlap
%
% OPTIMIZED VERSION with:
%   1. Pre-filtering: Skip pairs where centers are too far apart
%   2. Histogram overlap: Accurate, handles multi-modal distributions
%   3. Vectorized operations: Fast histogram computation
%
% INPUTS:
%   labels - Cluster labels [N x 1]
%   pca    - PCA structure with .score field [N x nComponents]
%   cfg    - Configuration with:
%            .nPCAcomp - Number of PCA components to use (default: all)
%            .T_merge  - Merge threshold (default: 0.05)
%            .nBins    - Number of histogram bins (default: 50)
%
% OUTPUT:
%   merge_matrix - Logical [nClusters x nClusters] merge recommendation

    %% Setup and defaults
    if ~isfield(cfg, 'nPCAcomp'), cfg.nPCAcomp = size(pca.score, 2); end
    if ~isfield(cfg, 'T_merge'),  cfg.T_merge = 0.05; end
    if ~isfield(cfg, 'nBins'),    cfg.nBins = 50; end
    
    nBins = cfg.nBins;
    
    % Get PCA scores
    nComp = min(cfg.nPCAcomp, size(pca.score, 2));
    X = pca.score(:, 1:nComp);
    
    unique_labels = unique(labels);
    nClusters = numel(unique_labels);
    
    % Handle edge cases
    if nClusters < 2
        merge_matrix = false(nClusters);
        return;
    end
    
    %% Precompute cluster statistics
    clusterData = cell(nClusters, 1);
    clusterCenters = zeros(nClusters, nComp);
    clusterStd = zeros(nClusters, 1);  % Overall spread for pre-filtering
    clusterCounts = zeros(nClusters, 1);
    validCluster = false(nClusters, 1);
    
    for i = 1:nClusters
        idx = (labels == unique_labels(i));
        data = X(idx, :);
        clusterData{i} = data;
        clusterCounts(i) = size(data, 1);
        
        if clusterCounts(i) >= 2
            clusterCenters(i, :) = mean(data, 1);
            % Use average std across dimensions as spread estimate
            clusterStd(i) = mean(std(data, 0, 1));
            validCluster(i) = true;
        end
    end
    
    %% Precompute pairwise center distances for pre-filtering
    % Skip pairs where distance > k * (std_i + std_j)
    preFilterK = 4;  % Number of std devs for pre-filter (conservative)
    
    centerDist = zeros(nClusters);
    combinedStd = zeros(nClusters);
    
    for i = 1:nClusters
        for j = i+1:nClusters
            centerDist(i,j) = norm(clusterCenters(i,:) - clusterCenters(j,:));
            centerDist(j,i) = centerDist(i,j);
            combinedStd(i,j) = clusterStd(i) + clusterStd(j);
            combinedStd(j,i) = combinedStd(i,j);
        end
    end
    
    % Pre-filter matrix: true if pair might overlap
    mightOverlap = centerDist < preFilterK * combinedStd;
    
    %% Compute pairwise overlaps (only for candidate pairs)
    overlap_matrix = zeros(nClusters);
    
    for i = 1:nClusters-1
        if ~validCluster(i)
            continue;
        end
        
        data_i = clusterData{i};
        center_i = clusterCenters(i, :);
        
        for j = i+1:nClusters
            if ~validCluster(j)
                continue;
            end
            
            % Pre-filter: skip if centers are too far apart
            if ~mightOverlap(i, j)
                overlap_matrix(i, j) = 0;
                overlap_matrix(j, i) = 0;
                continue;
            end
            
            data_j = clusterData{j};
            center_j = clusterCenters(j, :);
            
            % Direction between cluster centers
            diff_center = center_j - center_i;
            norm_diff = norm(diff_center);
            
            if norm_diff < eps
                % Coincident centers
                overlap_matrix(i, j) = 1;
                overlap_matrix(j, i) = 1;
                continue;
            end
            
            % Unit direction
            d = diff_center / norm_diff;
            
            % Project onto connecting line
            proj_i = data_i * d';
            proj_j = data_j * d';
            
            % Compute histogram-based overlap
            overlap = computeHistogramOverlap(proj_i, proj_j, nBins);
            
            overlap_matrix(i, j) = overlap;
            overlap_matrix(j, i) = overlap;
        end
    end
    
    %% Handle noise label
    noiseIdx = find(unique_labels == -1);
    if ~isempty(noiseIdx)
        overlap_matrix(noiseIdx, :) = NaN;
        overlap_matrix(:, noiseIdx) = NaN;
    end
    
    %% Adaptive thresholding
    valid_overlaps = overlap_matrix(:);
    valid_overlaps = valid_overlaps(~isnan(valid_overlaps) & valid_overlaps > 0);
    
    if isempty(valid_overlaps)
        merge_matrix = false(nClusters);
        return;
    end
    
    meanOverlap = mean(valid_overlaps);
    stdOverlap = std(valid_overlaps);
    threshold = meanOverlap + 1.68 * stdOverlap;
    
    merge_matrix = overlap_matrix > threshold;
    merge_matrix(1:nClusters+1:end) = false;
end

%% ========================================================================
%  HISTOGRAM-BASED OVERLAP COMPUTATION
%  ========================================================================

function overlap = computeHistogramOverlap(proj_i, proj_j, nBins)
% COMPUTEHISTOGRAMOVERLAP - Compute overlap using histogram intersection
%
% Histogram intersection = sum(min(pdf_i, pdf_j))
% Returns 1 for identical distributions, 0 for non-overlapping.
% Accurately handles multi-modal distributions.

    n_i = numel(proj_i);
    n_j = numel(proj_j);
    
    % Combined range
    minVal = min(min(proj_i), min(proj_j));
    maxVal = max(max(proj_i), max(proj_j));
    
    % Handle constant data
    range = maxVal - minVal;
    if range < eps
        overlap = 1;
        return;
    end
    
    % Bin edges with small padding to avoid edge effects
    padding = range * 0.01;
    edges = linspace(minVal - padding, maxVal + padding, nBins + 1);
    binWidth = edges(2) - edges(1);
    
    % Compute normalized histograms (probability density)
    counts_i = histcounts(proj_i, edges);
    counts_j = histcounts(proj_j, edges);
    
    % Convert to density (integral = 1)
    pdf_i = counts_i / (n_i * binWidth);
    pdf_j = counts_j / (n_j * binWidth);
    
    % Histogram intersection (integral of min)
    overlap = sum(min(pdf_i, pdf_j)) * binWidth;
    
    % Clamp to [0, 1]
    overlap = min(1, max(0, overlap));
end