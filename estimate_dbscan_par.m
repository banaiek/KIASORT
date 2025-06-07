function [optimalEps, k] = estimate_dbscan_par(data, k)

    if nargin < 2
        k = size(data,2) + 1;
    end
    
    [~, distances] = knnsearch(data, data, 'K', k+1, 'distance','minkowski','p',1);
    
    kDistances = distances(:, end);
    
    sortedKDistances = sort(kDistances);
    n = length(sortedKDistances);
    
    x = (1:n)';
    y = sortedKDistances;
    
    lineVec = [n-1, y(end) - y(1)];
    lineVecNorm = norm(lineVec);
    
    distancesToLine = abs((x - 1) * lineVec(2) - (y - y(1)) * lineVec(1)) / lineVecNorm;
    
    [~, kneeIdx] = max(distancesToLine);
    optimalEps = sortedKDistances(kneeIdx);
end