function [newLabels, changeType, timeLagChanged] = kiaSort_merge_cluster(merge1, merge2, maxXcorrLag, initLabels, counts, isiViolation)
% KIASORT_MERGE_CLUSTER - Sequentially merge clusters using two merge matrices
%
% This function applies hierarchical cluster merging using clique-based and
% connectivity-based merging strategies.
%
% INPUTS:
%   merge1       - First merge matrix (strict, clique-based)
%   merge2       - Second merge matrix (relaxed, connectivity-based)
%   maxXcorrLag  - Maximum cross-correlation lag matrix [N x N]
%   initLabels   - Initial cluster labels [N x 1]
%   counts       - Sample counts per cluster [N x 1]
%   isiViolation - ISI violation matrix [N x N]
%
% OUTPUTS:
%   newLabels      - Updated cluster labels
%   changeType     - Type of change for each cluster (0=none, 1=merged)
%   timeLagChanged - Time lag applied to merged clusters

n = length(initLabels);
currLabels = initLabels;

changeType = zeros(n, 1);
timeLagChanged = zeros(n, 1);

% First pass: strict clique-based merging
currLabels = applyMergeRule(merge1, maxXcorrLag, initLabels, counts, isiViolation, currLabels, true);

% Second pass: relaxed connectivity-based merging
newLabels = applyMergeRule(merge2, maxXcorrLag, initLabels, counts, isiViolation, currLabels, false);

% Propagate labels to ensure consistency
newLabels = propagate_labels(initLabels, newLabels);

% Record changes
for i = 1:n
    if newLabels(i) ~= initLabels(i)
        changeType(i) = 1;
        % Find matching original label
        matchIdx = find(initLabels == newLabels(i));
        if ~isempty(matchIdx)
            timeLagChanged(i) = maxXcorrLag(i, matchIdx(1));
        else
            timeLagChanged(i) = 0;
        end
    else
        changeType(i) = 0;
        timeLagChanged(i) = 0;
    end
end

end

%% ========================================================================
%  MERGE RULE APPLICATION
%  ========================================================================

function newLabels = applyMergeRule(mergeMat, maxXcorrLag, initLabels, counts, isiViolation, currLabels, useCliques)
% Apply merging rules to clusters

newLabels = currLabels;

% Build graph from merge matrix
G = graph(mergeMat);
comps = conncomp(G);

uniqueComps = unique(comps);

for comp = uniqueComps
    compIdx = find(comps == comp);
    
    % Only consider clusters with valid initial labels
    validComp = compIdx(initLabels(compIdx) ~= -1);
    if numel(validComp) < 2
        continue;
    end
    
    if useCliques
        % Clique-based merging: find maximal cliques and merge greedily
        newLabels = mergeUsingCliques(mergeMat, validComp, isiViolation, counts, currLabels, newLabels);
    else
        % Simple connectivity-based merging: merge entire component
        [~, repRel] = max(counts(validComp));
        rep = validComp(repRel);
        repLabel = currLabels(rep);
        
        for j = 1:length(validComp)
            idx = validComp(j);
            if newLabels(idx) ~= repLabel
                newLabels(idx) = repLabel;
            end
        end
    end
end

end

%% ========================================================================
%  CLIQUE-BASED MERGING
%  ========================================================================

function newLabels = mergeUsingCliques(mergeMat, validComp, isiViolation, counts, currLabels, newLabels)
% Merge clusters using maximal cliques with ISI violation optimization

% Extract subgraph adjacency
subAdj = mergeMat(validComp, validComp);

% Find all maximal cliques
candidateCliques = findMaximalCliques(subAdj);

% Compute average ISI violation for each clique
candidateGroups = {};
candidateAvgViol = [];

for i = 1:length(candidateCliques)
    clique = candidateCliques{i};
    if numel(clique) < 2
        continue;
    end
    
    groupGlobal = validComp(clique);
    candidateGroups{end+1} = groupGlobal; %#ok<AGROW>
    
    % Compute average ISI violation for this clique
    violMat = isiViolation(groupGlobal, groupGlobal);
    triuIdx = triu(ones(length(groupGlobal)), 1) == 1;
    candidateAvgViol(end+1) = mean(violMat(triuIdx)); %#ok<AGROW>
end

if isempty(candidateGroups)
    return;
end

% Greedy clique selection (prefer low ISI violation)
U = validComp;  % Unassigned clusters

while ~isempty(U) && numel(U) >= 2
    % Find candidate cliques that are subsets of unassigned clusters
    validCandIdx = [];
    validCandViol = [];
    
    for i = 1:length(candidateGroups)
        if all(ismember(candidateGroups{i}, U))
            validCandIdx(end+1) = i; %#ok<AGROW>
            validCandViol(end+1) = candidateAvgViol(i); %#ok<AGROW>
        end
    end
    
    if isempty(validCandIdx)
        break;
    end
    
    % Select clique with lowest ISI violation
    [~, bestLocalIdx] = min(validCandViol);
    bestIdx = validCandIdx(bestLocalIdx);
    chosenGroup = candidateGroups{bestIdx};
    
    if numel(chosenGroup) >= 2
        % Select representative (highest count)
        [~, repRel] = max(counts(chosenGroup));
        rep = chosenGroup(repRel);
        repLabel = currLabels(rep);
        
        % Merge all clusters in clique to representative
        for j = 1:length(chosenGroup)
            idx = chosenGroup(j);
            if newLabels(idx) ~= repLabel
                newLabels(idx) = repLabel;
            end
        end
    end
    
    % Remove merged clusters from unassigned set
    U = setdiff(U, chosenGroup);
end

end

%% ========================================================================
%  MAXIMAL CLIQUE FINDING (Bron-Kerbosch Algorithm)
%  ========================================================================

function cliques = findMaximalCliques(adj)
% Find all maximal cliques using Bron-Kerbosch algorithm

n = size(adj, 1);
cliques = {};
R = [];
P = 1:n;
X = [];

cliques = bronKerbosch(R, P, X, adj, cliques);
end

function cliques = bronKerbosch(R, P, X, adj, cliques)
% Bron-Kerbosch algorithm with pivoting

if isempty(P) && isempty(X)
    % Found a maximal clique
    cliques{end+1} = R;
    return;
end

% Choose pivot (vertex with most neighbors in P ∪ X)
PuX = union(P, X);
if isempty(PuX)
    return;
end

neighborCounts = sum(adj(PuX, P), 2);
[~, pivotLocalIdx] = max(neighborCounts);
pivot = PuX(pivotLocalIdx);

% Iterate over P \ N(pivot)
pivotNeighbors = find(adj(pivot, :));
toProcess = setdiff(P, pivotNeighbors);

for v = toProcess
    newR = [R, v];
    neighbors_v = find(adj(v, :));
    newP = intersect(P, neighbors_v);
    newX = intersect(X, neighbors_v);
    
    cliques = bronKerbosch(newR, newP, newX, adj, cliques);
    
    P = setdiff(P, v);
    X = union(X, v);
end

end