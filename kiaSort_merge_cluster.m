function [newLabels, changeType, timeLagChanged] = kiaSort_merge_cluster(merge1, merge2, maxXcorrLag, initLabels, counts, isiViolation)
% Sequentially merge clusters using merge1 and merge2.

currLabels = initLabels;
n = length(initLabels);

changeType = zeros(n,1);
timeLagChanged = zeros(n,1);
currLabels = apply_merge_rule(merge1, maxXcorrLag, initLabels, counts, isiViolation, currLabels, true);

% Second merge using merge2 with relaxed connectivity merging.
newLabels = apply_merge_rule(merge2, maxXcorrLag, initLabels, counts, isiViolation, currLabels, false);

newLabels = propagate_labels(initLabels, newLabels);

% Post-process: compare final labels with original labels.
for i = 1:n
    if newLabels(i) ~= initLabels(i)
        changeType(i) = 1;
        timeLagChanged(i) = maxXcorrLag(i, (initLabels==newLabels(i)));
    else
        changeType(i) = 0;
        timeLagChanged(i) = 0;
    end
end
end

function newLabels = apply_merge_rule(mergeMat, maxXcorrLag, initLabels, counts, isiViolation, currLabels, useCliques)

newLabels = currLabels;

G = graph(mergeMat);
comps = conncomp(G);

for comp = unique(comps)
    compIdx = find(comps == comp);
    % Only consider clusters with valid initial labels.
    validComp = compIdx(initLabels(compIdx) ~= -1);
    if numel(validComp) < 2
        continue;
    end
    
    if useCliques

        subAdj = mergeMat(validComp, validComp);
        candidateCliques = findMaximalCliques(subAdj); % indices relative to validComp
        candidateGroups = {};
        candidateAvgViol = [];
        for i = 1:length(candidateCliques)
            clique = candidateCliques{i};
            if numel(clique) < 2
                continue;
            end
            groupGlobal = validComp(clique);
            candidateGroups{end+1} = groupGlobal;
            violMat = isiViolation(groupGlobal, groupGlobal);
            triuIdx = find(triu(ones(length(groupGlobal)),1));
            candidateAvgViol(end+1) = mean(violMat(triuIdx));
        end
        
        U = validComp;
        while true
            candIdx = [];
            candViol = [];
            for i = 1:length(candidateGroups)
                group = candidateGroups{i};
                if all(ismember(group, U))
                    candIdx(end+1) = i; 
                    candViol(end+1) = candidateAvgViol(i); 
                end
            end
            if isempty(candIdx)
                break;
            end
            [~, bestCandLocalIdx] = min(candViol);
            bestIdx = candIdx(bestCandLocalIdx);
            chosenGroup = candidateGroups{bestIdx};
            if numel(chosenGroup) >= 2
                % Select the representative as the one with the highest count.
                [~, repRel] = max(counts(chosenGroup));
                rep = chosenGroup(repRel);
                repLabel = currLabels(rep);  
                for j = 1:length(chosenGroup)
                    idx = chosenGroup(j);
                    if newLabels(idx) ~= repLabel
                        newLabels(idx) = repLabel;
                    end
                end
            end
            U = setdiff(U, chosenGroup);
            if numel(U) < 2
                break;
            end
        end
        
    else
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

function cliques = findMaximalCliques(adj)

n = size(adj,1);
cliques = {};
R = [];
P = 1:n;
X = [];
cliques = bronKerbosch(R, P, X, adj, cliques);
end

function cliques = bronKerbosch(R, P, X, adj, cliques)
if isempty(P) && isempty(X)
    cliques{end+1} = R;
else
    for v = P
        newR = [R, v];
        newP = intersect(P, find(adj(v,:)));
        newX = intersect(X, find(adj(v,:)));
        cliques = bronKerbosch(newR, newP, newX, adj, cliques);
        P = setdiff(P, v);
        X = union(X, v);
    end
end
end
