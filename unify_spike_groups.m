function [out, unified] = unify_spike_groups(sortedSamples, sampleFeatures, jitter_gap)

%       find fraction of i's spikes that map to each cluster on channel j.
%       and export statistics for that.
% Then unify IDs for all good clusters across channels.

N = numel(sortedSamples);
fractionsCell    = cell(N,1);
maxClusLabelCell = cell(N,1);
fraction2Cell    = cell(N,1);
numMatchCell     = cell(N,1);

for i = 1:N
    if isempty(sortedSamples{i}), continue; end
    
    updLabels_i = sampleFeatures{i}.updatedLabels;
    uLabs_i     = sortedSamples{i}.clusteringInfo.clusterRelabeling.newUniqueLabels;
    keepVec_i   = sortedSamples{i}.clusteringInfo.clusterSelection.keep;
    maxChIdx_i  = sortedSamples{i}.clusteringInfo.clusterSelection.max_channelIdx;
    spk_idx_i   = sampleFeatures{i}.spk_idx;
    
    badMask_i   = ~keepVec_i;
    badLabs_i   = uLabs_i(badMask_i);
    badMCh_i    = maxChIdx_i(badMask_i);
    
    fractionsCell{i}    = cell(sum(badMask_i),1); 
    maxClusLabelCell{i} = cell(sum(badMask_i),1);    
    numMatchCell{i}     = cell(sum(badMask_i),1);
    
    for bc = 1:numel(badLabs_i)
        thisLab = badLabs_i(bc);
        tgtChan = badMCh_i(bc);
        
        inClustMask = (updLabels_i == thisLab);
        spkInClst_i = spk_idx_i(inClustMask);
        nSpk       = numel(spkInClst_i);
        
        if nSpk==0 || tgtChan<1 || tgtChan> N || isempty(sortedSamples{tgtChan})
            fractionsCell{i}{bc}    = [];
            maxClusLabelCell{i}{bc} = [];
            fraction2Cell{i}{bc}    = [];
            numMatchCell{i}{bc}     = [];
            continue;
        end
        
        [fracVec_j, ~, labelMax_j, numMatchedLabel] = localComputeFraction( ...
            spkInClst_i, ...
            sampleFeatures{tgtChan}.spk_idx, ...
            sampleFeatures{tgtChan}.updatedLabels, ...
            sortedSamples{tgtChan}.clusteringInfo.clusterRelabeling.newUniqueLabels, ...
            jitter_gap);
        
        fractionsCell{i}{bc}    = fracVec_j;
        numMatchCell{i}{bc}     = numMatchedLabel;
        maxClusLabelCell{i}{bc} = labelMax_j;
             
    end
end

idCount = 0;
for i = 1:N
    if isempty(sortedSamples{i}), continue; end
    
    uLabs_i       = sortedSamples{i}.clusteringInfo.clusterRelabeling.newUniqueLabels;
    waveform_i    = sortedSamples{i}.clusteringInfo.clusterRelabeling.newMeanWaveforms;
    keepVec_i     = sortedSamples{i}.clusteringInfo.clusterSelection.keep;
    detectblity_i = sortedSamples{i}.clusteringInfo.clusterSelection.detectblity_val;
    mainNegativePolarity_i = sortedSamples{i}.clusteringInfo.clusterSelection.mainNegativePolarity;
    sideNegativePolarity_i = sortedSamples{i}.clusteringInfo.clusterSelection.sideNegativePolarity;

    idxKept    = find(keepVec_i==1 & uLabs_i~=-1);
    keptClus   = uLabs_i(idxKept);
    keptWaveforms = waveform_i(idxKept, :, :);
    keptDetect = detectblity_i(idxKept);
    mainNegPolDetect = mainNegativePolarity_i(idxKept);
    sideNegPolDetect = sideNegativePolarity_i(idxKept);

    for kc = 1:numel(keptClus)
        idCount = idCount + 1;
        unified.label(idCount,1)       = idCount;
        unified.channelID(idCount,1)   = i;
        unified.labelInChannel(idCount,1)= keptClus(kc);
        unified.meanWaveforms(idCount,:,:)= keptWaveforms(kc, :, :);
        unified.detectblity(idCount,1) = keptDetect(kc);
        unified.mainNegativePolarity(idCount,1) = mainNegPolDetect(kc);
        unified.sideNegativePolarity(idCount,1) = sideNegPolDetect(kc);
    end
end

out.matchedFraction  = fractionsCell;
out.bestMatchedLabel = maxClusLabelCell;
out.numMatchedLabels = numMatchCell;      

end

function [fracVec, idxMax, labelMax, numMatchedLabel] = localComputeFraction( ...
    spk_idx_i, spk_idx_j, updLabels_j, uniqueLabels_j, jitter_gap)

nI = numel(spk_idx_i);
nJ = numel(spk_idx_j);
nClusters_j = numel(uniqueLabels_j);

if nI==0 || nJ==0 || nClusters_j==0
    fracVec  = [];
    idxMax   = [];
    labelMax = [];
    numMatchedLabel = [];
    return;
end

hitsCount_j = zeros(nClusters_j,1);
for j = 1:nClusters_j
    thisjCluster = uniqueLabels_j(j);
    inClusterjMask = (updLabels_j == thisjCluster);
    sIdx = spk_idx_j(inClusterjMask);
    [d1, d2] = nearest_distances(spk_idx_i, sIdx);
    hitsCount_j(j) = sum(d1<jitter_gap | d2<jitter_gap);
end

fracVec = hitsCount_j / nI;
[maxVal, iMax] = max(fracVec);
if maxVal > 0
    idxMax   = iMax;
    labelMax = uniqueLabels_j(iMax);
else
    idxMax   = [];
    labelMax = [];
end
numMatchedLabel = hitsCount_j;
end