function [clusterRelabeling] = realign_merge_Waveforms(meanClusterWaveform, meanClusterWaveformFull, clusterSampleCounts, clusterRelabeling)

if nargin < 4
    clusterRelabeling      = clusterSampleCounts;
    clusterSampleCounts    = meanClusterWaveformFull;
    meanClusterWaveformFull = [];
end

newUniqueLabels = unique(clusterRelabeling.newLabels);
numUniqueLabels = length(newUniqueLabels);

newMeanWaveforms = shiftAndMerge(meanClusterWaveform, clusterRelabeling, newUniqueLabels);
if ~isempty(meanClusterWaveformFull)
    clusterRelabeling.newMeanWaveformsFull = ...
        shiftAndMerge(meanClusterWaveformFull, clusterRelabeling, newUniqueLabels);
end

newSampleCounts        = zeros(numUniqueLabels, 1);
newClusterSpikeDensity = zeros(numUniqueLabels, size(clusterRelabeling.clusterSpikeDensity, 2));

for iLabel = 1:numUniqueLabels
    groupMembers = find(clusterRelabeling.newLabels == newUniqueLabels(iLabel));
    newSampleCounts(iLabel) = sum(clusterSampleCounts(groupMembers));

    for j = 1:length(groupMembers)
        newClusterSpikeDensity(iLabel, :) = newClusterSpikeDensity(iLabel, :) + ...
            clusterRelabeling.clusterSpikeDensity(groupMembers(j), :) .* ...
            clusterSampleCounts(groupMembers(j)) ./ newSampleCounts(iLabel);
    end
end

clusterRelabeling.newUniqueLabels        = newUniqueLabels;
clusterRelabeling.newMeanWaveforms       = newMeanWaveforms;
clusterRelabeling.newSampleCounts        = newSampleCounts;
clusterRelabeling.newClusterSpikeDensity = newClusterSpikeDensity;
end


function merged = shiftAndMerge(mw, clusterRelabeling, newUniqueLabels)
% timeLagChanged is in samples, so the same shift applies to either window
% length; reshape rather than squeeze so a single-channel footprint keeps
% its orientation.
[nClusters, nChannels, nSamples] = size(mw);
shifted = zeros(size(mw));

for i = 1:nClusters
    lag = -clusterRelabeling.timeLagChanged(i);
    wf  = reshape(mw(i, :, :), nChannels, nSamples);
    if lag > 0
        shifted(i, :, lag+1:end) = wf(:, 1:nSamples-lag);
    elseif lag < 0
        shifted(i, :, 1:nSamples+lag) = wf(:, -lag+1:end);
    else
        shifted(i, :, :) = wf;
    end
end

merged = zeros(numel(newUniqueLabels), nChannels, nSamples);
for iLabel = 1:numel(newUniqueLabels)
    mainIdx = find(clusterRelabeling.originalLabels == newUniqueLabels(iLabel));
    merged(iLabel, :, :) = shifted(mainIdx, :, :);
end
end
