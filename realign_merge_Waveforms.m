function [clusterRelabeling] = realign_merge_Waveforms(meanClusterWaveform, clusterSampleCounts, clusterRelabeling)

[nClusters, nChannels, nSamples] = size(meanClusterWaveform);
meanClusterWaveform_shifted = zeros(size(meanClusterWaveform));

% --- Shift with zero-padding instead of circshift ---
for i = 1:nClusters
    lag = -clusterRelabeling.timeLagChanged(i);  % shift amount along dim 3
    wf = squeeze(meanClusterWaveform(i, :, :));   % nChannels x nSamples
    if lag > 0
        % shift right: pad zeros on the left
        meanClusterWaveform_shifted(i, :, lag+1:end) = wf(:, 1:nSamples-lag);
    elseif lag < 0
        % shift left: pad zeros on the right
        meanClusterWaveform_shifted(i, :, 1:nSamples+lag) = wf(:, -lag+1:end);
    else
        meanClusterWaveform_shifted(i, :, :) = wf;
    end
end

% --- Merge: use the main class waveform, not weighted average ---
newUniqueLabels = unique(clusterRelabeling.newLabels);
numUniqueLabels = length(newUniqueLabels);
newMeanWaveforms = zeros(numUniqueLabels, nChannels, nSamples);
newSampleCounts  = zeros(numUniqueLabels, 1);
newClusterSpikeDensity = zeros(numUniqueLabels, size(clusterRelabeling.clusterSpikeDensity, 2));

for iLabel = 1:numUniqueLabels
    groupMembers = find(clusterRelabeling.newLabels == newUniqueLabels(iLabel));
    newSampleCounts(iLabel) = sum(clusterSampleCounts(groupMembers));

    % Find the "main" cluster: the one whose original label matches this unique label
    mainIdx = find(clusterRelabeling.originalLabels == newUniqueLabels(iLabel));
    newMeanWaveforms(iLabel, :, :) = meanClusterWaveform_shifted(mainIdx, :, :);

    % Spike density: still weighted average across all group members
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