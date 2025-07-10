function [clusterRelabeling] = realign_merge_Waveforms(meanClusterWaveform, clusterSampleCounts, clusterRelabeling)


meanClusterWaveform_shifted = zeros(size(meanClusterWaveform));

for i = 1 :length(clusterSampleCounts)
    meanClusterWaveform_shifted(i,:,:) = circshift(meanClusterWaveform(i,:,:), -clusterRelabeling.timeLagChanged(i),3);
end

newUniqueLabels = unique(clusterRelabeling.newLabels);
numUniqueLabels = length(newUniqueLabels);
newMeanWaveforms = zeros(numUniqueLabels, size(meanClusterWaveform,2), size(meanClusterWaveform,3));
newSampleCounts = zeros(numUniqueLabels,1);
newClusterSpikeDensity = zeros(numUniqueLabels, size(clusterRelabeling.clusterSpikeDensity,2));

for iLabel = 1: numUniqueLabels
    groupMembers = find(clusterRelabeling.newLabels == newUniqueLabels(iLabel));
    newSampleCounts(iLabel) = sum(clusterSampleCounts(groupMembers));
    % clusterID = find(clusterRelabeling.originalLabels == newUniqueLabels(iLabel));
    %  newMeanWaveforms (iLabel, :, :) = (meanClusterWaveform_shifted(clusterID, : ,:));
    for j = 1 : length(groupMembers)
    newMeanWaveforms (iLabel, :, :) = newMeanWaveforms (iLabel, :, :) + ...
        (meanClusterWaveform_shifted(groupMembers(j), : ,:) .* clusterSampleCounts(groupMembers(j))./newSampleCounts(iLabel));  
    newClusterSpikeDensity(iLabel, :) = clusterRelabeling.clusterSpikeDensity(iLabel, :) + ...
        (clusterRelabeling.clusterSpikeDensity(groupMembers(j), :) .* clusterSampleCounts(groupMembers(j))./newSampleCounts(iLabel));  
    end
end

clusterRelabeling.newUniqueLabels      = newUniqueLabels;
clusterRelabeling.newMeanWaveforms     = newMeanWaveforms;
clusterRelabeling.newSampleCounts      = newSampleCounts;
clusterRelabeling.newClusterSpikeDensity    = newClusterSpikeDensity;
end

