function [updatedLabels, realigned_spk_idx, realigned_waveform, changeType] = realignSpikes(labels, waveform, spk_idx, clusterRelabeling, cfg)

fs = cfg.samplingFrequency;
midPoint = floor(cfg.spikeDuration * fs/(2*1000))+1;
spike_length = floor(cfg.clusteringSpikeDuration * fs/(2*1000));
realigned_waveform = waveform;

uniqueLabels = clusterRelabeling.originalLabels(:);
numLabels    = length(uniqueLabels);
newLabels    = clusterRelabeling.newLabels(:);
timeLags     = clusterRelabeling.timeLagChanged(:);
changeID     = clusterRelabeling.changeType(:);
% Find the index of each label in uniqueLabels

updatedLabels      = labels;
realigned_spk_idx  = spk_idx;
realigned_waveform = waveform;
changeType         = zeros(size(labels));

[isMapped, loc] = ismember(labels, uniqueLabels);
mappedIdx = find(isMapped);

if ~isempty(mappedIdx)
    updatedLabels(mappedIdx)     = newLabels(loc(mappedIdx));
    realigned_spk_idx(mappedIdx) = spk_idx(mappedIdx) + timeLags(loc(mappedIdx));
    changeType(mappedIdx)        = changeID(loc(mappedIdx));
end


for i = 1:numLabels
    idx = find(loc == i);
    lag = - timeLags(i);
    if ~isempty(idx)
        realigned_waveform(idx, :, :) = circshift(waveform(idx, :, :), lag, 3);
    end
end

realigned_waveform = realigned_waveform(:,:,midPoint-spike_length:midPoint+spike_length);

end
