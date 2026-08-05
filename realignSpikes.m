function [updatedLabels, realigned_spk_idx, realigned_waveform, changeType] = realignSpikes(labels, waveform, spk_idx, clusterRelabeling, cfg)

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

end
