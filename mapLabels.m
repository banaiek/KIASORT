function globalLabels = mapLabels(labels, i, unified_labeling)
% mapLocalToGlobal maps local labels on channel i to global labels.
%
% Inputs:
%   labels        - local labels detected on channel i
%   i             - Channel number
%   label         - global labels
%   channelID     - channel for each global label
%   labelInChannel- local labels corresponding to global labels
if ~isempty(labels)
    idx = (unified_labeling.channelID == i);
    localLabels = unified_labeling.labelInChannel(idx);
    channelGlobalLabels = unified_labeling.label(idx);
    [~, loc] = ismember(labels, localLabels);
    globalLabels = channelGlobalLabels(loc);
else
    globalLabels = [];
end

end