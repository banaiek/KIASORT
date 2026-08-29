function [unif, lbl_all, nRemoved, remap] = kiaSort_compact_unit_table(unif, lbl_all)
%KIASORT_COMPACT_UNIT_TABLE  Drop dead unit rows and restore label == row index.
%
%   [unif, lbl_all, nRemoved, remap] = kiaSort_compact_unit_table(unif, lbl_all)
%
%   unif    crossChannelStats.unified_labels
%   lbl_all the /unifiedLabels column, same length as the other H5 outputs
%
%   Returns the table with every row whose label carries no spike removed,
%   the labels renumbered 1..N, and lbl_all remapped to match. remap is an
%   old-label -> new-label containers.Map. Both outputs must be written
%   together or the pair is left disagreeing.
% Drop unit-table rows whose label no longer carries any spike, then
% renumber so label == row index again.
%
% Phase 1 and Phase 2 reassign spikes away from absorbed and dropped units
% but leave their rows behind. unify_spike_groups establishes label == row
% index and the curation GUI depends on it -- numGroups is numel(label) and
% plotCCG / plotISI / plotDensity index groupList and the per-pair caches by
% label -- so a stale row is both an empty unit in the list and an offset
% into every cache after it. Renumbering means the spike labels have to be
% remapped in the same breath, which is why this runs on lbl_all too.
remap    = containers.Map('KeyType','double','ValueType','double');
nRemoved = 0;
if ~isfield(unif, 'label') || isempty(unif.label)
    return;
end

oldLab = unif.label(:);
nU     = numel(oldLab);
alive  = ismember(oldLab, unique(lbl_all(lbl_all > 0)));
if all(alive) || ~any(alive)
    % Nothing stale, or nothing survived -- in the second case the run has
    % bigger problems than the table and dropping every row would hide them.
    return;
end

% Per-unit field: a vector of length nU, or an array whose FIRST dimension
% indexes units (meanWaveforms is nUnits x nChan x nSample, or nUnits x nSample).
fn = fieldnames(unif);
for i = 1:numel(fn)
    v = unif.(fn{i});
    if isvector(v) && numel(v) == nU
        v = v(:);
        unif.(fn{i}) = v(alive);
    elseif ~isvector(v) && size(v, 1) == nU
        if ndims(v) == 3
            unif.(fn{i}) = v(alive, :, :);
        else
            unif.(fn{i}) = v(alive, :);
        end
    end
end

keptLab  = oldLab(alive);
nKept    = numel(keptLab);
nRemoved = nU - nKept;

lut = nan(max([keptLab(:); lbl_all(:); 1]), 1);
lut(keptLab) = (1:nKept)';
sel  = lbl_all > 0 & lbl_all <= numel(lut);
m    = lut(lbl_all(sel));
good = ~isnan(m);
idx  = find(sel);
lbl_all(idx(good)) = m(good);

unif.label = (1:nKept)';
for i = 1:nKept
    remap(keptLab(i)) = i;
end
end
