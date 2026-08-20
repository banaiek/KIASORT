function report = kiaSort_repair_label_index(outputPath, verbose)
%KIASORT_REPAIR_LABEL_INDEX  Restore label == row index in an existing sort.
%
%   report = kiaSort_repair_label_index(outputPath)
%
%   unify_spike_groups assigns unified_labels.label(i) = i, and the curation
%   GUI depends on it: plotCCG / plotISI / plotDensity take a unit's LABEL and
%   use it to index groupList, mergedFlag and the per-pair caches. A unit
%   inserted mid-list breaks the invariant, which shows up as panels drawing a
%   unit in another unit's colour.
%
%   Relabels every unit to its row position and remaps the spike labels to
%   match. No spike changes unit, no unit changes row -- only the numbering.
%   Safe to run repeatedly; a tree that already satisfies the invariant is
%   left untouched.

if nargin < 2, verbose = true; end
report = struct('nUnits', 0, 'nRelabelled', 0, 'changed', false, 'ok', false);

outputPath = char(outputPath);
uH5  = fullfile(outputPath, 'RES_Sorted',     'unifiedLabels.h5');
ssP  = fullfile(outputPath, 'Sorted_Samples', 'sorted_samples.mat');
if ~exist(uH5, 'file') || ~exist(ssP, 'file')
    if verbose, fprintf('Repair: required files missing.\n'); end
    return;
end

ssData = load(ssP, 'crossChannelStats');
if ~isfield(ssData, 'crossChannelStats') || ...
        ~isfield(ssData.crossChannelStats, 'unified_labels')
    return;
end
unif = ssData.crossChannelStats.unified_labels;
if ~isfield(unif, 'label') || isempty(unif.label), return; end

oldLab = unif.label(:);
nU     = numel(oldLab);
report.nUnits = nU;
if isequal(oldLab, (1:nU)')
    if verbose, fprintf('Repair: label == row index already holds for %d units.\n', nU); end
    report.ok = true;
    return;
end

lbl = double(h5read(uH5, '/unifiedLabels'));
sz  = size(lbl);
lbl = lbl(:);

maxOld = max([oldLab; lbl(lbl > 0)]);
lut = nan(maxOld, 1);
lut(oldLab) = (1:nU)';
sel  = find(lbl > 0 & lbl <= maxOld);
m    = lut(lbl(sel));
good = ~isnan(m);
lbl(sel(good)) = m(good);
unif.label = (1:nU)';

report.nRelabelled = sum(oldLab ~= (1:nU)');

lbl = reshape(lbl, sz);

% Two files, one logical change: back up and roll back on failure. Without
% this a failed .mat save leaves permuted spike labels against the old unit
% table, and the idempotence guard above would then re-permute on a re-run.
bk = kiaSort_backup_results(outputPath, 'prerepair', {uH5, ssP});
if ~bk.ok
    if verbose, fprintf('Repair: backup failed, not writing.\n'); end
    return;
end

try
    delete(uH5);
    h5create(uH5, '/unifiedLabels', size(lbl), 'Datatype', 'double');
    h5write(uH5,  '/unifiedLabels', lbl);

    ssData.crossChannelStats.unified_labels = unif;
    crossChannelStats = ssData.crossChannelStats;
    save(ssP, 'crossChannelStats', '-append');
catch ME
    kiaSort_restore_results(bk);
    if verbose
        fprintf('Repair: write failed (%s); rolled back from %s.\n', ME.message, bk.dir);
    end
    return;
end

report.changed = true;
report.ok      = true;
if verbose
    fprintf('Repair: renumbered %d of %d units; %d spikes remapped.\n', ...
        report.nRelabelled, nU, sum(good));
end
end
