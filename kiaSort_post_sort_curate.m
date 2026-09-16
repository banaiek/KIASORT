function postSortReport = kiaSort_post_sort_curate(outputPath, varargin)
%KIASORT_POST_SORT_CURATE  Post-hoc overlap removal, CCG cleanup, merging.
%
%   postSortReport = kiaSort_post_sort_curate(outputPath, ...)
%
%   Phases (each toggleable):
%     0) Overlap removal (overlap_removal): tiered, applied to the
%        lower-SNR side of each pair within coincYUm in y. A neighbour
%        kj only counts as the 'dominant' side when its SNR is at
%        least overlapSnrRatio (1.5) times the candidate's -- so two
%        co-firing real units of similar quality are never paired up
%        for a drop, no matter how high their overlap.
%          Tier 1: overlap >= overlapHighFrac (0.50) AND snr <
%                  overlapHighSnr (2.0)  -> drop unit
%          Tier 2: overlap >= overlapMidFrac  (0.20) AND snr <
%                  overlapMidSnr  (1.5)  -> drop unit
%          Tier 3: overlapLowFrac (0.05) <= overlap < overlapMidFrac
%                  (0.20) AND snr < overlapMidSnr (1.5) -> strip only
%                  the coincident spikes from the lower-SNR side
%     1) Cross-unit de-duplication (ccg_cleaning): a pair that double-
%        detects the same physical spike shows a sharp zero-lag CCG peak
%        while both ACGs stay clean. The lower-SNR unit (contaminated)
%        loses its coincident copies; the higher-SNR unit (owner) keeps
%        them, so recall is preserved. Gates before any strip:
%          - Ownership: owner = higher SNR. Equal SNR -> skip.
%          - CCG peak-at-zero above baseline (ccgPeakRatio).
%          - Waveform-similarity confirm: shared-channel similarity must
%            be >= dupWaveSim (0.8) -- same waveform = same source. The
%            same max_half_corr also gives the best-lag offset between
%            the two units' detections (the peak-vs-trough offset).
%          - Tight lag at the waveform offset: median signed lag within
%            dupLagTightSamples (2) of that best-lag offset, and at least
%            cleanLagConsistencyFrac of the coincident lags within
%            dupLagTightSamples of the median.
%          - Backstop: skip if removal would zero out more than
%            cleanMaxStripFrac of the contaminated unit.
%     2) Merge (merging): XCorr similarity + amp-similarity + PC
%        distance + multi-gate merged-ISI; shift absorbed spike times by
%        best-lag before relabel.
%
%   Per-pair errors are caught (try/catch around each pair) and the run
%   continues. unifiedLabels.h5 / spike_idx.h5 are rewritten only when
%   something actually changed.
%
%   Name/Value:
%       'ccg_cleaning'    (logical, true)   run Phase 1
%       'merging'         (logical, true)   run Phase 2
%       'xcorrThreshold'  (scalar, 0.9)     similarity (XCorr) cutoff
%       'thrAmp'          (scalar, 0.15)    amplitude-similarity gate
%       'thrPC'           (scalar, 0.45)    PC distance gate
%       'thrIsiAbs'       (scalar, 1.0)     hard cap on merged ISI %%
%       'thrIsiBudget'    (scalar, 1.5)     weighted-mean ISI budget
%       'thresholdISIms'  (scalar, 1)       refractory window (ms) for ISI
%       'ccgPeakRatio'    (scalar, 2)       CCG peak/median ratio gate
%       'coincYUm'        (scalar, 100)     y neighbourhood (um)
%       'cleanLagTightMs'         (scalar, 0.2)   tight window around median lag (ms)
%       'cleanLagConsistencyFrac' (scalar, 0.75)  min frac of coincident lags inside that window
%       'verbose'         (logical, false)  print per-pair / summary lines
%
%   INPUT FILES (relative to outputPath, all must exist for the run
%   to do any work; if any is missing the function returns silently)
%       RES_Sorted/spike_idx.h5
%       RES_Sorted/unifiedLabels.h5
%       RES_Sorted/channelNum.h5
%       RES_Samples/channel_info.mat
%       Sorted_Samples/sorted_samples.mat
%
%   OPTIONAL
%       RES_Sorted/drift_merge_report.mat  -- read for labelRemap so
%       we can pull pre-drift mean waveforms by the post-drift label.
%       If absent we assume identity remap.
%
%   OUTPUT
%       postSortReport.nClean        CCG strip count
%       postSortReport.nMerge        pair-merge count
%       postSortReport.nOverlapDrop  whole-unit overlap drops
%       postSortReport.nOverlapStrip per-spike overlap strips (Tier 3)
%       postSortReport.changed       true iff h5 files were rewritten
%       postSortReport.ok            true on clean finish

p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addParameter('ccg_cleaning', true,  @(x) islogical(x) || isnumeric(x));
p.addParameter('merging',      true,  @(x) islogical(x) || isnumeric(x));
p.addParameter('xcorrThreshold', 0.9, @(x) isscalar(x) && isnumeric(x));   % thrSim (XCorr metric)
p.addParameter('thrAmp',         0.1, @(x) isscalar(x) && isnumeric(x));  % auto-curate default
p.addParameter('thrPC',          0.45, @(x) isscalar(x) && isnumeric(x));  % calibrated once the gate actually fires
p.addParameter('thrIsiAbs',      1.0,  @(x) isscalar(x) && isnumeric(x));  % max ISI %
p.addParameter('thrIsiBudget',   1.5,  @(x) isscalar(x) && isnumeric(x));  % weighted-mean budget
p.addParameter('thresholdISIms', 1,    @(x) isscalar(x) && isnumeric(x));  % refractory window (ms)
p.addParameter('ccgPeakRatio',   2,   @(x) isscalar(x) && isnumeric(x));
p.addParameter('coincYUm',       100, @(x) isscalar(x) && isnumeric(x));
p.addParameter('cleanLagTightMs',          0.2, @(x) isscalar(x) && isnumeric(x));
p.addParameter('cleanLagConsistencyFrac',  0.75, @(x) isscalar(x) && isnumeric(x));
p.addParameter('cleanSnrDiff',             0.3, @(x) isscalar(x) && isnumeric(x));   % legacy, unused
p.addParameter('cleanMaxStripFrac',        0.5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('dupLagTightSamples',       2,   @(x) isscalar(x) && isnumeric(x));   % +-N sample consistency window
p.addParameter('dupWaveSim',               0.8, @(x) isscalar(x) && isnumeric(x));   % shared-channel waveform sim (0 disables)
p.addParameter('overlapMergeSim',          0.9, @(x) isscalar(x) && isnumeric(x));   % footprint-wide sim gate for merge escalation
p.addParameter('looseCorr',                0.8, @(x) isscalar(x) && isnumeric(x));   % pass-1 prefilter on the clustering means
p.addParameter('spikeCapN',               5000, @(x) isscalar(x) && isnumeric(x));   % spikes drawn per unit for pass 2
p.addParameter('minSpikesForMerge',         50, @(x) isscalar(x) && isnumeric(x));   % below this, fall back to the means
p.addParameter('mergeMaxSeparability',    0.92, @(x) isscalar(x) && isnumeric(x));   % cloud-overlap gate (0.5 = chance)
p.addParameter('mergedIsiMax',            0.20, @(x) isscalar(x) && isnumeric(x));   % refractory cap on the MERGED train
p.addParameter('wfCacheMB',                512, @(x) isscalar(x) && isnumeric(x));   % cap on the per-unit waveform cache
p.addParameter('overlap_removal', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('overlapHighFrac', 0.50, @(x) isscalar(x) && isnumeric(x));
p.addParameter('overlapHighSnr',  2.0,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('overlapMidFrac',  0.20, @(x) isscalar(x) && isnumeric(x));
p.addParameter('overlapMidSnr',   1.5,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('overlapLowFrac',  0.05, @(x) isscalar(x) && isnumeric(x));
p.addParameter('overlapSnrRatio', 1.5,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('verbose',        false, @(x) islogical(x) || isnumeric(x));
p.parse(outputPath, varargin{:});
opt = p.Results;
opt.ccg_cleaning    = logical(opt.ccg_cleaning);
opt.merging         = logical(opt.merging);
opt.overlap_removal = logical(opt.overlap_removal);
opt.verbose         = logical(opt.verbose);

postSortReport = struct('nClean', 0, 'nMerge', 0, ...
                        'nOverlapDrop', 0, 'nOverlapStrip', 0, ...
                        'changed', false, 'ok', false);

% Bail early if no phase is enabled.
if ~opt.ccg_cleaning && ~opt.merging && ~opt.overlap_removal
    postSortReport.ok = true;
    return;
end

% ---- Resolve paths ------------------------------------------------------
outputPath        = char(outputPath);
resSortedFolder   = fullfile(outputPath, 'RES_Sorted');
sortedSampFolder  = fullfile(outputPath, 'Sorted_Samples');
samplesFolder     = fullfile(outputPath, 'RES_Samples');

spikeIdxH5        = fullfile(resSortedFolder, 'spike_idx.h5');
unifiedLabelsH5   = fullfile(resSortedFolder, 'unifiedLabels.h5');
channelNumH5      = fullfile(resSortedFolder, 'channelNum.h5');
sortedSamplesPath = fullfile(sortedSampFolder, 'sorted_samples.mat');
channelInfoPath   = fullfile(samplesFolder,    'channel_info.mat');
driftReportPath   = fullfile(resSortedFolder, 'drift_merge_report.mat');

requiredFiles = {spikeIdxH5, unifiedLabelsH5, channelNumH5, ...
                 sortedSamplesPath, channelInfoPath};
for i = 1:numel(requiredFiles)
    if ~exist(requiredFiles{i}, 'file')
        if opt.verbose
            fprintf('Post-sort curate: %s missing, skipping.\n', requiredFiles{i});
        end
        return;
    end
end

% ---- Load HDF5 outputs --------------------------------------------------
try
    spk_all = double(h5read(spikeIdxH5,      '/spike_idx'));
    lbl_all = double(h5read(unifiedLabelsH5, '/unifiedLabels'));
    chn_all = double(h5read(channelNumH5,    '/channelNum'));
catch ME
    if opt.verbose
        fprintf('Post-sort curate: H5 read failed (%s), skipping.\n', ME.message);
    end
    return;
end
spk_all = spk_all(:); lbl_all = lbl_all(:); chn_all = chn_all(:);
if numel(spk_all) ~= numel(lbl_all) || numel(spk_all) ~= numel(chn_all)
    if opt.verbose
        fprintf('Post-sort curate: H5 lengths inconsistent, skipping.\n');
    end
    return;
end
if isempty(spk_all)
    postSortReport.ok = true;
    return;
end

% ---- Load probe geometry -----------------------------------------------
try
    chInfo = load(channelInfoPath);
catch ME
    if opt.verbose, fprintf('Post-sort curate: channel_info load failed (%s).\n', ME.message); end
    return;
end

% Peeked at ahead of the full sortedSamples load below, because the geometry
% fallback needs duplicateSearchChannels to size its neighbourhood.
cfgRaw0 = struct();
try
    ssPeek = load(sortedSamplesPath, 'sortedSamples');
    for iP = 1:numel(ssPeek.sortedSamples)
        if ~isempty(ssPeek.sortedSamples{iP}) && isfield(ssPeek.sortedSamples{iP}, 'cfg')
            cfgRaw0 = ssPeek.sortedSamples{iP}.cfg; break;
        end
    end
    clear ssPeek
catch
end
% Geometry is used for one thing only: deciding which unit pairs are close
% enough to be worth comparing. Returning without it meant a run with no
% channel map got NO merging and NO overlap removal at all -- and since the
% split pass has no such dependency, the result was a silently over-split
% output. Fall back to pairing by channel index instead, which is the right
% neighbourhood anyway when num_channel_extract is 0 (one channel per unit).
geomFallback = ~isfield(chInfo, 'channel_locations') || isempty(chInfo.channel_locations) ...
        || size(chInfo.channel_locations, 2) < 2;
if geomFallback
    nChTot = numel(chInfo.channel_inclusion);
    if nChTot < 1, nChTot = max(chn_all(~isnan(chn_all))); end
    ylocs  = (1:nChTot)';
    % In fallback the units are channel indices, not micrometres, so the
    % micrometre radius cannot apply. duplicateSearchChannels is the
    % existing knob for exactly this case.
    yRadius = 2;
    if isfield(cfgRaw0, 'duplicateSearchChannels') && ~isempty(cfgRaw0.duplicateSearchChannels)
        yRadius = double(cfgRaw0.duplicateSearchChannels);
    end
    warning('kiaSort:postSortCurate:noGeometry', ...
        ['No channel_locations (run without a channel map); pairing units by ' ...
         'channel index within +-%g instead of %g um.'], yRadius, opt.coincYUm);
    if opt.verbose
        fprintf('Post-sort curate: no geometry, pairing by channel index (+-%g).\n', yRadius);
    end
else
    ylocs   = chInfo.channel_locations(:, 2);
    yRadius = opt.coincYUm;
end

% ---- Load sortedSamples (for cfg + mean waveforms) ----------------------
try
    ssData = load(sortedSamplesPath, 'sortedSamples', 'crossChannelStats');
catch ME
    if opt.verbose, fprintf('Post-sort curate: sorted_samples load failed (%s).\n', ME.message); end
    return;
end
if ~isfield(ssData, 'sortedSamples') || ~isfield(ssData, 'crossChannelStats')
    if opt.verbose, fprintf('Post-sort curate: sorted_samples fields missing.\n'); end
    return;
end
sortedSamples     = ssData.sortedSamples;
% Recovered here, not inside Phase 1: the merge pass needs it too, and Phase 1
% is skipped when ccg_cleaning is off.
cfgRaw = [];
for iCfg = 1:numel(sortedSamples)
    if ~isempty(sortedSamples{iCfg}) && isfield(sortedSamples{iCfg}, 'cfg')
        cfgRaw = sortedSamples{iCfg}.cfg; break;
    end
end
crossChannelStats = ssData.crossChannelStats;
if ~isfield(crossChannelStats, 'unified_labels')
    if opt.verbose, fprintf('Post-sort curate: unified_labels missing.\n'); end
    return;
end
unif = crossChannelStats.unified_labels;
nU   = numel(unif.label);
if nU == 0
    postSortReport.ok = true;
    return;
end

fs = 30000;
for i = 1:numel(sortedSamples)
    if ~isempty(sortedSamples{i}) && isfield(sortedSamples{i}, 'cfg') ...
            && isfield(sortedSamples{i}.cfg, 'samplingFrequency')
        fs = sortedSamples{i}.cfg.samplingFrequency;
        break;
    end
end

% Pre-drift -> post-drift label remap (identity if drift_merge_report missing).
labelRemap = containers.Map('KeyType', 'double', 'ValueType', 'double');
if exist(driftReportPath, 'file')
    try
        dr = load(driftReportPath, 'driftReport');
        if isfield(dr, 'driftReport') && isfield(dr.driftReport, 'labelRemap') ...
                && isa(dr.driftReport.labelRemap, 'containers.Map')
            labelRemap = dr.driftReport.labelRemap;
        end
    catch
    end
end

% Per-unit info keyed by post-drift label (one representative per label,
% highest detectability wins).
postLabelPerUnif = zeros(nU, 1);
for u = 1:nU
    preL = unif.label(u);
    if isKey(labelRemap, double(preL))
        postLabelPerUnif(u) = labelRemap(double(preL));
    else
        postLabelPerUnif(u) = preL;
    end
end

uniqueLabels = unique(postLabelPerUnif);
% Drop the special "removed" label -1 (and any other negatives) -- only
% real units participate.
uniqueLabels = uniqueLabels(uniqueLabels >= 0);
keepUnifMask = ismember(postLabelPerUnif, uniqueLabels);
nUnits = numel(uniqueLabels);
if nUnits < 2
    postSortReport.ok = true;
    return;
end

unitInfo = struct('label',         num2cell(uniqueLabels(:)), ...
                  'channel',       num2cell(nan(nUnits, 1)), ...
                  'meanWF',        cell(nUnits, 1), ...
                  'spkRows',       cell(nUnits, 1), ...
                  'spkTimes',      cell(nUnits, 1), ...
                  'isiViol',       num2cell(zeros(nUnits, 1)), ...
                  'detectability', num2cell(zeros(nUnits, 1)), ...
                  'ny',            num2cell(nan(nUnits, 1)));

% Index by label for fast lookup.
labelToK = containers.Map('KeyType', 'double', 'ValueType', 'double');
for k = 1:nUnits
    labelToK(uniqueLabels(k)) = k;
end

% Pick representative unif entry per post-drift label (highest detectability).
for u = 1:nU
    if ~keepUnifMask(u), continue; end
    lab = postLabelPerUnif(u);
    if ~isKey(labelToK, lab), continue; end
    k = labelToK(lab);
    det = NaN;
    if isfield(unif, 'detectblity')
        det = unif.detectblity(u);
    elseif isfield(unif, 'detectability')
        det = unif.detectability(u);
    end
    if isnan(det), det = 0; end
    if isnan(unitInfo(k).detectability) || det > unitInfo(k).detectability
        unitInfo(k).detectability = det;
        ch  = double(unif.channelID(u));
        lL  = unif.labelInChannel(u);
        unitInfo(k).channel = ch;
        % Pull mean waveform from sortedSamples (best-effort).
        % Prefer the PER-UNIT template: labelInChannel is a per-channel class
        % id, and a post-hoc split child copies its parent's, so both would
        % resolve to the same stored waveform.
        unitInfo(k).meanWF = local_unitTemplate(unif, u);
        try
            rel = sortedSamples{ch}.clusteringInfo.clusterRelabeling;
            keptIdx = find(rel.newUniqueLabels == lL, 1);
            if ~isempty(keptIdx) && isempty(unitInfo(k).meanWF)
                % The full-length mean gives max_half_corr a wider lag
                % search than the clustering-length one.
                if isfield(rel, 'newMeanWaveformsFull') && ~isempty(rel.newMeanWaveformsFull)
                    src = rel.newMeanWaveformsFull;
                else
                    src = rel.newMeanWaveforms;
                end
                % reshape, not squeeze: a single-channel footprint would
                % otherwise come back transposed.
                mw = reshape(src(keptIdx, :, :), size(src,2), []);
                if ~isempty(mw)
                    unitInfo(k).meanWF = mw;
                end
            end
        catch
            % skip on any sortedSamples shape mismatch
        end
    end
end

% Spike rows per unit (rows into spk_all / lbl_all / chn_all).
% Group by label using accumarray-friendly indexing for speed.
[lblSorted, sortOrder] = sort(lbl_all);
% Find run boundaries.
boundaries = [0; find(diff(lblSorted) ~= 0); numel(lblSorted)];
for b = 1:numel(boundaries)-1
    rng = (boundaries(b)+1):boundaries(b+1);
    lab = lblSorted(rng(1));
    if ~isKey(labelToK, lab), continue; end
    k = labelToK(lab);
    rows = sortOrder(rng);
    unitInfo(k).spkRows  = rows;
    unitInfo(k).spkTimes = spk_all(rows);
    if numel(unitInfo(k).spkTimes) >= 2
        try
            [~, ~, isiV] = getISIViolations(unitInfo(k).spkTimes, fs, 2);
            unitInfo(k).isiViol = isiV;
        catch
            unitInfo(k).isiViol = 0;
        end
    end
    % Fall back to most common channel from the spike rows when the
    % unif representative didn't supply one.
    if isnan(unitInfo(k).channel)
        chs = chn_all(rows);
        chs = chs(~isnan(chs));
        if ~isempty(chs)
            unitInfo(k).channel = mode(chs);
        end
    end
    if ~isnan(unitInfo(k).channel)
        ch = unitInfo(k).channel;
        if ch >= 1 && ch <= numel(ylocs)
            unitInfo(k).ny = ylocs(ch);
        end
    end
end

% ---- Build neighbour list (pairs within yRadius in y) ------------------
% We do this once, up front, so both phases share the same candidate set.
yArr = [unitInfo.ny];
yArr = yArr(:);
% Drop any unit with no resolved channel; it can't be paired meaningfully.
hasY = ~isnan(yArr);
candIdx = find(hasY);
if numel(candIdx) < 2
    postSortReport.ok = true;
    return;
end

% Build pair list: every (k1, k2) with k1 < k2 and |yk1 - yk2| <= yRadius.
% Vectorise the distance check so this stays fast even with hundreds of
% units. Pre-allocate to the upper bound (n choose 2) and trim once.
nC = numel(candIdx);
maxPairs = nC * (nC - 1) / 2;
pairList = zeros(maxPairs, 2);
nPairs = 0;
for ii = 1:nC-1
    ki  = candIdx(ii);
    yi  = yArr(ki);
    rest = candIdx(ii+1:end);
    dy   = abs(yArr(rest) - yi);
    near = rest(dy <= yRadius);
    if isempty(near), continue; end
    blk = numel(near);
    pairList(nPairs+1:nPairs+blk, 1) = ki;
    pairList(nPairs+1:nPairs+blk, 2) = near;
    nPairs = nPairs + blk;
end
pairList = pairList(1:nPairs, :);
if isempty(pairList)
    postSortReport.ok = true;
    return;
end

% ---- Phase 0: overlap removal (tiered) ---------------------------------
nOverlapDrop  = 0;
nOverlapStrip = 0;
droppedLabels = [];
if opt.overlap_removal
    coincSamples = round(0.5e-3 * fs);
    droppedMask  = false(nUnits, 1);
    for ki = 1:nUnits
        if isnan(unitInfo(ki).ny), continue; end
        liveRows_i = unitInfo(ki).spkRows(lbl_all(unitInfo(ki).spkRows) == unitInfo(ki).label);
        spk_i = spk_all(liveRows_i);
        if numel(spk_i) < 10, continue; end
        snr_i = 1 + unitInfo(ki).detectability;
        if isnan(snr_i), snr_i = 0; end

        maxFrac    = 0;
        triggerRows = [];
        for kj = 1:nUnits
            if kj == ki || droppedMask(kj), continue; end
            if isnan(unitInfo(kj).ny), continue; end
            if abs(unitInfo(kj).ny - unitInfo(ki).ny) > yRadius, continue; end
            snr_j = 1 + unitInfo(kj).detectability;
            if isnan(snr_j), snr_j = 0; end
            if snr_j < snr_i * opt.overlapSnrRatio, continue; end
            liveRows_j = unitInfo(kj).spkRows(lbl_all(unitInfo(kj).spkRows) == unitInfo(kj).label);
            spk_j = spk_all(liveRows_j);
            if isempty(spk_j), continue; end
            d_ij = local_nearest_distance(spk_i, spk_j);
            coincMask = d_ij <= coincSamples;
            frac = sum(coincMask) / numel(spk_i);
            if frac > maxFrac
                maxFrac     = frac;
                triggerRows = liveRows_i(coincMask);
            end
        end

        if maxFrac >= opt.overlapHighFrac && snr_i < opt.overlapHighSnr
            droppedMask(ki) = true;
            droppedLabels(end+1,1) = unitInfo(ki).label; %#ok<AGROW>
            lbl_all(unitInfo(ki).spkRows) = -1;
            nOverlapDrop = nOverlapDrop + 1;
            if opt.verbose
                fprintf('overlap-drop T1: lbl %d (frac %.2f, snr %.2f)\n', ...
                    unitInfo(ki).label, maxFrac, snr_i);
            end
        elseif maxFrac >= opt.overlapMidFrac && snr_i < opt.overlapMidSnr
            droppedMask(ki) = true;
            droppedLabels(end+1,1) = unitInfo(ki).label; %#ok<AGROW>
            lbl_all(unitInfo(ki).spkRows) = -1;
            nOverlapDrop = nOverlapDrop + 1;
            if opt.verbose
                fprintf('overlap-drop T2: lbl %d (frac %.2f, snr %.2f)\n', ...
                    unitInfo(ki).label, maxFrac, snr_i);
            end
        elseif maxFrac >= opt.overlapLowFrac && maxFrac < opt.overlapMidFrac ...
                && snr_i < opt.overlapMidSnr && ~isempty(triggerRows)
            lbl_all(triggerRows) = -1;
            nOverlapStrip = nOverlapStrip + numel(triggerRows);
            if opt.verbose
                fprintf('overlap-strip T3: lbl %d (frac %.2f, %d spikes)\n', ...
                    unitInfo(ki).label, maxFrac, numel(triggerRows));
            end
        end
    end
end

% ---- Phase 1: cross-unit de-duplication --------------------------------
% A pair that double-detects the same physical spike shows a sharp zero-
% lag CCG peak while both ACGs stay clean. Owner = higher SNR (ties broken
% by spike count then label so equal-SNR fragments still get an owner).
% Strip the contaminated unit's coincident copies; the owner keeps them,
% so recall is preserved. Guards before any strip: CCG peak above
% baseline, high shared-channel waveform similarity, and tight lag at the
% waveform offset (systematic re-detection, not biology).
%   Below cleanMaxStripFrac: strip the confirmed duplicate copies.
%   Above it (mostly a duplicate): escalate to a recall-neutral MERGE --
%   strip the coincident copies and relabel + lag-shift the survivors into
%   the owner -- but only when a footprint-wide similarity (overlapMergeSim)
%   AND a merged-train refractory veto both pass; otherwise leave it be.
nClean = 0;
% spike times can change here (Phase 1 merge escalation shifts relabelled
% survivors) as well as in Phase 2; declare the flag + clamp bound once.
spkChanged = false;
maxSampAll = max(spk_all);
if opt.ccg_cleaning
    coincSamples = round(0.5e-3 * fs);
    tightSamples = max(1, round(opt.dupLagTightSamples));
    ccgBin       = round(1e-3 * fs);   % 1 ms
    ccgMaxLag    = 100 * ccgBin;       % +-100 ms
    ccgZeroTol   = 1;                  % +-1 ms peak window

    % Map the raw file (path from the saved cfg) to build genuine 2ms
    % templates for the waveform xcorr; fall back to padding the stored
    % 1ms template if the raw file isn't reachable.
    rawMap = []; chanMap = []; nSampRaw = 0; rawOK = false;
    dupHalf2 = round(1e-3 * fs);   % +-1ms -> 2ms window
    dupCapN  = 300;                % spikes averaged per template
    tplCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    if ~isempty(cfgRaw) && isfield(cfgRaw, 'fullFilePath') ...
            && (ischar(cfgRaw.fullFilePath) || isstring(cfgRaw.fullFilePath)) ...
            && exist(char(cfgRaw.fullFilePath), 'file') ...
            && isfield(cfgRaw, 'numChannels') && isfield(cfgRaw, 'dataType')
        try
            cfgRaw.outputFolder = outputPath;   % keep the memmap log local
            rawMap = map_input_file(char(cfgRaw.fullFilePath), cfgRaw);
            if isfield(chInfo, 'channel_mapping') && ~isempty(chInfo.channel_mapping)
                chanMap = double(chInfo.channel_mapping(:));
            else
                chanMap = (1:cfgRaw.numChannels)';
            end
            if isfield(chInfo, 'num_samples') && ~isempty(chInfo.num_samples)
                nSampRaw = double(chInfo.num_samples);
            else
                nSampRaw = size(rawMap.Data.data, 2);
            end
            rawOK = ~isempty(rawMap) && nSampRaw > 0 && ~isempty(chanMap);
        catch ME
            if opt.verbose
                fprintf('Post-sort de-dup: raw map failed (%s), padding instead.\n', ME.message);
            end
            rawMap = []; rawOK = false;
        end
    end

    for p = 1:size(pairList, 1)
        ki = pairList(p, 1); kj = pairList(p, 2);
        try
            % Re-pull spikes after possible earlier strips.
            kiRows = unitInfo(ki).spkRows;
            kjRows = unitInfo(kj).spkRows;
            kiLive = kiRows(lbl_all(kiRows) == unitInfo(ki).label);
            kjLive = kjRows(lbl_all(kjRows) == unitInfo(kj).label);
            spk_i = spk_all(kiLive);
            spk_j = spk_all(kjLive);
            if numel(spk_i) < 10 || numel(spk_j) < 10, continue; end

            % Owner = higher SNR; the other side is contaminated. Ties are
            % broken by spike count (more spikes = better-sampled owner),
            % then by label for determinism -- so two split fragments of
            % identical SNR still get an owner instead of being skipped.
            % The contaminated side loses spikes only if the pair then
            % clears every duplicate-signature gate below.
            snr_i = 1 + unitInfo(ki).detectability; if isnan(snr_i), snr_i = 0; end
            snr_j = 1 + unitInfo(kj).detectability; if isnan(snr_j), snr_j = 0; end
            ownerIsI = snr_i > snr_j || (snr_i == snr_j && ...
                (numel(spk_i) > numel(spk_j) || ...
                 (numel(spk_i) == numel(spk_j) && unitInfo(ki).label > unitInfo(kj).label)));
            if ownerIsI
                kCont = kj; kOwn = ki; contRows = kjLive;
                spkCont = spk_j; spkOwn = spk_i;
            else
                kCont = ki; kOwn = kj; contRows = kiLive;
                spkCont = spk_i; spkOwn = spk_j;
            end

            % Coincidence above chance: sharp zero-lag CCG peak.
            [ccg, zb] = local_pairCCG(spk_i, spk_j, ccgMaxLag, ccgBin);
            if ~local_ccgPeakAtZero(ccg, zb, ccgZeroTol, opt.ccgPeakRatio)
                continue;
            end

            % Shared-channel waveform confirm + expected offset. A true
            % duplicate is the same waveform on the contaminated unit's
            % channel; co-active different neurons are not. max_half_corr
            % also returns the best-lag offset where the coincident lags
            % should sit -- this is what lets the same spike picked up at
            % different peak/trough features (a non-zero but fixed offset)
            % still register, rather than only same-feature (lag 0) doubles.
            simWF = NaN; waveLag = 0;
            chC = unitInfo(kCont).channel;
            if rawOK && ~isnan(chC)
                tC = local_getTpl2ms(tplCache, kCont, chC, spkCont, ...
                    rawMap, chanMap, dupHalf2, nSampRaw, dupCapN);
                tO = local_getTpl2ms(tplCache, kOwn, chC, spkOwn, ...
                    rawMap, chanMap, dupHalf2, nSampRaw, dupCapN);
                if ~isempty(tC) && ~isempty(tO) && numel(tC) == numel(tO) && numel(tC) > 4
                    M2t = numel(tC);
                    [simWF, waveLag] = max_half_corr(tC(:)', tO(:)', 1, M2t, ...
                        max(1, round(M2t/4)), 0);
                end
            end
            if ~isfinite(simWF)
                % Raw unavailable (or too few clean snippets): pad the
                % stored 1ms template to 2ms instead.
                [simWF, waveLag] = local_pairWaveSim(unitInfo(kCont), unitInfo(kOwn), coincSamples);
            end
            if opt.dupWaveSim > 0 && (~isfinite(simWF) || simWF < opt.dupWaveSim)
                continue;
            end
            expectedLag = 0;
            if isfinite(waveLag), expectedLag = waveLag; end

            % Coincident contaminated spikes within +-0.5 ms of an owner spike.
            d_cont    = local_nearest_distance(spkCont, spkOwn);
            signedLag = local_signed_nearest_lag(spkCont, spkOwn);
            contMask  = d_cont <= coincSamples;
            if ~any(contMask), continue; end

            % Systematic re-detection: the coincident spikes' signed lags
            % must cluster within +-tightSamples of their median, and that
            % median must match the waveform offset (expectedLag) -- not a
            % random or synaptic coincidence.
            coincIdx  = find(contMask);
            coincLags = signedLag(coincIdx);
            coincLags = coincLags(~isnan(coincLags));
            if numel(coincLags) < 5, continue; end
            medLag = median(coincLags);
            if abs(medLag - expectedLag) > tightSamples, continue; end
            withinTight = abs(coincLags - medLag) <= tightSamples;
            if sum(withinTight) / numel(coincLags) < opt.cleanLagConsistencyFrac
                continue;
            end
            keepCoinc           = false(numel(signedLag), 1);
            keepCoinc(coincIdx) = withinTight;
            contMask            = contMask & keepCoinc;
            if ~any(contMask), continue; end

            % Backstop. Below the strip cap: strip the confirmed duplicate
            % copies (recall-neutral -- owner keeps its copy). Above the
            % cap the unit is mostly a duplicate; rather than walk away
            % (which is what left overlapping groups behind), escalate to a
            % recall-neutral MERGE: strip the coincident copies AND relabel
            % the unit's remaining spikes into the owner. Two guards make
            % the merge safe: footprint-wide waveform similarity (so two
            % look-alike co-active cells on one channel aren't merged) and
            % a merged-train refractory veto (a genuine two-neuron pair
            % shows an ISI violation when combined -> merge refused).
            if sum(contMask) / numel(contMask) > opt.cleanMaxStripFrac
                keepRows = contRows(~contMask);
                % Post-shift survivor times (their actual timing once merged)
                % drive the refractory veto.
                spkKeep = spk_all(keepRows);
                if isfinite(expectedLag) && expectedLag ~= 0 && ~isempty(spkKeep)
                    spkKeep = min(max(spkKeep + round(expectedLag), 1), maxSampAll);
                end
                % Merged-train refractory veto: owner + survivors must not
                % exceed the ISI cap (a genuine two-neuron pair would).
                combinedM  = sort([spkOwn(:); spkKeep(:)]);
                mergedISIok = true;
                if numel(combinedM) >= 2
                    try
                        [~, ~, isiM] = getISIViolations(combinedM, fs, opt.thresholdISIms);
                        mergedISIok = isiM <= opt.thrIsiAbs;
                    catch
                        mergedISIok = false;
                    end
                end
                fpSim = local_footprintSim(unitInfo(kCont), unitInfo(kOwn), 2, coincSamples);
                mergeSafe = isfinite(fpSim) && fpSim >= opt.overlapMergeSim && mergedISIok;
                if ~mergeSafe
                    continue;   % guards fail -> leave untouched (conservative)
                end
                lbl_all(contRows(contMask)) = -1;         % strip duplicate copies
                if ~isempty(keepRows)
                    spk_all(keepRows) = spkKeep;                 % commit lag-shift
                    if isfinite(expectedLag) && expectedLag ~= 0
                        spkChanged = true;
                    end
                    lbl_all(keepRows) = unitInfo(kOwn).label;   % relabel survivors
                    unitInfo(kOwn).spkRows = [unitInfo(kOwn).spkRows(:); keepRows(:)];
                end
                unitInfo(kCont).spkRows = [];   % contaminated unit merged away
                nClean = nClean + 1;
                continue;
            end

            lbl_all(contRows(contMask)) = -1;
            nClean = nClean + 1;
            liveAfter = unitInfo(kCont).spkRows;
            liveAfter = liveAfter(lbl_all(liveAfter) == unitInfo(kCont).label);
            if numel(liveAfter) >= 2
                try
                    [~, ~, isiNew] = getISIViolations(spk_all(liveAfter), fs, 2);
                    unitInfo(kCont).isiViol = isiNew;
                catch
                    unitInfo(kCont).isiViol = 0;
                end
            else
                unitInfo(kCont).isiViol = 0;
            end
        catch ME
            if opt.verbose
                fprintf('Post-sort de-dup: pair (%d,%d) failed: %s\n', ...
                    unitInfo(ki).label, unitInfo(kj).label, ME.message);
            end
        end
    end
end

% ---- Phase 2: Merging --------------------------------------------------
% Gates per pair: XCorr similarity, amp similarity, PC distance, merged-
% ISI multi-gate. On pass, the absorbed side's spike times are shifted
% by max_half_corr's bestLag before relabel.
nMerge     = 0;
mergedTo   = containers.Map('KeyType', 'double', 'ValueType', 'double');

if opt.merging
    chanPCA = containers.Map('KeyType','double','ValueType','any');

    % One draw per unit, reused across every pair it appears in. Without it a
    % unit is re-read once per pair (measured ~12x redundancy on a 122-unit
    % Utah array, and far worse on a dense probe).
    wfCache  = containers.Map('KeyType', 'double', 'ValueType', 'any');
    wfCacheB = 0;
    wfSrc = struct('ok', false);
    if ~isempty(cfgRaw)
        try
            wfSrc = kiaSort_waveform_source(outputPath, cfgRaw, opt.verbose);
        catch
            wfSrc = struct('ok', false);
        end
    end
    if opt.verbose
        if wfSrc.ok
            fprintf('Post-sort merge: judging pairs on individual spikes (%s).\n', wfSrc.mode);
        else
            fprintf('Post-sort merge: no waveform source, falling back to clustering means.\n');
        end
    end

    for p = 1:size(pairList, 1)
        ki = pairList(p, 1); kj = pairList(p, 2);
        try
            labI = local_resolveLabel(unitInfo(ki).label, mergedTo);
            labJ = local_resolveLabel(unitInfo(kj).label, mergedTo);
            if labI == labJ, continue; end

            wfI = unitInfo(ki).meanWF;
            wfJ = unitInfo(kj).meanWF;
            if isempty(wfI) || isempty(wfJ), continue; end

            chI = unitInfo(ki).channel;
            chJ = unitInfo(kj).channel;
            if isnan(chI) || isnan(chJ), continue; end

            % ---- pass 1: cheap prefilter on the clustering means -------
            % Polarity, merged-train refractoriness and a LOOSE correlation.
            % Nothing here reads spike waveforms, so a pair that is obviously
            % unrelated costs nothing.
            mw1 = local_rowOnChannel(wfI, chI, chI);
            mw2 = local_rowOnChannel(wfJ, chJ, chI);
            if isempty(mw1) || isempty(mw2), continue; end
            M2 = numel(mw1);
            if M2 < 5 || numel(mw2) ~= M2, continue; end
            if local_peakPolarity(mw1) ~= local_peakPolarity(mw2), continue; end

            spkI_now = spk_all(lbl_all == labI);
            spkJ_now = spk_all(lbl_all == labJ);
            if numel(spkI_now) < 1 || numel(spkJ_now) < 1, continue; end
            ok = local_checkMergedISI(spkI_now, spkJ_now, ...
                fs, opt.thresholdISIms, opt.thrIsiAbs, opt.thrIsiBudget, opt.mergedIsiMax);
            if ~ok, continue; end

            % Refractory cap on the merged train. Two cells cannot share a
            % refractory period, so this separates a genuine over-split from
            % a look-alike pair far more sharply than waveform shape does --
            % measured on this data the two groups sit at <=0.19 and >=0.25.
            % Kept separate from thrIsiAbs, which also caps each PARENT and
            % would otherwise block merging a slightly contaminated unit.
            if isfinite(opt.mergedIsiMax)
                try
                    [~, ~, isiMerged] = getISIViolations(sort([spkI_now(:); spkJ_now(:)]), ...
                        fs, opt.thresholdISIms);
                catch
                    isiMerged = Inf;
                end
                if ~isfinite(isiMerged) || isiMerged > opt.mergedIsiMax, continue; end
            end

            maxLag = max(1, round(M2/4));
            [simLoose, bestLag] = max_half_corr(mw1(:)', mw2(:)', 1, M2, maxLag, 0);
            if ~isfinite(simLoose) || simLoose < opt.looseCorr, continue; end

            % ---- pass 2: decide on the unit's OWN spikes ----------------
            % The clustering mean is a sample-stage artefact keyed by
            % labelInChannel, so two units can share it (a post-hoc split
            % child inherits its parent's). Templates rebuilt from the
            % matched spikes are per-unit by construction.
            simScore = simLoose;
            sepWI = []; sepWJ = []; sepLag = 0;
            if wfSrc.ok
                rowsI = unitInfo(ki).spkRows(lbl_all(unitInfo(ki).spkRows) == labI);
                rowsJ = unitInfo(kj).spkRows(lbl_all(unitInfo(kj).spkRows) == labJ);
                [WI, wfCacheB] = local_cachedWaveforms(wfCache, wfCacheB, labI, ...
                    rowsI, wfSrc, spk_all, chn_all, opt.spikeCapN, opt.wfCacheMB);
                [WJ, wfCacheB] = local_cachedWaveforms(wfCache, wfCacheB, labJ, ...
                    rowsJ, wfSrc, spk_all, chn_all, opt.spikeCapN, opt.wfCacheMB);
                if size(WI,1) >= opt.minSpikesForMerge && size(WJ,1) >= opt.minSpikesForMerge ...
                        && size(WI,2) == size(WJ,2)
                    eI = double(mean(WI, 1)); eJ = double(mean(WJ, 1));
                    Me = numel(eI);
                    [simScore, bestLag] = max_half_corr(eI, eJ, 1, Me, max(1,round(Me/4)), 0);
                    if ~isfinite(simScore) || simScore < opt.xcorrThreshold, continue; end
                    ampDiff = local_ampSimilarity(eI, eJ);
                    if ~isfinite(ampDiff) || ampDiff > opt.thrAmp, continue; end
                    % Two clouds that a clusterer cannot tell apart are one
                    % neuron. Recovery near chance -> merge; well separated
                    % -> two neurons, refuse.
                    sepWI = WI; sepWJ = WJ; sepLag = bestLag;
                else
                    % not enough spikes to judge on waveforms -- fall back to
                    % the mean-waveform gates rather than merging blind
                    if simScore < opt.xcorrThreshold, continue; end
                    ampDiff = local_ampSimilarity(mw1(:)', mw2(:)');
                    if ~isfinite(ampDiff) || ampDiff > opt.thrAmp, continue; end
                end
            else
                if simScore < opt.xcorrThreshold, continue; end
                ampDiff = local_ampSimilarity(mw1(:)', mw2(:)');
                if ~isfinite(ampDiff) || ampDiff > opt.thrAmp, continue; end
            end

            ok = local_checkPCDistance(wfI, wfJ, chI, chJ, sortedSamples, ...
                chanPCA, opt.thrPC);
            if ~ok, continue; end

            % Last, and only for pairs that already agree on everything else.
            % The clouds are aligned on bestLag first: comparing them unshifted
            % separates them on alignment alone, which is not evidence of two
            % neurons (a 1-sample offset scores identical populations at 1.00,
            % and 39% of correlation-passing pairs here carry a lag).
            %
            % The cut sits at 0.92, not 0.85: these clusters were produced by a
            % clusterer, so another clusterer can usually re-find the boundary
            % it drew. Measured on a Utah recording, a pair that is one neuron
            % by every other measure (corr 0.958, amplitude difference 0.031,
            % merged-train ISI 0.006%) still scored 0.905 here and was refused.
            % Above ~0.92 the two clouds really are distinct populations.
            if ~isempty(sepWI)
                sepAcc = local_cloudSeparability(sepWI, sepWJ, sepLag);
                if isfinite(sepAcc) && sepAcc > opt.mergeMaxSeparability, continue; end
            end

            nI = numel(spkI_now);
            nJ = numel(spkJ_now);
            if nI >= nJ
                primaryLab  = labI;
                absorbedLab = labJ;
                shiftSign   = -1;     % absorbed = J; max_half_corr lag = "J leads I", shift J by -bestLag
            else
                primaryLab  = labJ;
                absorbedLab = labI;
                shiftSign   = +1;     % absorbed = I; I leads J by -bestLag, shift I by +bestLag
            end

            if isfinite(bestLag) && bestLag ~= 0
                absRows = (lbl_all == absorbedLab);
                if any(absRows)
                    shifted = spk_all(absRows) + shiftSign * bestLag;
                    nSamp = max(spk_all);
                    shifted(shifted < 1)     = 1;
                    shifted(shifted > nSamp) = nSamp;
                    spk_all(absRows) = shifted;
                    spkChanged = true;
                end
            end

            lbl_all(lbl_all == absorbedLab) = primaryLab;
            mergedTo(absorbedLab) = primaryLab;
            % The primary's spike set just changed (and absorbed times may
            % have been shifted), so its cached draw is stale.
            if primaryLab == labI, kPri = ki; kAbs = kj; else, kPri = kj; kAbs = ki; end
            unitInfo(kPri).spkRows = [unitInfo(kPri).spkRows(:); unitInfo(kAbs).spkRows(:)];
            unitInfo(kAbs).spkRows = [];
            if isKey(wfCache, primaryLab),  remove(wfCache, primaryLab);  end
            if isKey(wfCache, absorbedLab), remove(wfCache, absorbedLab); end
            nMerge = nMerge + 1;
        catch ME
            if opt.verbose
                fprintf('Post-sort merge: pair (%d,%d) failed: %s\n', ...
                    unitInfo(ki).label, unitInfo(kj).label, ME.message);
            end
        end
    end
end

% ---- Write back if anything changed ------------------------------------
changed = (nClean > 0) || (nMerge > 0) || spkChanged || ...
          (nOverlapDrop > 0) || (nOverlapStrip > 0);
nRemoved   = 0;
compactMap = containers.Map('KeyType','double','ValueType','double');
if changed
    % Merging and dropping move spikes off a unit but left its row in the
    % table, so the two outputs disagreed: the GUI sizes itself from
    % numel(label) and showed the leftovers as empty units. Compact here,
    % while both are still in hand, so they are written in agreement.
    [unif, lbl_all, nRemoved, compactMap] = kiaSort_compact_unit_table(unif, lbl_all);

    % The labels and the unit table are separate files. Back both up so a
    % failure between the two writes can be rolled back rather than leaving
    % spike labels that disagree with the table.
    bk = kiaSort_backup_results(outputPath, 'postcurate', ...
        {unifiedLabelsH5, spikeIdxH5, sortedSamplesPath});

    try
        if exist(unifiedLabelsH5, 'file')
            delete(unifiedLabelsH5);
        end
        h5create(unifiedLabelsH5, '/unifiedLabels', size(lbl_all), 'Datatype', 'double');
        h5write(unifiedLabelsH5,  '/unifiedLabels', lbl_all);

        if spkChanged
            if exist(spikeIdxH5, 'file')
                delete(spikeIdxH5);
            end
            h5create(spikeIdxH5, '/spike_idx', size(spk_all), 'Datatype', 'double');
            h5write(spikeIdxH5,  '/spike_idx', spk_all);
        end

        if nRemoved > 0
            crossChannelStats.unified_labels = unif;
            save(sortedSamplesPath, 'crossChannelStats', '-append');
        end
    catch ME
        if opt.verbose
            fprintf('Post-sort curate: write failed (%s).\n', ME.message);
        end
        if isstruct(bk) && isfield(bk, 'ok') && bk.ok
            kiaSort_restore_results(bk);
        end
        postSortReport.nClean        = nClean;
        postSortReport.nMerge        = nMerge;
        postSortReport.nOverlapDrop  = nOverlapDrop;
        postSortReport.nOverlapStrip = nOverlapStrip;
        postSortReport.changed       = false;
        postSortReport.ok            = false;
        return;
    end
end

postSortReport.nClean        = nClean;
postSortReport.nMerge        = nMerge;
postSortReport.nOverlapDrop  = nOverlapDrop;
postSortReport.nOverlapStrip = nOverlapStrip;
postSortReport.droppedLabels = droppedLabels;
postSortReport.nEmptyRemoved = nRemoved;
postSortReport.compactRemap  = compactMap;
postSortReport.changed       = changed;
postSortReport.ok            = true;

if opt.verbose
    fprintf(['Post-sort curate: %d overlap drops, %d overlap strips, %d CCG strips, ' ...
             '%d merges, %d empty table rows removed (changed=%d).\n'], ...
        nOverlapDrop, nOverlapStrip, nClean, nMerge, nRemoved, changed);
end

end


% =========================================================================
% Local helpers
% =========================================================================

function [ccg, zeroBinIdx] = local_pairCCG(spkA, spkB, maxLagSamples, binSamples)
% CCG of (spkB - spkA) inside +-maxLagSamples, binSamples bin width.
spkA = sort(spkA(:));
spkB = sort(spkB(:));
nA   = numel(spkA);
nB   = numel(spkB);
halfBins   = floor(maxLagSamples / binSamples);
nBins      = 2*halfBins + 1;
zeroBinIdx = halfBins + 1;
ccg        = zeros(nBins, 1);
if nA == 0 || nB == 0, return; end

b_lo = 1; b_hi = 0;
diffsCell = cell(nA, 1);
for i = 1:nA
    t = spkA(i);
    while b_lo <= nB && spkB(b_lo) < t - maxLagSamples
        b_lo = b_lo + 1;
    end
    while b_hi < nB && spkB(b_hi + 1) <= t + maxLagSamples
        b_hi = b_hi + 1;
    end
    if b_lo <= b_hi
        diffsCell{i} = spkB(b_lo:b_hi) - t;
    end
end
if all(cellfun(@isempty, diffsCell)), return; end
allDiffs = vertcat(diffsCell{:});
if isempty(allDiffs), return; end
binIdx = round(double(allDiffs) / binSamples) + zeroBinIdx;
keep   = binIdx >= 1 & binIdx <= nBins;
binIdx = binIdx(keep);
if isempty(binIdx), return; end
ccg = accumarray(binIdx(:), 1, [nBins, 1]);
end


function tf = local_ccgPeakAtZero(ccg, zeroBinIdx, zeroTolBins, ratio)
% True iff the lag-0 window holds the global max and that peak exceeds
% ratio * 90th-percentile of the off-zero bins (falls back to mean of
% non-zero off-bins when the percentile is 0).
tf = false;
if isempty(ccg) || all(ccg == 0), return; end
nB = numel(ccg);
peakRange = max(1, zeroBinIdx-zeroTolBins) : min(nB, zeroBinIdx+zeroTolBins);
if isempty(peakRange), return; end
peakValue = max(ccg(peakRange));
offRange  = ccg;
offRange(peakRange) = [];
if isempty(offRange), return; end
if any(offRange > peakValue), return; end
if ~any(offRange > 0), return; end
baseline = prctile(offRange, 90);
if ~isfinite(baseline) || baseline <= 0
    nz = offRange(offRange > 0);
    if isempty(nz), return; end
    baseline = mean(nz);
    if baseline <= 0, return; end
end
tf = peakValue > ratio * baseline;
end


function d = local_nearest_distance(A, B)
% For each A(i), distance to nearest B (samples). +Inf if B is empty.
A = A(:); B = B(:);
nA = numel(A); nB = numel(B);
d  = inf(nA, 1);
if nA == 0 || nB == 0, return; end
[Bs, ~] = sort(B);
[As, srtA] = sort(A);
j = 1;
for i = 1:nA
    while j < nB && Bs(j+1) < As(i)
        j = j + 1;
    end
    cand = abs(As(i) - Bs(j));
    if j+1 <= nB
        cand = min(cand, abs(As(i) - Bs(j+1)));
    end
    if j-1 >= 1
        cand = min(cand, abs(As(i) - Bs(j-1)));
    end
    d(srtA(i)) = cand;
end
end


function lag = local_signed_nearest_lag(A, B)
% lag(i) = B_nearest - A(i). Positive = B after A. NaN when B is empty.
A = A(:); B = B(:);
nA = numel(A); nB = numel(B);
lag = nan(nA, 1);
if nA == 0 || nB == 0, return; end
[Bs, ~]    = sort(B);
[As, srtA] = sort(A);
j = 1;
for i = 1:nA
    while j < nB && Bs(j+1) < As(i)
        j = j + 1;
    end
    bestDiff = Bs(j) - As(i);
    bestAbs  = abs(bestDiff);
    if j+1 <= nB
        d2 = Bs(j+1) - As(i);
        if abs(d2) < bestAbs
            bestDiff = d2; bestAbs = abs(d2);
        end
    end
    if j-1 >= 1
        d2 = Bs(j-1) - As(i);
        if abs(d2) < bestAbs
            bestDiff = d2; bestAbs = abs(d2);
        end
    end
    lag(srtA(i)) = bestDiff;
end
end


function lab = local_resolveLabel(lab0, mergedTo)
% Follow the merged-to chain so transitive merges collapse correctly.
lab = lab0;
while isKey(mergedTo, lab)
    next = mergedTo(lab);
    if next == lab, break; end
    lab = next;
end
end


function row = local_rowOnChannel(wf, homeCh, targetCh)
% Row of (nLoc x T) mean waveform corresponding to global targetCh,
% given the footprint is centred on homeCh. [] if targetCh out of range.
row = [];
if isempty(wf) || ~ismatrix(wf), return; end
nLoc  = size(wf, 1);
half  = floor((nLoc - 1) / 2);
% Home channel is the CENTRE row, so it is +half+1, not +half. Without the
% +1 a single-channel footprint resolves to row 0 and every merge gate bails.
localIdx = targetCh - homeCh + half + 1;
if localIdx < 1 || localIdx > nLoc, return; end
row = wf(localIdx, :);
end


function tpl = local_getTpl2ms(cache, g, ch, spk, rawMap, chanMap, half2, nSamp, capN)
% Cached 2ms mean template (1 x 2*half2+1) for unit g on global channel
% ch, built by averaging up to capN of the unit's spike snippets read
% from the raw memmap. [] when ch is out of range or too few clean
% snippets land fully inside the recording. Keyed by (g, ch) so the same
% unit-on-channel is read once per run.
tpl = [];
key = sprintf('%d_%d', g, ch);
if isKey(cache, key)
    tpl = cache(key);
    return;
end
if ch >= 1 && ch <= numel(chanMap)
    rawRow = chanMap(ch);
    if ~isnan(rawRow) && rawRow >= 1
        s = round(spk(:));
        s = s(s > half2 & s <= nSamp - half2);
        if numel(s) >= 5
            if numel(s) > capN
                s = s(round(linspace(1, numel(s), capN)));
            end
            W   = 2*half2 + 1;
            acc = zeros(1, W);
            cnt = 0;
            for ii = 1:numel(s)
                a = s(ii) - half2;
                acc = acc + double(rawMap.Data.data(rawRow, a:a+W-1));
                cnt = cnt + 1;
            end
            if cnt >= 5, tpl = acc / cnt; end
        end
    end
end
cache(key) = tpl;
end


function [s, bestLag] = local_pairWaveSim(uCont, uOwn, padN)
% Shared-channel waveform similarity and best-lag offset (max_half_corr)
% between two units, measured on the contaminated unit's main channel.
% bestLag is the expected cont->owner detection offset, so a duplicate
% picked up at different peak/trough features still registers. [NaN, 0]
% if either footprint can't supply a row on that channel.
%
% padN zero-pads each side of the ~1ms template (typically 0.5ms worth)
% -> 2ms. max_half_corr's validMask ignores the zero pad, so it still
% correlates the real core but can search +-padN of lag -- enough to
% lock onto a peak/trough offset that the un-padded 1ms template clips.
s = NaN; bestLag = 0;
if isempty(uCont.meanWF) || isempty(uOwn.meanWF), return; end
chC = uCont.channel; chO = uOwn.channel;
if isnan(chC) || isnan(chO), return; end
rC = local_rowOnChannel(uCont.meanWF, chC, chC);
rO = local_rowOnChannel(uOwn.meanWF,  chO, chC);
if isempty(rC) || isempty(rO) || numel(rC) ~= numel(rO) || numel(rC) < 5
    return;
end
if nargin >= 3 && padN > 0
    z  = zeros(padN, 1);
    rC = [z; rC(:); z];
    rO = [z; rO(:); z];
end
M2 = numel(rC);
[s, bestLag] = max_half_corr(rC(:)', rO(:)', 1, M2, max(1, round(M2/4)), 0);
end


function s = local_footprintSim(uCont, uOwn, nSide, padN)
% Mean waveform similarity across the contaminated unit's main channel
% and +-nSide neighbours (footprint-wide), so two co-active cells that
% merely look alike on a single shared channel do NOT clear the gate.
% Only channels where BOTH units carry real signal are counted. NaN if
% fewer than 2 usable channels.
s = NaN;
if isempty(uCont.meanWF) || isempty(uOwn.meanWF), return; end
chC = uCont.channel; chO = uOwn.channel;
if isnan(chC) || isnan(chO), return; end
sims = [];
for dc = -nSide:nSide
    ch = chC + dc;
    rC = local_rowOnChannel(uCont.meanWF, chC, ch);
    rO = local_rowOnChannel(uOwn.meanWF,  chO, ch);
    if isempty(rC) || isempty(rO) || numel(rC) ~= numel(rO) || numel(rC) < 5
        continue;
    end
    if max(abs(rC)) <= eps || max(abs(rO)) <= eps, continue; end
    if nargin >= 4 && padN > 0
        z  = zeros(padN, 1);
        rC = [z; rC(:); z];
        rO = [z; rO(:); z];
    end
    M2 = numel(rC);
    c  = max_half_corr(rC(:)', rO(:)', 1, M2, max(1, round(M2/4)), 0);
    if isfinite(c), sims(end+1) = c; end %#ok<AGROW>
end
if numel(sims) >= 2, s = mean(sims); end
end


function [W, bytesOut] = local_cachedWaveforms(cache, bytesIn, lab, rows, src, spk, chn, capN, capMB)
% containers.Map is a handle, so the store persists in the caller. Cleared
% wholesale once it passes the byte cap -- simpler than an eviction policy
% and the pairs for one channel are processed together, so locality is good.
bytesOut = bytesIn;
if isKey(cache, lab)
    W = cache(lab);
    return;
end
W = kiaSort_read_waveforms(src, rows, spk, chn, capN);
W = single(W);
if bytesOut > capMB * 1e6
    remove(cache, keys(cache));
    bytesOut = 0;
end
cache(lab) = W;
bytesOut = bytesOut + numel(W) * 4;
end


function pol = local_peakPolarity(w)
% Sign of the dominant excursion. Two units with opposite polarity are not
% the same neuron, and this costs nothing to check.
w = w(:)';
if isempty(w) || all(~isfinite(w)), pol = 0; return; end
if max(w) >= abs(min(w)), pol = 1; else, pol = -1; end
end


function acc = local_cloudSeparability(WA, WB, lag)
% How well a blind 2-means on the pooled spikes recovers which unit each
% spike came from. ~0.5 means the two clouds are one population; high means
% they are genuinely distinct. Deterministic: seeded from the group means,
% so no RNG and no dependence on spike order.
acc = NaN;
nA = size(WA,1); nB = size(WB,1);
if nA < 5 || nB < 5 || size(WA,2) ~= size(WB,2), return; end
if nargin >= 3 && isfinite(lag) && lag ~= 0
    T = size(WA,2); L = round(lag);
    if abs(L) < T - 4
        if L > 0
            WA = WA(:, 1+L:T);  WB = WB(:, 1:T-L);
        else
            WA = WA(:, 1:T+L);  WB = WB(:, 1-L:T);
        end
    end
end
X = double([WA; WB]);
X = X - mean(X, 1);
nComp = min(5, min(size(X)) - 1);
if nComp < 1, return; end
try
    [~, F] = pca(X, 'Algorithm', 'svd', 'NumComponents', nComp);
catch
    return;
end
if isempty(F) || size(F,1) < 10, return; end
truth = [true(nA,1); false(nB,1)];
C = [mean(F(truth,:), 1); mean(F(~truth,:), 1)];
if any(~isfinite(C(:))), return; end
try
    k = kmeans(F, 2, 'Start', C, 'MaxIter', 300);
catch
    return;
end
a = mean((k == 1) == truth);
acc = max(a, 1 - a);
end


function ampDiff = local_ampSimilarity(mw1, mw2)
% Average of |peak diff|/maxAbs and |trough diff|/maxAbs. NaN if flat.
mw1 = mw1(:)'; mw2 = mw2(:)';
both = [mw1; mw2];
maxAbs = max(abs(both), [], 2);
maxP   = abs(max(both, [], 2));
maxN   = abs(min(both, [], 2));
denom  = max([maxAbs(1), maxAbs(2)]);
if ~isfinite(denom) || denom <= 0
    ampDiff = NaN;
    return;
end
ampDiffP = (maxP(1) - maxP(2)) / denom;
ampDiffN = (maxN(1) - maxN(2)) / denom;
ampDiff  = (abs(ampDiffP) + abs(ampDiffN)) / 2;
if ampDiff == 0, ampDiff = NaN; end
end


function ok = local_checkPCDistance(wfI, wfJ, chI, chJ, sortedSamples, chanPCA, thr)
% PC distance between mwA and mwB on chI, normalised by the basis's
% score range. Returns true (pass) when the cached PCA is missing.
ok = true;
if isnan(chI), return; end
mwA = local_rowOnChannel(wfI, chI, chI);
mwB = local_rowOnChannel(wfJ, chJ, chI);
if isempty(mwA) || isempty(mwB), return; end
if numel(mwA) ~= numel(mwB), return; end

if isKey(chanPCA, chI)
    PC = chanPCA(chI);
else
    PC = struct('coeff', [], 'mu', [], 'max', []);
    try
        if numel(sortedSamples) >= chI && ~isempty(sortedSamples{chI}) ...
                && isfield(sortedSamples{chI}, 'clusteringInfo') ...
                && isfield(sortedSamples{chI}.clusteringInfo, 'PCA')
            tmp = sortedSamples{chI}.clusteringInfo.PCA;
            if isfield(tmp, 'coeff') && ~isempty(tmp.coeff) ...
                    && isfield(tmp, 'mu')   && ~isempty(tmp.mu)
                PC.coeff = tmp.coeff;
                PC.mu    = tmp.mu;
                if isfield(tmp, 'max') && ~isempty(tmp.max)
                    PC.max = tmp.max;
                elseif isfield(tmp, 'score') && ~isempty(tmp.score)
                    PC.max = max(abs(tmp.score), [], 1);
                else
                    PC.max = 1;
                end
            end
        end
    catch
    end
    chanPCA(chI) = PC;
end

if isempty(PC.coeff) || isempty(PC.mu), return; end
% The basis is fit on the CLUSTERING crop and, for a wide footprint, on the
% flattened footprint -- while meanWF spans spikeDuration. Bailing on the
% length mismatch left this gate permanently inert. Rebuild the vector the
% basis expects instead. Cross-channel pairs are skipped: the two footprints
% are centred on different channels, so flattening them is not comparable.
Lb = size(PC.coeff, 1);
if chI ~= chJ, return; end
mwA = local_pcVector(wfI, Lb);
mwB = local_pcVector(wfJ, Lb);
if isempty(mwA) || isempty(mwB), return; end
pcA = (mwA - PC.mu) * PC.coeff;
pcB = (mwB - PC.mu) * PC.coeff;
scale = max(PC.max(:));
if ~isfinite(scale) || scale <= 0, scale = 1; end
ok = norm(pcA - pcB) / scale <= thr;
end


function v = local_pcVector(wf, Lb)
% Vector matching the stored PCA basis: footprint centre-cropped in time to
% Lb/nChannels samples, then flattened channel-fastest -- the order
% reshape(waveform, N, C*T) produced when the basis was fit.
v = [];
if isempty(wf) || ~ismatrix(wf), return; end
C = size(wf, 1); T = size(wf, 2);
if C < 1 || T < 1 || mod(Lb, C) ~= 0, return; end
Tc = Lb / C;
if Tc > T, return; end
o = floor((T - Tc) / 2);
v = reshape(wf(:, o+1:o+Tc), 1, []);
end


function ok = local_checkMergedISI(spkA, spkB, fs, threshMs, isiAbs, isiBudget, isiFloor)
% Multi-gate merged-ISI test: hard cap, per-parent cap, parent deviation
% from size-weighted mean, and size-shrunk merged budget.
ok = true;
NA = numel(spkA); NB = numel(spkB);
if NA < 1 || NB < 1, return; end
spkM = sort([spkA(:); spkB(:)]);
if numel(spkM) < 2, return; end
try
    [~, ~, isiM] = getISIViolations(spkM, fs, threshMs);
catch
    return;
end
isiA = 0; isiB = 0;
if NA >= 2
    try, [~, ~, isiA] = getISIViolations(spkA, fs, threshMs); catch, isiA = 0; end
end
if NB >= 2
    try, [~, ~, isiB] = getISIViolations(spkB, fs, threshMs); catch, isiB = 0; end
end

if nargin < 7 || isempty(isiFloor), isiFloor = 0; end
parentCap   = isiAbs * max(isiBudget, 1);
denom       = max(NA + NB, 1);
weightedIsi = (NA * isiA + NB * isiB) / denom;

weightedDevCap = 0.5;
devA = abs(isiA - weightedIsi);
devB = abs(isiB - weightedIsi);

sizeFactor = 2 * min(NA, NB) / denom;
effBudget  = 1 + (isiBudget - 1) * sizeFactor;
% Absolute floor. budgetCap alone is purely relative, so two immaculate
% parents (ISI ~0.01%) give a cap near zero and a merged unit at 0.06% --
% clean by any standard -- is refused for exceeding 1.5x of almost nothing.
% The floor makes it "clean in absolute terms OR no worse than the parents".
budgetCap  = max([effBudget * weightedIsi, isiFloor, 1e-6]);

ok = (isiM <= isiAbs) && ...
     (isiA <= parentCap) && ...
     (isiB <= parentCap) && ...
     (devA <= weightedDevCap) && ...
     (devB <= weightedDevCap) && ...
     (isiM <= budgetCap);
end


function mw = local_unitTemplate(unif, u)
% Per-unit mean waveform as (nChannels x nSamples), [] when unavailable.
mw = [];
if ~isfield(unif, 'meanWaveforms') || isempty(unif.meanWaveforms), return; end
MW = unif.meanWaveforms;
if size(MW,1) < u, return; end
if ndims(MW) == 3
    mw = reshape(MW(u,:,:), size(MW,2), []);
else
    mw = reshape(MW(u,:), 1, []);
end
if all(~isfinite(mw(:))) || all(mw(:) == 0), mw = []; end
end
