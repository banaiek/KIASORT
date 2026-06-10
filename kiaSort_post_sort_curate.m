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
%       'thrPC'           (scalar, 0.30)    PC1/PC2 distance gate
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
p.addParameter('thrPC',          0.30, @(x) isscalar(x) && isnumeric(x));  % auto-curate default
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
if ~isfield(chInfo, 'channel_locations') || isempty(chInfo.channel_locations) ...
        || size(chInfo.channel_locations, 2) < 2
    if opt.verbose, fprintf('Post-sort curate: channel_locations missing.\n'); end
    return;
end
ylocs = chInfo.channel_locations(:, 2);

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
        try
            rel = sortedSamples{ch}.clusteringInfo.clusterRelabeling;
            keptIdx = find(rel.newUniqueLabels == lL, 1);
            if ~isempty(keptIdx)
                mw = squeeze(rel.newMeanWaveforms(keptIdx, :, :));
                if ~isempty(mw) && ismatrix(mw)
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

% ---- Build neighbour list (pairs within coincYUm in y) -----------------
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

% Build pair list: every (k1, k2) with k1 < k2 and |yk1 - yk2| <= coincYUm.
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
    near = rest(dy <= opt.coincYUm);
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
            if abs(unitInfo(kj).ny - unitInfo(ki).ny) > opt.coincYUm, continue; end
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
% lag CCG peak while both ACGs stay clean. Strip the lower-SNR unit's
% coincident copies; the higher-SNR unit (owner) keeps them, so recall is
% preserved -- only double-counted spikes are removed, never a unit's
% independent activity. Guards: CCG peak above baseline, high shared-
% channel waveform similarity (same waveform = same source), and tight &
% near-zero lag (systematic re-detection, not biology).
nClean = 0;
if opt.ccg_cleaning
    coincSamples = round(0.5e-3 * fs);
    tightSamples = max(1, round(opt.dupLagTightSamples));
    ccgBin       = round(1e-3 * fs);   % 1 ms
    ccgMaxLag    = 100 * ccgBin;       % +-100 ms
    ccgZeroTol   = 1;                  % +-1 ms peak window

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

            % Owner = higher SNR; contaminated = lower SNR. Equal SNR ->
            % no clear owner, skip.
            snr_i = 1 + unitInfo(ki).detectability; if isnan(snr_i), snr_i = 0; end
            snr_j = 1 + unitInfo(kj).detectability; if isnan(snr_j), snr_j = 0; end
            if snr_i == snr_j, continue; end
            if snr_i < snr_j
                kCont = ki; kOwn = kj; contRows = kiLive;
                spkCont = spk_i; spkOwn = spk_j;
            else
                kCont = kj; kOwn = ki; contRows = kjLive;
                spkCont = spk_j; spkOwn = spk_i;
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
            [simWF, waveLag] = local_pairWaveSim(unitInfo(kCont), unitInfo(kOwn));
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

            % Backstop: refuse to gut the unit on a detection misfire.
            if sum(contMask) / numel(contMask) > opt.cleanMaxStripFrac
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
spkChanged = false;
mergedTo   = containers.Map('KeyType', 'double', 'ValueType', 'double');

if opt.merging
    chanPCA = containers.Map('KeyType','double','ValueType','any');

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

            mw1 = local_rowOnChannel(wfI, chI, chI);
            mw2 = local_rowOnChannel(wfJ, chJ, chI);
            if isempty(mw1) || isempty(mw2), continue; end
            M2 = numel(mw1);
            if M2 < 5 || numel(mw2) ~= M2, continue; end
            maxLag = max(1, round(M2/4));
            [simScore, bestLag] = max_half_corr(mw1(:)', mw2(:)', 1, M2, maxLag, 0);
            if ~isfinite(simScore) || simScore < opt.xcorrThreshold, continue; end

            ampDiff = local_ampSimilarity(mw1(:)', mw2(:)');
            if ~isfinite(ampDiff) || ampDiff > opt.thrAmp, continue; end

            ok = local_checkPCDistance(wfI, wfJ, chI, chJ, sortedSamples, ...
                chanPCA, opt.thrPC);
            if ~ok, continue; end

            spkI_now = spk_all(lbl_all == labI);
            spkJ_now = spk_all(lbl_all == labJ);
            if numel(spkI_now) < 1 || numel(spkJ_now) < 1, continue; end
            ok = local_checkMergedISI(spkI_now, spkJ_now, ...
                fs, opt.thresholdISIms, opt.thrIsiAbs, opt.thrIsiBudget);
            if ~ok, continue; end

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
if changed
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
    catch ME
        if opt.verbose
            fprintf('Post-sort curate: H5 write failed (%s).\n', ME.message);
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
postSortReport.changed       = changed;
postSortReport.ok            = true;

if opt.verbose
    fprintf('Post-sort curate: %d overlap drops, %d overlap strips, %d CCG strips, %d merges (changed=%d).\n', ...
        nOverlapDrop, nOverlapStrip, nClean, nMerge, changed);
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
localIdx = targetCh - (homeCh - half);
if localIdx < 1 || localIdx > nLoc, return; end
row = wf(localIdx, :);
end


function [s, bestLag] = local_pairWaveSim(uCont, uOwn)
% Shared-channel waveform similarity and best-lag offset (max_half_corr)
% between two units, measured on the contaminated unit's main channel.
% bestLag is the expected cont->owner detection offset, so a duplicate
% picked up at different peak/trough features still registers. [NaN, 0]
% if either footprint can't supply a row on that channel.
s = NaN; bestLag = 0;
if isempty(uCont.meanWF) || isempty(uOwn.meanWF), return; end
chC = uCont.channel; chO = uOwn.channel;
if isnan(chC) || isnan(chO), return; end
rC = local_rowOnChannel(uCont.meanWF, chC, chC);
rO = local_rowOnChannel(uOwn.meanWF,  chO, chC);
if isempty(rC) || isempty(rO) || numel(rC) ~= numel(rO) || numel(rC) < 5
    return;
end
M2 = numel(rC);
[s, bestLag] = max_half_corr(rC(:)', rO(:)', 1, M2, max(1, round(M2/4)), 0);
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
mwA = mwA(:)'; mwB = mwB(:)';
if numel(mwA) ~= size(PC.coeff, 1) || numel(mwB) ~= size(PC.coeff, 1)
    return;
end
pcA = (mwA - PC.mu) * PC.coeff;
pcB = (mwB - PC.mu) * PC.coeff;
scale = max(PC.max(:));
if ~isfinite(scale) || scale <= 0, scale = 1; end
ok = norm(pcA - pcB) / scale <= thr;
end


function ok = local_checkMergedISI(spkA, spkB, fs, threshMs, isiAbs, isiBudget)
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

parentCap   = isiAbs * max(isiBudget, 1);
denom       = max(NA + NB, 1);
weightedIsi = (NA * isiA + NB * isiB) / denom;

weightedDevCap = 0.5;
devA = abs(isiA - weightedIsi);
devB = abs(isiB - weightedIsi);

sizeFactor = 2 * min(NA, NB) / denom;
effBudget  = 1 + (isiBudget - 1) * sizeFactor;
budgetCap  = max(effBudget * weightedIsi, 1e-6);

ok = (isiM <= isiAbs) && ...
     (isiA <= parentCap) && ...
     (isiB <= parentCap) && ...
     (devA <= weightedDevCap) && ...
     (devB <= weightedDevCap) && ...
     (isiM <= budgetCap);
end
