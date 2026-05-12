function postSortReport = kiaSort_post_sort_curate(outputPath, varargin)
%KIASORT_POST_SORT_CURATE  Lightweight CCG-cleaning + similarity-merge step.
%
%   postSortReport = kiaSort_post_sort_curate(outputPath, ...)
%
%   Runs after drift merging on the saved sort. Two phases, both gated
%   by name/value flags so the caller can switch either off:
%
%     1) CCG cleaning  (ccg_cleaning, default true)
%        For each pair of units whose representative channels sit
%        within coincYUm in y, build the cross-correlogram and test
%        for a true peak at lag zero (peak in the lag-0 window is the
%        global max of the CCG and exceeds ccgPeakRatio x median of
%        off-zero bins). When a pair triggers, check that the
%        coincident spikes on the contaminated side have signed lags
%        (relative to the nearest clean-side spike) clustered tightly
%        around their median: at least cleanLagConsistencyFrac of the
%        coincidences must lie within +-cleanLagTightMs of the median
%        lag. Random coincidences spread the lag distribution; a real
%        duplicate produces a sharp sub-ms peak. Only the consistent-
%        lag subset has its spikes rewritten to label -1 in
%        unifiedLabels.h5. Identity is preserved -- only the
%        duplicate spikes are stripped.
%
%     2) Merging       (merging, default true)
%        Mirrors the full gate sequence from the GUI's Auto-curation
%        merge phase. A pair is merged ONLY when ALL gates pass:
%          a) Similarity (XCorr) >= xcorrThreshold (default 0.9):
%             max_half_corr on the shared-channel mean waveforms.
%          b) Amplitude similarity <= thrAmp (default 0.15):
%             same metric the GUI's onSimilarityEstimation builds
%             (mean of |peak diff|/maxAbs and |trough diff|/maxAbs).
%          c) PC distance <= thrPC (default 0.30):
%             projects the shared-channel mean waveforms onto the
%             per-channel PCA stored in sortedSamples{ch}, compares
%             PC1/PC2 distance normalised by score range. Skipped
%             (treated as pass) when the cached basis is missing.
%          d) Merged-ISI multi-gate (full GUI checkMergedISI):
%             merged isiM <= thrIsiAbs; per-parent <= thrIsiAbs *
%             thrIsiBudget; each parent within 0.5 %% of the size-
%             weighted mean; merged <= effBudget * weighted mean.
%        The smaller-spike-count side is absorbed into the larger.
%
%   This step does NOT remove any units. Removal is intentionally
%   left to manual / Auto-curation in the GUI; here we only clean
%   and merge.
%
%   The function tolerates missing inputs and per-pair errors: every
%   block is wrapped in try/catch so a bad single pair never stops
%   the whole pass. unifiedLabels.h5 is written back ONLY when
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
%       postSortReport.nClean   number of CCG-cleanups (spike strips)
%       postSortReport.nMerge   number of pair merges
%       postSortReport.changed  true iff unifiedLabels.h5 was rewritten
%       postSortReport.ok       true on a clean finish

% ---- Parse arguments ----------------------------------------------------
% Default merge thresholds match the GUI Auto-curation panel; only the
% similarity threshold defaults to a tighter 0.9 here per the request
% (the GUI's auto run uses 0.85). Each one can be overridden if you
% want to mirror a specific GUI tuning.
p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addParameter('ccg_cleaning', true,  @(x) islogical(x) || isnumeric(x));
p.addParameter('merging',      true,  @(x) islogical(x) || isnumeric(x));
p.addParameter('xcorrThreshold', 0.9, @(x) isscalar(x) && isnumeric(x));   % thrSim (XCorr metric)
p.addParameter('thrAmp',         0.15, @(x) isscalar(x) && isnumeric(x));  % auto-curate default
p.addParameter('thrPC',          0.30, @(x) isscalar(x) && isnumeric(x));  % auto-curate default
p.addParameter('thrIsiAbs',      1.0,  @(x) isscalar(x) && isnumeric(x));  % max ISI %
p.addParameter('thrIsiBudget',   1.5,  @(x) isscalar(x) && isnumeric(x));  % weighted-mean budget
p.addParameter('thresholdISIms', 1,    @(x) isscalar(x) && isnumeric(x));  % refractory window (ms)
p.addParameter('ccgPeakRatio',   2,   @(x) isscalar(x) && isnumeric(x));
p.addParameter('coincYUm',       100, @(x) isscalar(x) && isnumeric(x));
% Timing-consistency gate for CCG cleaning. After a pair triggers the
% peak-at-zero test, the coincident spikes on the contaminated side
% must have signed lags (relative to the nearest clean-side spike) that
% cluster tightly around their median: at least cleanLagConsistencyFrac
% of them within +-cleanLagTightMs of the median. This blocks stripping
% when the coincidences look random (uniform lag distribution inside
% +-0.5 ms) and only strips spikes whose lag falls inside the tight
% window -- exactly the duplicate-with-consistent-sub-ms-offset case.
p.addParameter('cleanLagTightMs',          0.2, @(x) isscalar(x) && isnumeric(x));
p.addParameter('cleanLagConsistencyFrac',  0.75, @(x) isscalar(x) && isnumeric(x));
p.addParameter('verbose',        false, @(x) islogical(x) || isnumeric(x));
p.parse(outputPath, varargin{:});
opt = p.Results;
opt.ccg_cleaning = logical(opt.ccg_cleaning);
opt.merging      = logical(opt.merging);
opt.verbose      = logical(opt.verbose);

postSortReport = struct('nClean', 0, 'nMerge', 0, 'changed', false, 'ok', false);

% Bail early if neither phase is enabled.
if ~opt.ccg_cleaning && ~opt.merging
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

% Sampling frequency: pulled from the first non-empty sortedSamples cfg.
fs = 30000;
for i = 1:numel(sortedSamples)
    if ~isempty(sortedSamples{i}) && isfield(sortedSamples{i}, 'cfg') ...
            && isfield(sortedSamples{i}.cfg, 'samplingFrequency')
        fs = sortedSamples{i}.cfg.samplingFrequency;
        break;
    end
end

% ---- Load drift label remap if present ---------------------------------
% After drift merge, unifiedLabels.h5 carries POST-drift labels but the
% unif struct (from sorted_samples.mat) was built from PRE-drift labels.
% labelRemap maps pre-drift -> post-drift. If the report is absent we
% fall back to identity, which is correct when drift merge was skipped
% or made no changes.
labelRemap = containers.Map('KeyType', 'double', 'ValueType', 'double');
if exist(driftReportPath, 'file')
    try
        dr = load(driftReportPath, 'driftReport');
        if isfield(dr, 'driftReport') && isfield(dr.driftReport, 'labelRemap') ...
                && isa(dr.driftReport.labelRemap, 'containers.Map')
            labelRemap = dr.driftReport.labelRemap;
        end
    catch
        % silently fall back to identity
    end
end

% ---- Build per-unit info keyed by post-drift label ---------------------
% For each unif entry compute its post-drift label, then collapse to one
% representative per unique label (highest detectability wins). The
% representative gives us a (channel, localLabel) pair we can feed into
% sortedSamples to recover a mean waveform.
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

% ---- Phase 1: CCG cleaning ---------------------------------------------
% Per-pair CCG peak-at-zero test. When triggered, strip the contaminated
% (higher-ISI) side's coincident rows by setting their lbl_all entries
% to -1. lbl_all is mutated in place; subsequent pairs see the cleaned
% counts.
nClean = 0;
if opt.ccg_cleaning
    coincSamples = round(0.5e-3 * fs);
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

            [ccg, zb] = local_pairCCG(spk_i, spk_j, ccgMaxLag, ccgBin);
            if ~local_ccgPeakAtZero(ccg, zb, ccgZeroTol, opt.ccgPeakRatio)
                continue;
            end

            % Coincidence mask in i's frame.
            d_ij = local_nearest_distance(spk_i, spk_j);
            coinc_i = d_ij <= coincSamples;

            % Decide contaminated side: higher ISI is contaminated;
            % ties broken by smaller spike count.
            isi_i = unitInfo(ki).isiViol; if isnan(isi_i), isi_i = 0; end
            isi_j = unitInfo(kj).isiViol; if isnan(isi_j), isi_j = 0; end
            if isi_i > isi_j || (isi_i == isi_j && numel(spk_i) <= numel(spk_j))
                kCont   = ki;
                contRows = kiLive;
                contMask = coinc_i;
                signedLag = local_signed_nearest_lag(spk_i, spk_j);
            else
                kCont   = kj;
                contRows = kjLive;
                d_ji    = local_nearest_distance(spk_j, spk_i);
                contMask = d_ji <= coincSamples;
                signedLag = local_signed_nearest_lag(spk_j, spk_i);
            end

            if ~any(contMask), continue; end

            % --- Timing-consistency gate -----------------------------
            % The coincident spikes' signed lags to the nearest clean
            % spike must cluster tightly around their median. A real
            % duplicate (same physical spike picked up twice with a
            % small sorter-induced offset) gives a sharp peak; random
            % coincidence gives a broad / uniform distribution.
            tightSamples = max(1, round(opt.cleanLagTightMs * 1e-3 * fs));
            coincIdx     = find(contMask);
            coincLags    = signedLag(coincIdx);
            coincLags    = coincLags(~isnan(coincLags));
            if numel(coincLags) < 5
                % Too few to assess consistency reliably; skip strip.
                continue;
            end
            medLag         = median(coincLags);
            withinTight    = abs(coincLags - medLag) <= tightSamples;
            consistencyFrac = sum(withinTight) / numel(coincLags);
            if consistencyFrac < opt.cleanLagConsistencyFrac
                if opt.verbose
                    fprintf(['Post-sort CCG-clean: pair (%d,%d) skipped, ' ...
                        'lag consistency %.2f < %.2f.\n'], ...
                        unitInfo(ki).label, unitInfo(kj).label, ...
                        consistencyFrac, opt.cleanLagConsistencyFrac);
                end
                continue;
            end
            % Restrict the strip mask to the consistent-lag subset.
            keepCoinc            = false(numel(signedLag), 1);
            keepCoinc(coincIdx)  = withinTight;
            contMask             = contMask & keepCoinc;
            if ~any(contMask), continue; end

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
                fprintf('Post-sort CCG-clean: pair (%d,%d) failed: %s\n', ...
                    unitInfo(ki).label, unitInfo(kj).label, ME.message);
            end
        end
    end
end

% ---- Phase 2: Merging --------------------------------------------------
% Mirrors the GUI Auto-curation merge stage: a candidate pair must clear
% ALL of these gates before it is merged.
%
%   1) Similarity (XCorr).
%      max_half_corr on the shared-channel mean waveforms (channel of
%      unit i). Must reach >= opt.xcorrThreshold (default 0.9).
%
%   2) Amplitude similarity.
%      Same metric the GUI's onSimilarityEstimation builds: average of
%      |peak-to-peak diffs| normalised by the larger of the two
%      maxAbs values. Must be <= opt.thrAmp.
%
%   3) PC distance.
%      Project both mean waveforms onto the per-channel PCA stored in
%      sortedSamples{ch}.clusteringInfo.PCA, then compare PC1/PC2.
%      Must be <= opt.thrPC. When the cached PCA is missing or invalid
%      this gate is skipped (treated as pass) -- same conservative
%      fallback the GUI uses when the basis is unavailable.
%
%   4) Merged ISI gates (full local_checkMergedISI logic).
%      Runs against the CURRENT spike trains for each label so that
%      Phase 1 strips and any earlier Phase 2 merges feed forward.
%
% Spike-time alignment: when the similarity metric finds a non-zero
% best-lag, the absorbed unit's spike times are shifted by -bestLag so
% the merged ACG's lag-0 bin really sits at lag zero (matches GUI).
nMerge = 0;
% Working remap from "label as it was" -> "label as it is now after
% any merges already applied this run". Allows transitive collapse.
mergedTo = containers.Map('KeyType', 'double', 'ValueType', 'double');

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
            [simScore, ~] = max_half_corr(mw1(:)', mw2(:)', 1, M2, maxLag, 0);
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
            else
                primaryLab  = labJ;
                absorbedLab = labI;
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
changed = (nClean > 0) || (nMerge > 0);
if changed
    try
        if exist(unifiedLabelsH5, 'file')
            delete(unifiedLabelsH5);
        end
        h5create(unifiedLabelsH5, '/unifiedLabels', size(lbl_all), 'Datatype', 'double');
        h5write(unifiedLabelsH5,  '/unifiedLabels', lbl_all);
    catch ME
        if opt.verbose
            fprintf('Post-sort curate: H5 write failed (%s).\n', ME.message);
        end
        postSortReport.nClean  = nClean;
        postSortReport.nMerge  = nMerge;
        postSortReport.changed = false;
        postSortReport.ok      = false;
        return;
    end
end

postSortReport.nClean  = nClean;
postSortReport.nMerge  = nMerge;
postSortReport.changed = changed;
postSortReport.ok      = true;

if opt.verbose
    fprintf('Post-sort curate: %d CCG strips, %d merges (changed=%d).\n', ...
        nClean, nMerge, changed);
end

end


% =========================================================================
% Local helpers
% =========================================================================

function [ccg, zeroBinIdx] = local_pairCCG(spkA, spkB, maxLagSamples, binSamples)
% Two-pointer pairwise-difference histogram for the CCG of (spkB - spkA)
% inside +-maxLagSamples, binned at binSamples. Same logic as the curate
% GUI's pairCCG -- inlined here so this file has no GUI dependency.
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
% Lag-0 peak test, stricter than the curate GUI's median-based version.
% A pair triggers ONLY when both:
%   * The lag-0 window holds the global max of the CCG (so the
%     biggest bin really sits at zero, not somewhere else).
%   * That peak is > ratio x the 90th-percentile of the off-zero
%     bins. The 90th percentile means the peak has to clear what
%     90 % of the rest of the CCG looks like, so a single tall
%     non-zero bin can't push the baseline up the way it can with
%     a mean, and structural temporal coupling (real biology with
%     wide lag bands) gets a much fairer comparison than against a
%     median that's usually 0 or 1 in sparse CCGs.
tf = false;
if isempty(ccg) || all(ccg == 0), return; end
nB = numel(ccg);
peakRange = max(1, zeroBinIdx-zeroTolBins) : min(nB, zeroBinIdx+zeroTolBins);
if isempty(peakRange), return; end
peakValue = max(ccg(peakRange));
offRange  = ccg;
offRange(peakRange) = [];
if isempty(offRange), return; end
% Off-window must NOT contain a strictly larger bin -- if it does,
% the actual global max sits away from zero and this isn't a "true
% peak at zero" pair.
if any(offRange > peakValue), return; end
if ~any(offRange > 0), return; end
baseline = prctile(offRange, 90);
if ~isfinite(baseline) || baseline <= 0
    % Sparse CCG: 90 % of bins are zero -> percentile is 0. Fall back
    % to the mean of the non-zero off-window bins so the peak still
    % has something meaningful to clear (and we don't trivially pass
    % every non-zero peak).
    nz = offRange(offRange > 0);
    if isempty(nz), return; end
    baseline = mean(nz);
    if baseline <= 0, return; end
end
tf = peakValue > ratio * baseline;
end


function d = local_nearest_distance(A, B)
% For each element of A, distance (samples) to its nearest neighbour in
% B. Two-pointer sweep on sorted inputs -- O(nA + nB). Returns +inf for
% A entries with no neighbour in B (so the caller's "<= tol" check is
% safe).
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
% For each element of A, SIGNED lag (samples) to its nearest neighbour
% in B: lag(i) = B_nearest - A(i). Positive = B comes after A,
% negative = B comes before A. NaN when B is empty or no neighbour
% within finite range. Used by the cleaning timing-consistency gate
% so we can tell whether the coincident spikes cluster around a
% specific sub-ms offset (real duplicate) or are scattered (random).
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
% Walk the merged-to chain to find the current top-level label for lab0.
% Used so transitive merges (A->B, B->C ⇒ A->C) collapse cleanly within
% a single Phase 2 pass without rewriting the whole label vector each
% step.
lab = lab0;
while isKey(mergedTo, lab)
    next = mergedTo(lab);
    if next == lab, break; end
    lab = next;
end
end


function row = local_rowOnChannel(wf, homeCh, targetCh)
% Pull the local row of a (nLocalChannels x T) mean-waveform tensor
% that corresponds to GLOBAL channel targetCh, given the waveform is
% centered on homeCh. Returns [] when the target falls outside the
% unit's local footprint.
row = [];
if isempty(wf) || ~ismatrix(wf), return; end
nLoc  = size(wf, 1);
half  = floor((nLoc - 1) / 2);
% Local row index 1 corresponds to global channel (homeCh - half).
localIdx = targetCh - (homeCh - half);
if localIdx < 1 || localIdx > nLoc, return; end
row = wf(localIdx, :);
end


function ampDiff = local_ampSimilarity(mw1, mw2)
% Mirrors the amplitude-similarity metric used in the GUI's
% onSimilarityEstimation: average of |peak diff|/maxAbs and
% |trough diff|/maxAbs across the pair. Smaller is more similar.
% Returns NaN when either trace is flat.
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
% Shared-channel PC1/PC2 distance between mean waveforms, normalised
% by the basis's score range. Mirrors the GUI's checkPCDistance, with
% the difference that we use the per-channel PCA cached in
% sortedSamples{chI}.clusteringInfo.PCA instead of recomputing it
% from raw spikes (the post-sort step deliberately does not re-read
% the input data file). When the cached basis is missing or
% malformed we conservatively pass.
ok = true;
if isnan(chI), return; end
% Pull rows on chI for both units.
mwA = local_rowOnChannel(wfI, chI, chI);
mwB = local_rowOnChannel(wfJ, chJ, chI);
if isempty(mwA) || isempty(mwB), return; end
if numel(mwA) ~= numel(mwB), return; end

% Cached PCA: load on first hit, reuse thereafter.
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
        % leave PC empty; gate will pass
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
% Same multi-condition merged-ISI gate the GUI's checkMergedISI uses:
%   1) Hard cap on the merged train: isiM <= isiAbs.
%   2) Per-parent cap: isiA, isiB <= isiAbs * max(isiBudget, 1).
%   3) Each parent's ISI within 0.5 % of the size-weighted mean.
%   4) Merged ISI <= effBudget * weighted mean, where effBudget
%      shrinks toward 1.0 as size disparity grows.
%
% spkA, spkB are sample-frame spike-time vectors. fs is the sampling
% rate in Hz, threshMs is the refractory threshold in ms passed
% straight to getISIViolations.
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
