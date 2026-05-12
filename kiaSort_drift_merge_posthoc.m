function driftReport = kiaSort_drift_merge_posthoc(outputPath, varargin)
% Post-hoc drift-aware merging step for KiaSort outputs.
%
% Single-pass pipeline with hard AND gates:
%
% For each geometrically-close pair (same x-column, |dy|<=maxDriftUm):
%   1. At least one unit must have a transient rate dip.
%   2. Amplitude ratio within ampRatioRange.
%   3. Waveform correlation at best drift-aligned y-shift >= corrThreshold.
%   4. Smoothed-density correlation <= densityAntiCorr (restricted to the
%      active support so mutually-silent tails don't dilute the signal).
%   5. Merged-train ISI violation <= isiThresholdPct.
%
% A final cohort-consistency filter rejects pairs whose drift-shift
% disagrees with the sign or magnitude of their local neighbours.
%
% Surviving pairs form groups via connected components; up to triples
% are allowed by default.
%
% INPUT FILE LAYOUT (relative to outputPath)
%   RES_Samples/channel_info.mat
%   Sorted_Samples/sorted_samples.mat   (sortedSamples, crossChannelStats)
%   Sorted_Samples/sample_features.mat  (sampleFeatures)
%   RES_Sorted/spike_idx.h5
%   RES_Sorted/unifiedLabels.h5
%   RES_Sorted/channelNum.h5
%
% OUTPUT FILES (written into outputPath/RES_Sorted/)
%   unifiedLabels_merged.h5     remapped unified labels (one entry per spike)
%   drift_merge_report.mat      candidate pairs, validated groups, label map
%
% NAME/VALUE PAIRS
%   'cfg'               (default [])     override cfg struct (samplingFrequency,
%                                        spikeDistance). If empty, cfg is pulled
%                                        from the first non-empty sortedSamples{ch}.cfg.
%   'maxDriftUm'        (default 80)      max COM y-distance (um)
%   'sameColumnXTol'    (default 8)       max COM x-distance (um)
%   'corrThreshold'     (default 0.85)    min waveform corr at best y-shift
%   'ampRatioRange'     (default [0.1 5]) allowed primary-channel amp ratio
%   'shiftStepUm'       (default 5)       drift-search granularity
%   'densityBinSec'     (default 60)      bin width for density
%   'densitySmoothSec'  (default 180)     Gaussian sigma for density (s)
%   'densityAntiCorr'   (default -0.2)    required r between member densities
%   'densityActiveOnly' (default true)    restrict corr to active support
%   'densityActiveFrac' (default 0.1)     a bin is "active" if rate > frac*peak
%   'dipFrac'           (default 0.2)     a unit dips when rate < dipFrac*peak
%   'minDipSec'         (default 60)      ... for at least this many seconds
%   'requireDip'        (default true)    dip gate as hard filter
%   'cohortRadiusUm'    (default 50)      radius for neighbour lookup (um)
%   'cohortDyTolUm'     (default 25)      |dy - median(nb dy)| tolerance (um)
%   'minCohortNeighbors'(default 2)       neighbours needed to run the test
%   'requireCohort'     (default true)    cohort consistency as hard filter
%   'isiThresholdPct'   (default 1.0)     max % ISI < refractoryMs (merged)
%   'refractoryMs'      (default 2)       refractory window (ms)
%   'minOverlapCh'      (default 4)       min overlapping global channels
%   'allowTriples'      (default true)    allow 3-member merge groups
%   'maxGroupSize'      (default 3)       max members per merge group
%   'overwriteOriginal' (default false)   overwrite canonical unifiedLabels.h5
%   'dryRun'            (default false)   skip HDF5 write
%   'verbose'           (default true)    print counters
%
% NOTES
% * The post-hoc merge is conservative by design. Cross-column merges are
%   forbidden. Each candidate pair must satisfy ALL four criteria. Triples
%   require a monotonic y-trajectory and joint density consistency (sum of
%   member densities must be more stable than the mean of individual CVs).
%   Transitive chaining (A~B, B~C => {A,B,C} regardless of A~C) is not
%   permitted.
% * The script does not modify any of the channel-wise .mat files. Only
%   unifiedLabels_merged.h5 (and optionally unifiedLabels.h5) is written.

% ---- Parse arguments ----------------------------------------------------
p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addParameter('cfg', [], @(x) isempty(x) || isstruct(x));
p.addParameter('maxDriftUm', 80, @isscalar);
p.addParameter('sameColumnXTol', 25, @isscalar);
p.addParameter('corrThreshold', 0.8, @isscalar);
p.addParameter('ampRatioMax', 1.5, @isscalar);         % relaxed (drifted units can lose up to ~5x on the edge channel)
p.addParameter('isiThresholdPct', 1.0, @isscalar);
p.addParameter('refractoryMs', 2, @isscalar);
p.addParameter('densityAntiCorr', -0.5, @isscalar);     % SOFT gate; coverage-gain does the real work
p.addParameter('densityBinSec', 5, @isscalar);         % 5s bins
p.addParameter('densitySmoothSec', 90, @isscalar);     % Gaussian sigma for smoothing (s)
p.addParameter('densityActiveOnly', false, @islogical); % restrict corr to active support
p.addParameter('densityActiveFrac', 0.1, @isscalar);   % bin counted active if > frac*max
% Per-unit transient drop-out (dip) required for drift candidacy.
p.addParameter('dipFrac', 0.5, @isscalar);             % rate must fall below dipFrac * max
p.addParameter('minDipSec', 60, @isscalar);            % for at least this many seconds
p.addParameter('requireDip', true, @islogical);        % hard gate on the pair
% Waveform drift estimation (full y-shift scan, snap B channels to nearest real channel).
p.addParameter('shiftStepUm', 5, @isscalar);           % scan step in um
p.addParameter('driftGainEps', -0.2, @isscalar);       % shifted corr must beat zero-shift corr by this
% Coverage-gain (replaces CV-ratio stationarity test).
%   For drift: the merged train fires in a LARGER fraction of the
%   recording than either alone because A's silent epochs are B's active
%   epochs. cov(A|B) - max(cov(A), cov(B)) >= coverageGainMin.
p.addParameter('coverageGainMin', 0.25, @isscalar);
p.addParameter('requireCoverage', true, @islogical);
% ACG (autocorrelogram) similarity: same neuron observed twice should
% have near-identical ACGs. Disable by setting acgCorrThreshold <= -1.
p.addParameter('acgCorrThreshold', -1, @isscalar);     % min Pearson corr of ACG_A vs ACG_B
p.addParameter('acgWindowMs',     50, @isscalar);      % +-windowMs lag range for ACG
p.addParameter('acgBinMs',         1, @isscalar);      % ACG bin width (ms)
p.addParameter('acgMinSpikes',    20, @isscalar);      % skip gate when either unit fires fewer spikes
% Cohort drift-field (fits over waveform-estimated dy_est values).
p.addParameter('cohortRadiusUm', 50, @isscalar);       % radius for cohort lookup (um)
p.addParameter('cohortDyTolUm', 20, @isscalar);        % max |dy_est - dy_hat(y_mid)| (um)
p.addParameter('requireCohort', true, @islogical);     % enforce cohort consistency
p.addParameter('minCohortNeighbors', 2, @isscalar);
p.addParameter('minOverlapCh', 4, @isscalar);
p.addParameter('allowTriples', true, @islogical);
p.addParameter('maxGroupSize', 6, @isscalar);
p.addParameter('overwriteOriginal', false, @islogical);
p.addParameter('dryRun', false, @islogical);
p.addParameter('verbose', true, @islogical);
p.addParameter('debugFigs', true, @islogical);       % TEMP: write diagnostic figs
p.addParameter('figFolder', '', @(x) ischar(x) || isstring(x));
p.addParameter('nFigsPass', 10, @isscalar);          % # top-corr figures among PASS pairs
p.addParameter('nFigsReject', 10, @isscalar);        % # top-corr figures among REJECT pairs
p.parse(outputPath, varargin{:});
opt = p.Results;

outputPath        = char(outputPath);
samplesFolder     = fullfile(outputPath, 'RES_Samples');
sortedSampFolder  = fullfile(outputPath, 'Sorted_Samples');
resSortedFolder   = fullfile(outputPath, 'RES_Sorted');

channelInfoPath   = fullfile(samplesFolder,    'channel_info.mat');
sortedSamplesPath = fullfile(sortedSampFolder, 'sorted_samples.mat');
sampleFeatsPath   = fullfile(sortedSampFolder, 'sample_features.mat');
spikeIdxH5        = fullfile(resSortedFolder,  'spike_idx.h5');
unifiedLabelsH5   = fullfile(resSortedFolder,  'unifiedLabels.h5');
channelNumH5      = fullfile(resSortedFolder,  'channelNum.h5');

requiredFiles = {channelInfoPath, sortedSamplesPath, sampleFeatsPath, ...
                 spikeIdxH5, unifiedLabelsH5, channelNumH5};
for i = 1:numel(requiredFiles)
    if ~exist(requiredFiles{i}, 'file')
        error('kiaSort_drift_merge_posthoc:missingInput', ...
            'Required file not found: %s', requiredFiles{i});
    end
end

% ---- Load inputs --------------------------------------------------------
if opt.verbose, fprintf('Loading channel info...\n'); end
chInfo = load(channelInfoPath, 'channel_locations', 'channel_inclusion', 'num_samples');
channel_locations = chInfo.channel_locations;
if isempty(channel_locations) || size(channel_locations, 2) < 2
    error('channel_locations must be (numChannels x 3) with x in col 1, y in col 2.');
end
if isfield(chInfo, 'num_samples') && ~isempty(chInfo.num_samples)
    num_samples = double(chInfo.num_samples);
else
    num_samples = [];
end

if opt.verbose, fprintf('Loading sorted samples...\n'); end
S = load(sortedSamplesPath, 'sortedSamples', 'crossChannelStats');
sortedSamples     = S.sortedSamples;
crossChannelStats = S.crossChannelStats;

if opt.verbose, fprintf('Loading sample features...\n'); end
SF = load(sampleFeatsPath, 'sampleFeatures');
sampleFeatures = SF.sampleFeatures;

% Resolve cfg: explicit override wins, otherwise pull from the first
% non-empty sortedSamples{ch}.cfg (saved by kiaSort_cluster_classify_Temp).
cfg = opt.cfg;
if isempty(cfg)
    for ch = 1:numel(sortedSamples)
        if ~isempty(sortedSamples{ch}) && isfield(sortedSamples{ch}, 'cfg')
            cfg = sortedSamples{ch}.cfg;
            break;
        end
    end
end
if isempty(cfg) || ~isstruct(cfg)
    error('kiaSort_drift_merge_posthoc:noCfg', ...
        'cfg not found. Pass via ''cfg'' or ensure sortedSamples{ch}.cfg is populated.');
end
if ~isfield(cfg, 'samplingFrequency') || isempty(cfg.samplingFrequency)
    error('cfg.samplingFrequency is required.');
end
fs = cfg.samplingFrequency;

if isfield(cfg, 'spikeDistance') && ~isempty(cfg.spikeDistance)
    jitter_gap = cfg.spikeDistance * fs / 1000;  % samples
else
    jitter_gap = 0.375 * fs / 1000;              % fallback to default from configs
end

if opt.verbose, fprintf('Reading final HDF5 outputs...\n'); end
spk_idx_all     = h5read(spikeIdxH5,      '/spike_idx');
unifiedLab_all  = h5read(unifiedLabelsH5, '/unifiedLabels');
spk_idx_all     = double(spk_idx_all(:));
unifiedLab_all  = double(unifiedLab_all(:));

if numel(spk_idx_all) ~= numel(unifiedLab_all)
    error('spike_idx and unifiedLabels HDF5 lengths disagree (%d vs %d).', ...
          numel(spk_idx_all), numel(unifiedLab_all));
end

% ---- Build per-global-unit info from crossChannelStats.unified_labels --
unif = crossChannelStats.unified_labels;
nUnits = numel(unif.label);
if opt.verbose, fprintf('Indexing %d global units...\n', nUnits); end

if nUnits == 0
    driftReport = struct('opt', opt, 'pairs', [], ...
        'groups', {{}}, 'labelRemap', containers.Map('KeyType','double','ValueType','double'), ...
        'primarySummary', zeros(0,4), 'nUnitsBefore', 0, 'nUnitsAfter', 0);
    if opt.verbose, fprintf('No global units present; nothing to merge.\n'); end
    return;
end

unitInfo = struct( ...
    'globalLabel',    cell(1, nUnits), ...
    'channel',        cell(1, nUnits), ...
    'localLabel',     cell(1, nUnits), ...
    'meanWF',         cell(1, nUnits), ...
    'tSet',           cell(1, nUnits), ...
    'sampleCount',    cell(1, nUnits), ...
    'detectability',  cell(1, nUnits), ...
    'spkIdxFull',     cell(1, nUnits), ...
    'density',        cell(1, nUnits), ...
    'com',            cell(1, nUnits), ...    % [x, y] center of mass (um)
    'spatialExtent',  cell(1, nUnits), ...    % sqrt(Var(xy)) weighted
    'peakAmp',        cell(1, nUnits), ...    % max per-ch peak-to-peak
    'localChGlobal',  cell(1, nUnits), ...    % global ch index for each local row
    'perChanAmp',     cell(1, nUnits), ...
    'hasDip',         cell(1, nUnits), ...    % transient rate drop (drift candidacy)
    'dipDepth',       cell(1, nUnits), ...    % min rate / max rate
    'dipDurSec',      cell(1, nUnits), ...    % longest sub-threshold run (s)
    'activeFrac',     cell(1, nUnits), ...    % fraction of bins with rate > activeFrac*max
    'activeMask',     cell(1, nUnits), ...    % boolean mask for coverage computation
    'acg',            cell(1, nUnits));       % autocorrelogram (binned, zero lag masked)

% Find original-template index for each unit so we can pull templateWaveforms
for u = 1:nUnits
    ch  = unif.channelID(u);
    lbl = unif.labelInChannel(u);

    unitInfo(u).globalLabel   = unif.label(u);
    unitInfo(u).channel       = ch;
    unitInfo(u).localLabel    = lbl;
    unitInfo(u).detectability = unif.detectblity(u);

    if isempty(sortedSamples{ch}), continue; end
    rel = sortedSamples{ch}.clusteringInfo.clusterRelabeling;

    keptIdx = find(rel.newUniqueLabels == lbl, 1);
    if isempty(keptIdx), continue; end

    unitInfo(u).meanWF = squeeze(rel.newMeanWaveforms(keptIdx, :, :));
    if isfield(rel, 'newSampleCounts') && numel(rel.newSampleCounts) >= keptIdx
        unitInfo(u).sampleCount = rel.newSampleCounts(keptIdx);
    else
        unitInfo(u).sampleCount = sum(sampleFeatures{ch}.updatedLabels == lbl);
    end

    % Template variants from sort-sample phase (preferred)
    tWF = sortedSamples{ch}.waveformInfo.templateWaveforms;
    origIdx = find(rel.newLabels == lbl);
    t = [];
    for oi = 1:numel(origIdx)
        tj = squeeze(tWF(origIdx(oi), :, :, :));
        if ndims(tj) < 3 && size(tWF, 3) == 1
            tj = reshape(tj, size(tj,1), 1, size(tj,2));
        end
        t = cat(1, t, tj);
    end
    if isempty(t) || all(t(:) == 0)
        mw = unitInfo(u).meanWF;
        t = reshape(mw, 1, size(mw,1), size(mw,2));
    end
    unitInfo(u).tSet = t;
end

% ---- Compute spatial center-of-mass per unit ---------------------------
% For each unit, map each local waveform row to its global channel index,
% then compute the amplitude-weighted (squared) xy center of mass and
% spatial extent. This operates in physical (um) coordinates and is
% insensitive to the interleaved Neuropixels channel ordering.
numGlobalCh  = size(channel_locations, 1);
N_halfWindow = cfg.num_channel_extract;
for u = 1:nUnits
    wf = unitInfo(u).meanWF;
    if isempty(wf)
        unitInfo(u).com = [NaN NaN];
        continue;
    end
    nLoc   = size(wf, 1);
    homeCh = unitInfo(u).channel;

    % Local rows map to global channels (home - N : home + N), clipped
    % to the valid channel range.
    localToGlobal = (homeCh - N_halfWindow) : (homeCh + N_halfWindow);
    valid         = (localToGlobal >= 1) & (localToGlobal <= numGlobalCh);
    chG = nan(nLoc, 1);
    if numel(localToGlobal) == nLoc
        chG(valid) = localToGlobal(valid);
    else
        % Defensive: if extraction clipped near an edge, align from center.
        midLoc = floor(nLoc/2) + 1;
        for rr = 1:nLoc
            gi = homeCh + (rr - midLoc);
            if gi >= 1 && gi <= numGlobalCh, chG(rr) = gi; end
        end
    end

    % Per-channel peak-to-peak amplitude (from template variants when present)
    tMean = unitInfo(u).meanWF; %squeeze(mean(unitInfo(u).tSet, 1));
    if isvector(tMean), tMean = tMean(:)'; end
    a_uc = max(tMean, [], 2) - min(tMean, [], 2);
    a_uc = a_uc(:);
    a_uc(isnan(chG)) = 0;

    vIdx = find(a_uc > 0 & ~isnan(chG));
    if isempty(vIdx)
        unitInfo(u).com           = [NaN NaN];
        unitInfo(u).spatialExtent = NaN;
        unitInfo(u).peakAmp       = 0;
    else
        xy = channel_locations(chG(vIdx), 1:2);
        w  = a_uc(vIdx).^2;
        com = (w' * xy) / sum(w);
        dxy = xy - com;
        varXY = (w' * (dxy.^2)) / sum(w);
        unitInfo(u).com           = com;
        unitInfo(u).spatialExtent = sqrt(sum(varXY));
        unitInfo(u).peakAmp       = max(a_uc);
    end
    unitInfo(u).localChGlobal = chG;
    unitInfo(u).perChanAmp    = a_uc;
end

% Attach full spike trains from HDF5
if opt.verbose, fprintf('Attaching full spike trains...\n'); end
for u = 1:nUnits
    mask = (unifiedLab_all == unitInfo(u).globalLabel);
    unitInfo(u).spkIdxFull = sort(spk_idx_all(mask));
end

% Precompute per-unit ACG so the pair loop doesn't recompute for every
% partner (each unit appears in O(nUnits) pairs). Empty for units that
% have too few spikes to produce a meaningful ACG.
if opt.acgCorrThreshold > -1
    if opt.verbose
        fprintf('Computing per-unit autocorrelograms (+-%d ms, %d ms bins)...\n', ...
                opt.acgWindowMs, opt.acgBinMs);
    end
    for u = 1:nUnits
        spk = unitInfo(u).spkIdxFull;
        if numel(spk) < max(2, opt.acgMinSpikes)
            unitInfo(u).acg = [];
        else
            unitInfo(u).acg = computeACG(spk, fs, opt.acgWindowMs, opt.acgBinMs);
        end
    end
end

% Compute density vectors over the full recording: binned counts, then
% Gaussian-smoothed (drift is slow relative to firing-rate fluctuations).
binSamples = max(1, round(opt.densityBinSec * fs));
if ~isempty(num_samples) && num_samples > 0
    totalSamples = num_samples;
else
    totalSamples = max([1; spk_idx_all]);
end
nBins = max(1, ceil(totalSamples / binSamples));

% Gaussian kernel over bins; sigma = densitySmoothSec
sigmaBins = max(1, opt.densitySmoothSec / opt.densityBinSec);
halfLen = ceil(3 * sigmaBins);
tK = -halfLen:halfLen;
kernel = exp(-0.5 * (tK / sigmaBins).^2);
kernel = kernel / sum(kernel);

% Edge-normalization: divide by the kernel mass that fell within the
% recording, so rate estimates near t=0 and t=end are not artificially
% damped by the Gaussian running off the edge.
edgeNorm = conv(ones(1, nBins), kernel, 'same');
edgeNorm(edgeNorm < 1e-6) = 1e-6;

if opt.verbose
    fprintf('Computing per-unit density (%d bins of %d s, sigma %d s)...\n', ...
            nBins, opt.densityBinSec, opt.densitySmoothSec);
end
for u = 1:nUnits
    if isempty(unitInfo(u).spkIdxFull)
        unitInfo(u).density = zeros(1, nBins);
        continue;
    end
    bins = min(nBins, max(1, ceil(unitInfo(u).spkIdxFull / binSamples)));
    d = accumarray(bins(:), 1, [nBins, 1]);
    dS = conv(d(:)', kernel, 'same') ./ edgeNorm;   % spikes per bin, edge-corrected
    unitInfo(u).density = dS / opt.densityBinSec;    % convert to Hz
end

% ---- Per-unit transient dip detection ----------------------------------
% A drifting neuron must show a window where its rate drops below
% dipFrac * peak for at least minDipSec seconds. Units that fire
% steadily throughout are not drift candidates.
for u = 1:nUnits
    r = unitInfo(u).density(:)';
    if isempty(r) || all(r == 0)
        unitInfo(u).hasDip     = false;
        unitInfo(u).dipDepth   = NaN;
        unitInfo(u).dipDurSec  = 0;
        unitInfo(u).activeFrac = 0;
        unitInfo(u).activeMask = false(1, nBins);
        continue;
    end
    rMax = max(r);
    if rMax <= 0
        unitInfo(u).hasDip     = false;
        unitInfo(u).dipDepth   = NaN;
        unitInfo(u).dipDurSec  = 0;
        unitInfo(u).activeFrac = 0;
        unitInfo(u).activeMask = false(1, numel(r));
        continue;
    end
    below = r < (opt.dipFrac * rMax);
    above = r >= (opt.dipFrac * rMax);
    [runLen, ~] = longestRun(below);
    dipDurSec = runLen * opt.densityBinSec;
    unitInfo(u).hasDip     = dipDurSec >= opt.minDipSec;
    unitInfo(u).dipDepth   = min(r) / rMax;
    unitInfo(u).dipDurSec  = dipDurSec;
    unitInfo(u).activeMask = above;
    unitInfo(u).activeFrac = mean(above);
end
nDipOK = sum(arrayfun(@(u) unitInfo(u).hasDip, 1:nUnits));
if opt.verbose
    fprintf('Units with transient firing dip (>= %gs below %.0f%% peak): %d / %d\n', ...
            opt.minDipSec, 100 * opt.dipFrac, nDipOK, nUnits);
end

% ---- Single-pass pair evaluation with OR-logic temporal evidence ------
%
% For each geometrically-close pair:
%   1. Amplitude ratio <= ampRatioMax (very relaxed)
%   2. Waveform correlation at best shift >= corrThreshold  (main FP defence)
%   3. At least ONE temporal signature of drift is present:
%        (a) at least one unit dips
%        (b) coverage gain >= covGainMin
%        (c) activity-mask corr <= maskCorrMax
%        (d) smoothed-density corr <= densCorrMax
%   4. Merged-train ISI violation <= isiThresholdPct
%
% A final defensive cohort outlier removal drops any pair whose drift
% estimate is > cohortDyTolUm from the local median of kept pairs
% (only when the local neighbourhood is dense enough to be trusted).

if opt.verbose, fprintf('Scanning candidate pairs...\n'); end

pitchY = inferYPitch(channel_locations);
yStep  = inferYStep(channel_locations);
if opt.verbose
    fprintf('Inferred y-pitch (within column) = %.2f um, y-step (row-to-row) = %.2f um\n', ...
            pitchY, yStep);
end

coms     = nan(nUnits, 2);
peakAmps = zeros(nUnits, 1);
for u = 1:nUnits
    if ~isempty(unitInfo(u).com), coms(u, :) = unitInfo(u).com; end
    peakAmps(u) = unitInfo(u).peakAmp;
end

pairs = struct('uA', {}, 'uB', {}, 'corr', {}, 'corrRaw', {}, 'dyEst', {}, ...
               'dx', {}, 'dy', {}, 'ampRatio', {}, ...
               'densCorr', {}, 'isiViol', {}, 'nOverlap', {}, ...
               'acgCorr', {}, ...
               'covA', {}, 'covB', {}, 'covM', {}, 'covGain', {});

diag.nGeomPairs  = 0;
diag.nEvaluated  = 0;
diag.rejAmp      = 0;
diag.rejDip      = 0;
diag.rejOverlap  = 0;
diag.rejCorr     = 0;
diag.rejDens     = 0;
diag.rejCoverage = 0;
diag.rejISI      = 0;
diag.rejCohort   = 0;
diag.rejACG      = 0;
diag.evidence    = struct('uA',{},'uB',{},'corr',{},'corrRaw',{},'dyEst',{}, ...
                          'dx',{},'dy',{}, ...
                          'ampRatio',{},'densCorr',{},'isiViol',{}, ...
                          'acgCorr',{}, ...
                          'nOverlap',{},'reject',{}, ...
                          'covA',{},'covB',{},'covM',{},'covGain',{});

% Infer geometry once (used as snap tolerance inside corrWithYShiftUm).
pitchY = inferYPitch(channel_locations);
if ~(pitchY > 0), pitchY = 20; end
yStep  = inferYStep(channel_locations);
if ~(yStep > 0),  yStep  = pitchY; end

for uA0 = 1:nUnits
    if any(isnan(coms(uA0, :))), continue; end
    tA0 = unitInfo(uA0).tSet;
    gA0 = unitInfo(uA0).localChGlobal;
    if isempty(tA0) || all(isnan(gA0)), continue; end

    for uB0 = (uA0 + 1):nUnits
        if any(isnan(coms(uB0, :))), continue; end

        % ---- Canonicalize pair by y-position -----------------------
        % Unit indices are arbitrary sorter labels and carry no
        % y-information. To make the sign of dy (and of the COM-
        % aligning shift dyEst) reflect drift DIRECTION physically --
        % not the accidental unit ordering -- always define uA as the
        % lower-y unit and uB as the higher-y unit. After this step
        % dy <= 0 for every pair, and pairs that share a physical
        % drift frame agree in dyEst sign as well as magnitude, so the
        % cohort consistency check becomes meaningful.
        if coms(uA0, 2) <= coms(uB0, 2)
            uA = uA0; uB = uB0;
            tA = tA0; gA = gA0;
        else
            uA = uB0; uB = uA0;
            tA = unitInfo(uA).tSet;
            gA = unitInfo(uA).localChGlobal;
            if isempty(tA) || all(isnan(gA)), continue; end
        end

        dx = coms(uA, 1) - coms(uB, 1);
        dy = coms(uA, 2) - coms(uB, 2);    % <= 0 after canonicalization
        if abs(dx) > opt.sameColumnXTol || abs(dy) > opt.maxDriftUm
            continue;
        end
        diag.nGeomPairs = diag.nGeomPairs + 1;
        diag.nEvaluated = diag.nEvaluated + 1;
        rejectReason = '';

        % Initialise fields for safe struct assignment
        corrRaw = NaN; bestR = NaN; dyEst = 0; nOverlap = 0;
        densCorr = NaN; ratio = NaN; isiViol = NaN;
        covA = NaN; covB = NaN; covM = NaN; covGain = NaN;
        acgCorr = NaN;

        % --- Dip gate: at least one unit must transiently fall silent -
        if opt.requireDip && ~(unitInfo(uA).hasDip && unitInfo(uB).hasDip)
            diag.rejDip = diag.rejDip + 1;
            rejectReason = 'dip';
        end

        % --- Amplitude ratio (hard) ---------------------------------
        if isempty(rejectReason)
            if peakAmps(uA) > 0 && peakAmps(uB) > 0
                ratio = max(peakAmps(uA), peakAmps(uB)) / min(peakAmps(uA), peakAmps(uB));
                if ratio > opt.ampRatioMax
                    diag.rejAmp = diag.rejAmp + 1;
                    rejectReason = 'amp';
                end
            else
                diag.rejAmp = diag.rejAmp + 1;
                rejectReason = 'amp';
            end
        end

        % --- Waveform drift estimation (COM-aligning shift) ----------
        % The ONLY physically valid drift shift for a pair is the one
        % that brings A's COM onto B's COM: dyShiftUm ~= dy = y_A - y_B.
        % If a merge would pass at dyShiftUm = 0 while |dy| is non-
        % trivial, the two units sit at different physical locations
        % and merely produce similar waveforms -- that is two separate
        % neurons, not a single drifting unit. Zero shift is therefore
        % acceptable ONLY when dy itself is already near zero.
        %
        % Concretely, we test the channel-aligned shifts that bracket
        % dy: s = floor(dy/pitchY)*pitchY and s = ceil(dy/pitchY)*pitchY.
        % When dy lands on a pitch multiple they coincide and a single
        % shift is tested; otherwise both bracket shifts are tested and
        % the larger correlation wins. No wider scan is performed --
        % any shift far from dy would violate the COM-alignment
        % constraint.
        %
        % For each candidate s, corrWithYShiftUm moves A's COM by s
        % toward B's COM and reads out A's waveform on B's observed
        % channels: for each B channel at (x, y) it finds A's channel
        % at (x, y + s) via same-x snap (tolY = 0.25 * pitchY) and
        % correlates the matched (A, B) waveform pairs.
        if isempty(rejectReason)
            gB = unitInfo(uB).localChGlobal; tB = unitInfo(uB).tSet;
            [corrRaw, ~] = corrOnGlobalIntersection(tA, gA, tB, gB); %#ok<ASGLU>  diagnostic only

            % Channel-aligned shift candidates bracketing dy.
            % The shift granularity is yStep = the smallest row-to-row
            % y-gap anywhere on the probe (NOT the within-column
            % pitch). On a staggered probe this means valid shifts
            % include odd multiples that cross between x-column pairs
            % (e.g. shift = yStep moves channels from one column-pair
            % to the other), which is a physically real drift regime
            % we would otherwise miss.
            %
            % s = 0 is only permitted when A and B cover the SAME set
            % of channels (identical footprints). If the footprints
            % differ, the two units sit at physically different
            % locations and zero shift cannot reconcile them -- a
            % valid drift shift must move A by at least one yStep.
            % For sub-yStep dy with non-identical footprints we pick
            % the nearest non-zero yStep multiple in the direction of
            % dy; for supra-yStep dy both bracket multiples are
            % tested and the higher correlation wins.
            gA_set = sort(gA(~isnan(gA)));
            gB_set = sort(gB(~isnan(gB)));
            identicalFootprint = isequal(gA_set, gB_set);

            if identicalFootprint
                kCand = 0;
            else
                kRatio = dy / yStep;
                kCand  = unique([floor(kRatio), ceil(kRatio)]);
                kCand  = kCand(kCand ~= 0);     % no zero shift when footprints differ
                if isempty(kCand)
                    if dy > 0
                        kCand = 1;
                    elseif dy < 0
                        kCand = -1;
                    else
                        kCand = [];             % dy == 0 with different footprints
                    end                         % -> no plausible drift shift
                end
            end
            shifts = kCand(:)' * yStep;

            bestR = -Inf; dyEst = NaN; nOverlap = 0;
            for s = shifts
                [rC, nC] = corrWithYShiftUm(tA, gA, tB, gB, ...
                            channel_locations, s, yStep);
                if isnan(rC) || nC < opt.minOverlapCh, continue; end
                if rC > bestR
                    bestR = rC; dyEst = s; nOverlap = nC;
                end
            end

            if ~isfinite(bestR) || nOverlap < opt.minOverlapCh
                diag.rejOverlap = diag.rejOverlap + 1;
                rejectReason = 'overlap';
                bestR = NaN;
            elseif bestR < opt.corrThreshold
                diag.rejCorr = diag.rejCorr + 1;
                rejectReason = 'corr';
            end
        end

        % --- Density anti-correlation (SOFT sanity check) ------------
        dA = unitInfo(uA).density;
        dB = unitInfo(uB).density;
        densCorr = safeCorr(dA(:)', dB(:)');
        if isempty(rejectReason)
            if ~isnan(densCorr) && densCorr > opt.densityAntiCorr
                diag.rejDens = diag.rejDens + 1;
                rejectReason = 'dens';
            end
        end

        % --- Coverage-gain (HARD) ------------------------------------
        % Merged activity mask union must cover substantially more
        % time than either unit individually. For true drift: A fires
        % in its epoch, B in its complementary epoch, union covers
        % the whole recording while each individual covers only half.
        mA = unitInfo(uA).activeMask;
        mB = unitInfo(uB).activeMask;
        Lm = min(numel(mA), numel(mB));
        if Lm > 0
            mA = mA(1:Lm); mB = mB(1:Lm);
            covA = mean(mA); covB = mean(mB);
            covM = mean(mA | mB);
            covGain = covM - max(covA, covB);
        end
        if isempty(rejectReason) && opt.requireCoverage
            if isnan(covGain) || covGain < opt.coverageGainMin
                diag.rejCoverage = diag.rejCoverage + 1;
                rejectReason = 'cov';
            end
        end

        % --- Refractory compatibility (hard) ------------------------
        if isempty(rejectReason)
            isiViol = mergedISIviolation(unitInfo(uA).spkIdxFull, ...
                unitInfo(uB).spkIdxFull, fs, opt.refractoryMs, jitter_gap);
            if isiViol > opt.isiThresholdPct
                diag.rejISI = diag.rejISI + 1;
                rejectReason = 'isi';
            end
        end

        % --- ACG similarity (hard, opt-in) --------------------------
        % The same neuron observed twice should have near-identical
        % autocorrelograms: same refractory hole, same burst/regular
        % structure. A pair with dissimilar ACGs therefore cannot be
        % a drift pair even if the waveforms look alike. Gate is
        % skipped when either unit fires fewer than acgMinSpikes.
        if isempty(rejectReason) && opt.acgCorrThreshold > -1
            acgA = unitInfo(uA).acg;
            acgB = unitInfo(uB).acg;
            if ~isempty(acgA) && ~isempty(acgB) && ...
               numel(acgA) == numel(acgB) && ...
               std(acgA) > 0 && std(acgB) > 0
                acgCorr = corr(acgA(:), acgB(:));
                if acgCorr < opt.acgCorrThreshold
                    diag.rejACG = diag.rejACG + 1;
                    rejectReason = 'acg';
                end
            end
        end

        clear ev; %#ok<CLEV>
        ev.uA = uA; ev.uB = uB;
        ev.corr = bestR; ev.corrRaw = corrRaw; ev.dyEst = dyEst;
        ev.dx = dx; ev.dy = dy;
        ev.ampRatio = ratio; ev.densCorr = densCorr;
        ev.isiViol = isiViol; ev.nOverlap = nOverlap;
        ev.acgCorr = acgCorr;
        ev.reject = rejectReason;
        ev.covA = covA; ev.covB = covB; ev.covM = covM; ev.covGain = covGain;
        diag.evidence(end+1) = ev; %#ok<AGROW>

        if ~isempty(rejectReason), continue; end
        pairs(end+1) = rmfield(ev, 'reject'); %#ok<AGROW>
    end
end

if opt.verbose, fprintf('Candidate pairs passing individual gates: %d\n', numel(pairs)); end

% ---- Cohort drift-field filter -----------------------------------------
% Fit a smooth dy_hat(y) field from the WAVEFORM-ESTIMATED drifts dyEst
% (not the noisy COM-based dy). Pairs whose dyEst disagrees with the
% local field (sign or magnitude beyond cohortDyTolUm) are rejected.
driftField = struct('pairY', [], 'pairDy', [], 'fieldY', [], 'fieldDy', []);
if opt.requireCohort && numel(pairs) >= 2
    nP = numel(pairs);
    mids = zeros(nP, 2);
    dys  = zeros(nP, 1);
    for p = 1:nP
        mids(p, :) = (coms(pairs(p).uA, :) + coms(pairs(p).uB, :)) / 2;
        dys(p)     = pairs(p).dyEst;       % waveform-estimated shift (um)
    end

    % Locally-weighted median field at each pair location
    dyHatPair = nan(nP, 1);
    for p = 1:nP
        d = sqrt(sum((mids - mids(p, :)).^2, 2));
        nb = (d <= opt.cohortRadiusUm) & ((1:nP)' ~= p);
        if nnz(nb) < opt.minCohortNeighbors, continue; end
        w = exp(-(d(nb).^2) / (2 * (opt.cohortRadiusUm/2)^2));
        dyHatPair(p) = weightedMedian(dys(nb), w);
    end

    keep = true(nP, 1);
    for p = 1:nP
        if isnan(dyHatPair(p)), continue; end       % isolated: keep
        % Sign disagreement (only if the field itself is non-zero)
        if sign(dyHatPair(p)) ~= 0 && sign(dys(p)) ~= 0 && ...
           sign(dys(p)) ~= sign(dyHatPair(p))
            keep(p) = false;
            diag.rejCohort = diag.rejCohort + 1;
            continue;
        end
        if abs(dys(p) - dyHatPair(p)) > opt.cohortDyTolUm
            keep(p) = false;
            diag.rejCohort = diag.rejCohort + 1;
        end
    end
    for p = 1:nP
        if ~keep(p)
            for e = 1:numel(diag.evidence)
                if diag.evidence(e).uA == pairs(p).uA && ...
                   diag.evidence(e).uB == pairs(p).uB && ...
                   isempty(diag.evidence(e).reject)
                    diag.evidence(e).reject = 'cohort';
                    break;
                end
            end
        end
    end

    % Save scatter + smoothed field for diagnostics
    driftField.pairY  = mids(:, 2);
    driftField.pairDy = dys;
    survMids = mids(keep, :);
    survDys  = dys(keep);
    if ~isempty(survDys)
        yGrid = linspace(min(survMids(:,2)), max(survMids(:,2)), 64);
        dyGrid = nan(size(yGrid));
        for gi = 1:numel(yGrid)
            d = abs(survMids(:,2) - yGrid(gi));
            mask = d <= opt.cohortRadiusUm;
            if nnz(mask) >= opt.minCohortNeighbors
                w = exp(-(d(mask).^2) / (2 * (opt.cohortRadiusUm/2)^2));
                dyGrid(gi) = weightedMedian(survDys(mask), w);
            end
        end
        driftField.fieldY  = yGrid;
        driftField.fieldDy = dyGrid;
    end

    pairs = pairs(keep);
end
diag.driftField = driftField;

if opt.verbose, fprintf('Pairs surviving cohort filter: %d\n', numel(pairs)); end

if opt.verbose
    fprintf(['  Geometry-close pairs     : %d\n' ...
             '    rejected by dip        : %d\n' ...
             '    rejected by amplitude  : %d\n' ...
             '    rejected by overlap    : %d\n' ...
             '    rejected by waveform corr (thr %.2f): %d\n' ...
             '    rejected by density (soft) : %d\n' ...
             '    rejected by coverage gain (>=%.2f) : %d\n' ...
             '    rejected by ISI        : %d\n' ...
             '    rejected by ACG corr (thr %.2f): %d\n' ...
             '    rejected by cohort     : %d\n'], ...
        diag.nGeomPairs, diag.rejDip, diag.rejAmp, diag.rejOverlap, ...
        opt.corrThreshold, diag.rejCorr, diag.rejDens, ...
        opt.coverageGainMin, diag.rejCoverage, diag.rejISI, ...
        opt.acgCorrThreshold, diag.rejACG, ...
        diag.rejCohort);

    if ~isempty(diag.evidence)
        wc = [diag.evidence.corr];  wc = wc(~isnan(wc));
        dc = [diag.evidence.densCorr]; dc = dc(~isnan(dc));
        fprintf('  Evidence summary over %d evaluated pairs:\n', numel(diag.evidence));
        if ~isempty(wc)
            fprintf('    waveform corr : min %.3f, median %.3f, max %.3f\n', ...
                min(wc), median(wc), max(wc));
        end
        if ~isempty(dc)
            fprintf('    density corr  : min %.3f, median %.3f, max %.3f\n', ...
                min(dc), median(dc), max(dc));
        end
    end
end

% ---- Diagnostic figures (TEMP) -----------------------------------------
% Plot top-N candidate pairs ranked by waveform correlation (whether they
% passed or not). Each figure shows density traces + multichannel waveforms
% + overlaid overlap channels, with the rejection reason in the title.
if opt.debugFigs
    if isempty(opt.figFolder)
        figFolder = fullfile(resSortedFolder, 'drift_debug_figs');
    else
        figFolder = char(opt.figFolder);
    end
    if ~exist(figFolder, 'dir')
        [okMk, mkMsg] = mkdir(figFolder);
        if ~okMk, error('mkdir failed for %s: %s', figFolder, mkMsg); end
    end
    if opt.verbose, fprintf('Debug-figure folder: %s\n', figFolder); end

    if isempty(diag.evidence)
        if opt.verbose
            fprintf(['No pairs were evaluated at all. Check:\n' ...
                     '  - channel_locations has same-x column members\n' ...
                     '  - any two units in the same x-column fall within maxDriftUm (%g um)\n' ...
                     '  - unit template sets (tSet) are non-empty\n'], opt.maxDriftUm);
        end
    else
        evAll = diag.evidence;
        rejects = arrayfun(@(e) ~isempty(e.reject), evAll);
        passIdx = find(~rejects);
        rejIdx  = find(rejects);

        % Rank each pool by corr descending so we see the strongest-looking
        % examples on each side.
        if ~isempty(passIdx)
            [~, ord] = sort([evAll(passIdx).corr], 'descend');
            passIdx = passIdx(ord);
        end
        if ~isempty(rejIdx)
            [~, ord] = sort([evAll(rejIdx).corr], 'descend');
            rejIdx = rejIdx(ord);
        end

        nPass = min(opt.nFigsPass,   numel(passIdx));
        nRej  = min(opt.nFigsReject, numel(rejIdx));
        if opt.verbose
            fprintf('Writing %d pass + %d reject diagnostic figures to %s\n', ...
                    nPass, nRej, figFolder);
        end

        tAxis = ((1:numel(unitInfo(1).density)) - 0.5) * opt.densityBinSec;
        dField = [];
        if isfield(diag, 'driftField'), dField = diag.driftField; end
        for ii = 1:nPass
            ev = evAll(passIdx(ii));
            plotDriftPair(figFolder, ii, 'pass', ev, unitInfo, ...
                          channel_locations, tAxis, opt, dField, fs);
        end
        for ii = 1:nRej
            ev = evAll(rejIdx(ii));
            plotDriftPair(figFolder, ii, 'rej', ev, unitInfo, ...
                          channel_locations, tAxis, opt, dField, fs);
        end
        if exist(figFolder, 'dir')
            dlist = dir(fullfile(figFolder, '*.png'));
            if opt.verbose, fprintf('  Wrote %d PNG(s) to disk.\n', numel(dlist)); end
        end
    end
end

% ---- Group construction (connected components, N>=2 supported) --------
% A drifting neuron may split across 2, 3, 4, or more channels, so we
% build merge groups as connected components of the pass-pair graph and
% then validate each component jointly. Groups that fail the joint
% validation are shrunk by iteratively dropping the weakest member.

groupRej = struct('accepted', 0, 'validate', 0, 'sizeCap', 0);
validGroups = {};

if ~isempty(pairs)
    % Union-find over all units present in pass pairs
    parent = 1:nUnits;
    for pi = 1:numel(pairs)
        parent = ufUnion(parent, pairs(pi).uA, pairs(pi).uB);
    end
    roots = nan(1, nUnits);
    for u = 1:nUnits
        [parent, roots(u)] = ufFind(parent, u);
    end
    % Collect components that contain at least one pair
    involved = false(1, nUnits);
    for pi = 1:numel(pairs)
        involved(pairs(pi).uA) = true;
        involved(pairs(pi).uB) = true;
    end
    uniqRoots = unique(roots(involved));

    for r = uniqRoots
        members = find(roots == r & involved);
        if numel(members) < 2, continue; end

        if numel(members) > opt.maxGroupSize
            % Keep the maxGroupSize highest-amplitude members
            [~, ord] = sort(arrayfun(@(u) unitInfo(u).peakAmp, members), 'descend');
            members = sort(members(ord(1:opt.maxGroupSize)));
            groupRej.sizeCap = groupRej.sizeCap + 1;
        end

        % Joint validation; iteratively peel weakest member if it fails
        m = members;
        while numel(m) >= 2
            if validateGroup(m, unitInfo, pairs, opt, fs, jitter_gap)
                validGroups{end+1} = sort(m); %#ok<AGROW>
                groupRej.accepted = groupRej.accepted + 1;
                break;
            else
                % Drop the weakest: unit with lowest mean pair-corr among group
                meanR = zeros(1, numel(m));
                cnt   = zeros(1, numel(m));
                for i = 1:numel(m)
                    for p = 1:numel(pairs)
                        if (pairs(p).uA == m(i)) || (pairs(p).uB == m(i))
                            if ismember(pairs(p).uA, m) && ismember(pairs(p).uB, m)
                                meanR(i) = meanR(i) + pairs(p).corr;
                                cnt(i) = cnt(i) + 1;
                            end
                        end
                    end
                end
                meanR(cnt > 0) = meanR(cnt > 0) ./ cnt(cnt > 0);
                meanR(cnt == 0) = -Inf;
                [~, wk] = min(meanR);
                m(wk) = [];
            end
        end
        if numel(m) < 2
            groupRej.validate = groupRej.validate + 1;
        end
    end
end

if opt.verbose
    fprintf(['Group construction: %d groups accepted, %d components rejected, ' ...
             '%d trimmed by size cap\n'], ...
        groupRej.accepted, groupRej.validate, groupRej.sizeCap);
end
diag.groupRej = groupRej;
if opt.verbose, fprintf('Validated drift groups: %d\n', numel(validGroups)); end

% ---- Choose primary per group and build label remap --------------------
labelRemap = containers.Map('KeyType', 'double', 'ValueType', 'double');
primarySummary = zeros(0, 4); % [groupIdx, primaryGlobal, primaryChannel, nMembers]

for gi = 1:numel(validGroups)
    members = validGroups{gi};
    counts  = arrayfun(@(u) double(unitInfo(u).sampleCount), members);
    detect  = arrayfun(@(u) double(unitInfo(u).detectability), members);
    [~, ord] = sortrows([-counts(:), -detect(:)]);
    primary = members(ord(1));
    primaryLabel = unitInfo(primary).globalLabel;
    primaryChan  = unitInfo(primary).channel;

    for m = members
        if m == primary, continue; end
        labelRemap(unitInfo(m).globalLabel) = primaryLabel;
    end
    primarySummary(end+1, :) = [gi, primaryLabel, primaryChan, numel(members)]; %#ok<AGROW>
end

% Apply remap to per-spike unifiedLabels
unifiedLab_merged = unifiedLab_all;
keysToRemap = cell2mat(labelRemap.keys);
for k = 1:numel(keysToRemap)
    src = keysToRemap(k);
    dst = labelRemap(src);
    unifiedLab_merged(unifiedLab_all == src) = dst;
end

% ---- Write outputs ------------------------------------------------------
mergedH5 = fullfile(resSortedFolder, 'unifiedLabels_merged.h5');
if ~opt.dryRun
    if exist(mergedH5, 'file'), delete(mergedH5); end
    h5create(mergedH5, '/unifiedLabels', size(unifiedLab_merged), ...
             'Datatype', class(unifiedLab_merged));
    h5write(mergedH5,  '/unifiedLabels', unifiedLab_merged);
    if opt.verbose, fprintf('Wrote %s\n', mergedH5); end

    if opt.overwriteOriginal
        backup = fullfile(resSortedFolder, 'unifiedLabels_predrift.h5');
        if ~exist(backup, 'file')
            copyfile(unifiedLabelsH5, backup);
            if opt.verbose, fprintf('Backed up original to %s\n', backup); end
        end
        delete(unifiedLabelsH5);
        h5create(unifiedLabelsH5, '/unifiedLabels', size(unifiedLab_merged), ...
                 'Datatype', class(unifiedLab_merged));
        h5write(unifiedLabelsH5, '/unifiedLabels', unifiedLab_merged);
        if opt.verbose, fprintf('Overwrote %s\n', unifiedLabelsH5); end
    end
end

driftReport.opt              = opt;
driftReport.pairs            = pairs;
driftReport.groups           = validGroups;
driftReport.labelRemap       = labelRemap;
driftReport.primarySummary   = primarySummary;
driftReport.nUnitsBefore     = nUnits;
driftReport.nUnitsAfter      = nUnits - sum(cellfun(@(g) numel(g)-1, validGroups));
driftReport.diagnostics      = diag;

if ~opt.dryRun
    reportPath = fullfile(resSortedFolder, 'drift_merge_report.mat');
    save(reportPath, 'driftReport', '-v7.3');
    if opt.verbose, fprintf('Wrote %s\n', reportPath); end
end

if opt.verbose
    fprintf('Drift merging done: %d -> %d units (%d groups merged).\n', ...
            driftReport.nUnitsBefore, driftReport.nUnitsAfter, numel(validGroups));
end
end


% =========================================================================
% Local helpers
% =========================================================================
function [r, nOverlap] = corrOnGlobalIntersection(tA, gA, tB, gB)
% Max correlation of A and B template variants over the physical channel
% overlap (same global channel index).
[commonCh, iA, iB] = intersect(gA, gB);
keep = ~isnan(commonCh);
commonCh = commonCh(keep); iA = iA(keep); iB = iB(keep);
nOverlap = numel(commonCh);
if nOverlap < 1, r = NaN; return; end
r = maxTemplateCorr(tA, tB, iA, iB);
end


function [r, nOverlap] = corrWithYShiftUm(tA, gA, tB, gB, chLocs, dyShiftUm, yStep)
% Test the drift hypothesis "A = B drifted by dyShiftUm" by shifting A's
% COM toward B and reading out A's waveform on B's observed channels.
%
% Convention: dyShiftUm = y_A - y_B aligns A with B (the COM-derived
% dy). For each B channel at (x_B, y_B), the corresponding A channel
% under this drift is whichever of A's OBSERVED channels (gA) lies
% closest in 2D to the drift-corrected target (x_B, y_B + dyShiftUm).
% The snap is NOT restricted to the same x-column: on a staggered
% probe, a drift of one yStep moves a neuron from one column-pair to
% the other, so the nearest A channel at the target position is
% legitimately in a different x-column than gB(k). A channel pair is
% accepted when the target lies within snapTol = 1.5 * yStep of an
% A observed channel; otherwise that B row is dropped from the
% correlation.
%
% Example (linear-column probe): A on channels 10:20, B on 15:25,
% dyShiftUm = dy ~ -100 um. Target for B ch 15 is at (x_{15}, y_{15}-
% 100), the closest A channel is ch 10, and so on -- correlation is
% the 11 (A, B) pairs.
%
% Example (staggered probe, yStep = 20 um): A at y around 2400 with
% x in {11, 43}, B at y around 2420 with x in {27, 59}, dyShiftUm =
% -20 um. For each B channel at (27, y_B) the target (27, y_B-20)
% has its nearest A channel at (11, y_B-20) or (43, y_B-20) at
% horizontal distance 16 um -- which snapTol = 30 comfortably admits.
r = NaN; nOverlap = 0;
if ~(yStep > 0), return; end

validA = ~isnan(gA);
if ~any(validA), return; end

posA_all = chLocs(gA(validA), 1:2);    % (nValidA x 2)
rowA_idx = find(validA);               % indices into tA rows

nLocB   = numel(gB);
snapTol = 1.5 * yStep;                 % um, covers min-diagonal snaps

rowA_for_B = nan(nLocB, 1);
for k = 1:nLocB
    if isnan(gB(k)), continue; end
    xyB    = chLocs(gB(k), 1:2);
    target = [xyB(1), xyB(2) + dyShiftUm];

    d         = sqrt(sum((posA_all - target).^2, 2));
    [dMin, i] = min(d);
    if dMin > snapTol, continue; end
    rowA_for_B(k) = rowA_idx(i);
end

validK   = find(~isnan(rowA_for_B));
nOverlap = numel(validK);
if nOverlap < 1, r = NaN; return; end

r = maxTemplateCorr(tA, tB, rowA_for_B(validK), validK);
end


function df = fitDriftField(seeds, coms, opt)
% Fit a smooth y-drift field from a set of seed pairs.
% seeds is a struct array with fields uA, uB, dyEst.
% The field is sampled on a regular y-grid via locally-weighted median.
df = struct('pairY', [], 'pairDy', [], 'fieldY', [], 'fieldDy', []);
if isempty(seeds)
    df.fieldY  = [];
    df.fieldDy = [];
    return;
end
n = numel(seeds);
ys = zeros(n, 1); dys = zeros(n, 1);
for p = 1:n
    ys(p)  = (coms(seeds(p).uA, 2) + coms(seeds(p).uB, 2)) / 2;
    dys(p) = seeds(p).dyEst;
end
df.pairY  = ys;
df.pairDy = dys;

if n < opt.minCohortNeighbors
    df.fieldY  = [];
    df.fieldDy = [];
    return;
end

yMin = min(ys); yMax = max(ys);
if yMax <= yMin
    df.fieldY  = yMin;
    df.fieldDy = weightedMedian(dys, ones(size(dys)));
    return;
end
yGrid = linspace(yMin, yMax, 128);
dyGrid = nan(size(yGrid));
sigma = max(opt.cohortRadiusUm / 2, 10);
for gi = 1:numel(yGrid)
    d = abs(ys - yGrid(gi));
    mask = d <= opt.cohortRadiusUm;
    if nnz(mask) >= opt.minCohortNeighbors
        w = exp(-(d(mask).^2) / (2 * sigma^2));
        dyGrid(gi) = weightedMedian(dys(mask), w);
    end
end
df.fieldY  = yGrid;
df.fieldDy = dyGrid;
end


function dyHat = interpField(df, y)
% Interpolate the fitted field at y. Returns NaN if y is outside the
% sampled range OR no field sample exists nearby.
dyHat = NaN;
if isempty(df.fieldY) || isempty(df.fieldDy), return; end
if y < df.fieldY(1) || y > df.fieldY(end), return; end
vm = ~isnan(df.fieldDy);
if nnz(vm) < 2, return; end
dyHat = interp1(df.fieldY(vm), df.fieldDy(vm), y, 'linear', NaN);
end


function [bestR, bestDy, bestN] = bestCorrInWindow(tA, gA, tB, gB, chLocs, dyCenter, slackUm, stepUm, minOverlapCh)
% Correlation in a small window [dyCenter-slack, dyCenter+slack] around
% the cohort-predicted shift. Returns best r and the corresponding dy.
bestR = NaN; bestDy = dyCenter; bestN = 0;
pitchY = inferYPitch(chLocs);
if pitchY <= 0, pitchY = 20; end
shifts = (dyCenter - slackUm):stepUm:(dyCenter + slackUm);
if isempty(shifts), shifts = dyCenter; end
for s = shifts
    if abs(s) < 1e-6
        [r, n] = corrOnGlobalIntersection(tA, gA, tB, gB);
    else
        [r, n] = corrWithYShiftUm(tA, gA, tB, gB, chLocs, s, pitchY);
    end
    if isnan(r) || n < minOverlapCh, continue; end
    if isnan(bestR) || r > bestR
        bestR = r; bestDy = s; bestN = n;
    end
end
end


function [bestR, bestDy, bestN] = estimateDriftShift(tA, gA, tB, gB, chLocs, maxDriftUm, stepUm, minOverlapCh)
% Search y-shifts of B from -maxDriftUm to +maxDriftUm in stepUm increments,
% return the shift that yields the highest template correlation subject
% to a minimum overlap constraint. This is the canonical per-pair drift
% estimate: the shift that makes the two templates most similar.
bestR = NaN; bestDy = 0; bestN = 0;
pitchY = inferYPitch(chLocs);
if pitchY <= 0, pitchY = 20; end
shifts = -maxDriftUm : stepUm : maxDriftUm;
for s = shifts
    if s == 0
        [r, n] = corrOnGlobalIntersection(tA, gA, tB, gB);
    else
        [r, n] = corrWithYShiftUm(tA, gA, tB, gB, chLocs, s, pitchY);
    end
    if isnan(r) || n < minOverlapCh, continue; end
    if isnan(bestR) || r > bestR
        bestR  = r;
        bestDy = s;
        bestN  = n;
    end
end
end


function m = weightedMedian(x, w)
% Weighted median of vector x with non-negative weights w.
x = x(:); w = w(:);
valid = ~isnan(x) & ~isnan(w) & w > 0;
x = x(valid); w = w(valid);
if isempty(x), m = NaN; return; end
[xs, ord] = sort(x);
ws = w(ord);
cw = cumsum(ws);
half = cw(end) / 2;
idx = find(cw >= half, 1, 'first');
m = xs(idx);
end


function [parent, r] = ufFind(parent, x)
while parent(x) ~= x
    parent(x) = parent(parent(x));
    x = parent(x);
end
r = x;
end


function parent = ufUnion(parent, a, b)
[parent, ra] = ufFind(parent, a);
[parent, rb] = ufFind(parent, b);
if ra ~= rb, parent(rb) = ra; end
end


function p = inferYPitch(chLocs)
% Minimum positive same-x y-difference across columns.
p = 0;
xs = chLocs(:,1); ys = chLocs(:,2);
[uniqX, ~, xg] = unique(xs);
diffsAll = [];
for ix = 1:numel(uniqX)
    y = sort(ys(xg == ix));
    d = diff(y); d = d(d > 0);
    diffsAll = [diffsAll; d(:)]; %#ok<AGROW>
end
if ~isempty(diffsAll), p = median(diffsAll); end
if p <= 0, p = 1; end
end


function acg = computeACG(spkSamples, fs, windowMs, binMs)
% Autocorrelogram of a spike train, binned and zero-lag masked.
%
% Bins the spike train at binMs resolution, then takes the zero-centred
% cross-correlation of the binned counts with itself out to windowMs on
% each side via xcorr. The zero-lag bin is then zeroed so the ACG
% reflects only refractory/bursting structure, not the trivial spike-
% self coincidence. Returned as a length-(2*halfBins + 1) row vector.
halfBins = floor(windowMs / binMs);
nBins    = 2 * halfBins + 1;
spk      = double(spkSamples(:));
if numel(spk) < 2
    acg = zeros(1, nBins);
    return;
end

% Bin spike times at binMs.
binSec  = binMs / 1000;
binIdx  = max(1, floor(spk / fs / binSec) + 1);
nT      = max(binIdx);
counts  = accumarray(binIdx, 1, [nT, 1]);

% ACG = xcorr of binned counts with itself at lags [-halfBins, halfBins].
acg = xcorr(counts, halfBins);   % column vector, length nBins
acg = acg(:)';
acg(halfBins + 1) = 0;           % mask zero lag
end


function s = inferYStep(chLocs)
% Smallest positive y-difference between ANY two probe channels.
%
% On a staggered / zig-zag probe (adjacent rows live in different
% x-columns), this row-to-row step is SMALLER than the per-column
% y-pitch and defines the finest drift shift the probe can resolve.
% A drift of one yStep typically moves the neuron from one column-
% pair to the other; a drift of 2*yStep stays within the original
% column-pair; and so on.
%
% For a non-staggered probe (all rows in every x-column), this
% reduces to the within-column y-pitch returned by inferYPitch.
ys = sort(unique(chLocs(:, 2)));
if numel(ys) < 2
    s = inferYPitch(chLocs);
    return;
end
d = diff(ys);
d = d(d > 0);
if isempty(d)
    s = inferYPitch(chLocs);
else
    s = min(d);
end
if s <= 0, s = 1; end
end


function [len, startIdx] = longestRun(b)
% Longest run of true values in logical vector b.
b = logical(b(:))';
len = 0; startIdx = 0;
if isempty(b), return; end
edges = diff([false, b, false]);
starts = find(edges == 1);
stops  = find(edges == -1) - 1;
if isempty(starts), return; end
runLens = stops - starts + 1;
[len, iMax] = max(runLens);
startIdx = starts(iMax);
end


function r = maxTemplateCorr(tA, tB, idxA, idxB)
% Max Pearson r across all (variantA, variantB) template pairs over the
% given matched channel indices. Templates are (nVariants, nChannels, T).
% idxA and idxB must have the same length and refer to the rows of tA/tB
% that correspond to the SAME physical (global) channel.
nA = size(tA, 1);
nB = size(tB, 1);
nT = size(tA, 3);
nO = numel(idxA);
if nO < 1, r = NaN; return; end
FA = reshape(tA(:, idxA, :), nA, nO * nT);
FB = reshape(tB(:, idxB, :), nB, nO * nT);
FA = FA - mean(FA, 2);
FB = FB - mean(FB, 2);
nrA = sqrt(sum(FA.^2, 2)); nrA(nrA == 0) = inf;
nrB = sqrt(sum(FB.^2, 2)); nrB(nrB == 0) = inf;
C = (FA ./ nrA) * (FB ./ nrB)';
r = max(C(:));
if isempty(r), r = NaN; end
end


function [r, cvA, cvB, cvM, cvRatio] = densityMetrics(dA, dB, activeOnly, activeFrac)
% Pair metrics on two smoothed firing-rate traces.
%   r        Pearson correlation (restricted to active support if activeOnly)
%   cvA,B,M  CV = std/mean of each rate trace and of the merged (A+B) trace,
%            evaluated on the same active support.
%   cvRatio  cvM / min(cvA, cvB); values <1 mean the merge makes the rate
%            more stationary (consistent with a single drifting unit).
dA = dA(:)'; dB = dB(:)';
n = min(numel(dA), numel(dB));
dA = dA(1:n); dB = dB(1:n);
r = NaN; cvA = NaN; cvB = NaN; cvM = NaN; cvRatio = NaN;
if activeOnly
    mA = max(dA); mB = max(dB);
    if mA <= 0 || mB <= 0, return; end
    mask = (dA > activeFrac * mA) | (dB > activeFrac * mB);
    if nnz(mask) < 5, return; end
    dA = dA(mask); dB = dB(mask);
end
r = safeCorr(dA, dB);
if mean(dA) > 0, cvA = std(dA) / mean(dA); end
if mean(dB) > 0, cvB = std(dB) / mean(dB); end
dM = dA + dB;
if mean(dM) > 0, cvM = std(dM) / mean(dM); end
if ~isnan(cvA) && ~isnan(cvB) && ~isnan(cvM)
    denom = min(cvA, cvB);
    if denom > 0, cvRatio = cvM / denom; end
end
end


function r = safeCorr(x, y)
x = x(:); y = y(:);
if numel(x) ~= numel(y) || numel(x) < 3 || std(x) == 0 || std(y) == 0
    r = NaN; return;
end
r = corr(x, y);
end




function viol = mergedISIviolation(spkA, spkB, fs, refractoryMs, jitter_gap)
% Merge two spike trains and drop near-simultaneous co-detections (the same
% physical spike picked up on both channels during drift) before evaluating
% the refractory violation percentage.
spkA = spkA(:); spkB = spkB(:);
if ~isempty(spkA) && ~isempty(spkB)
    [d1, d2] = nearest_distances(spkA, spkB);
    coMask = (d1 < jitter_gap) | (d2 < jitter_gap);
    spkA = spkA(~coMask);
end
allSpk = sort([spkA; spkB]);
if numel(allSpk) < 2, viol = 0; return; end
[~, ~, viol] = getISIViolations(allSpk, fs, refractoryMs);
end


function ok = validateGroup(members, unitInfo, pairs, opt, fs, jitter_gap)
% Joint validation for a merge group of any size >= 2.
%
% Criteria (all must pass):
%   (1) Chain connectivity: after sorting members by y-COM, every consecutive
%       pair must appear in the pass-pair list. This enforces a valid drift
%       trajectory without requiring all O(n^2) pairs to be directly tested.
%   (2) Consecutive y-gaps all <= maxDriftUm.
%   (3) Joint refractory: merged train (with iterative co-detection removal)
%       has ISI violation < isiThresholdPct.
%   (4) Coverage-gain: the fraction of time bins where the merged unit
%       is active must exceed the individual max by covGainMin. This
%       is the "merge fills in silent epochs" signature of drift.
ok = false;
n = numel(members);
if n < 2, return; end

% Sort by y-position to get the drift trajectory order
ys = arrayfun(@(u) unitInfo(u).com(2), members);
[ys_sorted, ord] = sort(ys(:));
members = members(ord);

% (1) Chain connectivity: consecutive members must be pass-pairs
for i = 1:(n-1)
    ui = members(i); uj = members(i+1);
    found = false;
    for p = 1:numel(pairs)
        if (pairs(p).uA == ui && pairs(p).uB == uj) || ...
           (pairs(p).uA == uj && pairs(p).uB == ui)
            found = true; break;
        end
    end
    if ~found, return; end
end

% (2) y-trajectory gaps
gaps = diff(ys_sorted);
if any(gaps > opt.maxDriftUm), return; end

% (3) Joint refractory check with co-detection removal
memberTrains = cell(n, 1);
for k = 1:n
    memberTrains{k} = unitInfo(members(k)).spkIdxFull(:);
end
for k = 2:n
    if isempty(memberTrains{k}), continue; end
    for kk = 1:(k-1)
        if isempty(memberTrains{kk}), continue; end
        [d1, d2] = nearest_distances(memberTrains{k}, memberTrains{kk});
        keepMask = ~((d1 < jitter_gap) | (d2 < jitter_gap));
        memberTrains{k} = memberTrains{k}(keepMask);
    end
end
allSpikes = sort(vertcat(memberTrains{:}));
if numel(allSpikes) > 1
    [~, ~, isiViol] = getISIViolations(allSpikes, fs, opt.refractoryMs);
else
    isiViol = 0;
end
if isiViol > opt.isiThresholdPct, return; end

% (4) Coverage-gain on the member union: the merged unit must light up
% more of the recording than the best individual member. This enforces
% the drift signature "complementary active epochs".
if opt.requireCoverage
    unionMask = unitInfo(members(1)).activeMask(:)';
    maxIndivCov = mean(unionMask);
    for k = 2:n
        mk = unitInfo(members(k)).activeMask(:)';
        L = min(numel(unionMask), numel(mk));
        unionMask = unionMask(1:L) | mk(1:L);
        maxIndivCov = max(maxIndivCov, mean(mk(1:L)));
    end
    covM = mean(unionMask);
    if (covM - maxIndivCov) < opt.coverageGainMin, return; end
end

ok = true;
end


function plotDriftPair(figFolder, rank, tag, ev, unitInfo, channel_locations, tAxis, opt, driftField, fs)
if nargin < 9, driftField = []; end
if nargin < 10, fs = []; end
% Diagnostic figure for a candidate drift pair.
%
% Layout (3x4):
%   (1,1) Firing-rate density of A and B over the whole recording
%   (1,2) Active-mask coverage; title shows covA/covB/covM/gain
%   (1,3) Text panel: PASS/REJECT, COMs, all criteria
%   (1,4) Stacked overlay of A and B on their overlapping global channels
%   (2,1) A waveforms placed at their physical (x,y) positions
%   (2,2) B waveforms placed at their physical (x,y) positions
%   (2,3) Cohort drift-field scatter (red = this pair dy_est)
%   (2,4) Drift-aligned overlay (B translated by dy_est)
%   (3,1) ACG of A and B overlaid (peak-normalised), title shows corr
try
    uA = ev.uA; uB = ev.uB;
    infoA = unitInfo(uA);
    infoB = unitInfo(uB);
    wfA = infoA.meanWF;   % [nLocalCh x T]
    wfB = infoB.meanWF;
    if isempty(wfA) || isempty(wfB), return; end

    gA = infoA.localChGlobal;
    gB = infoB.localChGlobal;
    [commonCh, iA, iB] = intersect(gA, gB);
    validMask = ~isnan(commonCh);
    commonCh = commonCh(validMask);
    iA = iA(validMask);
    iB = iB(validMask);

    binSec = opt.densityBinSec;

    fig = figure('Visible', 'off', 'Position', [50 50 1800 1500]);

    % --- Density panels --------------------------------------------------
    subplot(3,4,1); hold on;
    tMin = tAxis / 60;
    plot(tMin, infoA.density, 'b', 'LineWidth', 1.2);
    plot(tMin, infoB.density, 'r', 'LineWidth', 1.2);
    xlabel('Time (min)'); ylabel('Rate (Hz)');
    title(sprintf('Density (bin %gs, \\sigma %gs)  r=%.3f', ...
          binSec, opt.densitySmoothSec, ev.densCorr));
    legend(sprintf('A: ch %d', infoA.channel), ...
           sprintf('B: ch %d', infoB.channel), 'Location', 'best');
    grid on;

    % --- Active-mask coverage panel -----------------------------------
    % Shows where A is active (blue), where B is active (red), and the
    % union (black). covGain = mean(union) - max(mean(A), mean(B)); a
    % true drift pair has complementary epochs so covGain >> 0.
    subplot(3,4,2); hold on;
    mA = infoA.activeMask(:)';
    mB = infoB.activeMask(:)';
    L = min(numel(mA), numel(mB));
    mA = mA(1:L); mB = mB(1:L);
    mU = mA | mB;
    tM = tMin(1:L);
    stairs(tM, double(mA) * 0.8 + 2.1, 'Color', [0.3 0.5 0.9], 'LineWidth', 1.2);
    stairs(tM, double(mB) * 0.8 + 1.05, 'Color', [0.9 0.4 0.3], 'LineWidth', 1.2);
    stairs(tM, double(mU) * 0.8 + 0.0, 'k', 'LineWidth', 1.4);
    ylim([-0.1 3.1]); set(gca, 'YTick', [0.4 1.45 2.5], 'YTickLabel', {'A|B', 'B', 'A'});
    xlabel('Time (min)');
    title(sprintf('Active mask  covA=%.2f  covB=%.2f  covM=%.2f  gain=%.2f', ...
                  fallback(ev, 'covA', mean(mA)), fallback(ev, 'covB', mean(mB)), ...
                  fallback(ev, 'covM', mean(mU)), fallback(ev, 'covGain', NaN)));
    grid on;

    % --- Text panel ----------------------------------------------------
    subplot(3,4,3);
    comA = infoA.com; comB = infoB.com;
    corrRawStr = 'n/a';
    if isfield(ev, 'corrRaw') && ~isnan(ev.corrRaw)
        corrRawStr = sprintf('%.3f', ev.corrRaw);
    end
    dyEstVal = fallback(ev, 'dyEst', 0);
    if dyEstVal ~= 0
        shiftStr = sprintf('%+.1f um', dyEstVal);
    else
        shiftStr = 'none (raw)';
    end
    dipAstr = sprintf('%s (%.0fs)', yesno(infoA.hasDip), infoA.dipDurSec);
    dipBstr = sprintf('%s (%.0fs)', yesno(infoB.hasDip), infoB.dipDurSec);
    covGainVal = fallback(ev, 'covGain', NaN);
    txt = sprintf(['%s  RANK %d\n' ...
                   'uA = %d  (ch %d, gl %d, n=%d)  dip: %s\n' ...
                   '  COM = (%.1f, %.1f) um\n' ...
                   'uB = %d  (ch %d, gl %d, n=%d)  dip: %s\n' ...
                   '  COM = (%.1f, %.1f) um\n\n' ...
                   'dx = %.2f, dy = %.2f um (COM)\n' ...
                   'dy_est (waveform): %+.1f um\n' ...
                   'overlap ch      : %d\n' ...
                   'waveform corr   : %.3f   (raw %s)\n' ...
                   'drift shift used: %s\n' ...
                   'amp ratio       : %.3f\n' ...
                   'density corr    : %.3f\n' ...
                   'coverage gain   : %.3f  (>= %.2f needed)\n' ...
                   'isi violation %% : %.3f\n\n' ...
                   'status          : %s'], ...
                   upper(tag), rank, ...
                   uA, infoA.channel, infoA.globalLabel, numel(infoA.spkIdxFull), dipAstr, ...
                   comA(1), comA(2), ...
                   uB, infoB.channel, infoB.globalLabel, numel(infoB.spkIdxFull), dipBstr, ...
                   comB(1), comB(2), ...
                   ev.dx, ev.dy, dyEstVal, ev.nOverlap, ev.corr, corrRawStr, shiftStr, ...
                   ev.ampRatio, ev.densCorr, covGainVal, opt.coverageGainMin, ev.isiViol, ...
                   ternary(isempty(ev.reject), 'PASSED', ['REJECT: ' ev.reject]));
    axis off;
    text(0.02, 0.5, txt, 'FontName', 'Courier', 'FontSize', 10, ...
        'VerticalAlignment', 'middle');

    % --- Stacked overlay on overlapping channels -------------------------
    subplot(3,4,4); hold on;
    spacing = max([max(abs(wfA(:))), max(abs(wfB(:))), eps]) * 1.6;
    nOvlp = numel(commonCh);
    for oi = 1:nOvlp
        off = (nOvlp - oi) * spacing;
        plot(wfA(iA(oi), :) + off, 'b', 'LineWidth', 1.1);
        plot(wfB(iB(oi), :) + off, 'r', 'LineWidth', 1.1);
        text(-5, off, sprintf('g%d', commonCh(oi)), ...
             'FontSize', 7, 'HorizontalAlignment', 'right');
    end
    hold off; set(gca, 'YTick', []); xlabel('Time samples');
    title(sprintf('Overlap: %d ch, r=%.3f', nOvlp, ev.corr));
    legend('A', 'B', 'Location', 'best');

    % --- Geometric layout: A at physical (x,y) --------------------------
    subplot(3,4,5);
    plotGeometricWF(wfA, gA, channel_locations, spacing, 'b');
    title(sprintf('A at probe positions (home %d, gl %d)', ...
          infoA.channel, infoA.globalLabel));

    % --- Geometric layout: B at physical (x,y) --------------------------
    subplot(3,4,6);
    plotGeometricWF(wfB, gB, channel_locations, spacing, 'r');
    title(sprintf('B at probe positions (home %d, gl %d)', ...
          infoB.channel, infoB.globalLabel));

    % --- Cohort drift-field panel --------------------------------------
    subplot(3,4,7); hold on;
    if ~isempty(driftField) && ~isempty(driftField.pairY)
        scatter(driftField.pairY, driftField.pairDy, 18, [0.6 0.6 0.6], 'filled');
        if ~isempty(driftField.fieldY)
            plot(driftField.fieldY, driftField.fieldDy, 'k-', 'LineWidth', 1.8);
        end
    end
    yMid = (infoA.com(2) + infoB.com(2)) / 2;
    plot(yMid, fallback(ev, 'dyEst', 0), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
    hold off; grid on; box on;
    xlabel('y (um)'); ylabel('dy_{est} (um)');
    title('Cohort drift field (red = this pair dy_{est})');

    % --- Geometric overlay on REAL positions ----------------------------
    %   Leave the original (2,4,7) slot -- already re-used above.
    %   The drift-aligned overlay is still at (2,4,8).
    %   Skip the redundant "A+B on real positions" panel.

    % --- Drift-aligned overlay: B shifted by the waveform-estimated dy_est
    % Canonical convention: dy_est is the y-shift applied to B's channels
    % (in corrWithYShiftUm) that maximised the template correlation.
    % Positive dy_est => B shifts up; negative => B shifts down.
    subplot(3,4,8); hold on;
    plotGeometricWF(wfA, gA, channel_locations, spacing, 'b');
    locShifted = channel_locations;
    dyUsed = fallback(ev, 'dyEst', 0);
    locShifted(:,2) = locShifted(:,2) + dyUsed;
    plotGeometricWF(wfB, gB, locShifted, spacing, 'r');
    hold off;
    if dyUsed ~= 0
        title(sprintf('Drift-aligned (dy_{est} = %+.1f um, COM dy = %+.1f)', ...
                      dyUsed, ev.dy));
    else
        title(sprintf('No shift applied (COM dy = %+.1f um)', ev.dy));
    end

    % --- ACG overlay row -------------------------------------------------
    % Use the precomputed per-unit ACG when available; otherwise compute
    % on-the-fly so the overlay is populated even when the ACG gate is
    % disabled (acgCorrThreshold <= -1).
    subplot(3,4,9); hold on;
    acgA = [];
    acgB = [];
    if isfield(infoA, 'acg') && ~isempty(infoA.acg), acgA = infoA.acg; end
    if isfield(infoB, 'acg') && ~isempty(infoB.acg), acgB = infoB.acg; end
    if (isempty(acgA) || isempty(acgB)) && ~isempty(fs)
        winMs = opt.acgWindowMs; binMs = opt.acgBinMs;
        if isempty(acgA) && numel(infoA.spkIdxFull) >= 2
            acgA = computeACG(infoA.spkIdxFull, fs, winMs, binMs);
        end
        if isempty(acgB) && numel(infoB.spkIdxFull) >= 2
            acgB = computeACG(infoB.spkIdxFull, fs, winMs, binMs);
        end
    end
    if ~isempty(acgA) && ~isempty(acgB) && numel(acgA) == numel(acgB)
        halfBins = floor(opt.acgWindowMs / opt.acgBinMs);
        lags_ms  = (-halfBins:halfBins) * opt.acgBinMs;
        nrmA = max(acgA); if nrmA <= 0, nrmA = 1; end
        nrmB = max(acgB); if nrmB <= 0, nrmB = 1; end
        plot(lags_ms, acgA ./ nrmA, 'b', 'LineWidth', 1.2);
        plot(lags_ms, acgB ./ nrmB, 'r', 'LineWidth', 1.2);
        acgR = fallback(ev, 'acgCorr', NaN);
        if ~isnan(acgR)
            title(sprintf('ACG overlay (peak-normalised)  r=%.3f', acgR));
        else
            title('ACG overlay (peak-normalised)');
        end
        xlabel('Lag (ms)'); ylabel('Norm. count');
        xlim([lags_ms(1) lags_ms(end)]);
        legend('A', 'B', 'Location', 'best');
        grid on;
    else
        text(0.5, 0.5, 'ACG unavailable (too few spikes)', ...
             'HorizontalAlignment', 'center');
        axis off;
        title('ACG overlay');
    end

    status = ternary(isempty(ev.reject), 'pass', ev.reject);
    fname = sprintf('%s_rank%03d_%s_u%d_u%d_ch%d_ch%d.png', ...
        tag, rank, status, uA, uB, infoA.channel, infoB.channel);
    saveas(fig, fullfile(figFolder, fname));
    close(fig);
catch ME
    try close(fig); catch, end
    if opt.verbose
        fprintf(2, 'plotDriftPair FAILED for %s rank %d: %s\n', tag, rank, ME.message);
        for s = 1:min(3, numel(ME.stack))
            fprintf(2, '  at %s:%d\n', ME.stack(s).name, ME.stack(s).line);
        end
    end
    rethrow(ME);
end
end


function plotGeometricWF(wf, gChan, chLocs, ampScale, clr)
% Plot multi-channel waveform at each channel's physical (x, y) position.
%   wf     [nLocalCh x T]
%   gChan  [nLocalCh x 1] global channel index (NaN if row invalid)
%   chLocs (numGlobalCh x 3) channel_locations
%   ampScale   peak-to-peak amplitude scale used for per-trace vertical extent
%   clr    colour for the traces
hold on;
[nCh, T] = size(wf);
tt = linspace(0, 1, T);

% Estimate typical inter-channel spacing to choose trace width in um
validMask = ~isnan(gChan);
if nnz(validMask) < 2
    xs = 0; ys = 0;
else
    xy = chLocs(gChan(validMask), 1:2);
    xs = xy(:,1); ys = xy(:,2);
end
if numel(xs) >= 2
    dx = max(max(xs) - min(xs), 1);
    dy = max(max(ys) - min(ys), 1);
else
    dx = 1; dy = 1;
end

% Each trace occupies 80% of horizontal channel pitch, amplitude scaled
% to ~80% of vertical pitch
chPitchX = estimatePitch(xs, 'min_positive_diff');
chPitchY = estimatePitch(ys, 'min_positive_diff');
if chPitchX == 0, chPitchX = dx / max(numel(unique(xs)), 1); end
if chPitchY == 0, chPitchY = dy / max(numel(unique(ys)), 1); end
traceXW = 0.8 * chPitchX;
traceYH = 0.8 * chPitchY;
if ampScale <= 0, ampScale = 1; end

for c = 1:nCh
    gi = gChan(c);
    if isnan(gi), continue; end
    xy = chLocs(gi, 1:2);
    x = xy(1) + (tt - 0.5) * traceXW;
    y = xy(2) + wf(c, :) / ampScale * traceYH;
    plot(x, y, 'Color', clr, 'LineWidth', 1.0);
end
xlabel('x (um)'); ylabel('y (um)');
axis equal; box on; grid on;
end


function p = estimatePitch(v, ~)
v = unique(v(:));
if numel(v) < 2, p = 0; return; end
d = diff(v);
d = d(d > 0);
if isempty(d), p = 0; else, p = min(d); end
end


function v = ternary(cond, a, b)
if cond, v = a; else, v = b; end
end


function s = yesno(b)
if b, s = 'yes'; else, s = 'no'; end
end


function v = fallback(s, fld, default)
if isfield(s, fld) && ~isempty(s.(fld)) && ~any(isnan(s.(fld)(:)))
    v = s.(fld);
else
    v = default;
end
end


