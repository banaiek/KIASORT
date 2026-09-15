function report = kiaSort_auto_curate_nogui(outputPath, cfg, varargin)
%KIASORT_AUTO_CURATE_NOGUI  Headless equivalent of the GUI Auto-curation.
%
%   report = kiaSort_auto_curate_nogui(outputPath, cfg, ...)
%
%   Runs the same pipeline as the "Auto" button in kiaSort_curate_results
%   and writes the same artefacts into <outputPath>/RES_Sorted:
%       *_curated.h5, curated_sample.mat, curated_metrics.csv
%
%   Phases:
%     1) drop zero-rate units
%     2) overlap handling: chance-corrected duplicate detection, then
%        strip / merge (see local_pairEvidence)
%     3) gated merges (XCorr + amplitude + footprint shape + merged ISI)
%     4) isolation classification (SUA+/SUA/MUA+/MUA)
%     5) stable-interval trimming
%
%   Note this runs AFTER kiaSort_post_sort_curate, which kiaSort_main_sortData
%   already applies in place when cfg.postHocProcessing is true. The raw
%   unifiedLabels.h5 / spike_idx.h5 read here are therefore already
%   de-duplicated once; this pass is what produces the curated outputs.
%
%   Name/Value:
%       'overlapFrac'      (0.10)  aggregate chance-corrected overlap to act on
%       'distUm'           (200)   neighbourhood radius in y (um)
%       'ccgRatio'         (2)     CCG lag-0 peak / baseline ratio
%       'waveSim'          (0.80)  shared-channel waveform similarity gate
%       'mergeSim'         (0.90)  footprint-wide similarity gate for merges
%       'lagTightSamples'  (2)     +-N samples around the waveform offset
%       'lagConsistency'   (0.75)  min fraction of coincident lags in that window
%       'maxStripFrac'     (0.35)  per-unit cumulative strip budget
%       'isiAbs'           (1.0)   merged-ISI hard cap (%)
%       'isiBudget'        (1.5)   merged-ISI weighted budget
%       'ampDiff'          (0.15)  amplitude-similarity gate for merges
%       'simThr'           (0.85)  XCorr similarity gate for merges
%       'dropRate'         (2)     rate-drop factor for stable-interval trim
%       'thresholdISIms'   (1)     refractory window (ms)
%       'extractWaveforms' (true)  write waveforms_curated.h5 (large)
%       'verbose'          (true)

p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addRequired('cfg',        @isstruct);
p.addParameter('overlapFrac',      0.10, @(x) isscalar(x) && isnumeric(x));
p.addParameter('distUm',           200,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('ccgRatio',         2,    @(x) isscalar(x) && isnumeric(x));
p.addParameter('waveSim',          0.80, @(x) isscalar(x) && isnumeric(x));
p.addParameter('mergeSim',         0.90, @(x) isscalar(x) && isnumeric(x));
p.addParameter('lagTightSamples',  2,    @(x) isscalar(x) && isnumeric(x));
p.addParameter('lagConsistency',   0.75, @(x) isscalar(x) && isnumeric(x));
p.addParameter('maxStripFrac',     0.35, @(x) isscalar(x) && isnumeric(x));
p.addParameter('isiAbs',           1.0,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('isiBudget',        1.5,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('ampDiff',          0.15, @(x) isscalar(x) && isnumeric(x));
p.addParameter('simThr',           0.85, @(x) isscalar(x) && isnumeric(x));
p.addParameter('dropRate',         2,    @(x) isscalar(x) && isnumeric(x));
p.addParameter('thresholdISIms',   1,    @(x) isscalar(x) && isnumeric(x));
p.addParameter('minConfirmSpikes', 50,   @(x) isscalar(x) && isnumeric(x));
p.addParameter('extractWaveforms', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',          true, @(x) islogical(x) || isnumeric(x));
p.parse(outputPath, cfg, varargin{:});
opt = p.Results;
opt.verbose = logical(opt.verbose);

report = struct('ok', false, 'nZeroDrops', 0, 'nOverlapStrips', 0, ...
                'nOverlapMerges', 0, 'nMerges', 0, 'nUnitsKept', 0);

outputPath  = char(outputPath);
resFolder   = fullfile(outputPath, 'RES_Sorted');
sampFolder  = fullfile(outputPath, 'RES_Samples');
sortSampDir = fullfile(outputPath, 'Sorted_Samples');

chInfoPath  = fullfile(sampFolder,  'channel_info.mat');
sortSampPath= fullfile(sortSampDir, 'sorted_samples.mat');

if ~exist(chInfoPath,'file') || ~exist(sortSampPath,'file')
    error('kiaSort_auto_curate_nogui:missingInputs', ...
        'channel_info.mat or sorted_samples.mat not found under %s', outputPath);
end

fs = cfg.samplingFrequency;

%% Load
chInfo = load(chInfoPath);
ss     = load(sortSampPath, 'crossChannelStats');
unif          = ss.crossChannelStats.unified_labels;

% A probe on which no spike was detected writes no spike h5 files at all, and
% reading them unconditionally turned "this probe is silent" into a hard error.
% That is a real outcome rather than a fault: Wotan 112 D yielded a single unit
% at 3 channels and none at 5, because 51 of its 128 channels carry stimulation
% artefact and what remains is nearly silent. Treat the absent files as zero
% spikes and let curation return an empty result.
spkFile = fullfile(resFolder,'spike_idx.h5');
lblFile = fullfile(resFolder,'unifiedLabels.h5');
chnFile = fullfile(resFolder,'channelNum.h5');
if ~exist(spkFile,'file') || ~exist(lblFile,'file') || ~exist(chnFile,'file')
    warning('kiaSort_auto_curate_nogui:noSpikes', ...
        'No spike files under %s - treating as zero detected spikes.', resFolder);
    spk_idx = zeros(0,1); lbl = zeros(0,1); chn = zeros(0,1);
else
    spk_idx = double(h5read(spkFile, '/spike_idx'));
    lbl     = double(h5read(lblFile, '/unifiedLabels'));
    chn     = double(h5read(chnFile, '/channelNum'));
end
spk_idx = spk_idx(:); lbl = lbl(:); chn = chn(:);
if numel(spk_idx) ~= numel(lbl) || numel(spk_idx) ~= numel(chn)
    error('kiaSort_auto_curate_nogui:h5Mismatch', 'spike h5 lengths disagree.');
end

groupList    = double(unif.label(:));
channelList  = double(unif.channelID(:));
sampleWave   = unif.meanWaveforms;
detectblity  = double(unif.detectblity(:));
mainPolarity = unif.mainNegativePolarity;
sidePolarity = unif.sideNegativePolarity;
numGroups    = numel(groupList);

originalLabels     = lbl;
originalChannelNum = chn;
originalWaveform   = sampleWave;

% sampleWave rows stay centred on the unit's original channel even after a
% merge rewrites channelList, so template lookups must use this copy.
homeChannel = channelList;

ylocs        = chInfo.channel_locations(:,2);
channel_map  = chInfo.channel_mapping(:);
num_Samples  = double(chInfo.num_samples);
trialLength  = num_Samples / fs;
numChannelPlot = cfg.num_channel_extract;

unitNotes     = repmat({''},   numGroups, 1);
stable_length = [zeros(numGroups,1), trialLength * ones(numGroups,1)];
mergedFlag    = false(numGroups,1);
snrAll        = 1 + detectblity;

coincSamples = round(0.5e-3 * fs);
tightSamples = max(1, round(opt.lagTightSamples));

firingRate = local_metrics(groupList, lbl, spk_idx, fs, trialLength, opt.thresholdISIms);

%% Phase 1: zero-rate drops
for k = 1:numGroups
    if isnan(groupList(k)), continue; end
    if firingRate(k) == 0
        groupList(k)   = NaN;
        channelList(k) = NaN;
        mergedFlag(k)  = true;
        report.nZeroDrops = report.nZeroDrops + 1;
    end
end

%% Phase 2: overlap handling
% Evidence is gathered for every neighbour pair first, then applied per
% unit. Gathering before acting keeps the outcome independent of pair
% order, which the GUI's single-pass sweep is not.
liveCount0 = zeros(numGroups,1);
for k = 1:numGroups
    if isnan(groupList(k)), continue; end
    liveCount0(k) = sum(lbl == groupList(k));
end

edges = struct('cont', {}, 'own', {}, 'rows', {}, 'lag', {}, 'excess', {});
for k = 1:numGroups
    if isnan(groupList(k)) || isnan(channelList(k)), continue; end
    for j = (k+1):numGroups
        if isnan(groupList(j)) || isnan(channelList(j)), continue; end
        if groupList(k) == groupList(j), continue; end
        if abs(ylocs(channelList(j)) - ylocs(channelList(k))) > opt.distUm, continue; end

        ev = local_pairEvidence(k, j, groupList, homeChannel, snrAll, ...
            lbl, spk_idx, sampleWave, numChannelPlot, ...
            fs, coincSamples, tightSamples, num_Samples, opt);
        if ~isempty(ev)
            edges(end+1) = ev; %#ok<AGROW>
        end
    end
end

% Aggregate the confirmed duplicate rows per contaminated unit. Union of
% rows (not a sum of fractions) so a spike coincident with two neighbours
% is not counted twice.
stripBudget = zeros(numGroups,1);
for k = 1:numGroups
    stripBudget(k) = floor(opt.maxStripFrac * max(liveCount0(k),0));
end

consumed = false(numGroups,1);
for k = 1:numGroups
    if isnan(groupList(k)) || consumed(k), continue; end
    mine = find([edges.cont] == k);
    if isempty(mine), continue; end

    unionRows = [];
    excessSum = 0;
    for e = mine
        unionRows = union(unionRows, edges(e).rows);
        excessSum = excessSum + edges(e).excess;
    end
    unionRows = unionRows(lbl(unionRows) == groupList(k));
    if isempty(unionRows), continue; end

    aggFrac = numel(unionRows) / max(liveCount0(k), 1);
    if aggFrac < opt.overlapFrac || excessSum <= 0, continue; end

    if numel(unionRows) <= stripBudget(k)
        lbl(unionRows) = -1;
        report.nOverlapStrips = report.nOverlapStrips + numel(unionRows);
    else
        % Mostly duplicate: merging is recall-neutral where stripping past
        % the budget is not. Requires footprint-wide similarity and a
        % merged-train refractory veto against the dominant owner.
        [~, best] = max([edges(mine).excess]);
        own  = edges(mine(best)).own;
        eLag = edges(mine(best)).lag;
        if isnan(groupList(own)) || groupList(own) == groupList(k), continue; end

        keepRows = setdiff(find(lbl == groupList(k)), unionRows);
        spkKeep  = spk_idx(keepRows);
        if isfinite(eLag) && eLag ~= 0 && ~isempty(spkKeep)
            spkKeep = min(max(spkKeep + round(eLag), 1), num_Samples);
        end
        spkOwn = spk_idx(lbl == groupList(own));

        fpSim = local_footprintSim(k, own, homeChannel(k), sampleWave, ...
            homeChannel, numChannelPlot, 2, coincSamples);
        okISI = local_mergedISIok(spkOwn, spkKeep, fs, opt.thresholdISIms, opt.isiAbs);

        if isfinite(fpSim) && fpSim >= opt.mergeSim && okISI
            oldLab  = groupList(k);
            members = find(groupList == oldLab);
            lbl(unionRows) = -1;
            report.nOverlapStrips = report.nOverlapStrips + numel(unionRows);
            if ~isempty(keepRows)
                spk_idx(keepRows) = spkKeep;
                lbl(keepRows)     = groupList(own);
                chn(keepRows)     = channelList(own);
            end
            groupList(members)   = groupList(own);
            channelList(members) = channelList(own);
            mergedFlag(members)  = true;
            mergedFlag(own)      = true;
            consumed(members)    = true;
            report.nOverlapMerges = report.nOverlapMerges + 1;
        end
    end
end

%% Phase 3: gated merges
for a = 1:numGroups
    if isnan(groupList(a)) || isnan(channelList(a)), continue; end
    for b = (a+1):numGroups
        if isnan(groupList(b)) || isnan(channelList(b)), continue; end
        if groupList(a) == groupList(b), continue; end
        if abs(ylocs(channelList(b)) - ylocs(channelList(a))) > opt.distUm, continue; end

        ch  = homeChannel(a);
        mw1 = local_meanWaveform(a, ch, sampleWave, homeChannel, numChannelPlot);
        mw2 = local_meanWaveform(b, ch, sampleWave, homeChannel, numChannelPlot);
        if isempty(mw1) || isempty(mw2) || numel(mw1) ~= numel(mw2) || numel(mw1) < 5
            continue;
        end
        M2 = numel(mw1);
        [simScore, bestLag] = max_half_corr(mw1(:)', mw2(:)', 1, M2, max(1,round(M2/4)), 0);
        if ~isfinite(simScore) || simScore < opt.simThr, continue; end
        if local_ampDiff(mw1, mw2) > opt.ampDiff, continue; end
        % Footprint-wide shape agreement stands in for the GUI's PC-distance
        % gate: the stored PCA basis is fit on flattened multi-channel
        % waveforms and cannot be applied to a single-channel template.
        % NaN means too few usable channels to judge, which does not block.
        fpMerge = local_footprintSim(a, b, ch, sampleWave, homeChannel, ...
            numChannelPlot, 2, coincSamples);
        if isfinite(fpMerge) && fpMerge < opt.mergeSim, continue; end

        spkA = spk_idx(lbl == groupList(a));
        spkB = spk_idx(lbl == groupList(b));
        if isempty(spkA) || isempty(spkB), continue; end
        if ~local_mergedISIfull(spkA, spkB, fs, opt.thresholdISIms, opt.isiAbs, opt.isiBudget)
            continue;
        end

        if numel(spkA) >= numel(spkB)
            primary = a; absorbed = b; shiftSign = -1;
        else
            primary = b; absorbed = a; shiftSign = +1;
        end
        absLab = groupList(absorbed);
        if isfinite(bestLag) && bestLag ~= 0
            rows = (lbl == absLab);
            if any(rows)
                spk_idx(rows) = min(max(spk_idx(rows) + shiftSign*bestLag, 1), num_Samples);
            end
        end
        members = find(groupList == absLab);
        chn(lbl == absLab) = channelList(primary);
        lbl(lbl == absLab) = groupList(primary);
        groupList(members)   = groupList(primary);
        channelList(members) = channelList(primary);
        mergedFlag(members)  = true;
        report.nMerges = report.nMerges + 1;
    end
end

%% Phase 4: isolation classification on the curated trains
[firingRate, isiViol] = local_metrics(groupList, lbl, spk_idx, fs, trialLength, opt.thresholdISIms);
unitIsolation = local_classify(groupList, lbl, spk_idx, fs, isiViol, snrAll, numGroups);

%% Phase 5: stable-interval trimming
smoothN  = 0;
ccgLag   = 100;
binSize  = max(1, round(trialLength/(2*ccgLag)));
for k = 1:numGroups
    if isnan(groupList(k)), continue; end
    stable_length(k,:) = [0, trialLength];
    spk_k = unique(spk_idx(lbl == groupList(k)));
    if numel(spk_k) < 10, continue; end
    try
        [ctr, cnt] = local_presenceRatio(spk_k, fs, trialLength, binSize, smoothN);
        if numel(cnt) >= 3
            [sI, eI] = local_stableInterval(cnt, opt.dropRate);
            span = max(ctr(end), eps);
            stable_length(k,1) = max(0, min(trialLength, trialLength*(ctr(sI)-ctr(1))/span));
            stable_length(k,2) = max(0, min(trialLength, trialLength*ctr(eI)/span));
        end
    catch
    end
end

%% Outputs
if ~exist(resFolder,'dir'), mkdir(resFolder); end

% stable_length is in seconds; spike indices are samples.
validSpikes = zeros(size(spk_idx));
for k = 1:numGroups
    if isnan(groupList(k)), continue; end
    idx = find(lbl == groupList(k));
    if isempty(idx), continue; end
    lo = stable_length(k,1) * fs;
    hi = stable_length(k,2) * fs;
    validSpikes(idx(spk_idx(idx) >= lo & spk_idx(idx) <= hi)) = 1;
end

sorted_out.unifiedLabels_curated{1,1} = lbl;
sorted_out.channelNum_curated{1,1}    = chn;
sorted_out.spike_idx_curated{1,1}     = spk_idx;
sorted_out.inclusion_curated{1,1}     = validSpikes;
if logical(opt.extractWaveforms)
    sorted_out.waveforms_curated{1,1} = local_extractWaveforms(spk_idx, chn, ...
        channel_map, cfg, fs, num_Samples);
end
local_saveCuratedH5(resFolder, sorted_out);

[uniqLabels, repIdx] = unique(groupList);
keep       = ~isnan(uniqLabels);
uniqLabels = uniqLabels(keep);
repIdx     = repIdx(keep);

curatedSamples.unifiedLabels = uniqLabels;
curatedSamples.channelNum    = channelList(repIdx);
curatedSamples.waveform      = sampleWave(repIdx,:,:);
curatedSamples.unitIsolation = unitIsolation(repIdx);
curatedSamples.validInterval = stable_length(repIdx,:);
curatedSamples.spikePolarity = mainPolarity(repIdx,:);
curatedSamples.notes         = unitNotes(repIdx);

curatedSamples.original.unifiedLabels = originalLabels;
curatedSamples.original.channelNum    = originalChannelNum;
curatedSamples.original.waveform      = originalWaveform;
curatedSamples.original.unitIsolation = unitIsolation;

curatedSamples.session.groupList     = groupList;
curatedSamples.session.channelList   = channelList;
curatedSamples.session.unitIsolation = unitIsolation;
curatedSamples.session.stable_length = stable_length;
curatedSamples.session.mainPolarity  = mainPolarity;
curatedSamples.session.sidePolarity  = sidePolarity;
curatedSamples.session.mergedFlag    = mergedFlag;
curatedSamples.session.unitNotes     = unitNotes;
curatedSamples.session.detectblity   = detectblity;
curatedSamples.session.numGroups     = numGroups;
curatedSamples.session.spikeLabels   = lbl;
curatedSamples.session.spikeIdx      = spk_idx;
curatedSamples.session.spikeChannels = chn;

save(fullfile(resFolder,'curated_sample.mat'), 'curatedSamples', '-v7.3');

try
    T = table(uniqLabels(:), channelList(repIdx), firingRate(repIdx), ...
        isiViol(repIdx), snrAll(repIdx), unitIsolation(repIdx), ...
        mainPolarity(repIdx), stable_length(repIdx,1), stable_length(repIdx,2), ...
        unitNotes(repIdx), ...
        'VariableNames', {'Label','Channel','FiringRate_Hz','ISI_violation_pct', ...
                          'SNR','Isolation','Polarity','StableStart','StableEnd','Notes'});
    writetable(T, fullfile(resFolder,'curated_metrics.csv'));
catch ME
    warning('CSV export failed: %s', ME.message);
end

report.nUnitsKept = numel(uniqLabels);
report.ok = true;

if opt.verbose
    fprintf(['Auto curation (no-GUI): %d zero-rate drops, %d overlap strips, ' ...
        '%d overlap merges, %d merges, %d units kept.\n'], ...
        report.nZeroDrops, report.nOverlapStrips, report.nOverlapMerges, ...
        report.nMerges, report.nUnitsKept);
end

end


% =========================================================================

function [rate, isiv] = local_metrics(groupList, lbl, spk, fs, trialLength, thrMs)
n = numel(groupList);
rate = zeros(n,1);
isiv = zeros(n,1);
for i = 1:n
    if isnan(groupList(i)), continue; end
    s = spk(lbl == groupList(i));
    if isempty(s), continue; end
    rate(i) = numel(s) / trialLength;
    if numel(s) >= 2
        try
            [~,~,isiv(i)] = getISIViolations(s, fs, thrMs);
        catch
        end
    end
end
end


function iso = local_classify(groupList, lbl, spk, fs, isiViol, snrAll, n)
iso = repmat({'NA'}, n, 1);
for k = 1:n
    if isnan(groupList(k)), continue; end
    isi = isiViol(k);   if isnan(isi), isi = 2.0; end
    sr  = snrAll(k);    if isnan(sr),  sr  = 1.0; end
    acg = local_acgRatio(spk(lbl == groupList(k)), fs);
    if isnan(acg), acg = 1; end

    sIsi = max(0, 1 - isi/2.0);
    sSnr = max(0, min(1, (sr-1)/4));
    sAcg = max(0, 1 - acg/0.10);
    score = 0.4*sIsi + 0.4*sAcg + 0.2*sSnr;

    if score >= 0.85 && isi < 0.5 && acg < 0.05
        iso{k} = 'SUA+';
    elseif score >= 0.65 && isi < 1.0
        iso{k} = 'SUA';
    elseif score >= 0.45
        iso{k} = 'MUA+';
    else
        iso{k} = 'MUA';
    end
end
end


function r = local_acgRatio(spk, fs)
r = 1;
if numel(spk) < 50, return; end
d = diff(sort(double(spk))) / fs * 1000;
if isempty(d), return; end
baseline = sum(d >= 5 & d <= 25);
if baseline < 5, return; end
r = sum(d < 1.5) / baseline;
end


function ev = local_pairEvidence(k, j, groupList, channelList, snrAll, lbl, spk, ...
    sampleWave, nChanPlot, fs, coincSamples, tightSamples, nSamp, opt)
% Confirms a pair double-detects the same spikes, and returns the rows of
% the contaminated unit that are duplicate copies. Empty when unconfirmed.
%
% Confirmation requires all of: a lag-0 CCG peak above baseline, matching
% waveforms on the contaminated unit's channel, coincident lags clustered
% at the waveform offset, and a coincidence count above what independent
% firing would produce by chance.
ev = [];
chList = channelList;

rowsK = find(lbl == groupList(k));
rowsJ = find(lbl == groupList(j));
if numel(rowsK) < opt.minConfirmSpikes || numel(rowsJ) < opt.minConfirmSpikes
    return;
end
spkK = spk(rowsK);
spkJ = spk(rowsJ);

sk = snrAll(k); if isnan(sk), sk = 0; end
sj = snrAll(j); if isnan(sj), sj = 0; end
ownerIsK = sk > sj || (sk == sj && (numel(spkK) > numel(spkJ) || ...
    (numel(spkK) == numel(spkJ) && groupList(k) > groupList(j))));
if ownerIsK
    ci = j; oi = k; contRows = rowsJ; spkC = spkJ; spkO = spkK;
else
    ci = k; oi = j; contRows = rowsK; spkC = spkK; spkO = spkJ;
end

ccgBin = round(1e-3 * fs);
[ccg, zb] = local_ccg(spkC, spkO, 100*ccgBin, ccgBin);
if ~local_ccgPeak(ccg, zb, 1, opt.ccgRatio), return; end

chC = channelList(ci);
if isnan(chC), return; end
if nChanPlot == 0
    % No shared-channel view exists; compare each unit on its own channel.
    tC = local_meanWaveform(ci, chList(ci), sampleWave, chList, nChanPlot);
    tO = local_meanWaveform(oi, chList(oi), sampleWave, chList, nChanPlot);
else
    tC = local_meanWaveform(ci, chC, sampleWave, chList, nChanPlot);
    tO = local_meanWaveform(oi, chC, sampleWave, chList, nChanPlot);
end
if isempty(tC) || isempty(tO) || numel(tC) ~= numel(tO) || numel(tC) < 5
    return;
end
z  = zeros(coincSamples,1);
tC = [z; tC(:); z];
tO = [z; tO(:); z];
M2 = numel(tC);
[simWF, waveLag] = max_half_corr(tC(:)', tO(:)', 1, M2, max(1,round(M2/4)), 0);
if ~isfinite(simWF) || simWF < opt.waveSim, return; end
if ~isfinite(waveLag), waveLag = 0; end

[d1, d2] = nearest_distances(spkC, spkO);
dC = min(d1, d2);
coinc = dC <= coincSamples;
if ~any(coinc), return; end

lags = local_signedLag(spkC, spkO);
ci_idx = find(coinc);
cl = lags(ci_idx);
cl = cl(~isnan(cl));
if numel(cl) < 5, return; end
medLag = median(cl);
if abs(medLag - waveLag) > tightSamples, return; end
within = abs(cl - medLag) <= tightSamples;
if sum(within) / numel(cl) < opt.lagConsistency, return; end

tight = false(numel(lags),1);
tight(ci_idx) = within;
dupRows = contRows(coinc & tight);
if isempty(dupRows), return; end

% Chance coincidence for independent trains: the owner's rate times the
% tight window, applied to the contaminated unit's spike count. Without
% this a high-rate unit in a dense neighbourhood reads as contaminated.
rateO   = numel(spkO) / max(nSamp, 1);
chanceN = numel(spkC) * rateO * (2*tightSamples + 1);
excess  = numel(dupRows) - chanceN;
if excess <= 0 || excess < 4 * sqrt(max(chanceN, 1)), return; end

ev = struct('cont', ci, 'own', oi, 'rows', dupRows(:)', ...
            'lag', waveLag, 'excess', excess);
end


function [ccg, zeroBin] = local_ccg(A, B, maxLag, binSz)
A = sort(A(:)); B = sort(B(:));
half = floor(maxLag / binSz);
nB   = 2*half + 1;
zeroBin = half + 1;
ccg = zeros(nB,1);
if isempty(A) || isempty(B), return; end
lo = 1; hi = 0; nb = numel(B);
d = cell(numel(A),1);
for i = 1:numel(A)
    t = A(i);
    while lo <= nb && B(lo) < t - maxLag, lo = lo + 1; end
    while hi < nb && B(hi+1) <= t + maxLag, hi = hi + 1; end
    if lo <= hi, d{i} = B(lo:hi) - t; end
end
all_d = vertcat(d{:});
if isempty(all_d), return; end
bi = round(double(all_d)/binSz) + zeroBin;
bi = bi(bi >= 1 & bi <= nB);
if isempty(bi), return; end
ccg = accumarray(bi(:), 1, [nB, 1]);
end


function tf = local_ccgPeak(ccg, zeroBin, tolBins, ratio)
tf = false;
if isempty(ccg) || all(ccg == 0), return; end
n = numel(ccg);
pk = max(1, zeroBin-tolBins) : min(n, zeroBin+tolBins);
peak = max(ccg(pk));
off = ccg; off(pk) = [];
if isempty(off) || any(off > peak) || ~any(off > 0), return; end
if sum(off > 0) < 20, return; end
base = prctile(off, 90);
if ~isfinite(base) || base <= 0, return; end
tf = peak > ratio * base;
end


function lag = local_signedLag(A, B)
A = A(:); B = B(:);
lag = nan(numel(A),1);
if isempty(A) || isempty(B), return; end
Bs = sort(B);
[As, ord] = sort(A);
j = 1; nb = numel(Bs);
for i = 1:numel(As)
    while j < nb && Bs(j+1) < As(i), j = j + 1; end
    best = Bs(j) - As(i);
    for t = [j-1, j+1]
        if t >= 1 && t <= nb
            c = Bs(t) - As(i);
            if abs(c) < abs(best), best = c; end
        end
    end
    lag(ord(i)) = best;
end
end


function mw = local_meanWaveform(g, ch, sampleWave, chList, nChanPlot)
mw = [];
if isnan(chList(g)), return; end
li = ch - chList(g) + nChanPlot + 1;
if li < 1 || li > 2*nChanPlot + 1, return; end
if li > size(sampleWave, 2), return; end
mw = squeeze(sampleWave(g, li, :));
end


function s = local_footprintSim(gC, gO, ch0, sampleWave, chList, nChanPlot, nSide, padN)
% Similarity across the footprint rather than one channel: two co-active
% neurons can match on a shared channel but not across the footprint.
s = NaN;
if isnan(ch0), return; end
sims = [];
for dc = -nSide:nSide
    rC = local_meanWaveform(gC, ch0+dc, sampleWave, chList, nChanPlot);
    rO = local_meanWaveform(gO, ch0+dc, sampleWave, chList, nChanPlot);
    if isempty(rC) || isempty(rO) || numel(rC) ~= numel(rO) || numel(rC) < 5
        continue;
    end
    if max(abs(rC)) <= eps || max(abs(rO)) <= eps, continue; end
    z  = zeros(padN,1);
    a  = [z; rC(:); z];
    b  = [z; rO(:); z];
    M2 = numel(a);
    c  = max_half_corr(a(:)', b(:)', 1, M2, max(1,round(M2/4)), 0);
    if isfinite(c), sims(end+1) = c; end %#ok<AGROW>
end
if numel(sims) >= 2 || (nChanPlot == 0 && numel(sims) == 1), s = mean(sims); end
end


function d = local_ampDiff(mw1, mw2)
both = [mw1(:)'; mw2(:)'];
mx = max(abs(both), [], 2);
den = max(mx);
if ~isfinite(den) || den <= 0, d = Inf; return; end
p = abs(max(both, [], 2));
n = abs(min(both, [], 2));
d = (abs(p(1)-p(2)) + abs(n(1)-n(2))) / (2*den);
end


function ok = local_mergedISIok(spkA, spkB, fs, thrMs, isiAbs)
ok = true;
m = sort([spkA(:); spkB(:)]);
if numel(m) < 2, return; end
try
    [~,~,v] = getISIViolations(m, fs, thrMs);
    ok = v <= isiAbs;
catch
    ok = false;
end
end


function ok = local_mergedISIfull(spkA, spkB, fs, thrMs, isiAbs, isiBudget)
ok = true;
NA = numel(spkA); NB = numel(spkB);
if NA < 1 || NB < 1, return; end
m = sort([spkA(:); spkB(:)]);
if numel(m) < 2, return; end
try
    [~,~,isiM] = getISIViolations(m, fs, thrMs);
catch
    return;
end
isiA = 0; isiB = 0;
if NA >= 2
    try
        [~,~,isiA] = getISIViolations(spkA, fs, thrMs);
    catch
        isiA = 0;
    end
end
if NB >= 2
    try
        [~,~,isiB] = getISIViolations(spkB, fs, thrMs);
    catch
        isiB = 0;
    end
end

parentCap = isiAbs * max(isiBudget, 1);
den = max(NA + NB, 1);
wIsi = (NA*isiA + NB*isiB) / den;
sizeFactor = 2 * min(NA,NB) / den;
budgetCap = max((1 + (isiBudget-1)*sizeFactor) * wIsi, 1e-6);

ok = (isiM <= isiAbs) && (isiA <= parentCap) && (isiB <= parentCap) && ...
     (abs(isiA - wIsi) <= 0.5) && (abs(isiB - wIsi) <= 0.5) && (isiM <= budgetCap);
end


function [centers, counts] = local_presenceRatio(spk, fs, trialLength, binSize, smoothF)
t = double(spk) / fs;
overlap = binSize * smoothF;
maxStart = trialLength - binSize;
if maxStart < 0
    centers = []; counts = [];
    return;
end
starts  = (0:binSize:maxStart)';
centers = starts + binSize/2;
counts  = zeros(size(centers));
for i = 1:numel(centers)
    counts(i) = sum(t > starts(i)-overlap & t < starts(i)+binSize+overlap);
end
end


function [startIdx, endIdx] = local_stableInterval(d, dropFactor)
n = numel(d);
startIdx = 1; endIdx = n;
if n < 3, return; end
cp = findchangepts(d, 'Statistic', 'rms', 'MinThreshold', 5);
if isscalar(cp)
    if mean(d(1:cp)) < mean(d(cp+1:end)) / dropFactor
        startIdx = cp + 1;
    elseif mean(d(1:cp)) / dropFactor > mean(d(cp+1:end))
        endIdx = cp;
    end
elseif numel(cp) >= 2
    for ci = 1:numel(cp)
        if mean(d(1:cp(ci))) < mean(d(cp(ci)+1:end))/dropFactor && ...
           mean(d(1:cp(ci))) < mean(d(cp(ci)+1:min(2*cp(ci),n)))/dropFactor
            startIdx = cp(ci)+1;
        elseif mean(d(startIdx:cp(ci)))/dropFactor > mean(d(cp(ci)+1:end)) && ...
               mean(d(startIdx:cp(ci)))/dropFactor > mean(d(cp(ci)+1:min(2*cp(ci),n)))
            endIdx = cp(ci);
            break
        end
    end
end
end


function wf = local_extractWaveforms(spk, chn, channel_map, cfg, fs, nSamp)
half = round((cfg.spikeDuration/2) * fs / 1000);
W = 2*half + 1;
wf = zeros(numel(spk), W);
m = map_input_file(cfg.fullFilePath, cfg);
for i = 1:numel(spk)
    r = spk(i) + (-half:half);
    if r(1) < 1 || r(end) > nSamp, continue; end
    c = chn(i);
    if isnan(c) || c < 1 || c > numel(channel_map), continue; end
    wf(i,:) = double(m.Data.data(channel_map(c), r));
end
end


function local_saveCuratedH5(outputFolder, sorted_out)
% Only *curated*.h5 are cleared; the raw sort h5 in the same folder must
% survive (the shared saveh5SpikeData deletes every *.h5).
old = dir(fullfile(outputFolder, '*curated*.h5'));
for k = 1:numel(old)
    delete(fullfile(outputFolder, old(k).name));
end
flds = fieldnames(sorted_out);
for i = 1:numel(flds)
    fld  = flds{i};
    data = sorted_out.(fld){1,1};
    if isempty(data), continue; end
    f = fullfile(outputFolder, [fld '.h5']);
    sz = size(data);
    chunk = sz;
    chunk(1) = min(sz(1), 100);
    h5create(f, ['/' fld], sz, 'ChunkSize', chunk);
    h5write(f, ['/' fld], data);
end
end
