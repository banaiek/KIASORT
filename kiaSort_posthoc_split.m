function splitReport = kiaSort_posthoc_split(outputPath, varargin)
%KIASORT_POSTHOC_SPLIT  Split units whose amplitude distribution is bimodal.
%
%   splitReport = kiaSort_posthoc_split(outputPath, ...)
%
%   Template matching assigns a spike to a class on waveform shape within an
%   amplitude band, and unification merges classes on normalised cross-
%   correlation. Both are amplitude-blind, so two neurons with the same shape
%   at different amplitudes can end up as one unit. This pass finds those and
%   separates them.
%
%   Per unit:
%     1) Otsu the per-spike amplitude. Accept only a clear separation
%        (separability) -- no valley depth is assumed.
%     2) Veto drift: the two sides must overlap in time. A drifting single
%        unit is also bimodal but its halves are consecutive, not interleaved.
%     3) Re-cluster the waveforms of the two sides and keep the split only if
%        the waveforms agree with the amplitude cut. Child ISI is recorded
%        but not gated: splitting a cell out of noise leaves a dirty child
%        by design.
%     4) Require both children to be stable across the recording. A child
%        confined to part of the parent's span is an epoch, not a cell, and
%        the split is reverted.
%   Waveforms come from RES_Sorted/waveforms*.h5 when the sort saved them
%   (cfg.extractWaveform), otherwise they are read from the raw file named by
%   the stored cfg.fullFilePath and bandpass filtered.
%
%   Only unifiedLabels.h5 and the per-unit arrays in
%   crossChannelStats.unified_labels change. Nothing is written unless at
%   least one split is accepted.
%
%   Name/Value:
%       'separability'    (scalar, cfg.bimodalSeparability)    Otsu eta gate
%       'minSepDistance'  (scalar, 1.75)  alternative detection: mode distance
%                                         in pooled sd, for unbalanced splits
%       'minChildFrac'    (scalar, 0)     optional extra bar as a fraction of the unit
%       'minChildSpikes'  (scalar, 200)   absolute min child size
%       'minTimeOverlap'  (scalar, cfg.bimodalMinTimeOverlap)  drift veto
%       'maxMedianShift'  (scalar, cfg.bimodalMaxMedianShift)  drift veto
%       'waveSim'         (scalar, 0.6)   min agreement between the waveform
%                                         clustering and the amplitude cut
%       'minSpikes'       (scalar, 200)   skip units smaller than this
%       'maxSplits'       (scalar, Inf)   cap on accepted splits per run
%       'fitCap'          (scalar, 3000)  waveforms read when extracting raw
%       'assignMaxSpikes' (scalar, 2e5)   above this, fall back to the amplitude
%                                         cut rather than read every waveform
%       'ccgIndepMin'     (scalar, 0.5)   children must fire independently:
%                                         coincidence over chance in ccgBandMs
%       'ccgMinExpected'  (scalar, 20)    below this the test abstains -- with a
%                                         small minority "observed 0" means nothing
%       'ccgBandMs'       (1x2, [1 2])    lag band, above the detector dead time
%       'assignMinAgree'  (scalar, 0.75)  waveform assignment must agree with the
%                                         amplitude cut at least this well, else
%                                         the cut is kept
%       'minPresence'     (scalar, 0.8)   each child must appear in this
%                                         fraction of the parent's active bins
%       'nBins'           (scalar, 20)    bins used for that test
%       'minBinSpikes'    (scalar, 5)     parent spikes for a bin to count
%       'verbose'         (logical, false)

p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addParameter('separability',   [], @(x) isempty(x) || isscalar(x));
p.addParameter('minSepDistance', [], @(x) isempty(x) || isscalar(x));
p.addParameter('minChildFrac',    0, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minChildSpikes', 200, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minTimeOverlap', [], @(x) isempty(x) || isscalar(x));
p.addParameter('maxMedianShift', [], @(x) isempty(x) || isscalar(x));
p.addParameter('waveSim',      0.6,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('minSpikes',    200,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('maxSplits',    Inf,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('fitCap',       3000, @(x) isscalar(x) && isnumeric(x));
p.addParameter('assignMaxSpikes', 200000, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minPresence',  0.8,  @(x) isscalar(x) && isnumeric(x));
p.addParameter('ccgIndepMin',    [], @(x) isempty(x) || isscalar(x));
p.addParameter('ccgMinExpected', [], @(x) isempty(x) || isscalar(x));
p.addParameter('ccgBandMs', [1 2],   @(x) isnumeric(x) && numel(x)==2);
p.addParameter('assignMinAgree', [], @(x) isempty(x) || isscalar(x));
p.addParameter('nBins',        20,   @(x) isscalar(x) && isnumeric(x));
p.addParameter('minBinSpikes', 5,    @(x) isscalar(x) && isnumeric(x));
p.addParameter('verbose',      false, @(x) islogical(x) || isnumeric(x));
p.parse(outputPath, varargin{:});
opt = p.Results;
opt.verbose = logical(opt.verbose);

splitReport = struct('nSplit', 0, 'nTested', 0, 'newLabels', [], ...
                     'changed', false, 'ok', false, 'log', []);

outputPath       = char(outputPath);
resSortedFolder  = fullfile(outputPath, 'RES_Sorted');
sortedSampFolder = fullfile(outputPath, 'Sorted_Samples');
samplesFolder    = fullfile(outputPath, 'RES_Samples');

unifiedLabelsH5  = fullfile(resSortedFolder, 'unifiedLabels.h5');
spikeIdxH5       = fullfile(resSortedFolder, 'spike_idx.h5');
channelNumH5     = fullfile(resSortedFolder, 'channelNum.h5');
amplitudeH5      = fullfile(resSortedFolder, 'amplitude.h5');
sortedSamplesPath = fullfile(sortedSampFolder, 'sorted_samples.mat');
channelInfoPath   = fullfile(samplesFolder,   'channel_info.mat');

required = {unifiedLabelsH5, spikeIdxH5, channelNumH5, amplitudeH5, sortedSamplesPath};
for i = 1:numel(required)
    if ~exist(required{i}, 'file')
        if opt.verbose
            fprintf('Post-hoc split: %s missing, skipping.\n', required{i});
        end
        return;
    end
end

try
    lbl_all = double(h5read(unifiedLabelsH5, '/unifiedLabels'));
    spk_all = double(h5read(spikeIdxH5,      '/spike_idx'));
    chn_all = double(h5read(channelNumH5,    '/channelNum'));
    amp_all = double(h5read(amplitudeH5,     '/amplitude'));
catch ME
    if opt.verbose, fprintf('Post-hoc split: H5 read failed (%s).\n', ME.message); end
    return;
end
lbl_all = lbl_all(:); spk_all = spk_all(:); chn_all = chn_all(:); amp_all = amp_all(:);

n = numel(lbl_all);
if n == 0 || numel(spk_all) ~= n || numel(chn_all) ~= n || numel(amp_all) ~= n
    if opt.verbose, fprintf('Post-hoc split: H5 lengths inconsistent, skipping.\n'); end
    return;
end

try
    ssData = load(sortedSamplesPath, 'sortedSamples', 'crossChannelStats');
catch ME
    if opt.verbose, fprintf('Post-hoc split: sorted_samples load failed (%s).\n', ME.message); end
    return;
end
if ~isfield(ssData, 'crossChannelStats') || ~isfield(ssData.crossChannelStats, 'unified_labels')
    return;
end
sortedSamples = ssData.sortedSamples;
unif          = ssData.crossChannelStats.unified_labels;
if ~isfield(unif, 'label') || isempty(unif.label), return; end

cfgRaw = [];
for i = 1:numel(sortedSamples)
    if ~isempty(sortedSamples{i}) && isfield(sortedSamples{i}, 'cfg')
        cfgRaw = sortedSamples{i}.cfg; break;
    end
end
if isempty(cfgRaw) || ~isfield(cfgRaw, 'samplingFrequency'), return; end
fs = cfgRaw.samplingFrequency;

sepThr   = pick(opt.separability,   cfgRaw, 'bimodalSeparability',   0.75);
% Absolute count only. A fractional bar scales with the PARENT, so the better
% sampled a unit is the harder it becomes to find a small contaminant inside
% it -- backwards. Measured here it blocked three units whose waveform
% agreement was 1.000/0.998/1.000 while the waveform gate independently
% rejected the two genuinely unsupported candidates. Left as a knob at 0.
childFr  = opt.minChildFrac;
minOvlp  = pick(opt.minTimeOverlap, cfgRaw, 'bimodalMinTimeOverlap', 0.5);
maxShift = pick(opt.maxMedianShift, cfgRaw, 'bimodalMaxMedianShift', 0.3);
numBins  = pick([],                 cfgRaw, 'bimodalNumBins',        64);
sepMinDist = pick(opt.minSepDistance, cfgRaw, 'bimodalMinSepDistance', 1.75);
ccgMin   = pick(opt.ccgIndepMin,    cfgRaw, 'bimodalCcgIndepMin',    0.5);
ccgNeed  = pick(opt.ccgMinExpected, cfgRaw, 'bimodalCcgMinExpected',  20);
assignAgr = pick(opt.assignMinAgree, cfgRaw, 'bimodalAssignMinAgree', 0.75);

% Waveform source. The saved file is row-aligned with the other H5 outputs,
% so a unit's rows can be pulled straight out; otherwise fall back to the raw
% file, which costs one filtered read per spike and so is capped.
wf = localOpenWaveforms(resSortedFolder);
raw = struct('ok', false);
if ~wf.ok
    raw = localOpenRaw(cfgRaw, channelInfoPath, outputPath, opt.verbose);
end
if ~wf.ok && ~raw.ok
    if opt.verbose
        fprintf('Post-hoc split: no saved waveforms and raw file unreachable, skipping.\n');
    end
    return;
end

half = round((cfgRaw.spikeDuration/2) * fs / 1000);
recLen = double(max(spk_all)) ;
if isfield(cfgRaw, 'num_samples') && ~isempty(cfgRaw.num_samples)
    recLen = double(cfgRaw.num_samples);
end

labels    = unif.label(:);
newLbl    = lbl_all;
nextLabel = max([labels(:); lbl_all(:)]);
addRows   = struct('parent', {}, 'label', {}, 'meanWave', {});
splitLog  = [];
nSplit    = 0;
nTested   = 0;

for u = 1:numel(labels)
    if nSplit >= opt.maxSplits, break; end
    lab  = labels(u);
    rows = find(newLbl == lab);
    if numel(rows) < opt.minSpikes, continue; end

    a = abs(amp_all(rows));
    [thrA, sepA, ~, etaA] = kiaSort_otsu_split1d(a, numBins);
    rec = struct('label', lab, 'n', numel(rows), 'eta', etaA, ...
                 'timeOverlap', NaN, 'medShift', NaN, 'waveAgree', NaN, ...
                 'sep', sepA, 'isiChild', [NaN NaN], 'presence', [NaN NaN], ...
                 'ccgRatio', NaN, 'ccgExpected', NaN, 'assignAgree', NaN, ...
                 'assignSource', 'amplitude', 'split', false);
    nTested = nTested + 1;
    % Detection only -- the waveform, drift and stability gates below decide.
    % eta is mass-weighted, so a small but cleanly separated second population
    % scores low: measured here it admitted 2 of 163 units. sep (|dmean| over
    % the summed sd) is not mass-weighted and covers that case. The 1.8 cut is
    % where waveform agreement leaves chance (0.53 below it, 0.83 above), so a
    % nomination on sep alone still has to survive the same confirmation.
    if isnan(thrA) || (etaA < sepThr && sepA < sepMinDist)
        splitLog = appendRec(splitLog, rec); continue;
    end
    thrA = kiaSort_refine_1d_2means(a, thrA);

    lowSide = a <= thrA;
    minChild = max(opt.minChildSpikes, ceil(childFr * numel(rows)));
    if sum(lowSide) < minChild || sum(~lowSide) < minChild
        splitLog = appendRec(splitLog, rec); continue;
    end

    % Drift veto: co-active modes are two cells, consecutive ones are one
    % cell drifting.
    t1 = spk_all(rows(lowSide));  t2 = spk_all(rows(~lowSide));
    q1 = prctile(t1, [5 95]);     q2 = prctile(t2, [5 95]);
    span = min(q1(2)-q1(1), q2(2)-q2(1));
    if span <= 0, splitLog = appendRec(splitLog, rec); continue; end
    rec.timeOverlap = (min(q1(2),q2(2)) - max(q1(1),q2(1))) / span;
    rec.medShift    = abs(median(t1) - median(t2)) / max(1, recLen);
    if rec.timeOverlap < minOvlp || rec.medShift > maxShift
        splitLog = appendRec(splitLog, rec); continue;
    end

    % Confirm on the waveforms: cluster them into two and require that the
    % partition matches the amplitude cut.
    try
        [W, wRows] = localUnitWaveforms(rows, wf, raw, spk_all, chn_all, half, opt.fitCap);
    catch
        W = []; wRows = [];
    end
    if isempty(W) || size(W,1) < 2*20
        splitLog = appendRec(splitLog, rec); continue;
    end
    cl = localTwoCluster(W, lowSide(wRows));
    if isempty(cl), splitLog = appendRec(splitLog, rec); continue; end
    agree = mean(cl(:) == lowSide(wRows));
    rec.waveAgree = max(agree, 1-agree);
    if rec.waveAgree < opt.waveSim
        splitLog = appendRec(splitLog, rec); continue;
    end

    % Recorded, not gated. A unit that is one clean cell plus noise splits
    % into a clean child and a dirty one, and gating on the children's ISI
    % throws that split away -- which is the case worth keeping.
    [~,~,v1] = getISIViolations(t1, fs, 2);
    [~,~,v2] = getISIViolations(t2, fs, 2);
    rec.isiChild = [v1 v2];

    % Both children have to fire throughout the parent's span. A child that
    % only exists over part of it is an epoch of the same cell, not a second
    % cell, so the split is reverted.
    rec.presence = localPresence(spk_all(rows), t1, t2, opt.nBins, opt.minBinSpikes);
    if any(rec.presence < opt.minPresence)
        splitLog = appendRec(splitLog, rec); continue;
    end

    % Accepted. Assign on the WAVEFORMS, not the amplitude threshold: the
    % threshold cuts a straight line through the cloud while the clusterer
    % follows it. Measured on this data the amplitude cut left the two
    % children partly mutually exclusive (CCG 0.82 at 1-2 ms, i.e. still one
    % neuron) where the clustered assignment reached 2.33. Only units that
    % passed every gate pay for the full read -- kmeans on the leading PCs,
    % seeded from the amplitude cut so it stays deterministic.
    side = lowSide;
    if numel(rows) <= opt.assignMaxSpikes
        if numel(wRows) == numel(rows)
            Wall = W; iAll = wRows;                 % capped read already covered it
        else
            try
                [Wall, iAll] = localUnitWaveforms(rows, wf, raw, spk_all, chn_all, ...
                                                  half, opt.assignMaxSpikes);
            catch
                Wall = []; iAll = [];
            end
        end
        if ~isempty(Wall) && size(Wall,1) >= 2*opt.minChildSpikes
            clAll = localTwoCluster(Wall, lowSide(iAll));
            % The clusterer refines the amplitude cut; it must not overrule it.
            % On a low-SNR unit kmeans returns a near-balanced partition of the
            % noise, which agrees with the cut only at chance. Real refinements
            % here agree at 0.95-1.00, so anything near 0.5 is the clusterer
            % finding nothing and the amplitude cut stands.
            if ~isempty(clAll)
                ref = lowSide(iAll);
                if mean(clAll == ref) < 0.5, clAll = ~clAll; end  % keep the seed's orientation
                agr = mean(clAll == ref);
                rec.assignAgree = agr;
                if agr >= assignAgr
                    side(iAll) = clAll;              % spikes without a waveform keep the amplitude side
                    rec.assignSource = 'waveform';
                end
            end
        end
    end
    if sum(side) < opt.minChildSpikes || sum(~side) < opt.minChildSpikes
        side = lowSide;                              % clustering collapsed -- fall back
    end

    % Independence of the final children. Two cells can fire within 1-2 ms of
    % each other; one cell cut in two cannot, because the parent's refractory
    % period forbids it. Measured on this data the separation is clean where
    % the test has power. It ABSTAINS below ccgMinExpected: a small minority
    % gives an expected count near zero, where "observed 0" carries no
    % information and rejecting on it would throw away good splits.
    [ccgR, ccgE] = localPairIndependence(spk_all(rows(side)), spk_all(rows(~side)), ...
                                         fs, opt.ccgBandMs);
    rec.ccgRatio = ccgR; rec.ccgExpected = ccgE;
    if isfinite(ccgE) && ccgE >= ccgNeed && isfinite(ccgR) && ccgR < ccgMin
        splitLog = appendRec(splitLog, rec); continue;
    end

    % The larger side keeps the label so existing curation notes stay with
    % the dominant population.
    if sum(side) >= sum(~side)
        movedMask = ~side;
    else
        movedMask = side;
    end
    nextLabel = nextLabel + 1;
    newLbl(rows(movedMask)) = nextLabel;

    mw = [];
    if ~isempty(W)
        sel = movedMask(wRows);
        if any(sel), mw = mean(W(sel,:), 1, 'omitnan'); end
    end
    addRows(end+1) = struct('parent', u, 'label', nextLabel, 'meanWave', mw); %#ok<AGROW>

    rec.split = true;
    splitLog  = appendRec(splitLog, rec);
    nSplit    = nSplit + 1;
    if opt.verbose
        fprintf(['Post-hoc split: unit %g -> %g (eta %.3f, overlap %.2f, ' ...
                 'wave %.2f, presence %.2f/%.2f, n %d/%d)\n'], ...
            lab, nextLabel, etaA, rec.timeOverlap, rec.waveAgree, ...
            rec.presence(1), rec.presence(2), sum(~movedMask), sum(movedMask));
    end
end

splitReport.nTested   = nTested;
splitReport.nSplit    = nSplit;
splitReport.log       = splitLog;
splitReport.newLabels = [addRows.label];

if nSplit == 0
    splitReport.ok = true;
    return;
end

% Every per-unit array in unified_labels has to grow with it, or the curation
% GUI (numGroups = numel(label)) never shows the new unit. Insert before
% writing anything, so the renumbering below sees the final order.
try
    % Descending parent order: an insert shifts every later index, so
    % handling the deepest first keeps the remaining parent indices valid.
    [~, ordIns] = sort([addRows.parent], 'descend');
    for k = ordIns
        unif = localInsertUnit(unif, addRows(k).parent, addRows(k).label, addRows(k).meanWave);
    end
catch ME
    if opt.verbose
        fprintf('Post-hoc split: unified_labels update failed (%s).\n', ME.message);
    end
    return;
end

% unify_spike_groups assigns label == row index, and the curation GUI relies
% on it: plotCCG / plotISI / plotDensity take a label and use it to index
% groupList, mergedFlag and the per-pair caches. Inserting a row breaks that,
% so relabel to the new row order and remap the spike labels to match.
oldLab = unif.label(:);
nU     = numel(oldLab);
maxOld = max([oldLab; newLbl(newLbl > 0)]);
lut    = nan(maxOld, 1);
lut(oldLab) = (1:nU)';
sel  = find(newLbl > 0 & newLbl <= maxOld);
m    = lut(newLbl(sel));
good = ~isnan(m);
newLbl(sel(good)) = m(good);
unif.label = (1:nU)';
splitReport.newLabels = reshape(lut([addRows.label]), 1, []);

% The labels and the unit table are two files. Back them up first and roll
% back if the second write fails, so the pair is never left disagreeing.
bk = kiaSort_backup_results(outputPath, 'presplit', ...
    {unifiedLabelsH5, sortedSamplesPath});
if ~bk.ok
    if opt.verbose, fprintf('Post-hoc split: backup failed, not writing.\n'); end
    return;
end

try
    if exist(unifiedLabelsH5, 'file'), delete(unifiedLabelsH5); end
    h5create(unifiedLabelsH5, '/unifiedLabels', size(newLbl), 'Datatype', 'double');
    h5write(unifiedLabelsH5,  '/unifiedLabels', newLbl);

    ssData.crossChannelStats.unified_labels = unif;
    crossChannelStats = ssData.crossChannelStats;
    save(sortedSamplesPath, 'crossChannelStats', '-append');
catch ME
    kiaSort_restore_results(bk);
    if opt.verbose
        fprintf('Post-hoc split: write failed (%s); rolled back from %s.\n', ...
            ME.message, bk.dir);
    end
    return;
end

splitReport.changed = true;
splitReport.ok      = true;

if opt.verbose
    fprintf('Post-hoc split: %d of %d units split.\n', nSplit, nTested);
end

end


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function v = pick(override, cfg, name, dflt)
if ~isempty(override)
    v = override;
elseif isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = dflt;
end
end


function rec = appendRec(rec0, r)
if isempty(rec0), rec = r; else, rec = rec0; rec(end+1) = r; end
end


function wf = localOpenWaveforms(resSortedFolder)
% Saved per-spike waveforms, written either by the sorter (extractWaveform)
% or by kiaSort_export_waveforms. Both are row-aligned with spike_idx.h5.
wf = struct('ok', false, 'file', '', 'dset', '', 'is3d', false, 'nT', 0, 'sz', []);
d = dir(fullfile(resSortedFolder, 'waveforms*.h5'));
for i = 1:numel(d)
    f = fullfile(resSortedFolder, d(i).name);
    try
        info = h5info(f);
        if isempty(info.Datasets), continue; end
        sz = info.Datasets(1).Dataspace.Size;
        wf.file = f;
        wf.dset = ['/' info.Datasets(1).Name];
        wf.sz   = sz;
        wf.is3d = numel(sz) >= 3;
        wf.nT   = sz(end);
        wf.ok   = true;
        return;
    catch
    end
end
end


function raw = localOpenRaw(cfg, channelInfoPath, outputPath, verbose)
raw = struct('ok', false, 'map', [], 'chanMap', [], 'nSamp', 0, 'cfg', cfg);
if ~isfield(cfg, 'fullFilePath') || ~isfield(cfg, 'numChannels') || ~isfield(cfg, 'dataType')
    return;
end
fp = char(cfg.fullFilePath);
if isempty(fp) || ~exist(fp, 'file'), return; end
try
    cfg.outputFolder = outputPath;
    raw.map = map_input_file(fp, cfg);
    raw.nSamp = size(raw.map.Data.data, 2);
    raw.chanMap = (1:cfg.numChannels)';
    if exist(channelInfoPath, 'file')
        ci = load(channelInfoPath, 'channel_mapping');
        if isfield(ci, 'channel_mapping') && ~isempty(ci.channel_mapping)
            raw.chanMap = double(ci.channel_mapping(:));
        end
    end
    raw.cfg = cfg;
    raw.ok  = raw.nSamp > 0;
catch ME
    if verbose, fprintf('Post-hoc split: raw map failed (%s).\n', ME.message); end
end
end


function [W, rowsUsed] = localUnitWaveforms(rows, wf, raw, spk_all, chn_all, half, fitCap)
% Main-channel waveform per spike. Reading the saved file is cheap, so every
% spike is used; raw reads cost a filtered snippet each and are capped.
W = []; rowsUsed = [];
if wf.ok
    W = zeros(numel(rows), wf.nT);
    keep = true(numel(rows), 1);
    runs = localRuns(rows);
    at = 0;
    for r = 1:size(runs,1)
        s = runs(r,1); c = runs(r,2);
        if wf.is3d
            blk = h5read(wf.file, wf.dset, [s 1 1], [c wf.sz(2) wf.sz(3)]);
            mid = ceil(size(blk,2)/2);
            blk = reshape(blk(:,mid,:), c, []);
        else
            blk = h5read(wf.file, wf.dset, [s 1], [c wf.sz(2)]);
        end
        W(at+1:at+c, :) = double(blk);
        at = at + c;
    end
    bad = ~any(W, 2);
    W(bad,:) = []; keep(bad) = false;
    rowsUsed = find(keep);
    return;
end
if ~raw.ok, return; end

sel = (1:numel(rows))';
if numel(sel) > fitCap
    sel = round(linspace(1, numel(sel), fitCap))';
end
cfg = raw.cfg;
pad = max(4*half, 256);
W = nan(numel(sel), 2*half+1);
for i = 1:numel(sel)
    ri = rows(sel(i));
    ch = chn_all(ri);
    if ~isfinite(ch) || ch < 1 || ch > numel(raw.chanMap), continue; end
    row = raw.chanMap(ch);
    if ~isfinite(row) || row < 1, continue; end
    s = spk_all(ri);
    a = s - half - pad;
    b = s + half + pad;
    if a < 1 || b > raw.nSamp, continue; end
    seg = localBandpass(double(raw.map.Data.data(row, a:b)), cfg);
    c0  = half + pad + 1;
    W(i,:) = seg(c0-half : c0+half);
end
good = all(isfinite(W), 2);
W = W(good, :);
rowsUsed = sel(good);
end


function [ratio, expc] = localPairIndependence(tA, tB, fs, bandMs)
% Coincidence count between the two children in a lag band, over what chance
% would predict. The band starts above the detector's dead time (spikeDistance),
% where a same-channel pair is suppressed whatever its identity. Returns NaN
% when either side is too small to measure.
ratio = NaN; expc = NaN;
% Both sides deduplicated: interp1 below needs unique sample points, and the
% swap that follows can put either side in tB. Two merged units on different
% channels can carry the same sample index, so duplicates do occur.
tA = unique(sort(double(tA(:))));
tB = unique(sort(double(tB(:))));
if numel(tA) < 20 || numel(tB) < 20, return; end
if numel(tA) > numel(tB), tmp = tA; tA = tB; tB = tmp; end
span = max([tA; tB]) - min([tA; tB]);
if span <= 0, return; end
rB = numel(tB) / (span / fs);
lo = bandMs(1) * 1e-3 * fs;
hi = bandMs(2) * 1e-3 * fs;
j  = max(min(interp1(tB, 1:numel(tB), tA, 'nearest', 'extrap'), numel(tB)), 1);
cnt = 0;
for i = 1:numel(tA)
    i1 = max(1, j(i) - 60);
    i2 = min(numel(tB), j(i) + 60);
    d  = abs(tB(i1:i2) - tA(i));
    cnt = cnt + sum(d > lo & d <= hi);
end
expc  = numel(tA) * rB * 2 * (bandMs(2) - bandMs(1)) * 1e-3;
ratio = cnt / max(expc, 1e-9);
end


function pres = localPresence(tParent, t1, t2, nBins, minBinSpikes)
% Fraction of the parent's populated bins in which each child also fires.
pres = [0 0];
lo = min(tParent); hi = max(tParent);
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo, return; end
edges = linspace(lo, hi, max(2, round(nBins)) + 1);
hp = histcounts(tParent, edges);
use = hp >= minBinSpikes;
if ~any(use), return; end
h1 = histcounts(t1, edges);
h2 = histcounts(t2, edges);
pres = [mean(h1(use) > 0), mean(h2(use) > 0)];
end


function y = localBandpass(x, cfg)
% Same design as kiaSort_filter_signal, but the coefficients are cached and
% neither the pool nor the GPU is touched -- this runs once per spike.
persistent b a key
k = sprintf('%g_%g_%g', cfg.bandpass(1), cfg.bandpass(2), cfg.samplingFrequency);
if isempty(key) || ~strcmp(key, k)
    fsl = cfg.samplingFrequency;
    up  = min(1.05 * cfg.bandpass(2), 0.95 * (fsl/2));
    [b, a] = butter(4, [cfg.bandpass(1), up] / (fsl/2), 'bandpass');
    key = k;
end
y = filtfilt(b, a, x(:))';
end


function runs = localRuns(idx)
% Contiguous [start count] runs so scattered rows are read in few calls.
idx = idx(:);
brk = [1; find(diff(idx) ~= 1) + 1; numel(idx)+1];
runs = zeros(numel(brk)-1, 2);
for i = 1:numel(brk)-1
    runs(i,1) = idx(brk(i));
    runs(i,2) = brk(i+1) - brk(i);
end
end


function cl = localTwoCluster(W, seedLow)
% Two-means on the leading PCs, seeded from the amplitude cut so the result
% is deterministic and comparable with it.
cl = [];
W = W - mean(W, 1);
nComp = min(5, min(size(W)) - 1);
if nComp < 1, return; end
try
    [~, F] = pca(W, 'Algorithm', 'svd', 'NumComponents', nComp);
catch
    return;
end
if isempty(F) || size(F,1) < 4, return; end
seedLow = logical(seedLow(:));
if ~any(seedLow) || ~any(~seedLow), return; end
C = [mean(F(seedLow,:), 1); mean(F(~seedLow,:), 1)];
if any(~isfinite(C(:))), return; end
try
    k = kmeans(F, 2, 'Start', C, 'MaxIter', 300);
catch
    return;
end
cl = (k == 1);
end


function unif = localInsertUnit(unif, parentIdx, newLabel, meanWave)
% Placed directly after its parent, so a split unit stays with the rest of
% its channel instead of landing at the end of the list. Everything that
% describes where the unit lives is mirrored from the parent; only the label
% and the mean waveform differ.
vecFields = {'label', 'channelID', 'labelInChannel', 'detectblity', ...
             'mainNegativePolarity', 'sideNegativePolarity'};
for i = 1:numel(vecFields)
    f = vecFields{i};
    if isfield(unif, f) && isvector(unif.(f))
        unif.(f) = unif.(f)(:);
    end
end
at = parentIdx + 1;
for i = 1:numel(vecFields)
    f = vecFields{i};
    if ~isfield(unif, f) || numel(unif.(f)) < parentIdx, continue; end
    v = unif.(f);
    if strcmp(f, 'label')
        nv = newLabel;
    else
        nv = v(parentIdx);
    end
    unif.(f) = [v(1:at-1); cast(nv, 'like', v); v(at:end)];
end
if isfield(unif, 'meanWaveforms') && ~isempty(unif.meanWaveforms)
    MW = unif.meanWaveforms;
    if ndims(MW) == 3
        % The unit's own channel is the CENTRE of the footprint, not slot 1:
        % the GUI reads it as ch - channelList(g) + numChannelPlot + 1.
        row = MW(parentIdx, :, :);
        cCh = ceil(size(MW, 2) / 2);
        row(1, cCh, :) = reshape( ...
            localCentreFit(reshape(row(1, cCh, :), 1, []), meanWave), 1, 1, []);
        MW = cat(1, MW(1:at-1, :, :), row, MW(at:end, :, :));
    else
        row = localCentreFit(MW(parentIdx, :), meanWave);
        MW = cat(1, MW(1:at-1, :), row, MW(at:end, :));
    end
    unif.meanWaveforms = MW;
end
end


function w = localCentreFit(w, mw)
% meanWaveforms spans spikeDuration while the saved per-spike window spans
% clusteringSpikeDuration, so the two rarely have the same length. Put the
% child's mean in the centre, where the amplitude difference that justified
% the split actually lives, and leave the parent's flanks -- same channel,
% same neighbourhood.
if isempty(mw), return; end
mw = double(mw(:)');
n  = numel(w);
m  = numel(mw);
if m == n
    w = mw;
elseif m > n
    o = floor((m - n) / 2);
    w = mw(o+1 : o+n);
else
    o = floor((n - m) / 2);
    w(o+1 : o+m) = mw;
end
end
