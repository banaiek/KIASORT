function report = kiaSort_residual_units(outputPath, varargin)
%KIASORT_RESIDUAL_UNITS  Recover units hidden inside poorly-fitting spikes.
%
%   report = kiaSort_residual_units(outputPath, ...)
%
%   Template matching is argmin with no reject option: every detected spike
%   is handed to the least-bad template however badly it fits, and the only
%   filter is an amplitude band whose upper edge is open (ampBandHighFactor
%   Inf, because a finite cap truncated drifting units). So a neuron with no
%   template of its own -- one that appeared after the sample window, or was
%   dropped at the sample stage -- is absorbed spike by spike into whatever
%   unit sits nearest, and shows up as a second amplitude mode inside it.
%
%   This pass reads the match distance the sorter already computed
%   (features.h5), pools the spikes that fit badly, and asks what the pool
%   IS. The question is temporal, not morphological:
%
%     drift      the host stops as the pool appears -- the same spikes moved
%                across the amplitude cut, so the total rate is conserved.
%                Leave it alone; it is one cell and already one unit.
%     foreign    the host carries on unchanged while the pool fires. Two
%                cells coexist, so the pool becomes a unit of its own.
%     noise      incoherent, or too few, or refractory-violating. Leave it.
%
%   Shape cannot make this call. Drift moves the electrode, which changes
%   spatial sampling and so the waveform; measured here a genuine drifting
%   unit held shape correlation 0.995 across a 2.3x amplitude change, while
%   a foreign population sat at 0.939 of its host. Both near 1.
%
%   NOTHING IS EVER DELETED. The pass only moves flagged spikes from an
%   existing unit into a new one, so the worst case is an extra unit, never
%   lost data. It is off unless cfg.residualUnits is set.
%
%   Name/Value:
%       'zDefer'        (5)     fit-outlier cut, in robust sd of the unit's
%                               own match distances (see below)
%       'minPoolSpikes' (200)   pool smaller than this is never promoted
%       'rateKeepMin'   (0.6)   host rate inside the pool's window over
%                               outside; below this the host was depleted,
%                               i.e. drift, and the pool is left alone.
%                               Undefined when the pool fires in every bin,
%                               which is exactly what a foreign cell does --
%                               hence the regression below, which always works
%       'driftCorr'    (-0.5)   host-vs-pool rate correlation across bins
%       'driftSlope'   (-0.3)   normalised slope of that regression. Drift
%                               trades spikes one for one, so the host falls
%                               as the pool rises (measured -0.87 / -0.61 on a
%                               real drifting unit, +0.11 / +0.01 on a real
%                               foreign one). Both must fire to call drift
%       'minAmpRatio'   (2.0)   pool median amplitude over the host's
%       'minAmpGap'     (1.2)   pool 10th percentile over the host's 95th.
%                               The ratio alone is not enough: across 91 pools
%                               measured here the median ratio was 1.4 but the
%                               median gap 0.81, i.e. most pools are the upper
%                               tail of their host, not a separate cell. Only
%                               a pool that clears the host's distribution is
%                               promoted.
%       'maxPoolIsi'    (1.0)   pool's own ISI violation cap (%, at 2 ms)
%       'minCoherence'  (0.85)  mean correlation of pool spikes to the pool
%                               mean, when waveforms are available
%       'ccgIndepMin'   (0.5)   pool and host must fire independently
%       'ccgBandMs'     ([0.4 2]) lag band, above the detector dead time
%       'nBins'         (20)    bins for the coexistence test
%       'minBinSpikes'  (5)     host spikes for a bin to count
%       'maxPerUnit'    (2)     promotions per host unit, so a pool that
%                               splits cleanly can yield both halves
%       'spikeCap'      (4000)  waveforms read per pool
%       'rejectOutliers'(true)  unassign lone spikes far outside their unit's
%                               amplitude range once the pools are settled
%       'outlierFactor' (3)     multiple of the unit's 99th percentile. Only
%                               applied when there are fewer than
%                               minPoolSpikes of them, so a drifting
%                               population -- always many spikes -- can never
%                               be rejected this way. Measured here: 54 such
%                               spikes in 7.09M, median 1 per unit, never more
%                               than 2, i.e. isolated events rather than a
%                               population. They are set to -1 (unassigned),
%                               not removed: the row stays in every H5 output
%       'verbose'       (false)

p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addParameter('zDefer',          5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minPoolSpikes', 200, @(x) isscalar(x) && isnumeric(x));
p.addParameter('rateKeepMin',   0.6, @(x) isscalar(x) && isnumeric(x));
p.addParameter('driftCorr',    -0.5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('driftSlope',   -0.3, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minAmpRatio',   2.0, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minAmpGap',     1.2, @(x) isscalar(x) && isnumeric(x));
p.addParameter('maxPoolIsi',    1.0, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minCoherence', 0.85, @(x) isscalar(x) && isnumeric(x));
p.addParameter('ccgIndepMin',   0.5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('ccgBandMs', [0.4 2], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('nBins',          20, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minBinSpikes',    5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('maxPerUnit',      2, @(x) isscalar(x) && isnumeric(x));
p.addParameter('spikeCap',     4000, @(x) isscalar(x) && isnumeric(x));
p.addParameter('rejectOutliers', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('outlierFactor',    3, @(x) isscalar(x) && isnumeric(x));
p.addParameter('verbose',     false, @(x) islogical(x) || isnumeric(x));
p.parse(outputPath, varargin{:});
opt = p.Results;
opt.verbose = logical(opt.verbose);
opt.rejectOutliers = logical(opt.rejectOutliers);

report = struct('nPromoted', 0, 'nTested', 0, 'nOutliers', 0, 'newLabels', [], ...
                'changed', false, 'ok', false, 'log', []);

outputPath = char(outputPath);
resFolder  = fullfile(outputPath, 'RES_Sorted');
ssPath     = fullfile(outputPath, 'Sorted_Samples', 'sorted_samples.mat');
H = struct('lbl','unifiedLabels', 'spk','spike_idx', 'chn','channelNum', ...
           'amp','amplitude', 'ft','features');
f = fieldnames(H);
paths = struct();
for i = 1:numel(f)
    paths.(f{i}) = fullfile(resFolder, [H.(f{i}) '.h5']);
    if ~exist(paths.(f{i}), 'file')
        if opt.verbose, fprintf('Residual units: %s missing, skipping.\n', paths.(f{i})); end
        return;
    end
end
if ~exist(ssPath, 'file')
    if opt.verbose, fprintf('Residual units: sorted_samples.mat missing, skipping.\n'); end
    return;
end

try
    lbl = double(h5read(paths.lbl, ['/' H.lbl])); lbl = lbl(:);
    spk = double(h5read(paths.spk, ['/' H.spk])); spk = spk(:);
    chn = double(h5read(paths.chn, ['/' H.chn])); chn = chn(:);
    amp = double(h5read(paths.amp, ['/' H.amp])); amp = amp(:);
    ft  = double(h5read(paths.ft,  ['/' H.ft ])); ft  = ft(:);
catch ME
    if opt.verbose, fprintf('Residual units: H5 read failed (%s).\n', ME.message); end
    return;
end
n = numel(lbl);
if n == 0 || any([numel(spk) numel(chn) numel(amp) numel(ft)] ~= n)
    if opt.verbose, fprintf('Residual units: H5 lengths inconsistent, skipping.\n'); end
    return;
end

try
    ssData = load(ssPath, 'sortedSamples', 'crossChannelStats');
catch ME
    if opt.verbose, fprintf('Residual units: sorted_samples load failed (%s).\n', ME.message); end
    return;
end
if ~isfield(ssData, 'crossChannelStats') || ~isfield(ssData.crossChannelStats, 'unified_labels')
    return;
end
unif = ssData.crossChannelStats.unified_labels;
if ~isfield(unif, 'label') || isempty(unif.label), return; end

cfgRaw = [];
for i = 1:numel(ssData.sortedSamples)
    if ~isempty(ssData.sortedSamples{i}) && isfield(ssData.sortedSamples{i}, 'cfg')
        cfgRaw = ssData.sortedSamples{i}.cfg; break;
    end
end
if isempty(cfgRaw) || ~isfield(cfgRaw, 'samplingFrequency'), return; end
fs = cfgRaw.samplingFrequency;

wfSrc = struct('ok', false);
try
    wfSrc = kiaSort_waveform_source(outputPath, cfgRaw, opt.verbose);
catch
end

labels    = unif.label(:);
newLbl    = lbl;
nextLabel = max([labels(:); lbl(:)]);
addRows   = struct('parent', {}, 'label', {}, 'meanWave', {});
lg        = [];

for u = 1:numel(labels)
    lab  = labels(u);
    rows = find(newLbl == lab);
    if numel(rows) < 4*opt.minPoolSpikes, continue; end

    % ---- flag the fit outliers ------------------------------------------
    % Scale from the LEFT of the median only. A symmetric spread is inflated
    % by the very contamination this is looking for -- measured here a unit
    % holding a foreign population 2.9x its own amplitude scored 0.0% of its
    % spikes as outliers against its own symmetric MAD.
    fu  = ft(rows);
    md  = median(fu);
    low = fu(fu < md);
    if isempty(low), continue; end
    sd  = 1.4826 * median(md - low);
    if ~isfinite(sd) || sd <= 0, continue; end
    pool = (fu - md) / sd > opt.zDefer;
    if sum(pool) < opt.minPoolSpikes, continue; end

    report.nTested = report.nTested + 1;

    % A pool can hold more than one thing (a steady foreign cell plus a late
    % event). Offer the whole pool and, when its amplitude splits cleanly,
    % each half; the first candidate that qualifies is taken.
    cands = local_poolCandidates(abs(amp(rows(pool))), opt.minPoolSpikes);
    poolIdx = find(pool);
    nProm = 0;
    for ci = 1:numel(cands)
        if nProm >= opt.maxPerUnit, break; end
        sub  = poolIdx(cands{ci});
        pRow = rows(sub);
        % The comparison is against the unit's CORE, i.e. every pool spike
        % removed -- not merely this candidate. Leaving the other half of a
        % split pool in the host inflates its 95th percentile and the gap
        % test then rejects a candidate that is plainly separated from the
        % core, so only one half of a two-component pool ever got promoted.
        hRow = rows(~pool);
        [okProm, rec] = local_adjudicate(pRow, hRow, spk, amp, chn, fs, wfSrc, opt);
        rec.label  = lab;
        rec.nPool  = numel(pRow);
        rec.cand   = ci;
        rec.promoted = okProm;
        lg = local_append(lg, rec);
        if ~okProm, continue; end

        nextLabel = nextLabel + 1;
        newLbl(pRow) = nextLabel;
        mw = [];
        if wfSrc.ok
            try
                Wp = kiaSort_read_waveforms(wfSrc, pRow, spk, chn, min(opt.spikeCap, numel(pRow)), 'spread');
                if ~isempty(Wp), mw = mean(Wp, 1, 'omitnan'); end
            catch
            end
        end
        addRows(end+1) = struct('parent', u, 'label', nextLabel, 'meanWave', mw); %#ok<AGROW>
        nProm = nProm + 1;
        report.nPromoted = report.nPromoted + 1;
        if opt.verbose
            fprintf(['Residual units: %g -> %g  (n %d, amp x%.2f, rateKeep %s, ' ...
                     'gap %.2f, slope %+.2f, corr %+.2f, poolISI %.2f%%, coh %.2f, ccg %.2f/exp %.0f)\n'], ...
                lab, nextLabel, numel(pRow), rec.ampRatio, num2str(rec.rateKeep,'%.2f'), ...
                rec.ampGap, rec.rateSlope, rec.rateCorr, rec.poolIsi, rec.coherence, ...
                rec.ccgRatio, rec.ccgExpected);
        end
    end
end

% ---- lone amplitude outliers ------------------------------------------
% Template matching is argmin with no reject option, so a spike that belongs
% to nothing lands on the least-bad template however absurd the fit. Those
% survive the pool logic because there are only ever one or two of them per
% unit -- far below the count needed to argue they are a cell. Run AFTER the
% promotions, so the percentile they are judged against is the host's own
% once any real population has been lifted out.
outRows = [];
if opt.rejectOutliers
    for u = 1:numel(labels)
        rows = find(newLbl == labels(u));
        if numel(rows) < 20, continue; end
        a  = abs(amp(rows));
        cut = opt.outlierFactor * prctile(a, 99);
        bad = rows(a > cut);
        % Many of them would be a population, and a drifting unit is always a
        % population. Leave those to the pool logic above.
        if isempty(bad) || numel(bad) >= opt.minPoolSpikes, continue; end
        outRows = [outRows; bad(:)]; %#ok<AGROW>
    end
end
if ~isempty(outRows)
    newLbl(outRows) = -1;
    report.nOutliers = numel(outRows);
    if opt.verbose
        fprintf('Residual units: %d lone amplitude outliers unassigned.\n', numel(outRows));
    end
end

report.log = lg;
if report.nPromoted == 0 && report.nOutliers == 0
    report.ok = true;
    if opt.verbose, fprintf('Residual units: nothing promoted (%d pools tested).\n', report.nTested); end
    return;
end

% ---- write, atomically across the two files -----------------------------
bk = kiaSort_backup_results(outputPath, 'residual', {paths.lbl, ssPath});
if ~bk.ok
    if opt.verbose, fprintf('Residual units: backup failed, not writing.\n'); end
    return;
end
try
    % New units are appended; kiaSort_compact_unit_table renumbers so that
    % label == row index again, which the curation GUI relies on.
    for k = 1:numel(addRows)
        unif = local_appendUnit(unif, addRows(k).parent, addRows(k).label, addRows(k).meanWave);
    end
    if isempty(addRows), addRows = struct('parent',{},'label',{},'meanWave',{}); end
    [unif, newLbl, ~, remap] = kiaSort_compact_unit_table(unif, newLbl);
    if exist(paths.lbl, 'file'), delete(paths.lbl); end
    h5create(paths.lbl, ['/' H.lbl], size(newLbl), 'Datatype', 'double');
    h5write(paths.lbl,  ['/' H.lbl], newLbl);
    ssData.crossChannelStats.unified_labels = unif;
    crossChannelStats = ssData.crossChannelStats; %#ok<NASGU>
    save(ssPath, 'crossChannelStats', '-append');
    if isa(remap, 'containers.Map')
        nl = [addRows.label];
        out = nan(size(nl));
        for k = 1:numel(nl)
            if isKey(remap, nl(k)), out(k) = remap(nl(k)); end
        end
        report.newLabels = out;
    else
        report.newLabels = [addRows.label];
    end
catch ME
    kiaSort_restore_results(bk);
    if opt.verbose
        fprintf('Residual units: write failed (%s); rolled back from %s.\n', ME.message, bk.dir);
    end
    return;
end

report.changed = true;
report.ok      = true;
if opt.verbose
    fprintf('Residual units: %d promoted from %d pools, %d outliers unassigned.\n', ...
        report.nPromoted, report.nTested, report.nOutliers);
end
end


% =========================================================================
% Local functions
% =========================================================================

function c = local_poolCandidates(a, minN)
%LOCAL_POOLCANDIDATES  What to offer the adjudicator, in priority order.
%
% A pool can hold more than one thing -- a steady foreign cell plus a later
% amplitude event, say. When the pool's own amplitude splits cleanly, the
% HALVES are offered and the whole is not: promoting the whole first would
% hand back a unit as bimodal as the one the pass exists to clean up
% (measured: promoting intact gave a new unit at eta 0.92).
c = {(1:numel(a))'};
if numel(a) < 2*minN, return; end
try
    [thr, ~, ~, eta] = kiaSort_otsu_split1d(a, 64);
catch
    return;
end
if isnan(thr) || eta < 0.75, return; end
lo = find(a <= thr); hi = find(a > thr);
if numel(lo) >= minN && numel(hi) >= minN
    c = {lo, hi};                 % clean split: the halves replace the whole
elseif numel(lo) >= minN
    c = {lo, c{1}};
elseif numel(hi) >= minN
    c = {hi, c{1}};
end
end


function [ok, rec] = local_adjudicate(pRow, hRow, spk, amp, chn, fs, wfSrc, opt)
%LOCAL_ADJUDICATE  Is this pool a second cell, or the host having drifted?
rec = struct('ampRatio', NaN, 'ampGap', NaN, 'rateKeep', NaN, 'rateCorr', NaN, ...
             'rateSlope', NaN, 'poolIsi', NaN, ...
             'coherence', NaN, 'ccgRatio', NaN, 'ccgExpected', NaN, ...
             'why', "", 'label', NaN, 'nPool', NaN, 'cand', NaN, 'promoted', false);
ok = false;
if numel(pRow) < opt.minPoolSpikes || numel(hRow) < opt.minPoolSpikes
    rec.why = "too small"; return;
end

tp = spk(pRow); th = spk(hRow);
ap = abs(amp(pRow)); ah = abs(amp(hRow));
rec.ampRatio = median(ap) / max(median(ah), eps);

% ---- coexistence: the decisive test ------------------------------------
% Drift moves the host's own spikes into the pool, so the host thins out
% exactly where the pool appears. A second cell leaves the host untouched.
lo = min([tp; th]); hi = max([tp; th]);
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo, rec.why = "degenerate span"; return; end
edges = linspace(lo, hi, max(2, round(opt.nBins)) + 1);
hh = histcounts(th, edges);
hp = histcounts(tp, edges);
use = (hh + hp) >= opt.minBinSpikes;
if sum(use) < 4, rec.why = "too few usable bins"; return; end

% Direct form, when the pool leaves any bin empty.
act = hp > 0 & use;
oth = hp == 0 & use;
if any(act) && any(oth)
    rOut = sum(hh(oth)) / sum(oth);
    if rOut > 0, rec.rateKeep = (sum(hh(act)) / sum(act)) / rOut; end
end

% General form. A foreign cell that fires all session leaves no empty bin,
% so rateKeep is undefined exactly where it is most needed. Regressing the
% host's rate on the pool's works either way: drift trades spikes one for
% one, so the host falls as the pool rises.
x = hp(use); y = hh(use);
if std(x) > 0 && std(y) > 0
    rec.rateCorr = corr(double(x(:)), double(y(:)));
    b = [ones(numel(x),1) double(x(:))] \ double(y(:));
    rec.rateSlope = b(2) * mean(x) / max(mean(y), eps);
end

isDrift = (isfinite(rec.rateKeep) && rec.rateKeep < opt.rateKeepMin) || ...
          (isfinite(rec.rateCorr) && isfinite(rec.rateSlope) && ...
           rec.rateCorr < opt.driftCorr && rec.rateSlope < opt.driftSlope);
if isDrift
    rec.why = "host depleted -> drift"; return;    % leave it in the host
end

rec.ampGap = prctile(ap, 10) / max(prctile(ah, 95), eps);
if ~isfinite(rec.ampRatio) || rec.ampRatio < opt.minAmpRatio || ...
        ~isfinite(rec.ampGap) || rec.ampGap < opt.minAmpGap
    rec.why = "amplitude not separated"; return;
end

try
    [~, ~, rec.poolIsi] = getISIViolations(tp, fs, 2);
catch
    rec.poolIsi = Inf;
end
if ~isfinite(rec.poolIsi) || rec.poolIsi > opt.maxPoolIsi
    rec.why = "pool refractory violations"; return;
end

% ---- independence -------------------------------------------------------
[rec.ccgRatio, rec.ccgExpected] = local_coincidence(tp, th, fs, opt.ccgBandMs);
if isfinite(rec.ccgExpected) && rec.ccgExpected >= 20 && ...
        isfinite(rec.ccgRatio) && rec.ccgRatio < opt.ccgIndepMin
    rec.why = "not independent of host"; return;   % same cell, not a new one
end

% ---- the pool has to be one thing --------------------------------------
if wfSrc.ok
    try
        W = kiaSort_read_waveforms(wfSrc, pRow, spk, chn, min(opt.spikeCap, numel(pRow)), 'spread');
    catch
        W = [];
    end
    if ~isempty(W) && size(W,1) >= 20
        mw = mean(W, 1, 'omitnan');
        if any(isfinite(mw)) && norm(mw) > 0
            r = zeros(size(W,1), 1);
            for i = 1:size(W,1)
                v = W(i,:);
                if norm(v) > 0, r(i) = corr(v(:), mw(:)); end
            end
            rec.coherence = median(r, 'omitnan');
            if ~isfinite(rec.coherence) || rec.coherence < opt.minCoherence
                rec.why = "pool not coherent"; return;
            end
        end
    end
end

ok = true;
rec.why = "promoted";
end


function [ratio, expc] = local_coincidence(tA, tB, fs, bandMs)
%LOCAL_COINCIDENCE  Coincidences in a lag band over chance.
ratio = NaN; expc = NaN;
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


function unif = local_appendUnit(unif, parentIdx, newLabel, meanWave)
%LOCAL_APPENDUNIT  New unit at the end of the table, mirroring its parent.
vecFields = {'label', 'channelID', 'labelInChannel', 'detectblity', ...
             'mainNegativePolarity', 'sideNegativePolarity'};
for i = 1:numel(vecFields)
    f = vecFields{i};
    if isfield(unif, f) && isvector(unif.(f)), unif.(f) = unif.(f)(:); end
end
for i = 1:numel(vecFields)
    f = vecFields{i};
    if ~isfield(unif, f) || numel(unif.(f)) < parentIdx, continue; end
    v = unif.(f);
    if strcmp(f, 'label'), nv = newLabel; else, nv = v(parentIdx); end
    unif.(f) = [v; cast(nv, 'like', v)];
end
if isfield(unif, 'meanWaveforms') && ~isempty(unif.meanWaveforms)
    MW = unif.meanWaveforms;
    if ndims(MW) == 3
        row = MW(parentIdx, :, :);
        cCh = ceil(size(MW, 2) / 2);
        row(1, cCh, :) = reshape( ...
            local_centreFit(reshape(row(1, cCh, :), 1, []), meanWave), 1, 1, []);
        MW = cat(1, MW, row);
    else
        MW = cat(1, MW, local_centreFit(MW(parentIdx, :), meanWave));
    end
    unif.meanWaveforms = MW;
end
end


function w = local_centreFit(w, mw)
%LOCAL_CENTREFIT  Put the pool's mean in the centre of the parent's window.
if isempty(mw), return; end
mw = double(mw(:)');
n = numel(w); m = numel(mw);
if m == n
    w = mw;
elseif m > n
    o = floor((m - n) / 2); w = mw(o+1 : o+n);
else
    o = floor((n - m) / 2); w(o+1 : o+m) = mw;
end
end


function lg = local_append(lg, rec)
if isempty(lg), lg = rec; else, lg(end+1) = rec; end
end
