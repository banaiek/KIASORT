function report = kiaSort_echo_units(outputPath, varargin)
%KIASORT_ECHO_UNITS  Unassign units that are echoes of a neighbour's spikes.
%
%   report = kiaSort_echo_units(outputPath, ...)
%
%   A spike with a second threshold crossing 0.5-1.5 ms after its trough (a
%   slow repolarisation, or what is left once the first template has been
%   subtracted before the second detection round) is detected twice: at the
%   trough, and again at the later deflection, which lies outside the
%   detector dead time. The second detections have a shape of their own, so
%   sample sorting gives them a class and the sort a unit. That unit holds no
%   spikes of its own. Its mean waveform is dominated by the neighbour's
%   spike sitting at a fixed offset from the centre, and most of its spikes
%   follow a detected spike on the same or an adjacent channel by that
%   offset. Overlap removal cannot see it: the offset is beyond its +-0.5 ms
%   coincidence window.
%
%   Per unit:
%     1) Mean waveform over +-windowMs on the unit's own channel, read from
%        the raw file (the stored 1 ms window ends inside the offset).
%        Candidate when a deflection >= minOffsetMs from the centre is at
%        least minOffRatio of the centre deflection. The detection may sit on
%        the neighbour's repolarisation hump, which is as large as the trough
%        left in the mean once echoes of several parents are averaged, so
%        the off-centre deflection need not dominate -- but a real unit's
%        other phases are small next to the one it is aligned on.
%     2) Timing: the fraction of the unit's spikes that have a spike of a
%        live unit (channels within leaderChannels of the unit's own) at that
%        offset, +-lagTolMs. Confirmed when it is >= minEchoFrac and at least
%        minChanceRatio x the Poisson chance level.
%   Candidates are taken in descending order of the off-centre ratio and a
%   dropped unit leaves the leader pool, so when both phases of one spike
%   were detected the detection on the smaller phase is the one dropped, and
%   the larger one is not then matched against it.
%   A confirmed unit is unassigned whole (label -1). Measured on a Utah
%   array, the spikes of such a unit that had no matching leader were echoes
%   of unlabelled events with the same off-centre waveform, not a cell.
%
%   Only unifiedLabels.h5 and crossChannelStats.unified_labels change, and
%   only when at least one unit is dropped; the unit table is compacted so
%   that label == row index still holds.
%
%   Name/Value:
%       'windowMs'        (1.5)   half-window of the mean waveform (ms)
%       'minOffsetMs'     (0.4)   above the detector dead time (0.375 ms)
%       'minOffRatio'     (0.5)   off-centre deflection / centre deflection;
%                                 nomination only -- on a 189-unit reference
%                                 sort no unit reached it, and the timing
%                                 test is what confirms
%       'lagTolMs'        (0.3)   tolerance around the offset (ms)
%       'minEchoFrac'     (0.3)   min fraction of spikes with a leader
%       'minChanceRatio'  (3)     ...and at least this many times chance
%       'leaderChannels'  ([])    channel radius for leaders; default
%                                 cfg.duplicateSearchChannels, else 2
%       'capN'            (300)   waveforms read per unit (a mean; the
%                                 timing test uses every spike)
%       'minSpikes'       (50)    smaller units are left alone
%       'apply'           (true)  false = report only, nothing written
%       'verbose'         (false)
%
%   report.nTested / nDropped / dropped (labels) / log (per unit) /
%   changed / ok

p = inputParser;
p.addRequired('outputPath', @(x) ischar(x) || isstring(x));
p.addParameter('windowMs',       1.5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minOffsetMs',    0.4, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minOffRatio',    0.5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('lagTolMs',       0.3, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minEchoFrac',    0.3, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minChanceRatio',   3, @(x) isscalar(x) && isnumeric(x));
p.addParameter('leaderChannels',  [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
p.addParameter('capN',           300, @(x) isscalar(x) && isnumeric(x));
p.addParameter('minSpikes',       50, @(x) isscalar(x) && isnumeric(x));
p.addParameter('apply',         true, @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',      false, @(x) islogical(x) || isnumeric(x));
p.parse(outputPath, varargin{:});
opt = p.Results;
opt.apply   = logical(opt.apply);
opt.verbose = logical(opt.verbose);

report = struct('nTested', 0, 'nDropped', 0, 'dropped', [], ...
                'changed', false, 'ok', false, 'log', []);

outputPath = char(outputPath);
resFolder  = fullfile(outputPath, 'RES_Sorted');
ssPath     = fullfile(outputPath, 'Sorted_Samples', 'sorted_samples.mat');
H = struct('lbl','unifiedLabels', 'spk','spike_idx', 'chn','channelNum');
f = fieldnames(H);
paths = struct();
for i = 1:numel(f)
    paths.(f{i}) = fullfile(resFolder, [H.(f{i}) '.h5']);
    if ~exist(paths.(f{i}), 'file')
        if opt.verbose, fprintf('Echo units: %s missing, skipping.\n', paths.(f{i})); end
        return;
    end
end
if ~exist(ssPath, 'file')
    if opt.verbose, fprintf('Echo units: sorted_samples.mat missing, skipping.\n'); end
    return;
end

try
    lbl = double(h5read(paths.lbl, ['/' H.lbl])); lbl = lbl(:);
    spk = double(h5read(paths.spk, ['/' H.spk])); spk = spk(:);
    chn = double(h5read(paths.chn, ['/' H.chn])); chn = chn(:);
catch ME
    if opt.verbose, fprintf('Echo units: H5 read failed (%s).\n', ME.message); end
    return;
end
n = numel(lbl);
if n == 0 || numel(spk) ~= n || numel(chn) ~= n
    if opt.verbose, fprintf('Echo units: H5 lengths inconsistent, skipping.\n'); end
    return;
end

try
    ssData = load(ssPath, 'sortedSamples', 'crossChannelStats');
catch ME
    if opt.verbose, fprintf('Echo units: sorted_samples load failed (%s).\n', ME.message); end
    return;
end
if ~isfield(ssData, 'crossChannelStats') || ~isfield(ssData.crossChannelStats, 'unified_labels')
    return;
end
unif = ssData.crossChannelStats.unified_labels;
if ~isfield(unif, 'label') || isempty(unif.label), report.ok = true; return; end

cfgRaw = [];
for i = 1:numel(ssData.sortedSamples)
    if ~isempty(ssData.sortedSamples{i}) && isfield(ssData.sortedSamples{i}, 'cfg')
        cfgRaw = ssData.sortedSamples{i}.cfg; break;
    end
end
if isempty(cfgRaw) || ~isfield(cfgRaw, 'samplingFrequency'), return; end
fs = cfgRaw.samplingFrequency;

% Raw only: the saved window is the clustering window (1 ms), which ends
% inside the offsets this pass looks for.
src = struct('ok', false);
try
    src = kiaSort_waveform_source(outputPath, cfgRaw, opt.verbose, true);
catch
end
if ~src.ok
    if opt.verbose, fprintf('Echo units: raw file unreachable, skipping.\n'); end
    report.ok = true;
    return;
end
half     = round(opt.windowMs * fs / 1000);
src.half = half;
centre   = half + 1;
jit      = max(1, round(0.1e-3 * fs));      % alignment jitter allowed at the centre
tolS     = opt.lagTolMs * fs / 1000;

radius = opt.leaderChannels;
if isempty(radius)
    radius = 2;
    if isfield(cfgRaw, 'duplicateSearchChannels') && ~isempty(cfgRaw.duplicateSearchChannels)
        radius = double(cfgRaw.duplicateSearchChannels);
    end
end
durSec = double(max(spk)) / fs;

labels   = unif.label(:);
dropMask = false(n, 1);
dropped  = [];
lg = struct('label', {}, 'channel', {}, 'n', {}, 'peakAtMs', {}, ...
            'offRatio', {}, 'echoFrac', {}, 'chance', {}, 'dropped', {});

% ---- 1) where does each unit's energy sit? ------------------------------
offMask = abs((1:2*half+1) - centre) >= opt.minOffsetMs * fs / 1000;
cwin    = max(1, centre - jit):min(2*half+1, centre + jit);
for u = 1:numel(labels)
    lab  = labels(u);
    rows = find(lbl == lab);
    if numel(rows) < opt.minSpikes, continue; end
    ch  = mode(chn(rows));
    rec = struct('label', lab, 'channel', ch, 'n', numel(rows), 'peakAtMs', NaN, ...
                 'offRatio', NaN, 'echoFrac', NaN, 'chance', NaN, 'dropped', false);
    report.nTested = report.nTested + 1;
    try
        W = kiaSort_read_waveforms(src, rows, spk, chn, opt.capN, 'spread');
    catch
        W = [];
    end
    if ~isempty(W) && size(W, 1) >= 20 && size(W, 2) == 2*half + 1
        mw = mean(W, 1);
        cMag = max(abs(mw(cwin)));
        oAbs = abs(mw); oAbs(~offMask) = 0;
        [oMag, ix] = max(oAbs);
        if cMag > 0
            rec.peakAtMs = (ix - centre) / fs * 1000;
            rec.offRatio = oMag / cMag;
        end
    end
    lg(end+1) = rec; %#ok<AGROW>
end

% ---- 2) timing, largest off-centre ratio first ---------------------------
cand = find([lg.offRatio] >= opt.minOffRatio);
[~, order] = sort([lg(cand).offRatio], 'descend');
cand = cand(order);
for k = cand
    rec  = lg(k);
    rows = find(lbl == rec.label);
    live = lbl > 0 & ~dropMask;
    lead = spk(live & abs(chn - rec.channel) <= radius & lbl ~= rec.label);
    lead = unique(lead(:));
    if numel(lead) < 2, continue; end
    tB    = spk(rows) + rec.peakAtMs * fs / 1000;   % where the leader should be
    edges = [-Inf; lead; Inf];
    nHi   = discretize(tB + tolS, edges) - 1;        % leaders <= hi
    nLo   = discretize(tB - tolS, edges) - 1;        % leaders <  lo (times are integers)
    hit   = nHi > nLo;
    rec.echoFrac = mean(hit);
    rec.chance   = 1 - exp(-(numel(lead) / durSec) * (2 * opt.lagTolMs / 1000));
    if rec.echoFrac >= opt.minEchoFrac && rec.echoFrac >= opt.minChanceRatio * rec.chance
        rec.dropped    = true;
        dropMask(rows) = true;
        dropped(end+1) = rec.label; %#ok<AGROW>
        if opt.verbose
            fprintf(['Echo units: u%d (ch%d, n=%d) dropped: deflection %+.2f ms off centre ' ...
                     '(%.2f x centre), %.0f%% of spikes follow a leader (chance %.1f%%).\n'], ...
                rec.label, rec.channel, rec.n, rec.peakAtMs, rec.offRatio, ...
                100*rec.echoFrac, 100*rec.chance);
        end
    end
    lg(k) = rec;
end

report.log      = lg;
report.nDropped = numel(dropped);
report.dropped  = dropped;
if ~opt.apply || isempty(dropped)
    report.ok = true;
    if opt.verbose
        tail = '';
        if ~opt.apply, tail = ', not applied'; end
        fprintf('Echo units: %d of %d flagged%s.\n', numel(dropped), report.nTested, tail);
    end
    return;
end

% ---- write, atomically across the two files -----------------------------
bk = kiaSort_backup_results(outputPath, 'echo', {paths.lbl, ssPath});
if ~bk.ok
    if opt.verbose, fprintf('Echo units: backup failed, not writing.\n'); end
    return;
end
try
    newLbl = lbl;
    newLbl(dropMask) = -1;
    [unif, newLbl] = kiaSort_compact_unit_table(unif, newLbl);
    if exist(paths.lbl, 'file'), delete(paths.lbl); end
    h5create(paths.lbl, ['/' H.lbl], size(newLbl), 'Datatype', 'double');
    h5write(paths.lbl,  ['/' H.lbl], newLbl);
    ssData.crossChannelStats.unified_labels = unif;
    crossChannelStats = ssData.crossChannelStats; %#ok<NASGU>
    save(ssPath, 'crossChannelStats', '-append');
    report.changed = true;
    report.ok      = true;
catch ME
    kiaSort_restore_results(bk);
    if opt.verbose
        fprintf('Echo units: write failed (%s); rolled back from %s.\n', ME.message, bk.dir);
    end
    return;
end
if opt.verbose
    fprintf('Echo units: %d of %d units dropped.\n', numel(dropped), report.nTested);
end
end
