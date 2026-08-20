function snrReport = kiaSort_recompute_snr(outputPath, varargin)
%KIASORT_RECOMPUTE_SNR  Per-unit SNR from the sorted spikes and the detector threshold.
%
%   snrReport = kiaSort_recompute_snr(outputPath, ...)
%
%   The detectblity carried in unified_labels is measured at the sample stage
%   and keyed by labelInChannel, so two units that share a sample entry share
%   its value -- G60 and G61 on channel 26 both read 3.307 while their own
%   spikes give 3.39 and 5.82. This recomputes SNR after sorting, from each
%   unit's own amplitudes against its channel's detection threshold, and
%   writes it back as detectblity = SNR - 1 so every existing
%   "1 + detectblity" reader picks up the corrected value unchanged.
%
%   Options
%       'statistic' ('median')  'median' | 'mean' | 'p90' over |amplitude|
%       'write'     (true)      false computes and reports without saving
%       'verbose'   (false)
%
%   snrReport.snr       per-unit SNR, NaN where a unit has no spikes or its
%                       channel threshold is unavailable
%            .nUpdated  units whose stored value changed
%            .changed   true when sorted_samples.mat was rewritten
%            .ok        true when the pass completed

p = inputParser;
p.addParameter('statistic', 'median', @(x) ischar(x) || isstring(x));
p.addParameter('write',   true,  @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

snrReport = struct('snr', [], 'nUnits', 0, 'nUpdated', 0, ...
                   'changed', false, 'ok', false);

outputPath        = char(outputPath);
resSorted         = fullfile(outputPath, 'RES_Sorted');
samplesFolder     = fullfile(outputPath, 'RES_Samples');
sortedSampFolder  = fullfile(outputPath, 'Sorted_Samples');
unifiedLabelsH5   = fullfile(resSorted, 'unifiedLabels.h5');
amplitudeH5       = fullfile(resSorted, 'amplitude.h5');
sortedSamplesPath = fullfile(sortedSampFolder, 'sorted_samples.mat');

required = {unifiedLabelsH5, amplitudeH5, sortedSamplesPath};
for i = 1:numel(required)
    if ~exist(required{i}, 'file')
        if opt.verbose
            fprintf('Recompute SNR: %s missing, skipping.\n', required{i});
        end
        return;
    end
end

try
    lbl_all = double(h5read(unifiedLabelsH5, '/unifiedLabels'));
    amp_all = double(h5read(amplitudeH5,     '/amplitude'));
    ssData  = load(sortedSamplesPath, 'crossChannelStats');
catch ME
    if opt.verbose
        fprintf('Recompute SNR: read failed (%s).\n', ME.message);
    end
    return;
end
lbl_all = lbl_all(:);
amp_all = abs(amp_all(:));

if ~isfield(ssData, 'crossChannelStats') || ...
        ~isfield(ssData.crossChannelStats, 'unified_labels')
    return;
end
unif = ssData.crossChannelStats.unified_labels;
if ~isfield(unif, 'label') || ~isfield(unif, 'channelID'), return; end

thrCh = localChannelThresholds(samplesFolder, max(unif.channelID(:)));
if isempty(thrCh) || ~any(isfinite(thrCh))
    if opt.verbose
        fprintf('Recompute SNR: no channel thresholds found, skipping.\n');
    end
    return;
end

labels = unif.label(:);
nU     = numel(labels);
snr    = nan(nU, 1);
for u = 1:nU
    ch = unif.channelID(u);
    if ~isfinite(ch) || ch < 1 || ch > numel(thrCh), continue; end
    thr = thrCh(ch);
    if ~isfinite(thr) || thr <= 0, continue; end
    a = amp_all(lbl_all == labels(u));
    a = a(isfinite(a));
    if isempty(a), continue; end
    switch lower(char(opt.statistic))
        case 'mean', v = mean(a);
        case 'p90',  v = prctile(a, 90);
        otherwise,   v = median(a);
    end
    snr(u) = v / thr;
end

snrReport.snr    = snr;
snrReport.nUnits = nU;
snrReport.ok     = true;

old = nan(nU, 1);
if isfield(unif, 'detectblity') && numel(unif.detectblity) == nU
    old = double(unif.detectblity(:));
end
% Units the pass could not measure keep whatever they had.
newDet = old;
hit    = isfinite(snr);
newDet(hit) = snr(hit) - 1;
snrReport.nUpdated = sum(hit & (~isfinite(old) | abs(newDet - old) > 1e-9));

if opt.verbose
    fprintf('Recompute SNR: %d of %d units measured, %d changed.\n', ...
        sum(hit), nU, snrReport.nUpdated);
end
if ~opt.write || snrReport.nUpdated == 0, return; end

unif.detectblity = newDet;
unif.snr         = snr;

bk = kiaSort_backup_results(outputPath, 'snr', {sortedSamplesPath});
try
    ssData.crossChannelStats.unified_labels = unif;
    crossChannelStats = ssData.crossChannelStats;
    save(sortedSamplesPath, 'crossChannelStats', '-append');
catch ME
    kiaSort_restore_results(bk);
    snrReport.ok = false;
    if opt.verbose
        fprintf('Recompute SNR: write failed (%s); rolled back from %s.\n', ...
            ME.message, bk.dir);
    end
    return;
end
snrReport.changed = true;

end


function thr = localChannelThresholds(samplesFolder, nCh)
% mad_Thresh in channel_info.mat is the detector's per-channel threshold and
% covers every channel in one small file. The per-channel result files carry
% the same number split by polarity; they are only read when the shared file
% is unavailable, because each one also holds its waveforms.
thr = [];
if ~isfinite(nCh) || nCh < 1, return; end

infoPath = fullfile(samplesFolder, 'channel_info.mat');
if exist(infoPath, 'file')
    try
        info = load(infoPath, 'mad_Thresh');
        if isfield(info, 'mad_Thresh') && numel(info.mad_Thresh) >= nCh
            thr = double(info.mad_Thresh(:));
            return;
        end
    catch
    end
end

thr = nan(nCh, 1);
for c = 1:nCh
    f = fullfile(samplesFolder, sprintf('channel_%d_results.mat', c));
    if ~exist(f, 'file'), continue; end
    try
        S = load(f, 'out');
    catch
        continue;
    end
    if ~isfield(S, 'out'), continue; end
    v = [];
    if isfield(S.out, 'channel_thresholds_neg'), v = [v; double(S.out.channel_thresholds_neg(:))]; end %#ok<AGROW>
    if isfield(S.out, 'channel_thresholds_pos'), v = [v; double(S.out.channel_thresholds_pos(:))]; end %#ok<AGROW>
    v = v(isfinite(v) & v > 0);
    if ~isempty(v), thr(c) = max(v); end
end
end
