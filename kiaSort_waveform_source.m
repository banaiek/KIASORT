function src = kiaSort_waveform_source(outputPath, cfg, verbose)
%KIASORT_WAVEFORM_SOURCE  Resolve where per-spike waveforms can be read from.
%
%   src = kiaSort_waveform_source(outputPath, cfg, verbose)
%
%   Prefers RES_Sorted/waveforms*.h5, whose rows are aligned with the other
%   H5 outputs, so a unit's spikes can be pulled straight out. Falls back to
%   the raw file named by cfg.fullFilePath, which costs one filtered read per
%   spike and so should be capped by the caller.
%
%   src.ok    false when neither source is reachable
%   src.mode  'h5' | 'raw'
%   src.half  half-window in samples (from cfg.spikeDuration)
%
%   Pair with kiaSort_read_waveforms.

src = struct('ok', false, 'mode', '', 'half', 0, ...
             'file', '', 'dset', '', 'sz', [], 'is3d', false, 'nT', 0, ...
             'map', [], 'chanMap', [], 'nSamp', 0, 'cfg', cfg);
if nargin < 3, verbose = false; end
if ~isstruct(cfg) || ~isfield(cfg, 'samplingFrequency'), return; end

src.half = round((cfg.spikeDuration/2) * cfg.samplingFrequency / 1000);

resSorted = fullfile(char(outputPath), 'RES_Sorted');
d = dir(fullfile(resSorted, 'waveforms*.h5'));
for i = 1:numel(d)
    f = fullfile(resSorted, d(i).name);
    try
        info = h5info(f);
        if isempty(info.Datasets), continue; end
        sz = info.Datasets(1).Dataspace.Size;
        src.file = f;
        src.dset = ['/' info.Datasets(1).Name];
        src.sz   = sz;
        src.is3d = numel(sz) >= 3;
        src.nT   = sz(end);
        src.mode = 'h5';
        src.ok   = true;
        return;
    catch
    end
end

if ~isfield(cfg, 'fullFilePath') || ~isfield(cfg, 'numChannels') || ~isfield(cfg, 'dataType')
    return;
end
fp = char(cfg.fullFilePath);
if isempty(fp) || ~exist(fp, 'file'), return; end
try
    cfg.outputFolder = char(outputPath);
    src.map   = map_input_file(fp, cfg);
    src.nSamp = size(src.map.Data.data, 2);
    src.chanMap = (1:cfg.numChannels)';
    ciPath = fullfile(char(outputPath), 'RES_Samples', 'channel_info.mat');
    if exist(ciPath, 'file')
        ci = load(ciPath, 'channel_mapping');
        if isfield(ci, 'channel_mapping') && ~isempty(ci.channel_mapping)
            src.chanMap = double(ci.channel_mapping(:));
        end
    end
    src.cfg  = cfg;
    src.mode = 'raw';
    src.ok   = src.nSamp > 0;
catch ME
    if verbose, fprintf('Waveform source: raw map failed (%s).\n', ME.message); end
end
end
