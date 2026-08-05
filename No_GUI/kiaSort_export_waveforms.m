function outFile = kiaSort_export_waveforms(dataFilePath, outputPath, waveformDurMs, varargin)
%KIASORT_EXPORT_WAVEFORMS  Re-extract per-spike waveforms at any duration.
%
%   outFile = kiaSort_export_waveforms(dataFilePath, outputPath, waveformDurMs, ...)
%
%   Reads the sorted spike times from <outputPath>/RES_Sorted and pulls a
%   waveformDurMs window around each spike from the raw file, on that
%   spike's own channel. Writes <outputPath>/RES_Sorted/waveforms_export.h5
%   in the same [nSpikes x nSamples] layout as waveforms_curated.h5, so row
%   i lines up with row i of spike_idx.h5 / unifiedLabels.h5 / channelNum.h5.
%
%   Samples are raw (unfiltered), matching what the pipeline writes.
%   Windows running past either end of the recording are left as zeros.
%
%   samplingFrequency, numChannels and dataType are taken from the cfg
%   stored in Sorted_Samples/sorted_samples.mat; the name/value overrides
%   below are for recordings whose cfg is missing or wrong.
%
%   Name/Value:
%       'numChannels'       raw-file channel count
%       'samplingFrequency' Hz
%       'dataType'          e.g. 'int16'
%       'channelMapping'    raw row per channel number
%       'outputName'        default 'waveforms_export'
%       'batchSize'         spikes per write (default 20000)
%       'verbose'           default true

p = inputParser;
p.addRequired('dataFilePath',  @(x) ischar(x) || isstring(x));
p.addRequired('outputPath',    @(x) ischar(x) || isstring(x));
p.addRequired('waveformDurMs', @(x) isscalar(x) && isnumeric(x) && x > 0);
p.addParameter('numChannels',       [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
p.addParameter('samplingFrequency', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
p.addParameter('dataType',          '', @(x) ischar(x) || isstring(x));
p.addParameter('channelMapping',    [], @(x) isempty(x) || isnumeric(x));
p.addParameter('outputName', 'waveforms_export', @(x) ischar(x) || isstring(x));
p.addParameter('batchSize',  20000, @(x) isscalar(x) && isnumeric(x) && x > 0);
p.addParameter('verbose',    true,  @(x) islogical(x) || isnumeric(x));
p.parse(dataFilePath, outputPath, waveformDurMs, varargin{:});
opt = p.Results;
opt.verbose = logical(opt.verbose);

dataFilePath = char(dataFilePath);
outputPath   = char(outputPath);
outName      = char(opt.outputName);

if ~exist(dataFilePath, 'file')
    error('kiaSort_export_waveforms:noData', 'Data file not found: %s', dataFilePath);
end
resFolder = fullfile(outputPath, 'RES_Sorted');
if ~exist(resFolder, 'dir')
    error('kiaSort_export_waveforms:noResults', 'RES_Sorted not found under %s', outputPath);
end

%% cfg
cfg = struct();
ssPath = fullfile(outputPath, 'Sorted_Samples', 'sorted_samples.mat');
if exist(ssPath, 'file')
    try
        ss = load(ssPath, 'sortedSamples');
        for i = 1:numel(ss.sortedSamples)
            if ~isempty(ss.sortedSamples{i}) && isfield(ss.sortedSamples{i}, 'cfg')
                cfg = ss.sortedSamples{i}.cfg;
                break;
            end
        end
    catch
    end
end

if ~isempty(opt.samplingFrequency), cfg.samplingFrequency = opt.samplingFrequency; end
if ~isempty(opt.numChannels),       cfg.numChannels       = opt.numChannels;       end
if ~isempty(opt.dataType),          cfg.dataType          = char(opt.dataType);    end

missing = {};
if ~isfield(cfg,'samplingFrequency') || isempty(cfg.samplingFrequency), missing{end+1} = 'samplingFrequency'; end
if ~isfield(cfg,'numChannels')       || isempty(cfg.numChannels),       missing{end+1} = 'numChannels';       end
if ~isfield(cfg,'dataType')          || isempty(cfg.dataType),          missing{end+1} = 'dataType';          end
if ~isempty(missing)
    error('kiaSort_export_waveforms:missingCfg', ...
        ['Could not resolve %s from the stored cfg. Pass it as a name/value ' ...
         'argument.'], strjoin(missing, ', '));
end
cfg.outputFolder = outputPath;

fs = cfg.samplingFrequency;

%% spike table
spkFile = fullfile(resFolder, 'spike_idx.h5');
chnFile = fullfile(resFolder, 'channelNum.h5');
if ~exist(spkFile,'file') || ~exist(chnFile,'file')
    error('kiaSort_export_waveforms:noSpikes', ...
        'spike_idx.h5 / channelNum.h5 not found in %s', resFolder);
end
spk = double(h5read(spkFile, '/spike_idx'));
chn = double(h5read(chnFile, '/channelNum'));
spk = spk(:); chn = chn(:);
if numel(spk) ~= numel(chn)
    error('kiaSort_export_waveforms:h5Mismatch', 'spike_idx and channelNum lengths disagree.');
end
nSpikes = numel(spk);
if nSpikes == 0
    error('kiaSort_export_waveforms:noSpikes', 'No spikes found in %s', resFolder);
end

%% channel mapping
chanMap = opt.channelMapping;
if isempty(chanMap)
    ciPath = fullfile(outputPath, 'RES_Samples', 'channel_info.mat');
    if exist(ciPath, 'file')
        try
            ci = load(ciPath, 'channel_mapping');
            if isfield(ci, 'channel_mapping'), chanMap = ci.channel_mapping; end
        catch
        end
    end
end
if isempty(chanMap)
    chanMap = (1:cfg.numChannels)';
end
chanMap = double(chanMap(:));

%% window
half = round((waveformDurMs/2) * fs / 1000);
W    = 2*half + 1;

m = map_input_file(dataFilePath, cfg);
nSamp = size(m.Data.data, 2);

outFile = fullfile(resFolder, [outName '.h5']);
if exist(outFile, 'file')
    delete(outFile);
end
h5create(outFile, ['/' outName], [nSpikes W], ...
    'Datatype', 'double', 'ChunkSize', [min(nSpikes, 1024) W]);

batch   = round(opt.batchSize);
nSkipped = 0;

for s = 1:batch:nSpikes
    e   = min(s + batch - 1, nSpikes);
    idx = s:e;
    buf = zeros(numel(idx), W);

    for k = 1:numel(idx)
        i = idx(k);
        c = chn(i);
        if ~isfinite(c) || c < 1 || c > numel(chanMap)
            nSkipped = nSkipped + 1;
            continue;
        end
        row = chanMap(c);
        if ~isfinite(row) || row < 1 || row > cfg.numChannels
            nSkipped = nSkipped + 1;
            continue;
        end
        a = spk(i) - half;
        b = spk(i) + half;
        if a < 1 || b > nSamp
            nSkipped = nSkipped + 1;
            continue;
        end
        buf(k, :) = double(m.Data.data(row, a:b));
    end

    h5write(outFile, ['/' outName], buf, [s 1], [numel(idx) W]);
end

if opt.verbose
    fprintf('Exported %d waveforms (%.2f ms, %d samples) to %s\n', ...
        nSpikes, waveformDurMs, W, outFile);
    if nSkipped > 0
        fprintf('  %d left zero (edge of recording or unmapped channel)\n', nSkipped);
    end
end

end
