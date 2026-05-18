function out = kiaSort_load_results(resultsPath, curated, individual_waveform)
% Load kiaSort outputs.
%   resultsPath          - the outputFolder used during sorting
%   curated              - true => prefer RES_Sorted/<field>_curated.h5;
%                          fall back to non-curated with a warning if missing
%   individual_waveform  - true => also load per-spike waveforms when saved

if nargin < 2 || isempty(curated),             curated = false; end
if nargin < 3 || isempty(individual_waveform), individual_waveform = false; end

resFolder = fullfile(resultsPath, 'RES_Sorted');
if ~exist(resFolder, 'dir')
    error('RES_Sorted folder not found at %s', resFolder);
end

suffix = '';
if curated
    if exist(fullfile(resFolder, 'spike_idx_curated.h5'),     'file') && ...
       exist(fullfile(resFolder, 'unifiedLabels_curated.h5'), 'file')
        suffix = '_curated';
    else
        warning('kiaSort_load_results:noCurated', ...
            'Curated outputs missing in %s. Loading non-curated.', resFolder);
    end
end

out = struct();
out.source       = ['RES_Sorted' suffix];
out.spike_idx    = readH5(resFolder, ['spike_idx'     suffix]);
out.channelNum   = readH5(resFolder, ['channelNum'    suffix]);
out.unifiedLabels = readH5(resFolder, ['unifiedLabels' suffix]);

labelsFile = fullfile(resFolder, 'labels.h5');
if exist(labelsFile, 'file')
    out.labelOnChannel = readH5(resFolder, 'labels');
end

if individual_waveform
    wfFile = fullfile(resFolder, ['waveforms' suffix '.h5']);
    if exist(wfFile, 'file')
        out.waveforms = h5read(wfFile, ['/waveforms' suffix]);
    else
        warning('kiaSort_load_results:noWaveforms', ...
            'Waveform file %s not found.', wfFile);
    end
end

out.units = loadUnits(resultsPath, resFolder, suffix);
end


function v = readH5(folder, name)
v = h5read(fullfile(folder, [name '.h5']), ['/' name]);
v = v(:);
end


function units = loadUnits(resultsPath, resFolder, suffix)
units = struct('id', {}, 'channel', {}, 'meanWaveform', {}, 'isolation', {});

curatedMat = fullfile(resFolder, 'curated_sample.mat');
sortedSamp = fullfile(resultsPath, 'Sorted_Samples', 'sorted_samples.mat');

if strcmp(suffix, '_curated') && exist(curatedMat, 'file')
    s = load(curatedMat);
    if isfield(s, 'curatedSamples')
        c = s.curatedSamples;
        ids = c.unifiedLabels(:);
        for k = 1:numel(ids)
            units(k).id      = ids(k);
            if isfield(c, 'channelNum')   && numel(c.channelNum)   >= k
                units(k).channel = c.channelNum(k);
            end
            if isfield(c, 'waveform') && size(c.waveform, 1) >= k
                units(k).meanWaveform = squeeze(c.waveform(k, :, :));
            end
            if isfield(c, 'unitIsolation') && numel(c.unitIsolation) >= k
                units(k).isolation = c.unitIsolation{k};
            end
        end
        return;
    end
end

if exist(sortedSamp, 'file')
    s = load(sortedSamp);
    if isfield(s, 'crossChannelStats') && isfield(s.crossChannelStats, 'unified_labels')
        u = s.crossChannelStats.unified_labels;
        ids = u.label(:);
        for k = 1:numel(ids)
            units(k).id = ids(k);
            if isfield(u, 'channelID')     && numel(u.channelID)     >= k
                units(k).channel = u.channelID(k);
            end
            if isfield(u, 'meanWaveforms') && size(u.meanWaveforms, 1) >= k
                units(k).meanWaveform = squeeze(u.meanWaveforms(k, :, :));
            end
        end
    end
end
end
