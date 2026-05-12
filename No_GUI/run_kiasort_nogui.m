function run_kiasort_nogui(dataFilePath, outputFolder, channelMapFile, cfg_overrides)
% Run kiaSort end-to-end without the GUI.
%
%   dataFilePath   - .dat / .bin file
%   outputFolder   - directory for results
%   channelMapFile - optional .mat with channel map fields (empty to skip)
%   cfg_overrides  - optional struct of cfg fields to override after defaults

if nargin < 2
    error('dataFilePath and outputFolder are required');
end
if nargin < 3, channelMapFile = []; end
if nargin < 4 || isempty(cfg_overrides), cfg_overrides = struct(); end

if ~exist(dataFilePath, 'file')
    error('Input data file not found: %s', dataFilePath);
end
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

cfg = kiaSort_main_configs();
cfg = kiaSort_extended_configs(cfg);
cfg = kiaSort_hidden_configs(cfg);

ovrFields = fieldnames(cfg_overrides);
for i = 1:numel(ovrFields)
    cfg.(ovrFields{i}) = cfg_overrides.(ovrFields{i});
end

[cfg.inputFolder, ~, ~] = fileparts(dataFilePath);
cfg.fullFilePath = dataFilePath;
cfg.outputFolder = outputFolder;

hp = sorting_hyperparameters_in();

channel_mapping = [];
channel_inclusion = [];
channel_locations = [];
if ~isempty(channelMapFile) && exist(channelMapFile, 'file')
    cfg.channel_info = channelMapFile;
    [channel_mapping, channel_locations, channel_inclusion] = ...
        load_channel_map(channelMapFile, cfg);
elseif isfield(cfg, 'numChannels')
    channel_mapping   = 1:cfg.numChannels;
    channel_inclusion = true(cfg.numChannels, 1);
    channel_locations = [];
end

cfg.num_channel_extract = derive_num_channel_extract(channel_locations, ...
    cfg.waveform_radius, cfg.num_channel_extract);

fprintf('Input: %s\nOutput: %s\nChannels: %d, fs: %d Hz, dtype: %s, BP: [%d %d]\n', ...
    dataFilePath, outputFolder, cfg.numChannels, cfg.samplingFrequency, ...
    cfg.dataType, cfg.bandpass(1), cfg.bandpass(2));

kiaSort_main_extract_sample_data(cfg.fullFilePath, cfg.outputFolder, cfg, ...
    'channel_mapping',   channel_mapping, ...
    'channel_inclusion', channel_inclusion, ...
    'channel_locations', channel_locations);

if ~cfg.sort_only
    kiaSort_main_sort_samples(cfg.outputFolder, cfg, hp);
end

kiaSort_main_sortData(cfg.fullFilePath, cfg.outputFolder, cfg);

fprintf('Done. Results in %s\n', outputFolder);

end
