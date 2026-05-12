%% kiaSort_main_noGUI_run
% Edit paths, channel map, and config overrides below, then run.
% Defaults come from kiaSort_main_configs / kiaSort_extended_configs /
% kiaSort_hidden_configs; anything set in cfg_overrides replaces them.

dataFilePath   = 'C:/data/recording.dat';
outputFolder   = 'C:/data/results';
channelMapFile = '';                              % '' to skip

cfg_overrides = struct();
cfg_overrides.dataType            = 'int16';
cfg_overrides.numChannels         = 384;
cfg_overrides.samplingFrequency   = 3e4;
cfg_overrides.bandpass            = [500 6000];
cfg_overrides.commonRef           = 'median';       % 'median' | 'mean' | 'none'
cfg_overrides.denoising           = true;
cfg_overrides.extremeNoise        = false;
cfg_overrides.min_threshold       = 5.5;
cfg_overrides.waveform_radius     = 150;
cfg_overrides.num_channel_extract = 5;
cfg_overrides.numSampleChunks      = 300;
cfg_overrides.sortingChunkDuration = 120;
cfg_overrides.sort_only            = false;
cfg_overrides.useGPU               = true;
cfg_overrides.parallelProcessing   = false;
cfg_overrides.extractWaveform      = false;

run_kiasort_nogui(dataFilePath, outputFolder, channelMapFile, cfg_overrides);
