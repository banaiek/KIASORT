function out = load_sorted_samples_gui(outputPath)

logFile = fullfile(outputPath, 'KIASort_GUI_log.txt');

fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');

sampleFolder = fullfile(outputPath, 'RES_Samples');
sortedSamplessFolder = fullfile(outputPath, 'Sorted_Samples');
sortedSamplessPath = fullfile(sortedSamplessFolder, 'sorted_samples.mat');
sampleFeaturesPath = fullfile(sortedSamplessFolder, 'sample_features.mat');
channelInfoPath = fullfile(sampleFolder, 'channel_info.mat');

    try
    out.channel_info = load(channelInfoPath);
    fprintf(fid, 'Sorted Samples are loaded. %s\n', datetime);
    catch ME
    fprintf(fid, 'ERROR: Failed to load channel_info.mat from %s: %s.\n', channelInfoPath, ME.message);
    sprintf('ERROR: Failed to load channel channel_info.mat');
    end


    try
    load(sortedSamplessPath,'sortedSamples', 'crossChannelStats');
    load(sampleFeaturesPath,'sampleFeatures');
    out.sortedSamples     = sortedSamples;
    out.crossChannelStats = crossChannelStats;
    out.sampleFeatures     = sampleFeatures;
    catch ME
    fprintf(fid, 'ERROR: Failed to load sorted_samples.mat from %s: %s.\n', sortedSamplessPath, ME.message);
    sprintf('ERROR: Failed to load channel sorted_samples.mat');
    end

end