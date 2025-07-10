function kiaSort_main_sort_samples(outputPath, cfg, hp, varargin)

progressFcn = [];

for i = 1:2:length(varargin)
    name = lower(varargin{i});
    value = varargin{i+1};
    switch name
        case 'progressfcn'
            if isa(value, 'function_handle')
                progressFcn = value;
            end
        otherwise
            warning('Unknown parameter: %s', varargin{i});
    end
end

inputFolder = fullfile(outputPath, 'RES_Samples');

jitter_gap = cfg.spikeDistance * cfg.samplingFrequency / 1000;

outputFolder = fullfile(outputPath, 'Sorted_Samples');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% log file
logFile = fullfile(outputPath, 'KIASort_log.txt');
fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');
fprintf(fid, 'Sample sorting started at %s\n', datetime);

% Load channel info file
channelInfoPath = fullfile(inputFolder, 'channel_info.mat');
try
    channel_info = load(channelInfoPath);
    fprintf(fid, 'channel_info.mat loaded from %s.\n', channelInfoPath);
catch ME
    fullError = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(fid, 'ERROR: Failed to load %s:\n%s\n', channelInfoPath, fullError);
    fclose(fid);
    error(ME.message);
end

filePattern = fullfile(inputFolder, 'channel_*_results.mat');
files = dir(filePattern);

if isempty(files)
    fprintf(fid, 'ERROR: No files found matching the pattern %s\n', filePattern);
    fprintf(fid, 'Processing stopped due to errors.\n');
    fclose(fid);
    error('No files found matching the pattern %s', filePattern);
end

% Extract the channel number
channelNumbers = [];
for k = 1:length(files)
    tokens = regexp(files(k).name, 'channel_(\d+)_results\.mat', 'tokens');
    if ~isempty(tokens)
        channelNum = str2double(tokens{1}{1});
        if ~isnan(channelNum)
            channelNumbers(end+1) = channelNum;
        end
    end
end

numChannels = min(cfg.numChannels, length(channel_info.channel_mapping));

% Remove duplicates and sort the channel numbers
channelNumbers = unique(channelNumbers);
if ~all(ismember(1:numChannels, channelNumbers))
    fprintf(fid, 'ERROR: Not all channels are saved properly.\n');
    fprintf(fid, 'Processing stopped due to errors.\n');
    fclose(fid);
    error('Not all channels are saved properly.');
end

sortedSamples     = cell(numChannels,1);
sampleFeatures    = cell(numChannels,1);
crossChannelStats = [];
tic

logMessages = cell(numChannels,1);
hasError = false(numChannels,1);
channel_inclusion = channel_info.channel_inclusion;

if cfg.parallelProcessing
    setupParallel(cfg);
    parfor ch = 1:numChannels
        % Create temporary local log messages for each channel, since parfor cannot write to file directly
        msgs = {};
        data = [];

        inputFilename = fullfile(inputFolder, sprintf('channel_%d_results.mat', ch));

        if ~isfile(inputFilename)
            mainError = sprintf('File %s does not exist.', inputFilename);
            msgs{end+1} = sprintf('ERROR: %s\n', mainError);
            hasError(ch) = true;
        else
            try
                data = load(inputFilename, 'out').out;
            catch ME
                fullError = getReport(ME, 'extended', 'hyperlinks', 'off');
                msgs{end+1} = sprintf('ERROR: Failed to access %s:\n%s\n', inputFilename, fullError);
                hasError(ch) = true;
            end

            if ~isempty(data) && channel_inclusion(ch)
                try
                    [out, out_sampleFeatures] = kiaSort_cluster_classify(data, cfg, hp);
                    sortedSamples{ch,1} = out;
                    sampleFeatures{ch,1} = out_sampleFeatures;
                catch ME
                    fullError = getReport(ME, 'extended', 'hyperlinks', 'off');
                    msgs{end+1} = sprintf('ERROR: Error processing channel %d:\n%s\n', ch, fullError);
                    hasError(ch) = true;
                end
            else
                sortedSamples{ch,1} = [];
                sampleFeatures{ch,1} = [];
            end
        end
        logMessages{ch} = msgs;
    end
    % After parfor, if any error occurred, write log and throw error
    if any(hasError)
        for ch = 1:numChannels
            if ~isempty(logMessages{ch})
                for msgIdx = 1:length(logMessages{ch})
                    fprintf(fid, '%s', logMessages{ch}{msgIdx});
                end
            end
        end
        fclose(fid);
        error('Processing stopped due to errors. Check the log file: %s', logFile);
    end
else
    for ch = 1:numChannels
        msgs = {};
        data = [];

        inputFilename = fullfile(inputFolder, sprintf('channel_%d_results.mat', ch));

        if ~isfile(inputFilename)
            mainError = sprintf('File %s does not exist.', inputFilename);
            msgs{end+1} = sprintf('ERROR: %s\n', mainError);
            fprintf(fid, '%s', msgs{end});
            fclose(fid);
            error(mainError);
        else
            try
                data = load(inputFilename, 'out').out;
            catch ME
                fullError = getReport(ME, 'extended', 'hyperlinks', 'off');
                msgs{end+1} = sprintf('ERROR: Failed to access %s:\n%s\n', inputFilename, fullError);
                fprintf(fid, '%s', msgs{end});
                fclose(fid);
                error(ME.message);
            end

            if ~isempty(data) && channel_inclusion(ch)
                try
                    [out, out_sampleFeatures] = kiaSort_cluster_classify(data, cfg, hp);
                    sortedSamples{ch,1} = out;
                    sampleFeatures{ch,1} = out_sampleFeatures;
                catch ME
                    fullError = getReport(ME, 'extended', 'hyperlinks', 'off');
                    msgs{end+1} = sprintf('ERROR: Error processing channel %d:\n%s\n', ch, fullError);
                    fprintf(fid, '%s', msgs{end});
                    fclose(fid);
                    error(ME.message);
                end
            else
                sortedSamples{ch,1} = [];
                sampleFeatures{ch,1} = [];
            end
        end
        logMessages{ch} = msgs;
        if ~isempty(progressFcn)
            pct = ch / numChannels;
            msg = sprintf('Sorting samples channel %d of %d', ch, numChannels);
            progressFcn(pct, msg);
        end
    end
end

% log messages from loop (in non-parallel branch, errors would have already thrown)
for ch = 1:numChannels
    if ~isempty(logMessages{ch})
        for msgIdx = 1:length(logMessages{ch})
            fprintf(fid, '%s', logMessages{ch}{msgIdx});
        end
    end
end

% Unify detected clusters, transfered clusters, and labels accross channels

[sampleFeatures, sortedSamples] = unify_multiple_detections(cfg, sampleFeatures, sortedSamples);
[sampleFeatures, sortedSamples] = unify_cluster_relabling(cfg, sampleFeatures, sortedSamples);
[stats, unified_labels] = unify_spike_groups(sortedSamples, sampleFeatures, jitter_gap);
crossChannelStats.unified_labels = unified_labels;
crossChannelStats.stats = stats;
toc

outputFilename = fullfile(outputFolder, 'sorted_samples.mat');
featureFilename = fullfile(outputFolder, 'sample_features.mat');

% Save sorted samples
try
    save(outputFilename, 'sortedSamples', "crossChannelStats", '-v7.3');
catch ME
    fullError = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(fid, 'ERROR: Failed to save to %s:\n%s\n', outputFilename, fullError);
    fclose(fid);
    error(ME.message);
end

% Save sample features
try
    save(featureFilename, 'sampleFeatures', '-v7.3');
catch ME
    fullError = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(fid, 'ERROR: Failed to save to %s:\n%s\n', featureFilename, fullError);
    fclose(fid);
    error(ME.message);
end

fprintf(fid, '\nsamples were successfully sorted and results are saved at %s\n', datetime);
fclose(fid);
fprintf('All channels are processed and samples are sorted.\n');

end