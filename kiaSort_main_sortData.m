function kiaSort_main_sortData(inputPath, outputPath, cfg, varargin)

% main_sort_data gets inputPath and cfg and sort the data using sorted
% samples
% output of this function is HDF5 files, for sorted spikes. A list of
% sorted files are:
%                   spike_idx   :   indicies of spikes across all channels
%                   labels      :   class labels of spikes on their
%                   detected channel
%                   updatedLabels   : labels of spikes on their
%                   amplitude      : amplitude of the main channel spike
%                   unifiedLabels   :   lables of spikes across all channels
%                   features        :   low dimensional features of spikes
%                   (3comp)
%                   waveforms    :   waveforms for all detected spikes
%                   channelNum   :   channel where the spike was detected



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

% log file
logFile = fullfile(outputPath, 'KIASort_log.txt');
fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');
fprintf(fid, 'Main sorting started at %s\n', datetime);

warning('on','all');
lastwarn('');


try

    samplesFolder = fullfile(outputPath, 'RES_Samples');
    sortedSamplesFolder = fullfile(outputPath, 'Sorted_Samples');
    sortedSamplesPath = fullfile(sortedSamplesFolder, 'sorted_samples.mat');
    channelInfoPath = fullfile(samplesFolder, 'channel_info.mat');

    if exist(channelInfoPath, 'file')
        channel_info = load(channelInfoPath);
    else
        error('channel_info.mat not found in the specified inputPath.');
    end

    if exist(sortedSamplesPath, 'file')
        load(sortedSamplesPath,'sortedSamples', 'crossChannelStats');
    else
        error('sorted_samples.mat not found in Sorted_Samples folder.');
    end



    outputFolder = fullfile(outputPath, 'RES_Sorted');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    elseif isfield(cfg, 'sort_only')
        if cfg.sort_only
            baseName   = 'RES_Sorted';
            listing = dir(fullfile(outputPath, [baseName '_V*']));
            versionNums = [];

            for i = 1:numel(listing)
                name = listing(i).name;
                tok = regexp(name, [baseName '_V(\d+)$'], 'tokens');
                if ~isempty(tok)
                    versionNums(end+1) = str2double(tok{1}{1});  
                end
            end

            if isempty(versionNums)
                newVer = 1;
            else
                newVer = max(versionNums) + 1;
            end

            outputFolder = fullfile(outputPath, sprintf('%s_V%d', baseName, newVer));
            mkdir(outputFolder);
        end
    end




    channel_mapping = channel_info.channel_mapping;
    channel_inclusion = channel_info.channel_inclusion;

    num_channels = min(cfg.numChannels, length(channel_mapping));
    fs = cfg.samplingFrequency;

    % memory map
    m = map_input_file(inputPath, cfg);
    num_samples = size(m.Data.data,2);

    window_size = cfg.num_channel_extract * 2 + 1; % Number of channels in the window
    half_window = floor(window_size / 2);
    search_window = min(round(1.5 * half_window), round(num_channels/2)-3);
    release_window = 2 * search_window;

    % Parameters
    chunk_duration      = cfg.sortingChunkDuration;
    num_chunk_pts       = chunk_duration * fs;
    num_chunks          = ceil(num_samples / num_chunk_pts);
    marginOffset_pts    = cfg.borderMargin * fs / 1000;
    spikeDistance       = cfg.spikeDistance * fs / 1000;
    batch_ch_size       = min(cfg.batch_ch_size,num_channels);

    % Generate indices for batches
    channel_idx = [1:batch_ch_size:(num_channels-1),num_channels];

    % Maps to store intermediate results for this chunk
    mapNames = {
        'bp_channels', 'spk_idx_channels', ...
        'spk_Val_channels', 'spk_ID_channels', ...
        'waveform_channels', 'relabling_channels', 'keep_id_channels', ...
        'output_waveform', 'output_spk_idx', ...
        'output_labels', 'output_updatedLabels', 'output_features',...
        'output_amplitude','unevaluated_spk_idx'
        };

    for i = 1:length(mapNames)
        eval([mapNames{i} ' = containers.Map(''KeyType'', ''double'', ''ValueType'', ''any'');']);
    end

    channel_thresholds      = cell(num_channels, 1);
    sorted_out              = [];

    for chunk_i = 1:num_chunks

        chunk_channel_inclusion = true(num_channels,1);

        % Log start time for the chunk
        chunkStartTime = datetime;
        fprintf(fid, '\nChunk %d/%d started at %s\n', chunk_i, num_chunks, chunkStartTime);

        start_sample = (chunk_i - 1) * num_chunk_pts + 1;
        end_sample = min(chunk_i * num_chunk_pts + 2 * marginOffset_pts, num_samples);
        current_chunk_length = end_sample - start_sample + 1;
        if current_chunk_length <= 2 * marginOffset_pts
            continue;
        end
        
        chunk_limits = [start_sample,end_sample];
        

        discarded_spk_idx       = cell(num_channels, 1);
        unmatched_spk_idx       = cell(num_channels, 1);

        % Process each channel
        for ch_idx = 1:num_channels
            tic
            ch_indices_mapped  = [max(1, ch_idx - half_window) : min(num_channels, ch_idx + half_window)]';
            ch_indices_cluster = [max(1, ch_idx - release_window) : min(num_channels, ch_idx + release_window)]';


            last_batch_channel = min(ch_indices_mapped(end),num_channels);

            if ch_idx == 1
                batch_idx = 1;
                start_ch = channel_idx(batch_idx);
                end_ch = batch_idx * batch_ch_size;

                try
                    selected_data = batch_extract(m, chunk_limits, start_ch, end_ch, channel_mapping, channel_inclusion, cfg);
                catch ME
                    fprintf(fid, 'ERROR: Failed to extract batch %d: %s\n', batch_idx, ME.message);
                    if ~isempty(ME.stack)
                        fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
                    end
                    errorFlag = true;
                    error('Failed to extract data batch %d. Check the log file for details.', batch_idx);
                end

                if isfield(cfg, 'commonRef')
                    switch cfg.commonRef
                        case 'median'
                            selected_data(channel_inclusion(start_ch:end_ch),:) = selected_data(channel_inclusion(start_ch:end_ch),:) - median(selected_data(channel_inclusion(start_ch:end_ch),:),'omitmissing');
                            fprintf(fid, 'Samples median referenced.\n');
                        case 'mean'
                            selected_data(channel_inclusion(start_ch:end_ch),:) = selected_data(channel_inclusion(start_ch:end_ch),:) - mean(selected_data(channel_inclusion(start_ch:end_ch),:),'omitmissing');
                            fprintf(fid, 'Samples median referenced.\n');
                        case 'none'
                            % skip
                    end
                end

            elseif any(channel_idx == last_batch_channel) && last_batch_channel~=num_channels
                batch_idx = find(channel_idx == last_batch_channel);
                start_ch = channel_idx(batch_idx);
                end_ch = batch_idx * batch_ch_size;

                try
                    selected_data = batch_extract(m, chunk_limits, start_ch, end_ch, channel_mapping, channel_inclusion, cfg);
                catch ME
                    fprintf(fid, 'ERROR: Failed to extract batch %d: %s\n', batch_idx, ME.message);
                    if ~isempty(ME.stack)
                        fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
                    end
                    errorFlag = true;
                    error('Failed to extract data batch %d. Check the log file for details.', batch_idx);
                end

                if isfield(cfg, 'commonRef')
                    switch cfg.commonRef
                        case 'median'
                            selected_data(channel_inclusion(start_ch:end_ch),:) = selected_data(channel_inclusion(start_ch:end_ch),:) - median(selected_data(channel_inclusion(start_ch:end_ch),:),'omitmissing');
                            fprintf(fid, 'Samples median referenced.\n');
                        case 'mean'
                            selected_data(channel_inclusion(start_ch:end_ch),:) = selected_data(channel_inclusion(start_ch:end_ch),:) - mean(selected_data(channel_inclusion(start_ch:end_ch),:),'omitmissing');
                            fprintf(fid, 'Samples median referenced.\n');
                        case 'none'
                            % skip
                    end
                end

            end



            % filter and detect spikes for channels in the window
            for idx = 1:length(ch_indices_mapped)
                mapped_ch_idx = ch_indices_mapped(idx);

                if ~isKey(spk_idx_channels, mapped_ch_idx)

                    if ~channel_inclusion(mapped_ch_idx)
                        [waveOut, out_filtered, out_detected]   = setWaveformNull;

                        for j = 1:length(mapNames)
                            currentMap = eval(mapNames{j});
                            currentMap(mapped_ch_idx) = [];
                            eval([mapNames{j} ' = currentMap;']);
                        end

                        continue;
                    end

                    
                    batch_mapped_idx = mod(mapped_ch_idx,batch_ch_size);
                    if batch_mapped_idx == 0
                        batch_mapped_idx = batch_ch_size;
                    end
                    channel_data = double(selected_data(batch_mapped_idx, :));

                    out_filtered    = kiaSort_filter_signal(channel_data, cfg);
                    out_filtered.rms_bands = channel_info.rms_bands_channels(mapped_ch_idx,:);
                    out_filtered.scale_factor = channel_info.scale_factor(mapped_ch_idx);
                    out_filtered.adj_distant = channel_info.adj_distant(mapped_ch_idx);
                    out_filtered.Ns_seq     = channel_info.Ns_seq;
                    out_filtered.band_pairs = channel_info.band_pairs;
                    out_filtered.mad_Thresh = channel_info.mad_Thresh(mapped_ch_idx);
                    out_detected = kiaSort_detect_spike(out_filtered, cfg, 0);

                    spk_idx_channels(mapped_ch_idx)         = out_detected.spk_idx;
                    spk_Val_channels(mapped_ch_idx)         = out_detected.spk_Val;
                    spk_ID_channels(mapped_ch_idx)          = out_detected.spk_ID;
                    bp_channels(mapped_ch_idx)              = out_filtered.bandpass_signal;
                    channel_thresholds{mapped_ch_idx, 1}    = out_detected.rms_values;
                end
            end

            % Prepare data for waveform extraction
            inputWF.bandpass_data_window         = cell(length(ch_indices_mapped), 1);
            inputWF.spk_idx_all_channels         = cell(length(ch_indices_mapped), 1);
            inputWF.inclusion                    = zeros(length(ch_indices_mapped),1);

            for idx = 1:length(ch_indices_mapped)
                inputWF.main_channel_idx         = ch_idx;
                mapped_ch_idx                    = ch_indices_mapped(idx);
                inputWF.all_channel_idx          = ch_indices_mapped;
                inputWF.Ns_seq                   = channel_info.Ns_seq;
                inputWF.bandpass_data_window{idx, 1}   = bp_channels(mapped_ch_idx);
                inputWF.spk_idx_all_channels{idx, 1}   = spk_idx_channels(mapped_ch_idx);
                inputWF.spk_ID_all_channels{idx, 1}    = spk_ID_channels(mapped_ch_idx);
                inputWF.inclusion(idx, 1)              = channel_inclusion(mapped_ch_idx);
                inputWF.channel_thresholds{idx, 1}     = channel_thresholds{mapped_ch_idx, 1};
            end

            if channel_inclusion(ch_idx)

                % Extract waveforms if channel is included
                if inputWF.inclusion(ch_indices_mapped == ch_idx)
                    cancel_overlap = 1;
                    sample = 0;
                    [waveOut] = kiaSort_waveform_extraction(cfg, inputWF, cancel_overlap, sample);
                end

                waveform_channels(ch_idx) = waveOut;
                [out_bestChannel]     = kiaSort_best_channel_detection(waveOut, 100, cfg);
                waveOut.waveform      = waveOut.waveform(out_bestChannel.keep, :, :);
                waveOut.waveformInfo  = sortedSamples{ch_idx,1}.waveformInfo;
                tmp_idx               = spk_idx_channels(ch_idx);
                sorting_spk_idx       = tmp_idx(out_bestChannel.keep);
                amplitude             = spk_Val_channels(ch_idx) .* spk_ID_channels(ch_idx);
                amplitude             = amplitude(out_bestChannel.keep);
                chunk_channel_inclusion(ch_idx) = ~isempty(sorting_spk_idx);
                altChannels = out_bestChannel.max_altChannel(out_bestChannel.keep);
                keep_id_channels(ch_idx) = out_bestChannel;
                out_bestChannel       = [];
                unevaluated_spk_idx(ch_idx) = sorting_spk_idx;

                if chunk_channel_inclusion(ch_idx)

                    % Waveform preparation for classification
                    data_sorting = kiaSort_preprocess_waveforms(waveOut, cfg);
                    data_sorting.classifierInfo = sortedSamples{ch_idx, 1}.classifierInfo;
                    data_sorting.PCA = sortedSamples{ch_idx, 1}.clusteringInfo.PCA;
                    data_sorting.amplitude     = amplitude;

                    [predLabels, features, validKeep] = kiaSort_predict_spike_labels(data_sorting,cfg);
                    if iscategorical(predLabels)
                        predLabels = str2double(cellstr(predLabels));
                    end
                    clusterSelection = sortedSamples{ch_idx, 1}.clusteringInfo.clusterSelection;
                    clusterRelabeling = sortedSamples{ch_idx, 1}.clusteringInfo.clusterRelabeling;
                    [updatedLabels, realigned_spk_idx, realigned_waveform, ~] = ...
                        realignSpikes(predLabels, waveOut.waveform, sorting_spk_idx, clusterRelabeling, cfg);

                    relabling_out = kiaSort_evaluate_labels(updatedLabels, clusterSelection, altChannels, validKeep);
                    relabling_channels (ch_idx) = relabling_out;

                    % keep waveforms that pass label evaluation
                    output_labels(ch_idx)              = predLabels(relabling_out.keep);
                    output_updatedLabels(ch_idx)       = updatedLabels(relabling_out.keep);
                    output_features(ch_idx)            = features(relabling_out.keep, :);
                    output_spk_idx(ch_idx)             = realigned_spk_idx(relabling_out.keep);
                    output_waveform(ch_idx)            = realigned_waveform(relabling_out.keep,:,:);
                    output_amplitude(ch_idx)          = amplitude(relabling_out.keep);
                else

                    relabling_channels(ch_idx)         = [];
                    output_labels(ch_idx)              = [];
                    output_updatedLabels(ch_idx)       = [];
                    output_features(ch_idx)            = [];
                    output_spk_idx(ch_idx)             = [];
                    output_waveform(ch_idx)            = [];
                    output_amplitude(ch_idx)          = [];
                end

            end

            if ch_idx > search_window

                if ch_idx < num_channels
                    iCh = (ch_idx - search_window);
                    if channel_inclusion(iCh) && chunk_channel_inclusion(iCh)

                        relabling_in = relabling_channels (iCh);
                        spike_idx_in = unevaluated_spk_idx(iCh);
                        num_relabling  = relabling_in.nRelabling;

                        for iRlbl = 1 : num_relabling

                            relabeled_channel = relabling_in.lookUpChannel(iRlbl);
                            not_kept_idx = relabling_in.not_kept_idx{iRlbl, 1};
                            relabled_idx = spike_idx_in(not_kept_idx);

                            if abs(relabeled_channel - iCh) <= search_window && relabeled_channel ~= iCh

                                temp_spike_idx = spk_idx_channels(relabeled_channel);
                                temp_spike_val = spk_Val_channels(relabeled_channel) .* spk_ID_channels(relabeled_channel);
                                if isempty(temp_spike_idx)
                                    discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; relabled_idx];
                                    continue,
                                end
                                temp_waveform = waveform_channels(relabeled_channel).waveform;
                                keep_alt_idx = keep_id_channels(relabeled_channel).keep;
                                idx_keep = find(~keep_alt_idx);
                                [d, alt_idx_subset] = nearest_index(relabled_idx, temp_spike_idx(idx_keep));
                                valid_replacement = (d <= 1.5 * spikeDistance);
                                replaced_id = idx_keep(alt_idx_subset(valid_replacement));
                                replaced_idx = temp_spike_idx(replaced_id);
                                amplitude    = temp_spike_val(replaced_id);
                                non_replaced_idx = relabled_idx(d > 1.5 * spikeDistance);
                                discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; non_replaced_idx];

                                if any(valid_replacement)
                                    clusterSelection = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterSelection;
                                    clusterRelabeling  = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterRelabeling;
                                    data_sorting.waveform      = temp_waveform(replaced_id,:,:);
                                    data_sorting.waveformInfo  = sortedSamples{relabeled_channel,1}.waveformInfo;
                                    data_sorting     = kiaSort_preprocess_waveforms(data_sorting, cfg);
                                    data_sorting.PCA = sortedSamples{relabeled_channel, 1}.clusteringInfo.PCA;
                                    data_sorting.classifierInfo = sortedSamples{relabeled_channel, 1}.classifierInfo;
                                    data_sorting.amplitude     = amplitude;
                                    [predLabels, features, validKeep]      = kiaSort_predict_spike_labels(data_sorting, cfg);
                                    if iscategorical(predLabels)
                                        predLabels = str2double(cellstr(predLabels));
                                    end
                                    [updatedLabels, realigned_spk_idx, realigned_waveform, ~] = ...
                                        realignSpikes(predLabels, temp_waveform(replaced_id,:,:), replaced_idx, clusterRelabeling, cfg);

                                    relabling_out    = evaluate_labels(updatedLabels, clusterSelection, validKeep);
                                    unmatched_idx    = replaced_id(relabling_out.keep ~= 1);

                                    kept_idx = realigned_spk_idx(relabling_out.keep == 1, :);
                                    predLabels = predLabels(relabling_out.keep == 1, :);
                                    updatedLabels = updatedLabels(relabling_out.keep == 1, :);
                                    features = features(relabling_out.keep == 1, :);
                                    waveforms = realigned_waveform(relabling_out.keep == 1,:,:);
                                    amplitude = amplitude(relabling_out.keep == 1);

                                    output_waveform(relabeled_channel)          = [output_waveform(relabeled_channel); waveforms];
                                    output_spk_idx(relabeled_channel)           = [output_spk_idx(relabeled_channel); kept_idx];
                                    output_updatedLabels(relabeled_channel)     = [output_updatedLabels(relabeled_channel); updatedLabels];
                                    output_labels(relabeled_channel)            = [output_labels(relabeled_channel); predLabels];
                                    output_features(relabeled_channel)          = [output_features(relabeled_channel); features];
                                    output_amplitude(relabeled_channel)        = [output_amplitude(relabeled_channel); amplitude];
                                    unmatched_spk_idx{relabeled_channel, 1}     = [unmatched_spk_idx{relabeled_channel, 1}; unmatched_idx];

                                end

                            else
                                discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; relabled_idx];
                            end
                        end
                    end


                else

                    for iCh = (ch_idx - search_window) : num_channels
                        if channel_inclusion(iCh)  && chunk_channel_inclusion(iCh)
                            relabling_in = relabling_channels (iCh);
                            spike_idx_in = unevaluated_spk_idx(iCh);
                            num_relabling  = relabling_in.nRelabling;

                            for iRlbl = 1 : num_relabling

                                relabeled_channel = relabling_in.lookUpChannel(iRlbl);
                                not_kept_idx = relabling_in.not_kept_idx{iRlbl, 1};
                                relabled_idx = spike_idx_in(not_kept_idx);

                                if abs(relabeled_channel - iCh) <= search_window && relabeled_channel ~= iCh

                                    temp_spike_idx = spk_idx_channels(relabeled_channel);
                                    temp_spike_val = spk_Val_channels(relabeled_channel) .* spk_ID_channels(relabeled_channel);

                                    if isempty(temp_spike_idx)
                                        discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; relabled_idx];
                                        continue,
                                    end
                                    temp_waveform = waveform_channels(relabeled_channel).waveform;
                                    keep_alt_idx = keep_id_channels(relabeled_channel).keep;
                                    idx_keep = find(~keep_alt_idx);
                                    [d, alt_idx_subset] = nearest_index(relabled_idx, temp_spike_idx(idx_keep));
                                    valid_replacement = (d <= 1.5 * spikeDistance);
                                    replaced_id = idx_keep(alt_idx_subset(valid_replacement));
                                    replaced_idx = temp_spike_idx(replaced_id);
                                    amplitude    = temp_spike_val(replaced_id);
                                    non_replaced_idx = relabled_idx(d > 1.5 * spikeDistance);
                                    discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; non_replaced_idx];

                                    if any(valid_replacement)
                                        clusterSelection = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterSelection;
                                        clusterRelabeling  = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterRelabeling;
                                        data_sorting.waveform      = temp_waveform(replaced_id,:,:);
                                        data_sorting.waveformInfo  = sortedSamples{relabeled_channel,1}.waveformInfo;
                                        data_sorting     = kiaSort_preprocess_waveforms(data_sorting, cfg);
                                        data_sorting.PCA = sortedSamples{relabeled_channel, 1}.clusteringInfo.PCA;
                                        data_sorting.classifierInfo = sortedSamples{relabeled_channel, 1}.classifierInfo;
                                        data_sorting.amplitude     = amplitude;
                                        [predLabels, features, validKeep]      = kiaSort_predict_spike_labels(data_sorting, cfg);
                                        if iscategorical(predLabels)
                                            predLabels = str2double(cellstr(predLabels));
                                        end
                                        [updatedLabels, realigned_spk_idx, realigned_waveform, ~] = ...
                                            realignSpikes(predLabels, temp_waveform(replaced_id,:,:), replaced_idx, clusterRelabeling, cfg);

                                        relabling_out    = evaluate_labels(updatedLabels, clusterSelection, validKeep);
                                        unmatched_idx    = replaced_id(relabling_out.keep ~= 1);

                                        kept_idx = realigned_spk_idx(relabling_out.keep == 1, :);
                                        predLabels = predLabels(relabling_out.keep == 1, :);
                                        updatedLabels = updatedLabels(relabling_out.keep == 1, :);
                                        features = features(relabling_out.keep == 1, :);
                                        waveforms = realigned_waveform(relabling_out.keep == 1,:,:);
                                        amplitude = amplitude(relabling_out.keep == 1);

                                        output_waveform(relabeled_channel)          = [output_waveform(relabeled_channel); waveforms];
                                        output_spk_idx(relabeled_channel)           = [output_spk_idx(relabeled_channel); kept_idx];
                                        output_updatedLabels(relabeled_channel)     = [output_updatedLabels(relabeled_channel); updatedLabels];
                                        output_labels(relabeled_channel)            = [output_labels(relabeled_channel); predLabels];
                                        output_features(relabeled_channel)          = [output_features(relabeled_channel); features];
                                        output_amplitude(relabeled_channel)        = [output_amplitude(relabeled_channel); amplitude];
                                        unmatched_spk_idx{relabeled_channel, 1}     = [unmatched_spk_idx{relabeled_channel, 1}; unmatched_idx];

                                    end

                                else
                                    discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; relabled_idx];
                                end
                            end
                        end
                    end
                end
            end



            % reoder the sorted spikes the labels
            if ch_idx > release_window
                if ch_idx < num_channels

                    iCh = ch_idx - release_window;
                    if channel_inclusion(iCh) && chunk_channel_inclusion(iCh)
                        [output_spk_idx(iCh),...
                            output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                            output_amplitude(iCh)] = ...
                            reorder_sorted_data(output_spk_idx(iCh),...
                            output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                            output_amplitude(iCh));

                        unified_labeling = crossChannelStats.unified_labels;
                        unifiedLabels  = mapLabels(output_updatedLabels(iCh), iCh, unified_labeling);

                        if cfg.extractWaveform
                            sorted_out.waveforms{iCh, 1}        = output_waveform(iCh);
                        end
                        sorted_out.labels{iCh, 1}           = output_labels(iCh);
                        sorted_out.upadatedLabels{iCh, 1}   = output_updatedLabels(iCh);
                        sorted_out.features{iCh, 1}         = output_features(iCh);
                        sorted_out.amplitude{iCh, 1}        = output_amplitude(iCh);
                        sorted_out.unifiedLabels{iCh, 1}    = unifiedLabels;
                        sorted_out.spike_idx{iCh, 1}        = output_spk_idx(iCh) + start_sample - 1;
                        sorted_out.channelNum{iCh, 1}       = iCh * ones(size(output_spk_idx(iCh)));
                        % sorted_out.unmatched_spike_idx{iCh, 1}    = unmatched_spk_idx{iCh, 1} + start_sample - 1;
                        % sorted_out.discarded_spike_idx{iCh, 1}    = discarded_spk_idx{iCh, 1} + start_sample - 1;
                    end
                else
                    for iCh = ch_idx - release_window : num_channels
                        if channel_inclusion(iCh) && chunk_channel_inclusion(iCh)
                            [output_spk_idx(iCh),...
                                output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                                output_amplitude(iCh)] = ...
                                reorder_sorted_data(output_spk_idx(iCh),...
                                output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                                output_amplitude(iCh));

                            unified_labeling = crossChannelStats.unified_labels;
                            unifiedLabels  = mapLabels(output_updatedLabels(iCh), iCh, unified_labeling);

                            if cfg.extractWaveform
                                sorted_out.waveforms{iCh, 1}        = output_waveform(iCh);
                            end
                            sorted_out.labels{iCh, 1}           = output_labels(iCh);
                            sorted_out.upadatedLabels{iCh, 1}   = output_updatedLabels(iCh);
                            sorted_out.features{iCh, 1}         = output_features(iCh);
                            sorted_out.amplitude{iCh, 1}        = output_amplitude(iCh);
                            sorted_out.unifiedLabels{iCh, 1}    = unifiedLabels;
                            sorted_out.spike_idx{iCh, 1}        = output_spk_idx(iCh) + start_sample - 1;
                            sorted_out.channelNum{iCh, 1}       = iCh * ones(size(output_spk_idx(iCh)));
                            % sorted_out.unmatched_spike_idx{iCh, 1}    = unmatched_spk_idx{iCh, 1} + start_sample - 1;
                            % sorted_out.discarded_spike_idx{iCh, 1}    = discarded_spk_idx{iCh, 1} + start_sample - 1;
                        end
                    end
                end
            end

            fprintf('Processed channel %d/%d\n', ...
                ch_idx, num_channels);
            toc

            % remove channels that are no longer needed
            toRemoveChannels = setdiff(cell2mat(keys(bp_channels)), ch_indices_cluster);
            if ch_idx == num_channels
                toRemoveChannels = cell2mat(keys(bp_channels));
            end
            keys_to_remove = num2cell(toRemoveChannels);

            for j = 1:length(mapNames)
                currentMap = eval(mapNames{j});
                remove(currentMap, keys_to_remove);
                eval([mapNames{j} ' = currentMap;']);
            end

        end

        if ~exist('prevSize','var')
            fields = fieldnames(sorted_out);
            prevSize = zeros(size(fields));
        end

        currentSize = saveh5SpikeData(outputFolder, sorted_out, prevSize);
        prevSize = prevSize + currentSize;
        sorted_out = [];

        % Log end time for the chunk
        chunkEndTime = datestr(now);
        fprintf(fid, 'Chunk %d/%d ended at %s\n', chunk_i, num_chunks, chunkEndTime);

        if ~isempty(progressFcn)
            pct = chunk_i / num_chunks;
            msg = sprintf('Sorting data chunk %d of %d', chunk_i, num_chunks);
            progressFcn(pct, msg);
        end

    end

    disp('Sample chunk processed. Spike detection and extraction completed for all channels.');

    % Log final completion time of all chunks
    fprintf(fid, '\nAll chunks processed at %s\n', datestr(now));
    [wmsg,wid] = lastwarn;
    if ~isempty(wmsg)
        fprintf(fid, 'WARNING: %s (ID: %s)\n', wmsg, wid);
    end
    fclose(fid);

catch ME
    % Log errors
    fprintf(fid, '\nERROR OCCURRED:\n');
    fprintf(fid, 'Message: %s\n', ME.message);
    for s = 1:length(ME.stack)
        fprintf(fid, 'File: %s\nFunction: %s\nLine: %d\n', ...
            ME.stack(s).file, ME.stack(s).name, ME.stack(s).line);
    end
    fclose(fid);
    rethrow(ME);
end
end


function [waveOut, out_filter, out] = setWaveformNull()
waveOut.waveform         = [];
waveOut.wavformChanelIdx = [];
waveOut.channel_inclusion       = [];
waveOut.channel_thresholds_neg  = [];
waveOut.channel_thresholds_pos  = [];

out_filter.bandpass_signal = [];
out_filter.rms_bands = [];
out_filter.band_pairs = [];
out_filter.Ns_seq = [];
out_filter.scale_factor = [];

out.spk_idx      = [];
out.spk_Val      = [];
out.spk_ID       = [];
out.scale_factor = [];
out.rms_values   = [];
out.SNR          = [];
out.inclusion    = 0;
end

function selected_data = batch_extract(m, chunk_limits, start_ch, end_ch, channel_mapping, channel_inclusion, cfg)

selected_data = double(m.Data.data(channel_mapping(start_ch:end_ch), chunk_limits(1):chunk_limits(2)));

if cfg.extremeNoise
    selected_data(channel_inclusion(start_ch:end_ch),:) = remove_res_shared_noise(selected_data(channel_inclusion(start_ch:end_ch),:), cfg);
end

if cfg.denoising
    selected_data(channel_inclusion(start_ch:end_ch),:) = remove_shared_noise(selected_data(channel_inclusion(start_ch:end_ch),:), cfg);
end

end