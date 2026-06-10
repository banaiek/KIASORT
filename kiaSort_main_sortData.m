function kiaSort_main_sortData(inputPath, outputPath, cfg, varargin)

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

if isfield(cfg, 'modelType')
    if strcmp(cfg.modelType,'template')
        useTemplate = true;
    end
else
    useTemplate = false;
end

if isfield(cfg, 'residual_template')
    useResidualTemplate = cfg.residual_template;
else
    useResidualTemplate = true;
end

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

    if isfield(channel_info, 'final_bandpass') && ~isempty(channel_info.final_bandpass)
        cfg.bandpass = channel_info.final_bandpass;
    end

    num_channels = min(cfg.numChannels, length(channel_mapping));
    fs = cfg.samplingFrequency;

    m = map_input_file(inputPath, cfg);
    num_samples = size(m.Data.data,2);

    window_size = cfg.num_channel_extract * 2 + 1;
    half_window = floor(window_size / 2);
    search_window = min(round(1.5 * half_window), round(num_channels/2)-3);
    release_window = 2 * search_window;

    chunk_duration      = cfg.sortingChunkDuration;
    num_chunk_pts       = chunk_duration * fs;
    num_chunks          = ceil(num_samples / num_chunk_pts);
    marginOffset_pts    = cfg.borderMargin * fs / 1000;
    spikeDistance        = cfg.spikeDistance * fs / 1000;
    batch_ch_size       = min(cfg.batch_ch_size,num_channels);

    channel_idx = [1:batch_ch_size:(num_channels-1),num_channels];

    mapNames = {
        'bp_channels', 'spk_idx_channels', ...
        'spk_Val_channels', 'spk_ID_channels', ...
        'waveform_channels', 'relabling_channels', 'keep_id_channels', ...
        'output_waveform', 'output_spk_idx', ...
        'output_labels', 'output_updatedLabels', 'output_features',...
        'output_amplitude','unevaluated_spk_idx', 'out_filtered_channels', ...
        'output_bestVariant', ...
        'reserve_spk_idx_channels', 'reserve_spk_Val_channels', 'reserve_spk_ID_channels'
        };

    for i = 1:length(mapNames)
        eval([mapNames{i} ' = containers.Map(''KeyType'', ''double'', ''ValueType'', ''any'');']);
    end

    % Per-channel cache for prepareTemplateInfo to skip redundant rebuilds.
    tmplInfoCache = containers.Map('KeyType','double','ValueType','any');

    channel_thresholds      = cell(num_channels, 1);
    sorted_out              = [];

    for chunk_i = 1:num_chunks

        chunk_channel_inclusion = true(num_channels,1);

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

        for ch_idx = 1:num_channels
            tic
            ch_indices_mapped  = [max(1, ch_idx - half_window) : min(num_channels, ch_idx + half_window)]';
            ch_indices_cluster = [max(1, ch_idx - release_window - half_window) : min(num_channels, ch_idx + release_window)]';

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
                    end
                end

            end

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
                    out_filtered_channels(mapped_ch_idx)    = out_filtered;
                    reserve_spk_idx_channels(mapped_ch_idx) = out_detected.reserve_spk_idx;
                    reserve_spk_Val_channels(mapped_ch_idx) = out_detected.reserve_spk_Val;
                    reserve_spk_ID_channels(mapped_ch_idx)  = out_detected.reserve_spk_ID;
                end
            end

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

                    if useTemplate
                        data_sorting.waveform = waveOut.waveform;
                        if isKey(tmplInfoCache, ch_idx)
                            data_sorting.templateInfo = tmplInfoCache(ch_idx);
                        else
                            data_sorting.templateInfo = prepareTemplateInfo(sortedSamples{ch_idx, 1});
                            tmplInfoCache(ch_idx) = data_sorting.templateInfo;
                        end
                        data_sorting.amplitude = amplitude;
                        [predLabels, features, validKeep, bestVarIdx] = kiaSort_template_matching(data_sorting, cfg);
                    else
                        data_sorting = kiaSort_preprocess_waveforms(waveOut, cfg);
                        data_sorting.classifierInfo = sortedSamples{ch_idx, 1}.classifierInfo;
                        data_sorting.PCA = sortedSamples{ch_idx, 1}.clusteringInfo.PCA;
                        data_sorting.amplitude = amplitude;
                        [predLabels, features, validKeep] = kiaSort_predict_spike_labels(data_sorting, cfg);
                        bestVarIdx = ones(size(predLabels));
                    end

                    if iscategorical(predLabels)
                        predLabels = str2double(cellstr(predLabels));
                    end
                    clusterSelection = sortedSamples{ch_idx, 1}.clusteringInfo.clusterSelection;
                    clusterRelabeling = sortedSamples{ch_idx, 1}.clusteringInfo.clusterRelabeling;
                    [updatedLabels, realigned_spk_idx, realigned_waveform, ~] = ...
                        realignSpikes(predLabels, waveOut.waveform, sorting_spk_idx, clusterRelabeling, cfg);

                    relabling_out = kiaSort_evaluate_labels(updatedLabels, clusterSelection, altChannels, validKeep);
                    relabling_channels(ch_idx) = relabling_out;

                    output_labels(ch_idx)              = predLabels(relabling_out.keep);
                    output_updatedLabels(ch_idx)       = updatedLabels(relabling_out.keep);
                    output_features(ch_idx)            = features(relabling_out.keep, :);
                    output_spk_idx(ch_idx)             = realigned_spk_idx(relabling_out.keep);
                    output_waveform(ch_idx)            = realigned_waveform(relabling_out.keep,:,:);
                    output_amplitude(ch_idx)           = amplitude(relabling_out.keep);
                    output_bestVariant(ch_idx)         = bestVarIdx(relabling_out.keep);

                else

                    relabling_channels(ch_idx)         = [];
                    output_labels(ch_idx)              = [];
                    output_updatedLabels(ch_idx)       = [];
                    output_features(ch_idx)            = [];
                    output_spk_idx(ch_idx)             = [];
                    output_waveform(ch_idx)            = [];
                    output_amplitude(ch_idx)           = [];
                    output_bestVariant(ch_idx)         = [];
                end

            end

            if ch_idx > search_window

                if ch_idx < num_channels
                    iCh = (ch_idx - search_window);
                    if channel_inclusion(iCh) && chunk_channel_inclusion(iCh)

                        relabling_in = relabling_channels(iCh);
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
                                    searchReserveSpikes(relabled_idx, relabeled_channel, ...
                                        reserve_spk_idx_channels, reserve_spk_Val_channels, reserve_spk_ID_channels, ...
                                        bp_channels, spk_idx_channels, spk_ID_channels, ...
                                        channel_thresholds, channel_inclusion, channel_info, ...
                                        sortedSamples, num_channels, spikeDistance, useTemplate, cfg, ...
                                        output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
                                        output_features, output_amplitude, output_bestVariant);
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

                                if ~isempty(non_replaced_idx)
                                    searchReserveSpikes(non_replaced_idx, relabeled_channel, ...
                                        reserve_spk_idx_channels, reserve_spk_Val_channels, reserve_spk_ID_channels, ...
                                        bp_channels, spk_idx_channels, spk_ID_channels, ...
                                        channel_thresholds, channel_inclusion, channel_info, ...
                                        sortedSamples, num_channels, spikeDistance, useTemplate, cfg, ...
                                        output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
                                        output_features, output_amplitude, output_bestVariant);
                                end

                                if any(valid_replacement)
                                    clusterSelection = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterSelection;
                                    clusterRelabeling  = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterRelabeling;

                                    if useTemplate
                                        data_sorting.waveform = temp_waveform(replaced_id,:,:);
                                        if isKey(tmplInfoCache, relabeled_channel)
                                            data_sorting.templateInfo = tmplInfoCache(relabeled_channel);
                                        else
                                            data_sorting.templateInfo = prepareTemplateInfo(sortedSamples{relabeled_channel, 1});
                                            tmplInfoCache(relabeled_channel) = data_sorting.templateInfo;
                                        end
                                        data_sorting.amplitude = amplitude;
                                        [predLabels, features, validKeep, bestVarIdx] = kiaSort_template_matching(data_sorting, cfg);
                                    else
                                        data_sorting.waveform = temp_waveform(replaced_id,:,:);
                                        data_sorting.waveformInfo = sortedSamples{relabeled_channel,1}.waveformInfo;
                                        data_sorting = kiaSort_preprocess_waveforms(data_sorting, cfg);
                                        data_sorting.PCA = sortedSamples{relabeled_channel, 1}.clusteringInfo.PCA;
                                        data_sorting.classifierInfo = sortedSamples{relabeled_channel, 1}.classifierInfo;
                                        data_sorting.amplitude = amplitude;
                                        [predLabels, features, validKeep] = kiaSort_predict_spike_labels(data_sorting, cfg);
                                        bestVarIdx = ones(size(predLabels));
                                    end

                                    if iscategorical(predLabels)
                                        predLabels = str2double(cellstr(predLabels));
                                    end
                                    [updatedLabels, realigned_spk_idx, realigned_waveform, ~] = ...
                                        realignSpikes(predLabels, temp_waveform(replaced_id,:,:), replaced_idx, clusterRelabeling, cfg);

                                    relabling_out    = evaluate_labels(updatedLabels, clusterSelection, validKeep);

                                    kept_idx = realigned_spk_idx(relabling_out.keep == 1, :);
                                    predLabels = predLabels(relabling_out.keep == 1, :);
                                    updatedLabels = updatedLabels(relabling_out.keep == 1, :);
                                    features = features(relabling_out.keep == 1, :);
                                    waveforms = realigned_waveform(relabling_out.keep == 1,:,:);
                                    amplitude = amplitude(relabling_out.keep == 1);
                                    bestVarIdx = bestVarIdx(relabling_out.keep == 1);

                                    output_waveform(relabeled_channel)          = [output_waveform(relabeled_channel); waveforms];
                                    output_spk_idx(relabeled_channel)           = [output_spk_idx(relabeled_channel); kept_idx];
                                    output_updatedLabels(relabeled_channel)     = [output_updatedLabels(relabeled_channel); updatedLabels];
                                    output_labels(relabeled_channel)            = [output_labels(relabeled_channel); predLabels];
                                    output_features(relabeled_channel)          = [output_features(relabeled_channel); features];
                                    output_amplitude(relabeled_channel)         = [output_amplitude(relabeled_channel); amplitude];
                                    output_bestVariant(relabeled_channel)       = [output_bestVariant(relabeled_channel); bestVarIdx];

                                    secondaryTransferStep(relabling_out, replaced_idx, relabeled_channel, ...
                                        spk_idx_channels, spk_Val_channels, spk_ID_channels, waveform_channels, ...
                                        keep_id_channels, sortedSamples, spikeDistance, useTemplate, cfg, ...
                                        output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
                                        output_features, output_amplitude, output_bestVariant);

                                end

                            else
                                discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; relabled_idx];
                            end
                        end
                    end


                else

                    for iCh = (ch_idx - search_window) : num_channels
                        if channel_inclusion(iCh)  && chunk_channel_inclusion(iCh)
                            relabling_in = relabling_channels(iCh);
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
                                        searchReserveSpikes(relabled_idx, relabeled_channel, ...
                                            reserve_spk_idx_channels, reserve_spk_Val_channels, reserve_spk_ID_channels, ...
                                            bp_channels, spk_idx_channels, spk_ID_channels, ...
                                            channel_thresholds, channel_inclusion, channel_info, ...
                                            sortedSamples, num_channels, spikeDistance, useTemplate, cfg, ...
                                            output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
                                            output_features, output_amplitude, output_bestVariant);
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

                                    if ~isempty(non_replaced_idx)
                                        searchReserveSpikes(non_replaced_idx, relabeled_channel, ...
                                            reserve_spk_idx_channels, reserve_spk_Val_channels, reserve_spk_ID_channels, ...
                                            bp_channels, spk_idx_channels, spk_ID_channels, ...
                                            channel_thresholds, channel_inclusion, channel_info, ...
                                            sortedSamples, num_channels, spikeDistance, useTemplate, cfg, ...
                                            output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
                                            output_features, output_amplitude, output_bestVariant);
                                    end

                                    if any(valid_replacement)
                                        clusterSelection = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterSelection;
                                        clusterRelabeling  = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterRelabeling;

                                        if useTemplate
                                            data_sorting.waveform = temp_waveform(replaced_id,:,:);
                                            if isKey(tmplInfoCache, relabeled_channel)
                                                data_sorting.templateInfo = tmplInfoCache(relabeled_channel);
                                            else
                                                data_sorting.templateInfo = prepareTemplateInfo(sortedSamples{relabeled_channel, 1});
                                                tmplInfoCache(relabeled_channel) = data_sorting.templateInfo;
                                            end
                                            data_sorting.amplitude = amplitude;
                                            [predLabels, features, validKeep, bestVarIdx] = kiaSort_template_matching(data_sorting, cfg);
                                        else
                                            data_sorting.waveform = temp_waveform(replaced_id,:,:);
                                            data_sorting.waveformInfo = sortedSamples{relabeled_channel,1}.waveformInfo;
                                            data_sorting = kiaSort_preprocess_waveforms(data_sorting, cfg);
                                            data_sorting.PCA = sortedSamples{relabeled_channel, 1}.clusteringInfo.PCA;
                                            data_sorting.classifierInfo = sortedSamples{relabeled_channel, 1}.classifierInfo;
                                            data_sorting.amplitude = amplitude;
                                            [predLabels, features, validKeep] = kiaSort_predict_spike_labels(data_sorting, cfg);
                                            bestVarIdx = ones(size(predLabels));
                                        end

                                        if iscategorical(predLabels)
                                            predLabels = str2double(cellstr(predLabels));
                                        end
                                        [updatedLabels, realigned_spk_idx, realigned_waveform, ~] = ...
                                            realignSpikes(predLabels, temp_waveform(replaced_id,:,:), replaced_idx, clusterRelabeling, cfg);

                                        relabling_out    = evaluate_labels(updatedLabels, clusterSelection, validKeep);

                                        kept_idx = realigned_spk_idx(relabling_out.keep == 1, :);
                                        predLabels = predLabels(relabling_out.keep == 1, :);
                                        updatedLabels = updatedLabels(relabling_out.keep == 1, :);
                                        features = features(relabling_out.keep == 1, :);
                                        waveforms = realigned_waveform(relabling_out.keep == 1,:,:);
                                        amplitude = amplitude(relabling_out.keep == 1);
                                        bestVarIdx = bestVarIdx(relabling_out.keep == 1);

                                        output_waveform(relabeled_channel)          = [output_waveform(relabeled_channel); waveforms];
                                        output_spk_idx(relabeled_channel)           = [output_spk_idx(relabeled_channel); kept_idx];
                                        output_updatedLabels(relabeled_channel)     = [output_updatedLabels(relabeled_channel); updatedLabels];
                                        output_labels(relabeled_channel)            = [output_labels(relabeled_channel); predLabels];
                                        output_features(relabeled_channel)          = [output_features(relabeled_channel); features];
                                        output_amplitude(relabeled_channel)         = [output_amplitude(relabeled_channel); amplitude];
                                        output_bestVariant(relabeled_channel)       = [output_bestVariant(relabeled_channel); bestVarIdx];

                                        secondaryTransferStep(relabling_out, replaced_idx, relabeled_channel, ...
                                            spk_idx_channels, spk_Val_channels, spk_ID_channels, waveform_channels, ...
                                            keep_id_channels, sortedSamples, spikeDistance, useTemplate, cfg, ...
                                            output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
                                            output_features, output_amplitude, output_bestVariant);

                                    end

                                else
                                    discarded_spk_idx{iCh, 1} = [discarded_spk_idx{iCh, 1}; relabled_idx];
                                end
                            end
                        end
                    end
                end
            end

            if ch_idx > release_window
                if ch_idx < num_channels

                    iCh = ch_idx - release_window;
                    if channel_inclusion(iCh) && chunk_channel_inclusion(iCh)
                        
                        if useResidualTemplate
                            [output_spk_idx, output_waveform, output_features, output_labels, ...
                             output_updatedLabels, output_amplitude, output_bestVariant] = ...
                                performSecondRound(iCh, bp_channels, output_spk_idx, output_waveform, ...
                                output_features, output_labels, output_updatedLabels, output_amplitude, ...
                                output_bestVariant, ...
                                out_filtered_channels, spk_idx_channels, channel_thresholds, channel_inclusion, ...
                                sortedSamples, num_channels, search_window, half_window, useTemplate, cfg);
                        end

                        [output_spk_idx(iCh),...
                            output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                            output_amplitude(iCh), output_bestVariant(iCh)] = ...
                            reorder_sorted_data(output_spk_idx(iCh),...
                            output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                            output_amplitude(iCh), output_bestVariant(iCh));

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
                    end
                else
                    for iCh = ch_idx - release_window : num_channels
                        if channel_inclusion(iCh) && chunk_channel_inclusion(iCh)
                            
                            if useResidualTemplate
                                [output_spk_idx, output_waveform, output_features, output_labels, ...
                                 output_updatedLabels, output_amplitude, output_bestVariant] = ...
                                    performSecondRound(iCh, bp_channels, output_spk_idx, output_waveform, ...
                                    output_features, output_labels, output_updatedLabels, output_amplitude, ...
                                    output_bestVariant, ...
                                    out_filtered_channels, spk_idx_channels, channel_thresholds, channel_inclusion, ...
                                    sortedSamples, num_channels, search_window, half_window, useTemplate, cfg);
                            end
                            
                            [output_spk_idx(iCh),...
                                output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                                output_amplitude(iCh), output_bestVariant(iCh)] = ...
                                reorder_sorted_data(output_spk_idx(iCh),...
                                output_waveform(iCh),  output_features(iCh), output_labels(iCh), output_updatedLabels(iCh),...
                                output_amplitude(iCh), output_bestVariant(iCh));

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
                        end
                    end
                end
            end

            fprintf('Processed channel %d/%d\n', ch_idx, num_channels);
            toc

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

        chunkEndTime = datestr(now);
        fprintf(fid, 'Chunk %d/%d ended at %s\n', chunk_i, num_chunks, chunkEndTime);

        if ~isempty(progressFcn)
            pct = chunk_i / num_chunks;
            msg = sprintf('Sorting data chunk %d of %d', chunk_i, num_chunks);
            progressFcn(pct, msg);
        end

    end

    disp('Sample chunk processed. Spike detection and extraction completed for all channels.');

    fprintf(fid, '\nAll chunks processed at %s\n', datestr(now));
    [wmsg,wid] = lastwarn;
    if ~isempty(wmsg)
        fprintf(fid, 'WARNING: %s (ID: %s)\n', wmsg, wid);
    end
    fclose(fid);

catch ME
    fprintf(fid, '\nERROR OCCURRED:\n');
    fprintf(fid, 'Message: %s\n', ME.message);
    for s = 1:length(ME.stack)
        fprintf(fid, 'File: %s\nFunction: %s\nLine: %d\n', ...
            ME.stack(s).file, ME.stack(s).name, ME.stack(s).line);
    end
    fclose(fid);
    rethrow(ME);
end

if isfield(cfg,'postHocProcessing')
    if cfg.postHocProcessing
        if ~cfg.sort_only
            kiaSort_drift_merge_posthoc_iterative(outputPath, ...
                'overwrite', true, ...
                'verbose',   false, ...
                'mainArgs',  {'debugFigs', false});
            try
                kiaSort_post_sort_curate(outputPath, ...
                    'ccg_cleaning',   true, ...
                    'merging',        true, ...
                    'xcorrThreshold', 0.9, ...
                    'verbose',        false);
            catch

            end
        end
    end
end

end


%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function secondaryTransferStep(relabling_out, replaced_idx, relabeled_channel, ...
    spk_idx_channels, spk_Val_channels, spk_ID_channels, waveform_channels, ...
    keep_id_channels, sortedSamples, spikeDistance, useTemplate, cfg, ...
    output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
    output_features, output_amplitude, output_bestVariant)

if relabling_out.nRelabling < 1
    return;
end

for iRlbl2 = 1:relabling_out.nRelabling
    sec_channel = relabling_out.lookUpChannel(iRlbl2);
    sec_not_kept = relabling_out.not_kept_idx{iRlbl2};
    sec_spk_times = replaced_idx(sec_not_kept);

    if sec_channel == relabeled_channel
        continue;
    end
    if ~isKey(spk_idx_channels, sec_channel) || ~isKey(waveform_channels, sec_channel)
        continue;
    end

    sec_spike_idx = spk_idx_channels(sec_channel);
    if isempty(sec_spike_idx), continue; end

    sec_spike_val = spk_Val_channels(sec_channel) .* spk_ID_channels(sec_channel);
    sec_wf = waveform_channels(sec_channel).waveform;
    sec_keep = keep_id_channels(sec_channel).keep;
    sec_pool = find(~sec_keep);
    if isempty(sec_pool), continue; end

    [d2, idx2] = nearest_index(sec_spk_times, sec_spike_idx(sec_pool));
    valid2 = (d2 <= 1.5 * spikeDistance);
    if ~any(valid2), continue; end

    sec_replaced_id  = sec_pool(idx2(valid2));
    sec_replaced_idx = sec_spike_idx(sec_replaced_id);
    sec_amplitude    = sec_spike_val(sec_replaced_id);

    sec_clusterSelection  = sortedSamples{sec_channel, 1}.clusteringInfo.clusterSelection;
    sec_clusterRelabeling = sortedSamples{sec_channel, 1}.clusteringInfo.clusterRelabeling;

    if useTemplate
        ds.waveform     = sec_wf(sec_replaced_id,:,:);
        ds.templateInfo = prepareTemplateInfo(sortedSamples{sec_channel, 1});
        ds.amplitude    = sec_amplitude;
        [pL, ft, vK, bVI] = kiaSort_template_matching(ds, cfg);
    else
        ds.waveform     = sec_wf(sec_replaced_id,:,:);
        ds.waveformInfo = sortedSamples{sec_channel,1}.waveformInfo;
        ds = kiaSort_preprocess_waveforms(ds, cfg);
        ds.PCA             = sortedSamples{sec_channel, 1}.clusteringInfo.PCA;
        ds.classifierInfo  = sortedSamples{sec_channel, 1}.classifierInfo;
        ds.amplitude       = sec_amplitude;
        [pL, ft, vK] = kiaSort_predict_spike_labels(ds, cfg);
        bVI = ones(size(pL));
    end

    if iscategorical(pL), pL = str2double(cellstr(pL)); end

    [uL, rIdx, rWf, ~] = realignSpikes(pL, sec_wf(sec_replaced_id,:,:), sec_replaced_idx, sec_clusterRelabeling, cfg);
    rOut = evaluate_labels(uL, sec_clusterSelection, vK);

    kept = rOut.keep == 1;
    if any(kept)
        output_waveform(sec_channel)       = [output_waveform(sec_channel); rWf(kept,:,:)];
        output_spk_idx(sec_channel)        = [output_spk_idx(sec_channel); rIdx(kept)];
        output_updatedLabels(sec_channel)  = [output_updatedLabels(sec_channel); uL(kept)];
        output_labels(sec_channel)         = [output_labels(sec_channel); pL(kept)];
        output_features(sec_channel)       = [output_features(sec_channel); ft(kept,:)];
        output_amplitude(sec_channel)      = [output_amplitude(sec_channel); sec_amplitude(kept)];
        output_bestVariant(sec_channel)    = [output_bestVariant(sec_channel); bVI(kept)];
    end
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

% Whitening uses an overlapping batch: pull a margin of extra channels on
% each side, whiten the extended slice, return only the central rows so
% cross-batch noise components actually get decorrelated.
needWhitening = (isfield(cfg,'extremeNoise') && cfg.extremeNoise) ...
             || (isfield(cfg,'denoising')    && cfg.denoising);

if ~needWhitening
    selected_data = double(m.Data.data(channel_mapping(start_ch:end_ch), ...
                                       chunk_limits(1):chunk_limits(2)));
    return;
end

% Default margin = 0 -> overlap path is opt-in via cfg.denoising_margin_ch.
% With margin = 0 the fallback below runs, which is byte-identical to the
% original non-overlapping behaviour.
if isfield(cfg, 'denoising_margin_ch') && ~isempty(cfg.denoising_margin_ch)
    margin = max(0, round(cfg.denoising_margin_ch));
else
    margin = 0;
end

maxCh     = min(numel(channel_mapping), numel(channel_inclusion));
ext_start = max(1, start_ch - margin);
ext_end   = min(maxCh, end_ch + margin);

if ext_start >= start_ch && ext_end <= end_ch
    selected_data = double(m.Data.data(channel_mapping(start_ch:end_ch), ...
                                       chunk_limits(1):chunk_limits(2)));
    if cfg.extremeNoise
        selected_data(channel_inclusion(start_ch:end_ch), :) = ...
            remove_res_shared_noise(selected_data(channel_inclusion(start_ch:end_ch), :), cfg);
    end
    if cfg.denoising
        selected_data(channel_inclusion(start_ch:end_ch), :) = ...
            remove_shared_noise(selected_data(channel_inclusion(start_ch:end_ch), :), cfg);
    end
    return;
end

ext_data = double(m.Data.data(channel_mapping(ext_start:ext_end), ...
                              chunk_limits(1):chunk_limits(2)));

ext_inclusion = channel_inclusion(ext_start:ext_end);
incl_idx = find(ext_inclusion);
if ~isempty(incl_idx)
    sub = ext_data(incl_idx, :);
    if cfg.extremeNoise
        sub = remove_res_shared_noise(sub, cfg);
    end
    if cfg.denoising
        sub = remove_shared_noise(sub, cfg);
    end
    ext_data(incl_idx, :) = sub;
end

central_offset = start_ch - ext_start + 1;
central_len    = end_ch - start_ch + 1;
selected_data  = ext_data(central_offset : central_offset + central_len - 1, :);

end

function templateInfo = prepareTemplateInfo(sortedSample)

if isfield(sortedSample.waveformInfo, 'templateWaveforms')
    templateInfo.templates = sortedSample.waveformInfo.templateWaveforms;
    if isfield(sortedSample.waveformInfo, 'templateWeights')
        templateInfo.templateWeights = sortedSample.waveformInfo.templateWeights;
    end
else
    meanWF = sortedSample.waveformInfo.meanWaveform;
    [nClasses, C, T] = size(meanWF);
    templateInfo.templates = reshape(meanWF, [nClasses, 1, C, T]);
end

templateInfo.classLabels = sortedSample.clusteringInfo.classLabels(:);
templateInfo.clusterStatus = sortedSample.clusteringInfo.clusterRelabeling.changeType(:);

if isfield(sortedSample, 'classifierInfo')
    templateInfo.polarity   = sortedSample.classifierInfo.class_polarity(:);
    templateInfo.lowAmpThr  = sortedSample.classifierInfo.lowAmpThr(:);
    templateInfo.highAmpThr = sortedSample.classifierInfo.highAmpThr(:);
else
    nClasses = length(templateInfo.classLabels);
    templateInfo.polarity = ones(nClasses, 1);
    for i = 1:nClasses
        if templateInfo.classLabels(i) == -1
            continue;
        end
        if ndims(templateInfo.templates) == 4
            template_i = squeeze(templateInfo.templates(i, 1, :, :));
        else
            template_i = squeeze(templateInfo.templates(i, :, :));
        end
        if abs(min(template_i(:))) > abs(max(template_i(:)))
            templateInfo.polarity(i) = -1;
        end
    end
    templateInfo.lowAmpThr  = -inf(nClasses, 1);
    templateInfo.highAmpThr = inf(nClasses, 1);
end

end


function searchReserveSpikes(transfer_idx, relabeled_channel, ...
    reserve_spk_idx_channels, reserve_spk_Val_channels, reserve_spk_ID_channels, ...
    bp_channels, spk_idx_channels, spk_ID_channels, ...
    channel_thresholds, channel_inclusion, channel_info, ...
    sortedSamples, num_channels, spikeDistance, useTemplate, cfg, ...
    output_waveform, output_spk_idx, output_updatedLabels, output_labels, ...
    output_features, output_amplitude, output_bestVariant)
%SEARCHRESERVESPIKES  Search sub-threshold reserve spikes for transfer matches.
%
%   For transferred spikes that could not be matched to above-threshold
%   detections on the reference channel, searches the reserve pool (spikes
%   detected at 0.5x threshold). Matched reserve spikes are extracted,
%   classified, and appended to output.

if isempty(transfer_idx) || ~isKey(reserve_spk_idx_channels, relabeled_channel)
    return;
end

reserve_idx = reserve_spk_idx_channels(relabeled_channel);
if isempty(reserve_idx)
    return;
end

if ~isKey(bp_channels, relabeled_channel)
    return;
end

[d, res_subset] = nearest_index(transfer_idx, reserve_idx);
valid = (d <= 1.5 * spikeDistance);
if ~any(valid)
    return;
end

matched_idx = reserve_idx(res_subset(valid));

fs  = cfg.samplingFrequency;
num_channel_extract = cfg.num_channel_extract;

fe_ch_range = max(1, relabeled_channel - num_channel_extract) : ...
              min(num_channels, relabeled_channel + num_channel_extract);
n_fe_ch = length(fe_ch_range);

inputWF_fe.main_channel_idx     = relabeled_channel;
inputWF_fe.all_channel_idx      = fe_ch_range(:);
inputWF_fe.Ns_seq               = channel_info.Ns_seq;
inputWF_fe.bandpass_data_window = cell(n_fe_ch, 1);
inputWF_fe.spk_idx_all_channels = cell(n_fe_ch, 1);
inputWF_fe.spk_ID_all_channels  = cell(n_fe_ch, 1);
inputWF_fe.inclusion            = zeros(n_fe_ch, 1);
inputWF_fe.channel_thresholds   = cell(n_fe_ch, 1);

for idx = 1:n_fe_ch
    ch = fe_ch_range(idx);
    if isKey(bp_channels, ch)
        inputWF_fe.bandpass_data_window{idx, 1} = bp_channels(ch);
    else
        inputWF_fe.bandpass_data_window{idx, 1} = [];
    end
    if ch == relabeled_channel
        inputWF_fe.spk_idx_all_channels{idx, 1} = matched_idx(:);
        inputWF_fe.spk_ID_all_channels{idx, 1}  = [];
    elseif isKey(spk_idx_channels, ch)
        inputWF_fe.spk_idx_all_channels{idx, 1} = spk_idx_channels(ch);
        inputWF_fe.spk_ID_all_channels{idx, 1}  = [];
    else
        inputWF_fe.spk_idx_all_channels{idx, 1} = [];
        inputWF_fe.spk_ID_all_channels{idx, 1}  = [];
    end
    inputWF_fe.inclusion(idx, 1) = channel_inclusion(ch);
    if ~isempty(channel_thresholds{ch, 1})
        inputWF_fe.channel_thresholds{idx, 1} = channel_thresholds{ch, 1};
    else
        inputWF_fe.channel_thresholds{idx, 1} = struct('bandpass_min_threshold', 0);
    end
end

[waveOut_fe] = kiaSort_waveform_extraction(cfg, inputWF_fe, 1, 0);
if isempty(waveOut_fe.waveform) || size(waveOut_fe.waveform, 1) == 0
    return;
end

mid_ch_wf = ceil(size(waveOut_fe.waveform, 2) / 2);
mid_t_wf  = ceil(size(waveOut_fe.waveform, 3) / 2);
fe_amplitude = squeeze(waveOut_fe.waveform(:, mid_ch_wf, mid_t_wf));

waveOut_fe.waveformInfo = sortedSamples{relabeled_channel, 1}.waveformInfo;

clusterSelection  = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterSelection;
clusterRelabeling = sortedSamples{relabeled_channel, 1}.clusteringInfo.clusterRelabeling;

if useTemplate
    data_fe.waveform     = waveOut_fe.waveform;
    data_fe.templateInfo = prepareTemplateInfo(sortedSamples{relabeled_channel, 1});
    data_fe.amplitude    = fe_amplitude;
    [predLabels_fe, features_fe, validKeep_fe, bestVarIdx_fe] = kiaSort_template_matching(data_fe, cfg);
else
    data_fe = struct();
    data_fe.waveform     = waveOut_fe.waveform;
    data_fe.waveformInfo = sortedSamples{relabeled_channel, 1}.waveformInfo;
    data_fe = kiaSort_preprocess_waveforms(data_fe, cfg);
    data_fe.classifierInfo = sortedSamples{relabeled_channel, 1}.classifierInfo;
    data_fe.PCA            = sortedSamples{relabeled_channel, 1}.clusteringInfo.PCA;
    data_fe.amplitude      = fe_amplitude;
    [predLabels_fe, features_fe, validKeep_fe] = kiaSort_predict_spike_labels(data_fe, cfg);
    bestVarIdx_fe = ones(size(predLabels_fe));
end

if iscategorical(predLabels_fe)
    predLabels_fe = str2double(cellstr(predLabels_fe));
end

[updatedLabels_fe, realigned_spk_idx_fe, realigned_waveform_fe, ~] = ...
    realignSpikes(predLabels_fe, waveOut_fe.waveform, matched_idx(:), clusterRelabeling, cfg);

relabling_out_fe = evaluate_labels(updatedLabels_fe, clusterSelection, validKeep_fe);

kept = relabling_out_fe.keep == 1;
if any(kept)
    output_waveform(relabeled_channel)      = [output_waveform(relabeled_channel); realigned_waveform_fe(kept,:,:)];
    output_spk_idx(relabeled_channel)       = [output_spk_idx(relabeled_channel); realigned_spk_idx_fe(kept)];
    output_updatedLabels(relabeled_channel) = [output_updatedLabels(relabeled_channel); updatedLabels_fe(kept)];
    output_labels(relabeled_channel)        = [output_labels(relabeled_channel); predLabels_fe(kept)];
    output_features(relabeled_channel)      = [output_features(relabeled_channel); features_fe(kept,:)];
    output_amplitude(relabeled_channel)     = [output_amplitude(relabeled_channel); fe_amplitude(kept)];
    output_bestVariant(relabeled_channel)   = [output_bestVariant(relabeled_channel); bestVarIdx_fe(kept)];
end

end


function [output_spk_idx, output_waveform, output_features, output_labels, ...
          output_updatedLabels, output_amplitude, output_bestVariant] = ...
    performSecondRound(iCh, bp_channels, output_spk_idx, output_waveform, ...
    output_features, output_labels, output_updatedLabels, output_amplitude, ...
    output_bestVariant, ...
    out_filtered_channels, spk_idx_channels, channel_thresholds, channel_inclusion, ...
    sortedSamples, num_channels, search_window, half_window, useTemplate, cfg)
%PERFORMSECONDROUND  Template-matched subtraction + re-detection on residual.
%
%   Supports both template and classifier paths:
%   - Template path: uses best-variant index stored in output_bestVariant
%     for per-spike subtraction; recovered spikes classified via template matching.
%   - Classifier path: subtraction uses meanWaveform (variant 1 for all spikes);
%     recovered spikes classified via kiaSort_predict_spike_labels.

fs  = cfg.samplingFrequency;
num_channel_extract = cfg.num_channel_extract;
middle_channel      = num_channel_extract + 1;

spikeDuration_pts = round(cfg.spikeDuration * fs / 1000) + 1;
midPoint          = floor(spikeDuration_pts / 2) + 1;
T = spikeDuration_pts;

spike_length_cl = floor(cfg.clusteringSpikeDuration * fs / (2*1000));
cl_start = max(1, midPoint - spike_length_cl);
cl_end   = min(T, midPoint + spike_length_cl);

wf_range          = max(1, iCh - half_window) : min(num_channels, iCh + half_window);
subtraction_range = max(1, iCh - search_window) : min(num_channels, iCh + search_window);
wf_range_start    = wf_range(1);
n_wf_ch           = length(wf_range);
spkDist           = cfg.spikeDistance * fs / 1000;

if ~isKey(bp_channels, iCh), return; end

% =========================================================================
% PHASE 1: Build local bp_matrix
% =========================================================================
sig_len   = length(bp_channels(iCh));
bp_matrix = zeros(n_wf_ch, sig_len);
for i = 1:n_wf_ch
    ch = wf_range(i);
    if isKey(bp_channels, ch)
        tmp = bp_channels(ch);
        if numel(tmp) == sig_len
            bp_matrix(i, :) = tmp(:)';
        end
    end
end

% =========================================================================
% PHASE 2: Collect all first-round sorted spikes in subtraction range
%          Now also collects the stored best-variant index per spike.
% =========================================================================
total_count = 0;
for sub_ch = subtraction_range
    if isKey(output_spk_idx, sub_ch)
        s = output_spk_idx(sub_ch);
        if ~isempty(s), total_count = total_count + numel(s); end
    end
end
if total_count == 0, return; end

all_spk_times     = zeros(total_count, 1);
all_spk_labels    = zeros(total_count, 1);
all_spk_variants  = zeros(total_count, 1);
all_spk_amps      = zeros(total_count, 1);
all_spk_source_ch = zeros(total_count, 1);
pos = 0;
for sub_ch = subtraction_range
    if ~isKey(output_spk_idx, sub_ch), continue; end
    idx_ch = output_spk_idx(sub_ch);
    if isempty(idx_ch), continue; end
    n = numel(idx_ch);
    rng = pos+1 : pos+n;
    all_spk_times(rng)     = idx_ch(:);
    all_spk_labels(rng)    = output_labels(sub_ch);
    all_spk_amps(rng)      = output_amplitude(sub_ch);
    all_spk_source_ch(rng) = sub_ch;
    if isKey(output_bestVariant, sub_ch) && ~isempty(output_bestVariant(sub_ch))
        all_spk_variants(rng) = output_bestVariant(sub_ch);
    else
        all_spk_variants(rng) = 1;
    end
    pos = pos + n;
end
all_spk_times     = all_spk_times(1:pos);
all_spk_labels    = all_spk_labels(1:pos);
all_spk_variants  = all_spk_variants(1:pos);
all_spk_amps      = all_spk_amps(1:pos);
all_spk_source_ch = all_spk_source_ch(1:pos);
N_total = pos;
if N_total == 0, return; end

t_starts_all = all_spk_times - midPoint + 1;
t_ends_all   = t_starts_all + T - 1;
valid_spk = (all_spk_labels > 0) & (t_starts_all >= 1) & (t_ends_all <= sig_len);
if ~any(valid_spk), return; end

% =========================================================================
% PHASE 3a: Group by (source_ch, label, VARIANT)
%
%   Because each spike already knows its variant, each group maps to
%   exactly ONE template. No variant loop needed — same structure as
%   single-template subtraction, fully vectorised alpha per group.
%
%   prepareTemplateInfo is cached per source_ch to avoid redundant calls.
% =========================================================================

% Encode group key as (source_ch, label, variant) triple
max_lbl = max(all_spk_labels) + 1;
max_var = max(all_spk_variants) + 1;
grp_key = all_spk_source_ch * (max_lbl * max_var) + all_spk_labels * max_var + all_spk_variants;
grp_key(~valid_spk) = -1;

[unique_keys, ~, grp_idx] = unique(grp_key);
invalid_grp = find(unique_keys == -1);
n_groups    = length(unique_keys);

% Time weight + edge taper (computed once)
peak_hw = max(3, floor(T / 8));
w_time  = 0.3 * ones(1, T);
w_time(max(1, midPoint - peak_hw) : min(T, midPoint + peak_hw)) = 1.0;

taper_samples = max(3, floor(T / 10));
edge_taper = ones(1, T);
edge_taper(1:taper_samples) = 0.5 * (1 - cos(pi * (0:taper_samples-1) / taper_samples));
edge_taper(end-taper_samples+1:end) = 0.5 * (1 + cos(pi * (1:taper_samples) / taper_samples));

% Cache prepareTemplateInfo per source_ch
unique_source_chs = unique(all_spk_source_ch(valid_spk));
tmplInfo_cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
for ci = 1:length(unique_source_chs)
    sch = unique_source_chs(ci);
    tmplInfo_cache(sch) = prepareTemplateInfo(sortedSamples{sch, 1});
end

% Per-group: ONE template each (K x T arrays, not cell-of-cells)
grp_tmpl       = cell(n_groups, 1);   % K x T
grp_tw         = cell(n_groups, 1);   % K x T
grp_W          = cell(n_groups, 1);   % K x T
grp_bp_rows    = cell(n_groups, 1);   % K x 1
grp_denom      = zeros(n_groups, 1);
grp_tmpl_eng   = zeros(n_groups, 1);
grp_ok         = false(n_groups, 1);

for g = 1:n_groups
    if unique_keys(g) == -1, continue; end

    first_in  = find(grp_idx == g, 1);
    source_ch = all_spk_source_ch(first_in);
    label     = all_spk_labels(first_in);
    variant   = all_spk_variants(first_in);

    templateInfo_g = tmplInfo_cache(source_ch);
    all_templates  = templateInfo_g.templates;          % [nClasses, nTPC, C_wf, T_tmpl]
    classLabels_g  = templateInfo_g.classLabels(:);

    classIdx = find(classLabels_g == label, 1);
    if isempty(classIdx), continue; end

    [~, nTPC, C_wf, T_tmpl] = size(all_templates);

    % Clamp variant index to valid range
    v = min(variant, nTPC);
    if v < 1, v = 1; end

    % Extract the SPECIFIC variant (not looping over all)
    tmpl_v = squeeze(all_templates(classIdx, v, :, :));
    if isvector(tmpl_v)
        if C_wf > 1, tmpl_v = reshape(tmpl_v, C_wf, []);
        else,        tmpl_v = tmpl_v(:)'; end
    end
    if all(tmpl_v(:) == 0), continue; end

    % Centre-align to spikeDuration_pts
    [~, Tv] = size(tmpl_v);
    if Tv > T
        off = floor((Tv - T) / 2);
        tmpl_v = tmpl_v(:, off+1 : off+T);
    elseif Tv < T
        pad = zeros(C_wf, T);
        off = floor((T - Tv) / 2);
        pad(:, off+1 : off+Tv) = tmpl_v;
        tmpl_v = pad;
    end

    % Cosine edge taper
    tmpl_v = tmpl_v .* edge_taper;

    % Channel mapping
    ch_map = zeros(C_wf, 1);
    for c = 1:C_wf
        local_idx = (source_ch + c - middle_channel) - wf_range_start + 1;
        if local_idx >= 1 && local_idx <= n_wf_ch
            ch_map(c) = local_idx;
        end
    end
    valid_ch = ch_map > 0;
    if ~any(valid_ch), continue; end
    bp_rows = ch_map(valid_ch);
    K       = length(bp_rows);
    tv      = tmpl_v(valid_ch, :);       % K x T

    % Channel weighting
    ch_peak = max(abs(tv), [], 2);
    mid_ch_local  = min(middle_channel, C_wf);
    valid_indices = find(valid_ch);
    mid_in_valid  = find(valid_indices == mid_ch_local, 1);
    if ~isempty(mid_in_valid)
        ch_peak(mid_in_valid) = ch_peak(mid_in_valid) * 2;
    end
    ch_w = ch_peak / max(ch_peak(:) + eps);
    W    = ch_w .* w_time;
    tw   = W .* tv;
    den  = sum(tw(:) .* tv(:));
    if den < eps, continue; end

    grp_tmpl{g}     = tv;
    grp_tw{g}       = tw;
    grp_W{g}        = W;
    grp_bp_rows{g}  = bp_rows;
    grp_denom(g)    = den;
    grp_tmpl_eng(g) = sum(tv(:).^2);
    grp_ok(g)       = true;
end

% Mark spikes with no valid template
for s = 1:N_total
    if valid_spk(s) && ~grp_ok(grp_idx(s))
        valid_spk(s) = false;
    end
end
valid_idx = find(valid_spk);
n_valid   = length(valid_idx);
if n_valid == 0, return; end

% =========================================================================
% PHASE 3b: Vectorised alpha — ONE template per group, no variant loop
% =========================================================================
alphas = zeros(N_total, 1);

active_groups = setdiff(find(grp_ok), invalid_grp);
for gi = 1:length(active_groups)
    g = active_groups(gi);
    mask = valid_spk & (grp_idx == g);
    if ~any(mask), continue; end

    grp_times  = all_spk_times(mask);
    N_g        = sum(mask);
    bp_rows    = grp_bp_rows{g};
    tv         = grp_tmpl{g};
    tw         = grp_tw{g};
    W          = grp_W{g};
    den        = grp_denom(g);
    K          = length(bp_rows);
    grp_indices = find(mask);

    % Vectorised segment extraction
    t_s     = grp_times - midPoint + 1;
    idx_mat = t_s + (0:T-1);                        % N_g x T
    col_3d  = reshape(idx_mat, N_g, 1, T);
    row_3d  = reshape(bp_rows, 1, K, 1);
    lin_3d  = row_3d + (col_3d - 1) * n_wf_ch;
    segments = bp_matrix(lin_3d);                    % N_g x K x T

    % Vectorised alpha (single template per group)
    tw_3d = reshape(tw, 1, K, T);
    W_3d  = reshape(W,  1, K, T);
    num   = sum(sum(segments .* tw_3d, 3), 2);
    ga    = num / den;

    % Quality gate
    predicted = ga .* reshape(tv, 1, K, T);
    residual  = segments - predicted;
    res_e     = sum(sum(residual.^2 .* W_3d, 3), 2);
    obs_e     = sum(sum(segments.^2  .* W_3d, 3), 2);

    bad = (ga < 0.2) | (obs_e < eps) | (res_e >= 0.80 * obs_e);
    ga(bad) = 0;
    good = ga > 0;
    ga(good) = min(ga(good), 2.5);

    alphas(grp_indices) = ga;
end

% =========================================================================
% PHASE 3c: Subtract in amplitude order
% =========================================================================
sub_priority = abs(alphas(valid_idx)) .* grp_tmpl_eng(grp_idx(valid_idx));
[~, sub_order] = sort(sub_priority, 'descend');
n_subtracted = 0;

for oi = 1:n_valid
    si = sub_order(oi);
    s  = valid_idx(si);
    a  = alphas(s);
    if a == 0, continue; end

    g   = grp_idx(s);
    t_s = t_starts_all(s);
    t_e = t_ends_all(s);

    bp_matrix(grp_bp_rows{g}, t_s:t_e) = ...
        bp_matrix(grp_bp_rows{g}, t_s:t_e) - a * grp_tmpl{g};
    n_subtracted = n_subtracted + 1;
end

if n_subtracted == 0, return; end

% =========================================================================
% PHASE 4: Re-detect on residual (original threshold)
% =========================================================================
if ~isKey(out_filtered_channels, iCh), return; end
main_wf_idx = iCh - wf_range_start + 1;
if main_wf_idx < 1 || main_wf_idx > n_wf_ch, return; end

out_filtered_residual = out_filtered_channels(iCh);
original_bp = out_filtered_residual.bandpass_signal;

residual_bp = bp_matrix(main_wf_idx, :);
if iscolumn(original_bp)
    out_filtered_residual.bandpass_signal = residual_bp(:);
else
    out_filtered_residual.bandpass_signal = residual_bp;
end

out_detected_residual = kiaSort_detect_spike(out_filtered_residual, cfg, 0);
if isempty(out_detected_residual.spk_idx), return; end

% PHASE 5: drop residual re-detections within the two jitter zones --
% 0.3 * spkDist around any raw round-1 detection, 0.6 * spkDist around
% saved round-1 spike positions.
n_r = numel(out_detected_residual.spk_idx);
is_new = true(n_r, 1);

tight_tol = round(0.3 * spkDist);
if isKey(spk_idx_channels, iCh)
    raw_first = spk_idx_channels(iCh);
    if ~isempty(raw_first)
        d_raw = nearest_index(out_detected_residual.spk_idx(:), raw_first(:));
        is_new = is_new & (d_raw > tight_tol);
    end
end

loose_tol = round(0.6 * spkDist);
if isKey(output_spk_idx, iCh)
    saved_first = output_spk_idx(iCh);
    if ~isempty(saved_first)
        d_saved = nearest_index(out_detected_residual.spk_idx(:), saved_first(:));
        is_new = is_new & (d_saved > loose_tol);
    end
end

if ~any(is_new), return; end
out_detected_residual.spk_idx = out_detected_residual.spk_idx(is_new);
out_detected_residual.spk_Val = out_detected_residual.spk_Val(is_new);
out_detected_residual.spk_ID  = out_detected_residual.spk_ID(is_new);

if numel(out_detected_residual.spk_idx) < 2, return; end

% =========================================================================
% PHASE 6: Waveform extraction on residual
% =========================================================================
inputWF_r.main_channel_idx     = iCh;
inputWF_r.all_channel_idx      = wf_range(:);
inputWF_r.Ns_seq               = out_filtered_channels(iCh).Ns_seq;
inputWF_r.bandpass_data_window = cell(n_wf_ch, 1);
inputWF_r.spk_idx_all_channels = cell(n_wf_ch, 1);
inputWF_r.spk_ID_all_channels  = cell(n_wf_ch, 1);
inputWF_r.inclusion            = zeros(n_wf_ch, 1);
inputWF_r.channel_thresholds   = cell(n_wf_ch, 1);

is_col = iscolumn(original_bp);
for idx = 1:n_wf_ch
    mapped_ch = wf_range(idx);
    if is_col
        inputWF_r.bandpass_data_window{idx, 1} = bp_matrix(idx, :)';
    else
        inputWF_r.bandpass_data_window{idx, 1} = bp_matrix(idx, :);
    end
    if mapped_ch == iCh
        inputWF_r.spk_idx_all_channels{idx, 1} = out_detected_residual.spk_idx;
        inputWF_r.spk_ID_all_channels{idx, 1}  = out_detected_residual.spk_ID;
    elseif isKey(spk_idx_channels, mapped_ch)
        inputWF_r.spk_idx_all_channels{idx, 1} = spk_idx_channels(mapped_ch);
        inputWF_r.spk_ID_all_channels{idx, 1}  = [];
    else
        inputWF_r.spk_idx_all_channels{idx, 1} = [];
        inputWF_r.spk_ID_all_channels{idx, 1}  = [];
    end
    inputWF_r.inclusion(idx, 1) = channel_inclusion(mapped_ch);
    if ~isempty(channel_thresholds{mapped_ch, 1})
        inputWF_r.channel_thresholds{idx, 1} = channel_thresholds{mapped_ch, 1};
    else
        inputWF_r.channel_thresholds{idx, 1} = struct('bandpass_min_threshold', 0);
    end
end

[waveOut_r] = kiaSort_waveform_extraction(cfg, inputWF_r, 1, 0);
if isempty(waveOut_r.waveform) || size(waveOut_r.waveform, 1) == 0
    return;
end

% =========================================================================
% PHASE 7: Best-channel check (same strictness as first round)
% =========================================================================
[out_bc] = kiaSort_best_channel_detection(waveOut_r, 100, cfg);
if ~any(out_bc.keep), return; end

keep = out_bc.keep;
waveOut_r.waveform     = waveOut_r.waveform(keep, :, :);
waveOut_r.waveformInfo = sortedSamples{iCh, 1}.waveformInfo;

sorting_spk_idx_r = out_detected_residual.spk_idx(keep);
amplitude_r       = out_detected_residual.spk_Val(keep) .* ...
                    out_detected_residual.spk_ID(keep);
altChannels_r     = out_bc.max_altChannel(keep);

% =========================================================================
% PHASE 8: Classification of recovered spikes
% =========================================================================
if isKey(tmplInfo_cache, iCh)
    templateInfo = tmplInfo_cache(iCh);
else
    templateInfo = prepareTemplateInfo(sortedSamples{iCh, 1});
end

if useTemplate
    ds.waveform     = waveOut_r.waveform;
    ds.templateInfo = templateInfo;
    ds.amplitude    = amplitude_r;
    [pL, ft, vK, recovered_varIdx] = kiaSort_template_matching(ds, cfg);
else
    ds_c = struct();
    ds_c.waveform     = waveOut_r.waveform;
    ds_c.waveformInfo = sortedSamples{iCh, 1}.waveformInfo;
    ds_c = kiaSort_preprocess_waveforms(ds_c, cfg);
    ds_c.classifierInfo = sortedSamples{iCh, 1}.classifierInfo;
    ds_c.PCA            = sortedSamples{iCh, 1}.clusteringInfo.PCA;
    ds_c.amplitude      = amplitude_r;
    [pL, ft, vK] = kiaSort_predict_spike_labels(ds_c, cfg);
    recovered_varIdx = ones(size(pL));
end

if iscategorical(pL), pL = str2double(cellstr(pL)); end

% =========================================================================
% PHASE 8b: Vectorised multichannel correlation validation
%           Groups recovered spikes by predicted label, then batch-computes
%           weighted multichannel correlation using matrix operations.
% =========================================================================
N_recovered = numel(pL);
corrValid   = true(N_recovered, 1);

classLabels_ch = templateInfo.classLabels(:);
tmpl_data = [];
if isfield(templateInfo, 'templates')
    tmpl_data = templateInfo.templates;
end

if ~isempty(tmpl_data) && N_recovered > 0
    n_wf_ch_tmpl = size(waveOut_r.waveform, 2);
    T_full_wf    = size(waveOut_r.waveform, 3);
    cl_s = max(1, cl_start);
    cl_e = min(T_full_wf, cl_end);
    cl_len = cl_e - cl_s + 1;
    mid_ch_use = ceil(n_wf_ch_tmpl / 2);

    [~, ~, n_tmpl_ch_raw, T_tmpl_raw] = size(tmpl_data);
    n_ch_use = min(n_tmpl_ch_raw, n_wf_ch_tmpl);
    L_t = min(T_tmpl_raw, cl_len);

    % Group by predicted label (typically 2–10 unique labels)
    unique_pL = unique(pL(pL > 0));

    for li = 1:length(unique_pL)
        label_i = unique_pL(li);
        spike_mask = (pL == label_i);
        spike_indices = find(spike_mask);
        N_label = length(spike_indices);
        if N_label == 0, continue; end

        origIdx = find(classLabels_ch == label_i, 1);
        if isempty(origIdx)
            corrValid(spike_mask) = false;
            continue;
        end

        % Use highest-energy variant for correlation reference
        [~, nTPC_val, ~, ~] = size(tmpl_data);
        best_v = 1;
        best_e = -inf;
        for v = 1:nTPC_val
            tv_c = squeeze(tmpl_data(origIdx, v, :, :));
            e = sum(tv_c(:).^2);
            if e > best_e, best_e = e; best_v = v; end
        end

        tmpl_all_ch = squeeze(tmpl_data(origIdx, best_v, 1:n_ch_use, 1:L_t));
        if n_ch_use == 1, tmpl_all_ch = tmpl_all_ch(:)'; end

        % Channel weights (template peak amplitude, mid boosted 3x)
        ch_w = max(abs(tmpl_all_ch), [], 2);
        ch_w(mid_ch_use) = ch_w(mid_ch_use) * 3;
        total_w = sum(ch_w);
        if total_w < eps
            corrValid(spike_mask) = false;
            continue;
        end

        % Template stats (once per label)
        tmpl_means = mean(tmpl_all_ch(:, 1:L_t), 2);
        tmpl_centered = tmpl_all_ch(:, 1:L_t) - tmpl_means;
        tmpl_norms = sqrt(sum(tmpl_centered.^2, 2));
        tmpl_norms(tmpl_norms < eps) = eps;

        % Extract all spikes of this label: N_label x n_ch_use x cl_len
        spk_batch = waveOut_r.waveform(spike_indices, 1:n_ch_use, cl_s:cl_e);

        % Vectorised weighted multichannel correlation (channel-at-a-time)
        weighted_corr = zeros(N_label, 1);
        for c = 1:n_ch_use
            spk_c = squeeze(spk_batch(:, c, 1:L_t));
            if N_label == 1, spk_c = spk_c(:)'; end

            spk_means = mean(spk_c, 2);
            spk_centered = spk_c - spk_means;
            spk_norms = sqrt(sum(spk_centered.^2, 2));
            spk_norms(spk_norms < eps) = eps;

            dots = spk_centered * tmpl_centered(c, :)';
            ch_corr = dots ./ (tmpl_norms(c) * spk_norms);
            weighted_corr = weighted_corr + ch_w(c) * ch_corr;
        end
        weighted_corr = weighted_corr / total_w;

        corrValid(spike_indices(weighted_corr < 0.4)) = false;
    end
end

vK = vK & corrValid;

% =========================================================================
% PHASE 9: Label evaluation + realignment
% =========================================================================
clusterSelection  = sortedSamples{iCh, 1}.clusteringInfo.clusterSelection;
clusterRelabeling = sortedSamples{iCh, 1}.clusteringInfo.clusterRelabeling;

[uL, rIdx, rWf, ~] = realignSpikes(pL, waveOut_r.waveform, ...
                      sorting_spk_idx_r, clusterRelabeling, cfg);

rOut = kiaSort_evaluate_labels(uL, clusterSelection, altChannels_r, vK);
if ~any(rOut.keep), return; end

% =========================================================================
% PHASE 10: Label-aware deduplication
% =========================================================================
kept_mask = rOut.keep;

% (a) Same-channel
first_round_spk_ch = output_spk_idx(iCh);
first_round_lbl_ch = output_updatedLabels(iCh);

if ~isempty(first_round_spk_ch) && any(kept_mask)
    kept_indices = find(kept_mask);
    new_times    = rIdx(kept_indices);
    new_labels   = uL(kept_indices);

    [d_near, idx_near] = nearest_index(new_times(:), first_round_spk_ch(:));

    for ki = 1:length(kept_indices)
        if d_near(ki) <= spkDist && ...
                first_round_lbl_ch(idx_near(ki)) == new_labels(ki)
            kept_mask(kept_indices(ki)) = false;
        end
    end
end

% (b) Adjacent channels (tight jitter)
adj_jitter = round(0.4 * spkDist);
adj_channels = [iCh - 1, iCh + 1];
adj_channels = adj_channels(adj_channels >= 1 & adj_channels <= num_channels);

if any(kept_mask)
    kept_indices = find(kept_mask);
    new_times    = rIdx(kept_indices);
    new_labels   = uL(kept_indices);

    for ac = adj_channels
        if ~isKey(output_spk_idx, ac), continue; end
        adj_spk = output_spk_idx(ac);
        adj_lbl = output_updatedLabels(ac);
        if isempty(adj_spk), continue; end

        [d_adj, idx_adj] = nearest_index(new_times(:), adj_spk(:));

        for ki = 1:length(kept_indices)
            if ~kept_mask(kept_indices(ki)), continue; end
            if d_adj(ki) <= adj_jitter && ...
                    adj_lbl(idx_adj(ki)) == new_labels(ki)
                kept_mask(kept_indices(ki)) = false;
            end
        end
    end
end

if ~any(kept_mask), return; end

% =========================================================================
% PHASE 11: Append validated recovered spikes (including their variant idx)
% =========================================================================
k = kept_mask;
output_labels(iCh)        = [output_labels(iCh);        pL(k)];
output_updatedLabels(iCh) = [output_updatedLabels(iCh); uL(k)];
output_features(iCh)      = [output_features(iCh);      ft(k, :)];
output_spk_idx(iCh)       = [output_spk_idx(iCh);       rIdx(k)];
output_waveform(iCh)      = cat(1, output_waveform(iCh), rWf(k, :, :));
output_amplitude(iCh)     = [output_amplitude(iCh);      amplitude_r(k)];
output_bestVariant(iCh)   = [output_bestVariant(iCh);    recovered_varIdx(k)];

end