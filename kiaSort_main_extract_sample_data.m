function kiaSort_main_extract_sample_data(inputPath, outputPath, cfg, varargin )

progressFcn = [];

for i = 1:2:length(varargin)
    name = lower(varargin{i});
    value = varargin{i+1};
    switch name
        case 'channel_mapping'
            channel_mapping = value;
        case 'channel_inclusion'
            channel_inclusion = value;
            chan_wave_inclusion = value;
        case 'channel_locations'
            channel_locations = value;
        case 'progressfcn'
            if isa(value, 'function_handle')
                progressFcn = value;
            end
        otherwise
            warning('Unknown parameter: %s', varargin{i});
    end
end

% Main path
outputFolder = fullfile(outputPath, 'RES_Samples');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

logFile = fullfile(outputPath, 'KIASort_log.txt');
fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');
fprintf(fid, 'Sample extraction started at %s\n', datetime);

errorFlag = false; 

try
    % Load channel_map.mat
    if ~exist('channel_mapping', 'var')
        channelMapPath = fullfile(inputPath, 'channel_map.mat');
        if exist(channelMapPath, 'file')
            try
                [channel_mapping, channel_locations, channel_inclusion] = load_channel_map(channelMapPath, cfg);
                chan_wave_inclusion = channel_inclusion;
                fprintf(fid, 'Loaded channel_map.mat successfully.\n');
            catch ME
                fprintf(fid, 'ERROR: Failed to load channel_map.mat: %s\n', ME.message);
                if ~isempty(ME.stack)
                    fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
                end
                errorFlag = true;
                error('Failed to load channel_map.mat. Check the log file for details.');
            end
        else
            fprintf(fid, 'channel_map.mat not found in the specified inputPath.\n');
            fprintf(fid, 'Assuming all channels are in correct orders.\n');
            channel_mapping = [1 : cfg.numChannels]';
        end
    end

    if min(channel_mapping) == 0
        channel_mapping = channel_mapping + 1;
        fprintf(fid, 'Adjusted channel_mapping by adding 1 to zero-based indices.\n');
    end

    % Memory mapping
    m = map_input_file(inputPath, cfg);
    num_samples = size(m.Data.data,2);
    window_size = cfg.num_channel_extract * 2 + 1; % Number of channels in the window
    half_window = floor(window_size / 2);
    middleChannel_idx = half_window +1;
    num_channels = min(cfg.numChannels, length(channel_mapping));
    fs = cfg.samplingFrequency;

    num_sample_chunks   = cfg.numSampleChunks;
    chunk_duration      = cfg.sampleChunkDuration;
    maxSampleSpan       = cfg.maxSampleSpan * 60;
    total_duration      = num_samples / fs;
    sampleSpan          = min(maxSampleSpan, total_duration); % Total duration in seconds
    num_chunk_pts       = chunk_duration * fs; % Number of samples per chunk
    interval            = sampleSpan / num_sample_chunks; % Interval between each chunk in seconds
    interval_samples    = floor(interval * fs); % Interval between chunks in samples
    qualityCheckLength  = cfg.qualityCheckLength * fs;
    batch_ch_size       = min(cfg.batch_ch_size,num_channels);

    edge_length           = round(cfg.spikeDuration * fs / 2000);
    sample_spike_distance = round(cfg.spikeSampleDistance * fs / 1000); % in ms
    spk_Distance = round(cfg.spikeDistance * 3e4 / 1000);

    edge_exclusion = [];
    adj_distant = spk_Distance * ones(num_channels, 1);

    if ~exist("channel_inclusion",'var')
        channel_inclusion    = ones(num_channels, 1);
        chan_wave_inclusion  = ones(num_channels, 1);
    end

    if ~exist("channel_locations",'var')
        channel_locations      = zeros(num_channels, 3);
        channel_locations(:,2) = 1 : num_channels;
    end

    try
        qualityCheck_idx = randi(round(sampleSpan * fs/2),1) + round(sampleSpan * fs/4);
        qualityCheckSignal = double(m.Data.data(channel_mapping, qualityCheck_idx:qualityCheck_idx + qualityCheckLength-1));
        if isfield(cfg, 'commonRef')
            switch cfg.commonRef
                case 'median'
                    qualityCheckSignal = qualityCheckSignal - median(qualityCheckSignal,'omitmissing');
                    fprintf(fid, 'Samples median referenced.\n');
                case 'mean'
                    qualityCheckSignal = qualityCheckSignal - mean(qualityCheckSignal,'omitmissing');
                    fprintf(fid, 'Samples median referenced.\n');
                case 'none'
                    % skip
            end
        end
    catch ME
        fprintf(fid, 'ERROR: Failed to quality check the data: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
        errorFlag = true;
        error('Failed to quality check the data %d. Check the log file for details.', i);
    end

    [qualityOut] = kiaSort_signal_quality_check(qualityCheckSignal, cfg);
    scale_factor = qualityOut.scale_factor;
    thresh_MAD   = qualityOut.thresh_MAD;
    band_pairs = qualityOut.band_pairs;
    Ns_seq = qualityOut.Ns_seq;
    clear qualityCheckSignal

    % Generate random start indices for the chunks
    chunk_limits = zeros(num_sample_chunks,2);
    for i = 1:num_sample_chunks
        start_sample = round((i - 1) * interval_samples + randi(interval_samples - num_chunk_pts + 1));
        end_sample = start_sample + num_chunk_pts - 1;
        chunk_limits(i,:) = [start_sample, end_sample];
        edge_exclusion = [edge_exclusion, i*num_chunk_pts-edge_length:i*num_chunk_pts+edge_length];
    end

    % Generate indices for batches
    channel_idx = [1:batch_ch_size:(num_channels-1),num_channels];


    % Initialize maps and variables

    % Maps to store intermediate results for this chunk
    mapNames = {
        'bp_channels', 'spk_idx_channels', ...
        'spk_Val_channels', 'spk_ID_channels', ...
        'scale_factor_channels', 'rms_values_channels'
        };

    % Initialize each map as a local variable
    for i = 1:length(mapNames)
        eval([mapNames{i} ' = containers.Map(''KeyType'', ''double'', ''ValueType'', ''any'');']);
    end



    snr_status_channels     = cell(num_channels, 1);
    rms_bands_channels      = zeros(num_channels, cfg.numBands+1);
    channel_thresholds      = cell(num_channels, 1);
    mad_Thresh              = zeros(num_channels, 1);

    
    % Process each channel
    for ch_idx = 1:num_channels
        
        physical_ch_idx = channel_mapping(ch_idx); % Physical channel index mapping
        ch_indices_mapped = (max(1, ch_idx - half_window) : min(num_channels, ch_idx + half_window))';  % neighboring channels for windowed extraction

        last_batch_channel = min(ch_indices_mapped(end),num_channels);
        
        if ch_idx == 1
            batch_idx = 1;
            start_ch = channel_idx(batch_idx);
            end_ch = batch_idx * batch_ch_size;

            try
                selected_data = batch_extract(m, chunk_limits, start_ch, end_ch, channel_mapping, channel_inclusion, cfg);
                fprintf(fid, 'Batch %d successfully extracted! %s\n', i);
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
                fprintf(fid, 'Batch %d successfully extracted! %s\n', i);
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

        fprintf(fid, 'Processing channel %d/%d...\n', ch_idx, num_channels);
        try
            tic

            % detect spikes for channels in the window if not already done
            for idx = 1:length(ch_indices_mapped)
                mapped_ch_idx = ch_indices_mapped(idx);

                if ~isKey(bp_channels, mapped_ch_idx)
                    try
                        if ~channel_inclusion(mapped_ch_idx)

                           [waveOut, out_filtered, out_detected]   = setWaveformNull;

                            for j = 1:length(mapNames)
                                currentMap = eval(mapNames{j});
                                currentMap(mapped_ch_idx) = [];
                                eval([mapNames{j} ' = currentMap;']);
                            end

                            continue;
                        end
    
                        % Extract chunk data
                        batch_mapped_idx = mod(mapped_ch_idx,batch_ch_size);
                        if batch_mapped_idx == 0
                            batch_mapped_idx = batch_ch_size;
                        end
                        channel_data = double(selected_data(batch_mapped_idx, :));

                        % filter the data for the chunk
                        out_filtered = kiaSort_filter_signal(channel_data, cfg);
                        out_filtered.rms_bands = qualityOut.rms_bands(:, mapped_ch_idx);
                        out_filtered.scale_factor = qualityOut.scale_factor(mapped_ch_idx);
                        out_filtered.band_pairs = qualityOut.band_pairs;
                        out_filtered.Ns_seq = qualityOut.Ns_seq;

                        % Detect spikes
                        out_detected = kiaSort_detect_spike(out_filtered, cfg, 1);

                        % Store results in maps
                        spk_idx_channels(mapped_ch_idx)         = out_detected.spk_idx;
                        spk_Val_channels(mapped_ch_idx)         = out_detected.spk_Val;
                        spk_ID_channels(mapped_ch_idx)          = out_detected.spk_ID;
                        scale_factor_channels(mapped_ch_idx)    = out_detected.scale_factor;
                        rms_values_channels(mapped_ch_idx)      = out_detected.rms_values;
                        bp_channels(mapped_ch_idx)              = out_filtered.bandpass_signal;
                        snr_status_channels{mapped_ch_idx}      = out_detected.SNR;
                        rms_bands_channels(mapped_ch_idx, :)    = out_filtered.rms_bands;
                        channel_inclusion(mapped_ch_idx)        = out_detected.inclusion;
                        chan_wave_inclusion(mapped_ch_idx)      = out_detected.waveInclusion;
                        adj_distant(mapped_ch_idx)              = out_detected.adj_distant;
                        channel_thresholds{mapped_ch_idx, 1}    = out_detected.rms_values;
                        if isfield(out_detected.rms_values,"bandpass_min_threshold")
                            mad_Thresh(mapped_ch_idx)               = out_detected.rms_values.bandpass_min_threshold;
                        end
                    catch ME
                        fprintf(fid, 'ERROR: Failed to decompose/detect spikes for channel %d: %s\n', mapped_ch_idx, ME.message);
                        if ~isempty(ME.stack)
                            fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
                        end
                        errorFlag = true;
                        error('Failed to decompose/detect spikes for channel %d. Check the log file for details.', mapped_ch_idx);
                    end
                end
            end

            %  waveform extraction
            inputWF.bandpass_data_window         = cell(length(ch_indices_mapped), 1);
            inputWF.spk_idx_all_channels         = cell(length(ch_indices_mapped), 1);
            inputWF.inclusion                    = zeros(length(ch_indices_mapped),1);
            for idx = 1:length(ch_indices_mapped)
                inputWF.main_channel_idx         = ch_idx;
                mapped_ch_idx                    = ch_indices_mapped(idx);
                inputWF.all_channel_idx          = ch_indices_mapped;
                inputWF.bandpass_data_window{idx, 1}   = bp_channels(mapped_ch_idx);
                inputWF.spk_idx_all_channels{idx, 1}   = spk_idx_channels(mapped_ch_idx);
                inputWF.spk_ID_all_channels{idx, 1}    = spk_ID_channels(mapped_ch_idx);
                inputWF.inclusion(idx, 1)              = channel_inclusion(mapped_ch_idx);
                inputWF.channel_thresholds{idx, 1}     = channel_thresholds{mapped_ch_idx, 1};
            end

            current_channel_idx = find(ch_indices_mapped == ch_idx);
            inclusion_flag = any(inputWF.inclusion(current_channel_idx));
            if inclusion_flag
                try

                    waveOut = kiaSort_waveform_extraction(cfg, inputWF, 1, 1);
                catch ME
                    fprintf(fid, 'ERROR: Failed to extract waveforms for channel %d: %s\n', ch_idx, ME.message);
                    if ~isempty(ME.stack)
                        fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
                    end
                    errorFlag = true;
                    error('Failed to extract waveforms for channel %d. Check the log file for details.', ch_idx);
                end
                
                % extracting waveforms that are absent on the main channel
                side_ch = (floor(half_window/2));
                if side_ch
                    ch_m1 = find(inputWF.inclusion(side_ch+1:end), 1, 'first') + side_ch -1;
                    ch_m3 = find(inputWF.inclusion(1:end-side_ch), 1, 'last');
                    idx_m1 = inputWF.spk_idx_all_channels{ch_m1, 1};
                    idx_m3 = inputWF.spk_idx_all_channels{ch_m3, 1};
                    ch_start = find(inputWF.inclusion, 1, 'first');
                    ch_end = find(inputWF.inclusion, 1, 'last');
                    idx_start = inputWF.spk_idx_all_channels{ch_start, 1};
                    idx_end = inputWF.spk_idx_all_channels{ch_end, 1};
                    idx_middle = inputWF.spk_idx_all_channels{current_channel_idx, 1};
                    d_start = nearest_index(idx_start,[idx_end; idx_middle]);
                    d_end = nearest_index(idx_end,[idx_start; idx_middle]);
                    d_m1 = nearest_index(idx_m1,[idx_start; idx_middle; idx_end; idx_m3]);
                    d_m3 = nearest_index(idx_m3,[idx_start; idx_middle; idx_end; idx_m1]);
                    side_inds = [idx_start(d_start > sample_spike_distance); idx_end(d_end > sample_spike_distance); idx_m1(d_m1 > sample_spike_distance); idx_m3(d_m3 > sample_spike_distance) ];

                    if length(side_inds) > 1600
                        side_inds = datasample(side_inds,1600,'Replace',false);
                    end

                    if length(side_inds) > cfg.nPCAcomp
                        side_inputWF = inputWF;
                        side_inputWF.spk_idx_all_channels{current_channel_idx, 1} = side_inds;
                        side_waveOut = kiaSort_waveform_extraction(cfg, side_inputWF, 0, 1);
                    else
                        side_waveOut = [];
                    end
                else
                    side_waveOut = [];
                end
            else
                waveOut.waveform         = [];
                waveOut.wavformChanelIdx = [];
                waveOut.channel_inclusion       = [];
                waveOut.channel_thresholds_neg  = [];
                waveOut.channel_thresholds_pos  = [];
                side_waveOut = [];
            end


            % Filter valid spikes

            tmp_id  = spk_ID_channels(ch_idx);
            tmp_idx = spk_idx_channels(ch_idx);
            tmp_val = spk_Val_channels(ch_idx);

            
            
            valid_idx_edg  = ~ismember(tmp_idx, edge_exclusion);
            neighbor_Spk_Idx = cell2mat(inputWF.spk_idx_all_channels);
            [d1, d2] = nearest_distances_nz(tmp_idx, neighbor_Spk_Idx, 1 * spk_Distance);
            valid_idx      = valid_idx_edg & d1>sample_spike_distance & d2>sample_spike_distance ;
                        
            dist_temp       = tmp_idx;
            left_diff       = [Inf; diff(dist_temp)];
            right_diff      = [diff(dist_temp); Inf];
            
            valid_idx  = (left_diff > sample_spike_distance) & (right_diff > sample_spike_distance) & valid_idx;            

            out.spk_Val_full            = tmp_val(valid_idx);
            out.spk_ID_full             = tmp_id(valid_idx);
            out.spk_idx_full            = tmp_idx(valid_idx);
            out.wavformChanelIdx        = waveOut.wavformChanelIdx;
            out.channel_thresholds_neg  = waveOut.channel_thresholds_neg;
            out.channel_thresholds_pos  = waveOut.channel_thresholds_pos;
            out.channel_inclusion       = waveOut.channel_inclusion;
            out.scale_factor_full       = scale_factor_channels(ch_idx);
            out.rms_values_full         = rms_values_channels(ch_idx);
            out.Ns_seq                  = qualityOut.Ns_seq;
            out.side_waveforms          = side_waveOut;
            out.waveform_bp_full        = waveOut.waveform(valid_idx, :, :);
            waveOut.waveform        = waveOut.waveform(valid_idx, :, :);
            [out_bestChannel]     = best_channel_detection(waveOut, 5, cfg, 0);
            
            if sum(out_bestChannel.keep) > cfg.nPCAcomp
                out.spk_Val_full            = out.spk_Val_full(out_bestChannel.keep);
                out.spk_ID_full             = out.spk_ID_full(out_bestChannel.keep);
                out.spk_idx_full            = out.spk_idx_full(out_bestChannel.keep);
                out.wavformChanelIdx        = waveOut.wavformChanelIdx;
                out.channel_thresholds_neg  = waveOut.channel_thresholds_neg;
                out.channel_thresholds_pos  = waveOut.channel_thresholds_pos;
                out.channel_inclusion       = waveOut.channel_inclusion;
                out.scale_factor_full       = scale_factor_channels(ch_idx);
                out.rms_values_full         = rms_values_channels(ch_idx);
                out.Ns_seq                  = qualityOut.Ns_seq;
                out.waveform_bp_full        = out.waveform_bp_full(out_bestChannel.keep,:,:);
            end
           
            % Save extracted sample spikes for each channel
            channel_result_filename = fullfile(outputFolder, sprintf('channel_%d_results.mat', ch_idx));
            try
                save(channel_result_filename, 'out', 'ch_idx', '-v7.3');
                fprintf(fid, 'Saved results for channel %d to %s.\n', ch_idx, channel_result_filename);
            catch ME
                fprintf(fid, 'ERROR: Failed to save results for channel %d: %s\n', ch_idx, ME.message);
                if ~isempty(ME.stack)
                    fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
                end
                errorFlag = true;
                error('Failed to save results for channel %d. Check the log file for details.', ch_idx);
            end

            fprintf('Processed channel %d/%d (Physical Channel %d)\n', ch_idx, num_channels, physical_ch_idx);
            toc

            % remove channels that are out no longer needed
            try
                toRemoveChannels = setdiff(cell2mat(keys(bp_channels)), ch_indices_mapped);
                keys_to_remove = num2cell(toRemoveChannels);

                for j = 1:length(mapNames)
                    currentMap = eval(mapNames{j});
                    remove(currentMap, keys_to_remove);
                    eval([mapNames{j} ' = currentMap;']);
                end
            catch ME
                fprintf(fid, 'WARNING: Failed to remove obsolete channels: %s\n', ME.message);
                if ~isempty(ME.stack)
                    fprintf(fid, '         In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
                end
                % Continue processing even if removal fails
            end
        catch ME
            fprintf(fid, 'ERROR: Unexpected error processing channel %d: %s\n', ch_idx, ME.message);
            if ~isempty(ME.stack)
                fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
            end
            errorFlag = true;
            error('Unexpected error processing channel %d. Check the log file for details.', ch_idx);
        end

        if ~isempty(progressFcn)
            pct = ch_idx / num_channels;
            msg = sprintf('Extracting samples from channel %d of %d', ch_idx, num_channels);
            progressFcn(pct, msg);
        end
                
    end

    % Save channel info for main sorting
    try
        channelInfo_filename = fullfile(outputFolder, 'channel_info.mat');
        save(channelInfo_filename, 'channel_inclusion', 'num_samples', 'channel_locations', 'chan_wave_inclusion', 'channel_mapping',...
            'snr_status_channels', 'rms_bands_channels', 'scale_factor', 'adj_distant', 'thresh_MAD', 'mad_Thresh', 'Ns_seq','band_pairs', '-v7.3');
        fprintf(fid, 'Saved channel info to %s.\n', channelInfo_filename);
    catch ME
        fprintf(fid, 'ERROR: Failed to save channel info: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
        errorFlag = true;
        error('Failed to save channel info. Check the log file for details.');
    end

catch ME
    fprintf(fid, 'ERROR: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
    errorFlag = true;
    fclose(fid);
    error('An unexpected error occurred. Processing stopped. Check the log file for details.');
end

if ~errorFlag
    try
        fprintf(fid, 'All channels processed successfully.\n');
        fprintf(fid, 'Sample chunk processed. Spike detection and extraction completed for all channels at %s\n', datetime);
        fclose(fid);
        fprintf('All channels processed successfully. Summary written to %s.\n', logFile);
    catch ME
        fprintf('ERROR: Failed to write summary to log file: %s\n', ME.message);
    end
else
    try
        fprintf(fid, '\nSummary Report:\n');
        fprintf(fid, 'Processing stopped due to errors. Check the log file for details.\n');
        fclose(fid);
        error('Processing stopped due to errors. Check the log file: %s', logFile);
    catch ME
        fprintf('ERROR: Failed to write summary to log file: %s\n', ME.message);
        error('Processing stopped due to errors. Additionally, failed to write summary to log file.');
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

out.spk_idx      = [];
out.spk_Val      = [];
out.spk_ID       = [];
out.scale_factor = [];
out.rms_values   = [];
out.SNR          = [];
out.inclusion    = 0;
out.Ns_seq       = [];
end

function selected_data = batch_extract(m, chunk_limits, start_ch, end_ch, channel_mapping, channel_inclusion, cfg)

num_sample_chunks   = cfg.numSampleChunks;
chunk_duration      = cfg.sampleChunkDuration;
num_chunk_pts       = chunk_duration * cfg.samplingFrequency; % Number of samples per chunk

selected_data = zeros(end_ch-start_ch+1,num_chunk_pts*num_sample_chunks);
for i = 1:num_sample_chunks
    selected_data(:, (i-1)*num_chunk_pts + 1:i*num_chunk_pts) = double(m.Data.data(channel_mapping(start_ch:end_ch), chunk_limits(i,1):chunk_limits(i,2)));
end

if cfg.extremeNoise
    selected_data(channel_inclusion(start_ch:end_ch),:) = remove_res_shared_noise(selected_data(channel_inclusion(start_ch:end_ch),:), cfg);
end

if cfg.denoising
    selected_data(channel_inclusion(start_ch:end_ch),:) = remove_shared_noise(selected_data(channel_inclusion(start_ch:end_ch),:), cfg);
end



end



