function [out] = kiaSort_waveform_extraction(cfg, data, cancel_overlap, sample)

% Extract waveforms with GPU acceleration for band pass signal
% overlapping spikes are canceled out through a sigmoid window multiplication

if nargin < 4
    sample = 0;
end

if nargin < 3
    cancel_overlap = 1;
end

main_channel_idx       = data.main_channel_idx;
all_channel_idx        = data.all_channel_idx;
bandpass_data_window   = data.bandpass_data_window;
spk_idx_all_channels   = data.spk_idx_all_channels;
inclusion              = data.inclusion;
all_channel_thresholds = data.channel_thresholds;
clear data

% setup gpu arrays if available
useGPU = false;
if isfield(cfg, 'useGPU') && cfg.useGPU && gpuDeviceCount > 0
    useGPU = true;
    main_channel_idx = gpuArray(main_channel_idx);
    all_channel_idx  = gpuArray(all_channel_idx);
    inclusion        = gpuArray(inclusion);
    for i = 1:numel(bandpass_data_window)
        bandpass_data_window{i} = gpuArray(bandpass_data_window{i});
    end
    for i = 1:numel(spk_idx_all_channels)
        spk_idx_all_channels{i} = gpuArray(spk_idx_all_channels{i});
    end
    for i = 1:numel(all_channel_thresholds)
        if isfield(all_channel_thresholds{i}, 'bandpass_min_threshold')
            all_channel_thresholds{i}.bandpass_min_threshold = gpuArray(all_channel_thresholds{i}.bandpass_min_threshold);
        end
        if isfield(all_channel_thresholds{i}, 'bandpass_max_threshold')
            all_channel_thresholds{i}.bandpass_max_threshold = gpuArray(all_channel_thresholds{i}.bandpass_max_threshold);
        end
    end
end

fs                  = cfg.samplingFrequency;
num_channel_extract = cfg.num_channel_extract;
spikeDuration       = round(cfg.spikeDuration * fs / 1000) + 1;
jitter_gap          = round(1 * round(cfg.spikeDistance * fs / 1000));
midPoint            = floor(spikeDuration / 2) + 1;
middle_channel      = num_channel_extract + 1;

% sigmoid window library
mask_lib = zeros(midPoint, spikeDuration);
for i = 1:midPoint
    mask_lib(i,:) = smooth_square(spikeDuration, 2*i, 2);
end
mask_lib = mask_lib./max(mask_lib,[],2);
left_mask_lib  = mask_lib;
right_mask_lib = mask_lib;
left_mask_lib(:, midPoint:end)  = 1;
right_mask_lib(:, 1:midPoint)     = 1;
if useGPU
    left_mask_lib = gpuArray(left_mask_lib);
    right_mask_lib= gpuArray(right_mask_lib);
end

main_channel_pos = find(all_channel_idx == main_channel_idx);
channel_positions = (main_channel_pos - num_channel_extract) : (main_channel_pos + num_channel_extract);
num_channels = length(channel_positions);
main_spk_inds = spk_idx_all_channels{main_channel_pos};
num_spikes = length(main_spk_inds);

if useGPU
    wavformChanelIdx       = gpuArray(nan(num_channels, 1));
    channel_thresholds_neg = gpuArray(nan(num_channels, 1));
    channel_thresholds_pos = gpuArray(nan(num_channels, 1));
    channel_inclusion      = gpuArray(zeros(num_channels, 1));
else
    wavformChanelIdx       = nan(num_channels, 1);
    channel_thresholds_neg = nan(num_channels, 1);
    channel_thresholds_pos = nan(num_channels, 1);
    channel_inclusion      = zeros(num_channels, 1);
end


if useGPU
    waveform_bp = gpuArray.zeros(num_spikes, num_channels, spikeDuration);
else
    waveform_bp = zeros(num_spikes, num_channels, spikeDuration);
end


if num_spikes > 0
    extraction_idx = main_spk_inds + (1:spikeDuration) - midPoint;
    for c = 1:num_channels
        channel_pos = channel_positions(c);
        if channel_pos < 1 || channel_pos > length(all_channel_idx)
            continue;
        else
            wavformChanelIdx(c,1) = all_channel_idx(channel_pos,1);
            if inclusion(channel_pos)
                data_length = length(bandpass_data_window{channel_pos});
                valid_idx = extraction_idx >= 1 & extraction_idx <= data_length;
                
                if all(valid_idx(:))
                    channel_bandpass_data = bandpass_data_window{channel_pos}(extraction_idx);
                    waveform_bp(:, c, :) = reshape(channel_bandpass_data, [num_spikes, 1, spikeDuration]);
                    channel_inclusion(c,1) = inclusion(channel_pos);
                    channel_thresholds_neg(c,1) = all_channel_thresholds{channel_pos,1}.bandpass_min_threshold;
                    if isfield(all_channel_thresholds{channel_pos,1}, 'bandpass_max_threshold')
                        channel_thresholds_pos(c,1) = all_channel_thresholds{channel_pos,1}.bandpass_max_threshold;
                    else
                        channel_thresholds_pos(c,1) = all_channel_thresholds{channel_pos,1}.bandpass_min_threshold;
                    end
                end
            end
        end
    end


if cancel_overlap 
    % sigmoid window to the main channel
    spike_distance = diff(main_spk_inds);
    overlap_id = find(spike_distance < 0.75 * spikeDuration & spike_distance > jitter_gap);
    overlap_length = round(spike_distance(overlap_id)/2);
    right_spike_overlap = overlap_id;
    left_spike_overlap  = overlap_id + 1;
    right_mask = reshape(right_mask_lib(overlap_length, :), [length(overlap_id), 1, spikeDuration]);
    left_mask  = reshape(left_mask_lib(overlap_length, :), [length(overlap_id), 1, spikeDuration]);

    waveform_bp(right_spike_overlap, middle_channel, :) = waveform_bp(right_spike_overlap, middle_channel, :) .* right_mask;
    waveform_bp(left_spike_overlap, middle_channel, :)  = waveform_bp(left_spike_overlap, middle_channel, :)  .* left_mask;

    % sigmoid window to other than the main channels
    for c = 1:num_channels
        channel_pos = channel_positions(c);
        if channel_pos < 1 || channel_pos > length(all_channel_idx)
            continue;
        elseif inclusion(channel_pos) && channel_pos ~= main_channel_pos
            spk_inds_channel = spk_idx_all_channels{channel_pos};
            if useGPU
                [d1, d2] = nearest_distances_nz(gather(main_spk_inds), gather(spk_inds_channel), gather(jitter_gap));
                d1 = gpuArray(d1);
                d2 = gpuArray(d2);
            else
                [d1, d2] = nearest_distances_nz(main_spk_inds, spk_inds_channel, jitter_gap);
            end
            right_overlap_idx = find(d1 < 0.75 * spikeDuration );
            left_overlap_idx  = find(d2 < 0.75 * spikeDuration );
            right_overlap_length = round(d1(right_overlap_idx)/2);
            left_overlap_length = round(d2(left_overlap_idx)/2);
            right_mask = reshape(right_mask_lib(right_overlap_length, :), [length(right_overlap_idx), 1, spikeDuration]);
            left_mask  = reshape(left_mask_lib(left_overlap_length, :), [length(left_overlap_idx), 1, spikeDuration]);
            if sample
                if useGPU
                    randRight = (.5-gpuArray.rand(size(right_mask))) .* (1-right_mask);
                    randLeft = (.5-gpuArray.rand(size(left_mask))) .* (1-left_mask);
                else
                    randRight = (.5-rand(size(right_mask))) .* (1-right_mask);
                    randLeft = (.5-rand(size(left_mask))) .* (1-left_mask);
                end
                waveform_bp(right_overlap_idx, c, :) = waveform_bp(right_overlap_idx, c, :) .* right_mask + (randRight.*abs(waveform_bp(right_overlap_idx,c , midPoint)/5));
                waveform_bp(left_overlap_idx, c, :)  = waveform_bp(left_overlap_idx,c , :) .* left_mask + (randLeft.*abs(waveform_bp(left_overlap_idx,c , midPoint)/5));
            else
                waveform_bp(right_overlap_idx, c, :) = waveform_bp(right_overlap_idx, c, :) .* right_mask ;
                waveform_bp(left_overlap_idx, c, :)  = waveform_bp(left_overlap_idx,c , :) .* left_mask ;
            end
        end
    end
end

end


out.waveform      = waveform_bp;
out.wavformChanelIdx = wavformChanelIdx;
out.channel_thresholds_neg = channel_thresholds_neg;
out.channel_thresholds_pos = channel_thresholds_pos;
out.channel_inclusion      = channel_inclusion;

if useGPU
    out = gather_out(out);
end

end

function out = gather_out(out)
fn = fieldnames(out);
for i = 1:length(fn)
    if ~isempty(out.(fn{i})) && isa(out.(fn{i}), 'gpuArray')
        out.(fn{i}) = gather(out.(fn{i}));
    end
end
end