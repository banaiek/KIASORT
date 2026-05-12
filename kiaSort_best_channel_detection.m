function out = kiaSort_best_channel_detection(data, p, cfg, sample)

if nargin < 4
    sample = 0;
end

jitter_gap = ceil(cfg.spikeDistance * cfg.samplingFrequency / 1000);
threshold_pos = data.channel_thresholds_pos(:);
threshold_neg = data.channel_thresholds_neg(:);
waveforms = data.waveform;
channelIdx = data.wavformChanelIdx;
[n_clusters, n_channels, n_timepoints] = size(waveforms);
if n_clusters < 1
    out.keep = [];
    out.max_val = [];
    out.max_channel = [];
    out.max_channelIdx = [];
    out.max_altChannel = [];
    if sample
        out.rank = [];
    end
    return;
end
midpoint = ceil(n_timepoints/2);
window_start = max(1, midpoint - jitter_gap);
window_end = min(n_timepoints, midpoint + jitter_gap);
middle_channel = ceil(n_channels/2);
%
if sample
mask = gaussianMask(n_channels, n_timepoints, jitter_gap + 1, 3);
mMask(1,:,:) = mask;
waveforms = waveforms .* mMask;
end
waveform_window = waveforms(:, :, window_start:window_end);

mid_windowed_point = (window_end-window_start)/2+1;

T = size(waveform_window, 3);
if T >= 3 & sample
    d = diff(waveform_window, 1, 3);
    prodDiff = d(:,:,1:end-1) .* d(:,:,2:end);
    mask = false(size(waveform_window));
    mask(:,:,2:T-1) =  prodDiff < 0;
    waveform_window(~mask) = 0;
end
 



if sample
    max_w = max(waveform_window, [], 3);
    min_w = min(waveform_window, [], 3);
    pos_deviation = max_w - threshold_pos';
    neg_deviation = abs(min_w) - threshold_neg';
    max_detected = max(abs(waveform_window), [], 3);
    midChanDev = max_detected(:, middle_channel);
    [maxK_val, maxK_idx] = maxk(max_detected, 2, 2);

    detectability_ratio = max(pos_deviation ./ threshold_pos', neg_deviation ./ threshold_neg');
    out.detectblity_val = max(detectability_ratio, [], 2);
    out.mainNegativePolarity = (neg_deviation(:, middle_channel) - pos_deviation(:, middle_channel) > 0) & ...
        (neg_deviation(:, middle_channel) > 0);
    out.sideNegativePolarity = any((neg_deviation - pos_deviation > 0) & (neg_deviation > 0), 2);
    p_keep = prctile(max_detected, p, 2);
    keep = midChanDev >= p_keep;

    % check if main channel amplitude is low and other channels are even lower
    main_threshold = max(threshold_pos(middle_channel), threshold_neg(middle_channel));
    main_amp_low = midChanDev < 1.25 * main_threshold;

    % Check if other channels are less than 80% of main channel
    other_channels = setdiff(1:n_channels, middle_channel);
    other_amps_low = all(max_detected(:, other_channels) < 1 * midChanDev, 2);

    lowAmpNotKept = main_amp_low & other_amps_low;
    out.lowAmpNotKept = lowAmpNotKept;

    % Set keep to 0 for clusters meeting both conditions
    keep(lowAmpNotKept) = 0;

    
    if p > 99
        [~, sortedIdx] = sort(max_detected, 2, 'descend');
        rank_mid = zeros(n_clusters, 1);
        for i = 1:n_clusters
            rank_mid(i) = find(sortedIdx(i, :) == middle_channel, 1);
        end
        out.rank = rank_mid;
    end
else
    [max_val, max_idx] = max(abs(waveform_window(:,:,mid_windowed_point)),[], 2);
    max_detected = max(abs(waveform_window), [], 3);
    [maxK_val, maxK_idx] = maxk(max_detected, 2, 2);
    maxK_val(:,1) = max_val;
    maxK_idx(:,1) = max_idx;
    keep = middle_channel == maxK_idx(:,1);
end
out.keep = keep;
out.max_val = maxK_val(:,1);
out.max_channel = maxK_idx(:,1);
out.max_channelIdx = channelIdx(maxK_idx(:,1));
out.max_altChannel = channelIdx(maxK_idx(:,2));
out.max_altChannel((maxK_val(:,2)-threshold_pos(maxK_idx(:,2)))<0) = 0;
end

function mask = gaussianMask(n, m, jitterGap, sigma)
mask = zeros(n, m);
centerRow = ceil(n / 2);
centerCol = ceil(m / 2);
for i = 1:n
    d = abs(i - centerRow);
    halfWidth = ceil(jitterGap * exp(-(d ^ 2) / (2 * sigma ^ 2)));
    if halfWidth > 0
        colStart = max(1, centerCol - halfWidth) ;
        colEnd = min(m, centerCol + halfWidth) ;
        mask(i, colStart:colEnd) = 1;
    end
end
end