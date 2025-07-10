function out = best_channel_detection(data, p, cfg, sample)
if nargin < 4
    sample = 0;
end

cfg.bestChannelDetection = 'Amp';

if sample
    factor = 1;
elseif p > 50
    factor = .01;
else
    factor = 1;
end
jitter_gap = round(factor * cfg.spikeDistance * cfg.samplingFrequency / 1000);

threshold_pos = data.channel_thresholds_pos(:);
threshold_neg = data.channel_thresholds_neg(:);
waveforms     = data.waveform;
channelIdx    = data.wavformChanelIdx;

[n_clusters, n_channels, n_timepoints] = size(waveforms);
if n_clusters < 1
    out.keep = [];
    out.max_val = [];
    out.max_channel = [];
    out.max_channelIdx = [];
    if sample
        out.rank = [];
    end
    return;
end

midpoint = ceil(n_timepoints/2);
window_start = max(1, midpoint - jitter_gap);
window_end   = min(n_timepoints, midpoint + jitter_gap);
middle_channel = ceil(n_channels/2);

mask = gaussianMask(n_channels, n_timepoints, jitter_gap + 1, 3);
mMask(1,:,:) = mask;
waveforms = waveforms .* mMask;

waveform_window = waveforms(:, :, window_start:window_end);
midChannelVal   = waveforms(:, middle_channel, midpoint);

max_w = max(waveform_window, [], 3);
min_w = min(waveform_window, [], 3);

pos_deviation = max_w - threshold_pos';
neg_deviation = abs(min_w) - threshold_neg';

if strcmp(cfg.bestChannelDetection, 'Amp')
    max_detected = max(abs(waveform_window), [], 3);
    midChanDev   = max_detected(:, middle_channel);
    validMask    = (pos_deviation > 0) | (neg_deviation > 0);
else
    midChanDev   = min(abs(midChannelVal - threshold_pos(middle_channel)), ...
                       abs(midChannelVal + threshold_neg(middle_channel)));
    max_detected = max(pos_deviation, neg_deviation);
end

[max_val, max_channel] = max(max_detected, [], 2);

if sample
    detectability_ratio      = max(pos_deviation ./ threshold_pos', neg_deviation ./ threshold_neg');
    out.detectblity_val      = max(detectability_ratio, [], 2);
    out.mainNegativePolarity = (neg_deviation(:, middle_channel) - pos_deviation(:, middle_channel) > 0) & ...
                               (neg_deviation(:, middle_channel) > 0);
    out.sideNegativePolarity = any((neg_deviation - pos_deviation > 0) & (neg_deviation > 0), 2);
    p_keep = prctile(max_detected, p, 2);
    keep   = midChanDev >= p_keep;

    [~, sortedIdx] = sort(max_detected, 2, 'descend');
    rank_mid       = zeros(n_clusters, 1);
    for i = 1:n_clusters
        rank_mid(i) = find(sortedIdx(i, :) == middle_channel, 1);
    end
    out.rank = rank_mid;
else
    p_keep = prctile(max_detected(:, sum(max_detected, 1) > 0), p, 2);
    keep   = midChanDev >= p_keep;
end

out.keep            = keep;
out.max_val         = max_val;
out.max_channel     = max_channel;
out.max_channelIdx  = channelIdx(max_channel);
end

function mask = gaussianMask(n, m, jitterGap, sigma)
mask = zeros(n, m);
centerRow = ceil(n / 2);
centerCol = ceil(m / 2);
for i = 1:n
    d = abs(i - centerRow);
    halfWidth = ceil(jitterGap * exp(-(d ^ 2) / (2 * sigma ^ 2)));
    if halfWidth > 0
        colStart = max(1, centerCol - halfWidth) - 1;
        colEnd   = min(m, centerCol + halfWidth) + 1;
        mask(i, colStart:colEnd) = 1;
    end
end
end