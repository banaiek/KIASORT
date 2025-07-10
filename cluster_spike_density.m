function [spike_density, time_bins, groups] = cluster_spike_density(spk_idx, group_labels, cfg)
% spike density per group using a Gaussian kernel.
%
% Inputs:
%   spk_idx      - Spike indices.
%   group_labels - Group IDs for each spike.
%   cfg          - samplingFrequency, sampleChunkDuration.
%
% Outputs:
%   spike_density - Matrix (groups x time bins).
%   time_bins     - bin centers
%   groups        - Unique group IDs

    samplingFrequency = cfg.samplingFrequency;
    sampleChunkDuration = cfg.sampleChunkDuration/10;
    total_time = cfg.numSampleChunks * cfg.sampleChunkDuration; 
    spk_idx = spk_idx(:);
    group_labels = group_labels(:);
    spike_times = spk_idx / samplingFrequency;
    
    groups = unique(group_labels);
    n_groups = numel(groups);
    
    % total_time = ceil(max(spike_times) / sampleChunkDuration) * sampleChunkDuration;
    bin_edges = 0:sampleChunkDuration:total_time;
    time_bins = bin_edges(1:end-1) + sampleChunkDuration/2;
    n_bins = numel(time_bins);
    
    smoothness = 3;
    spike_density = zeros(n_groups, n_bins);
    kernel_width = 6 * smoothness;
    num_kernel_pts = ceil(kernel_width / sampleChunkDuration) + mod(ceil(kernel_width / sampleChunkDuration),2); % odd number
    half_len = floor(num_kernel_pts / 2);
    t_kernel = (-half_len:half_len) * sampleChunkDuration;
    kernel = exp(-0.5 * (t_kernel / smoothness).^2);
    kernel = kernel / sum(kernel); 
    
    normalization = conv(ones(1, n_bins), kernel, 'same');
    
    for i = 1:n_groups
        group = groups(i);
        group_spike_times = spike_times(group_labels == group);
        if isempty(group_spike_times)
            continue;
        end
        counts = histcounts(group_spike_times, bin_edges);
        density_raw = conv(counts, kernel, 'same');
        density_corrected = density_raw ./ normalization;
        spike_density(i, :) = density_corrected / sampleChunkDuration;
    end
end