
function [bin_centers, bin_counts, presence_ratio] = presenceRatio(spike_idx, fs, trial_length, bin_size, smoothness_factor)

    spike_times = double(spike_idx) / fs;

    overlap = bin_size * smoothness_factor ;
    
    % Generate bin start times
    max_start = trial_length - bin_size;
    if max_start < 0
        bin_starts = [];
    else
        bin_starts = (0:bin_size:max_start)';
        bin_starts = bin_starts(bin_starts <= max_start);
    end
        
    if isempty(bin_starts)
        bin_centers = [];
        bin_counts = [];
        presence_ratio = 0;
        return;
    end
    
    bin_centers = bin_starts + bin_size/2;
    bin_counts = zeros(size(bin_centers));

    for i =1:numel(bin_centers)
    bin_counts(i) = sum(spike_times > bin_starts(i) - overlap & spike_times < bin_starts(i) + bin_size + overlap);
    end

    
    presence_ratio = sum(bin_counts > 0) / numel(bin_starts);
end