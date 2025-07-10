function [counts, centers, violations] = getISIViolations(spike_idx, fs , threshold)    
    
    if nargin < 2
        fs = 1; 
    end

    if nargin < 3
        threshold = 2; 
    end
    
    spike_times = sort(spike_idx(:))/fs;
    isi = diff(spike_times);
    threshold = threshold/1000;

    edges = [0, threshold];
    max_isi = max(isi);
    if max_isi > threshold
        log_edges = logspace(log10(threshold), log10(max_isi), 50);
        edges = [edges, log_edges];
    end
    edges = unique(edges);
    
    counts = histcounts(isi, edges);
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    violations = 100*sum(isi < threshold)/numel(isi);

end