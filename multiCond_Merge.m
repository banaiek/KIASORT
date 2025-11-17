function merge_matrix = multiCond_Merge(A, options)
    % Robust comparison with temporal alignment (±2 samples)
    % A: N×Nch×T matrix (N groups, Nch channels, T time points)
    
    if nargin < 2
        options.corr_weight = 0.5;
        options.shape_weight = 0.5;
        options.threshold = 0.9;
        options.max_shift = 15;  % Check shifts from -2 to +2
    end
    
    N = size(A, 1);
    Nch = size(A, 2);
    T = size(A, 3);
    shifts = -options.max_shift:options.max_shift;  % [-2, -1, 0, 1, 2]
    
    % Initialize similarity matrices
    max_corr_sim = zeros(N, N);
    max_shape_sim = zeros(N, N);
    
    % Compare all pairs with temporal shifts
    for i = 1:N
        for j = i:N
            % Initialize temporary similarity scores for all shifts
            corr_scores = zeros(length(shifts), 1);
            shape_scores = zeros(length(shifts), 1);
            
            % Try all circular shifts
            for s_idx = 1:length(shifts)
                shift = shifts(s_idx);
                
                % Shift waveform j in time dimension
                A_j_shifted = circshift(A(j, :, :), shift, 3);
                
                % 1. Correlation similarity
                A_i_flat = reshape(A(i, :, :), 1, []);
                A_j_flat = reshape(A_j_shifted, 1, []);
                
                % Normalize
                A_i_norm = (A_i_flat - mean(A_i_flat)) / std(A_i_flat);
                A_j_norm = (A_j_flat - mean(A_j_flat)) / std(A_j_flat);
                
                % Correlation
                corr_scores(s_idx) = (A_i_norm * A_j_norm') / (Nch * T);
                
                % 2. Shape similarity (peak-to-peak ratios)
                peaks_i = squeeze(max(A(i,:,:), [], 3) - min(A(i,:,:), [], 3));
                peaks_j = squeeze(max(A_j_shifted, [], 3) - min(A_j_shifted, [], 3));
                
                if norm(peaks_i) + norm(peaks_j) > 0
                    shape_scores(s_idx) = 1 - norm(peaks_i - peaks_j) / (norm(peaks_i) + norm(peaks_j));
                else
                    shape_scores(s_idx) = 1;  % Both are zero
                end
            end
            
            % Take maximum similarity across all shifts
            max_corr_sim(i, j) = max(corr_scores);
            max_corr_sim(j, i) = max_corr_sim(i, j);
            
            max_shape_sim(i, j) = max(shape_scores);
            max_shape_sim(j, i) = max_shape_sim(i, j);
        end
    end
    
    % Combine metrics using weights
    combined_sim = options.corr_weight * max_corr_sim + ...
                   options.shape_weight * max_shape_sim;
    
    % Generate merge matrix
    merge_matrix = combined_sim > options.threshold;
    merge_matrix(logical(eye(N))) = 0;  % Don't merge with self
end