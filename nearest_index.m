function [d, idx] = nearest_index(A, B)
    
    %   [d, idx] = nearest_distance_and_index(A, B) returns two vectors:
    %   d   - Minimum absolute distance to the nearest element in B for each element in A
    %   idx - Index in B corresponding to the nearest element for each element in A

     A = A(:);
    B = B(:);
    
    % Sort B and track original indices
    [B_sorted, sort_idx_B] = sort(B, 'ascend');
    
    % Sort A and track original indices
    [A_sorted, sort_idx_A] = sort(A, 'ascend');
    [~, original_idx_A] = sort(sort_idx_A, 'ascend');
    
    % Define bin edges
    edges = [-inf; B_sorted; inf];
    
    % Assign each A to a bin
    idx_bin = discretize(A_sorted, edges);
    
    % Initialize distances and indices
    d1_sorted = Inf(size(A_sorted));
    d2_sorted = Inf(size(A_sorted));
    idx1_sorted = NaN(size(A_sorted));
    idx2_sorted = NaN(size(A_sorted));
    
    % Nearest after
    valid_d1 = idx_bin <= length(B_sorted);
    d1_sorted(valid_d1) = B_sorted(idx_bin(valid_d1)) - A_sorted(valid_d1);
    idx1_sorted(valid_d1) = sort_idx_B(idx_bin(valid_d1));
    
    % Nearest before
    valid_d2 = idx_bin > 1;
    d2_sorted(valid_d2) = A_sorted(valid_d2) - B_sorted(idx_bin(valid_d2) - 1);
    idx2_sorted(valid_d2) = sort_idx_B(idx_bin(valid_d2) - 1);
    
    % Compute absolute distances
    abs_d1 = abs(d1_sorted);
    abs_d2 = abs(d2_sorted);
    
    % Determine minimum distance and corresponding index
    [d_sorted, min_idx] = min([abs_d1, abs_d2], [], 2);
    idx_sorted = NaN(size(A_sorted));
    idx_sorted(min_idx == 1) = idx1_sorted(min_idx == 1);
    idx_sorted(min_idx == 2) = idx2_sorted(min_idx == 2);
    
    % Reorder to original A order
    d = d_sorted(original_idx_A);
    idx = idx_sorted(original_idx_A);
end