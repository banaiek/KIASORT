function [d1, d2] = nearest_distances_nz(A, B, jitter_gap)
    % nearest_distances: Finds the distances to the nearest points in B that are
    % beyond a specified jitter_gap relative to each element in A.
    A = A(:);
    B = sort(B(:));
    [A_sorted, sort_idx] = sort(A);
    [~, original_idx] = sort(sort_idx);


    A_right = A_sorted + jitter_gap;
    A_left  = A_sorted - jitter_gap;

    edges = [-inf; B; inf];

    idx_right = discretize(A_right, edges);
    d1_sorted = Inf(size(A_sorted));
    valid_right = idx_right <= length(B);

    d1_sorted(valid_right) = B(idx_right(valid_right)) - A_sorted(valid_right);


    idx_left = discretize(A_left, edges);
    d2_sorted = Inf(size(A_sorted));
    valid_left = idx_left > 1;
    d2_sorted(valid_left) = A_sorted(valid_left) - B(idx_left(valid_left) - 1);

    d1 = d1_sorted(original_idx);
    d2 = d2_sorted(original_idx);
end