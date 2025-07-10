function [d1, d2] = nearest_distances(A, B)
    %  distances to nearest points after and before each element in A from B
    %   d1 : Distance to nearest point after each element of A in B (B(j) >= A(i))
    %   d2 : Distance to nearest point before each element of A in B (B(j) <= A(i))

    A = A(:);
    B = sort(B(:));

    [A_sorted, sort_idx] = sort(A);
    [~, original_idx] = sort(sort_idx);

    edges = [-inf; B; inf];
    idx = discretize(A_sorted, edges);

    d1_sorted = Inf(size(A_sorted));
    d2_sorted = Inf(size(A_sorted));

    valid_d1 = idx <= length(B);
    d1_sorted(valid_d1) = B(idx(valid_d1)) - A_sorted(valid_d1);

    valid_d2 = idx > 1;
    d2_sorted(valid_d2) = A_sorted(valid_d2) - B(idx(valid_d2) - 1);

    d1 = d1_sorted(original_idx);
    d2 = d2_sorted(original_idx);
end