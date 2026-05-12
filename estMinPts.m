function pt = estMinPts(X)
% ESTMINPTS  Smallest dense cluster size via mutual k-NN graph.
    [N, d] = size(X);
    ns = min(N, 100000);
    S  = X(randperm(N, ns), :);
    k  = min(ceil(sqrt(d) + log(ns)), ns - 1);
    I  = knnsearch(S, S, 'K', k+1);  I = I(:, 2:end);
    r  = repmat((1:ns)', 1, k);
    A  = sparse(r(:), I(:), 1, ns, ns);
    A  = A & A';                                   % mutual edges only
    [~, ~, b] = dmperm(A | speye(ns));             % connected components
    sz = diff(b);  sz = sz(sz >= ceil(log(ns)));   % drop noise specks
    if isempty(sz), pt = max(ceil(log(N)), d+1); 
    else,           pt = max(round(min(sz)*N/ns), d+1); 
    end
end
