function thr = kiaSort_refine_1d_2means(v, thr0)
% Deterministic 1-D two-means (Lloyd) seeded from the Otsu cut. No RNG, so
% it cannot perturb any other random draw in the pipeline.
thr = thr0;
v   = v(:);
for it = 1:50
    a = v(v <= thr);  b = v(v > thr);
    if isempty(a) || isempty(b), return; end
    c1 = mean(a);     c2 = mean(b);
    newThr = (c1 + c2) / 2;
    if ~isfinite(newThr) || abs(newThr - thr) < eps(max(abs(thr),1))
        thr = newThr;
        return;
    end
    thr = newThr;
end
end
