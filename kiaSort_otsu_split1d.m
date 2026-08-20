function [thr, sep, valley, eta] = kiaSort_otsu_split1d(t, nb)
% eta is Otsu's separability: the between-class variance at the optimal cut
% over the total variance. It measures how well the best possible split
% distinguishes the two modes, so no valley depth has to be assumed.
thr = NaN; sep = 0; valley = 1; eta = 0;
tmin = min(t); tmax = max(t);
if ~isfinite(tmin) || ~isfinite(tmax) || (tmax - tmin) <= eps, return; end
edges = linspace(tmin, tmax, nb+1);
h     = histcounts(t, edges);
ctr   = (edges(1:end-1) + edges(2:end))/2;
p     = h / sum(h);
omega = cumsum(p);
muc   = cumsum(p .* ctr);
muT   = muc(end);
den   = omega .* (1 - omega);
sigmaB = zeros(size(den));
ok = den > eps;
sigmaB(ok) = (muT*omega(ok) - muc(ok)).^2 ./ den(ok);
[sigmaBmax, kb] = max(sigmaB);
thr = edges(kb+1);
varT = sum(p .* (ctr - muT).^2);
if varT > eps
    eta = sigmaBmax / varT;
end
a = t(t <= thr);  b = t(t > thr);
if numel(a) < 3 || numel(b) < 3, thr = NaN; return; end
sd = std(a) + std(b);
if sd <= eps, thr = NaN; return; end
sep = abs(mean(a) - mean(b)) / sd;
hs = movmean(h, 3);
[~, iL] = min(abs(ctr - mean(a)));
[~, iR] = min(abs(ctr - mean(b)));
if iL > iR, tmp = iL; iL = iR; iR = tmp; end
if iR - iL < 2, thr = NaN; return; end
valley = min(hs(iL:iR)) / max(eps, min(max(hs(1:iL)), max(hs(iR:end))));
end
