function out = overlapping_waveform(X, tol, counts, t, minLag, maxLag, spread, Z)
% Detect overlapping waveforms by checking if the residual (Ai - Bj)
% matches any shifted/spread template.
%
% OPTIMIZATION: Precompute ALL candidate waveforms (template × spread ×
% shift) into a matrix, then test each residual against all candidates in
% one matrix–vector multiply instead of a 5-deep nested loop.
%
% Test:  ‖R − C‖² / ‖R‖² < tol²
%  ⟹    (1−tol²)·‖R‖² − 2·R·Cᵀ + ‖C‖² < 0
%
% The dot products R·Cᵀ for all candidates are computed in a single call.

[n, ch, m] = size(X);
[r, ~,  ~] = size(Z);

globalFactor = max(max(abs(X(:))), max(abs(Z(:))));
if globalFactor == 0
    out = false(n, n);
    return
end
Xnorm = X / globalFactor;
Znorm = Z / globalFactor;

% ------------------------------------------------------------------
% 1. Precompute candidate matrix  (nCand × ch*m)
% ------------------------------------------------------------------
% Estimate maximum number of candidates to pre-allocate
maxCand = r * spread * (2*maxLag + 1);
candMat = zeros(maxCand, ch*m);
nCand   = 0;

for k = 1:r
    Ck = squeeze(Znorm(k,:,:));                    % ch × m
    for curSpread = 1:spread
        % Build spread version
        spreadCk = zeros(ch, m);
        for h = -curSpread:curSpread
            spreadCk = spreadCk + shiftSignal(Ck, h);
        end
        spreadCk = spreadCk / (2*curSpread + 1);
        % Apply every valid shift
        for s = -maxLag:maxLag
            if abs(s) < minLag
                continue
            end
            nCand = nCand + 1;
            candidate = shiftSignal(spreadCk, s);
            candMat(nCand,:) = candidate(:)';
        end
    end
end
candMat = candMat(1:nCand, :);                     % trim

candNormSq = sum(candMat.^2, 2);                   % nCand × 1
tol2 = tol^2;

% ------------------------------------------------------------------
% 2. Vectorised residual matching
% ------------------------------------------------------------------
out = false(n, n);

% Flatten cluster waveforms once
Xflat = zeros(n, ch*m);
for i = 1:n
    tmp = squeeze(Xnorm(i,:,:));
    Xflat(i,:) = tmp(:)';
end

for i = 1:n
    Ai = Xflat(i,:);
    for j = 1:n
        if i == j || counts(i) > t * counts(j)
            continue
        end
        R      = Ai - Xflat(j,:);
        Rnorm2 = R * R';                             % ‖R‖²
        if Rnorm2 < eps
            continue
        end

        % All dot products at once: nCand × 1
        dots = candMat * R';
        % Test: (1-tol²)·‖R‖² + ‖C‖² − 2·dot < 0
        test = (1 - tol2) * Rnorm2 + candNormSq - 2*dots;
        if any(test < 0)
            out(i,j) = true;
        end
    end
end
end

% ==================================================================
function Y = shiftSignal(X, s)
[ch, m] = size(X);
Y = zeros(ch, m);
if s > 0
    Y(:, s+1:end) = X(:, 1:end-s);
elseif s < 0
    s = abs(s);
    Y(:, 1:end-s) = X(:, s+1:end);
else
    Y = X;
end
end