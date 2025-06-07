function out = overlapping_waveform(X, tol, counts, t, minLag, maxLag, spread, Z)
[n,ch,m] = size(X);
[r,~,~] = size(Z);
globalFactor = max(max(abs(X(:))), max(abs(Z(:))));
if globalFactor == 0
    Xnorm = X;
    Znorm = Z;
else
    Xnorm = X / globalFactor;
    Znorm = Z / globalFactor;
end
out = false(n,n);
for i = 1:n
    Ai = squeeze(Xnorm(i,:,:));
    for j = 1:n
        if i == j || counts(i) > t * counts(j)
            continue;
        end
        Bj = squeeze(Xnorm(j,:,:));
        R = Ai - Bj;
        foundMatch = false;
        for k = 1:r
            Ck = squeeze(Znorm(k,:,:));
            for curSpread = 1:spread
                spreadCk = zeros(ch, m);
                for h = -curSpread:curSpread
                    spreadCk = spreadCk + shiftSignal(Ck, h);
                end
                spreadCk = spreadCk / (2*curSpread + 1);
                for s = -maxLag:maxLag
                    if abs(s) < minLag
                        continue;
                    end
                    candidate = shiftSignal(spreadCk, s);
                    if norm(R - candidate, 'fro') / norm(R, 'fro') < tol
                        out(i,j) = true;
                        foundMatch = true;
                        break;
                    end
                end
                if foundMatch
                    break;
                end
            end
            if foundMatch
                break;
            end
        end
    end
end
end

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