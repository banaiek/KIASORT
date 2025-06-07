function [bestCorr, bestLag] = max_half_corr(waveA, waveB, N, M2, maxLag, renorm)
% Find the max normalized correlation and corresponding lag for half-centered waveforms.
%
% Inputs:
%   waveA, waveB: 1 x (N*M2) flattened signals.
%   N: number of channels.
%   M2: total length.
%   maxLag: maximum lag (typically floor(M2/4)).
%
% Outputs:
%   bestCorr: maximum correlation.
%   bestLag: lag at which bestCorr occurs.

    bestCorr = -Inf;
    bestLag  = 0;

    waveA2D = reshape(waveA, [N, M2]);
    waveB2D = reshape(waveB, [N, M2]);

    if renorm
        z_A = waveA2D == 0;
        z_B = waveB2D == 0;
        waveA2D(z_A) = 1;
        waveB2D(z_B) = 1;
        waveA2D = waveA2D./max(abs(waveA2D),[],2);
        waveB2D = waveB2D./max(abs(waveB2D),[],2);
        waveA2D(z_A) = 0;
        waveB2D(z_B) = 0;
    end

    halfWidth = floor(M2/4);
    aStart = floor((M2/2) + (-halfWidth + 1));
    aEnd   = floor((M2/2) + halfWidth);
    snippetLen = aEnd - aStart + 1;
    if snippetLen <= 0
        bestCorr = 0;
        bestLag  = 0;
        return
    end

    fixedA = waveA2D(:, aStart:aEnd);
    fixedA = fixedA(:)';  

    for L = -maxLag : maxLag
        bStart = floor((M2/2) + (-halfWidth + 1 + L));
        bEnd   = floor((M2/2) + (halfWidth + L));
        if bStart < 1 || bEnd > M2
            continue
        end

        moveB = waveB2D(:, bStart:bEnd);
        moveB = moveB(:)';

        validMask = (fixedA ~= 0) & (moveB ~= 0);
        if ~any(validMask)
            continue
        end

        subA = fixedA(validMask);
        subB = moveB(validMask);
        c = corr(subA', subB', 'Type', 'Pearson');

        if c > bestCorr
            bestCorr = c;
            bestLag  = -L;
        end
    end

    if isinf(bestCorr)
        bestCorr = 0;
        bestLag  = 0;
    end
end