function isRefr = refractoryMatrix(spikeIdx, labels, fs, unitIDs)


BIN_MS        = 1;      % Bin width (ms)
WIN_MS        = 500;    % Half-width of correlogram (ms)
REF_MS        = 2;      % Central refractory half-width (ms)
ACG_RATIO_THR = 0.10;   % Refractory ratio threshold (ACG)
CCG_RATIO_THR = 0.25;   % Refractory ratio threshold (CCG)
ACG_P_THR     = 0.20;   % Poisson tail p-value threshold (ACG)
CCG_P_THR     = 0.05;   % Poisson tail p-value threshold (CCG)

binSamp  = round(BIN_MS * 1e-3 * fs);
winSamp  = round(WIN_MS * 1e-3 * fs);
refSamp  = round(REF_MS * 1e-3 * fs);
shoulderRange = [round(10e-3 * fs), round(50e-3 * fs)];  % 10-50 ms baseline

N = numel(unitIDs);
isRefr = false(N, N);


spkByUnit = cell(N, 1);
for iu = 1:N
    spkByUnit{iu} = sort(spikeIdx(labels == unitIDs(iu)));
end


for ii = 1:N
    s1 = spkByUnit{ii};
    if isempty(s1)
        continue;
    end
    
    for jj = ii:N
        s2 = spkByUnit{jj};
        if isempty(s2)
            continue;
        end
        
        % Compute coincidences using efficient sweep
        [centCnt, baseCnt, baseBins] = computeCoincidences(s1, s2, ...
            winSamp, binSamp, refSamp, shoulderRange);
        
        if isempty(baseBins) || baseBins == 0
            continue;
        end
        
        % Compute baseline rate (conservative: use maximum)
        R = max(baseCnt) / baseBins;
        if R == 0
            continue;
        end
        
        % Central bins to test
        nCentralBins = floor(refSamp / binSamp) + 1;
        kList = 0:(nCentralBins - 1);
        nk = centCnt(1:nCentralBins);
        
        % Expected counts under Poisson model
        E_nk = (2 * kList + 1) * R;
        E_nk = max(E_nk, eps);  % Avoid division by zero
        
        % Refractory ratio
        ratio = min(nk ./ E_nk);
        
        % Poisson cumulative probability
        pVals = poisscdf(nk, E_nk);
        qMin = min(pVals);
        
        % Apply appropriate threshold
        if ii == jj
            % Autocorrelogram test
            pass = (ratio < ACG_RATIO_THR) && (qMin < ACG_P_THR);
        else
            % Cross-correlogram test
            pass = (ratio < CCG_RATIO_THR) && (qMin < CCG_P_THR);
        end
        
        isRefr(ii, jj) = pass;
        isRefr(jj, ii) = pass;
    end
end

end


function [centCnt, baseCnt, baseBins] = computeCoincidences(a, b, winS, binS, refS, shRange)


if isempty(a) || isempty(b)
    centCnt = 0;
    baseCnt = [];
    baseBins = [];
    return;
end

if numel(a) > numel(b)
    [a, b] = deal(b, a);
end

nCentralBins = floor(refS / binS) + 1;
centCnt = zeros(1, nCentralBins);
leftCnt = 0;
rightCnt = 0;

ja = 1;
nb = numel(b);

for ia = 1:numel(a)
    tA = a(ia);
    
    while ja <= nb && b(ja) < tA - winS
        ja = ja + 1;
    end
    
    jb = ja;
    while jb <= nb && b(jb) <= tA + winS
        dt = b(jb) - tA;
        absDt = abs(dt);
        
        if absDt <= refS
            % Central region
            k = floor(absDt / binS) + 1;
            if k <= nCentralBins
                centCnt(k) = centCnt(k) + 1;
            end
        elseif dt < -shRange(1) && dt >= -shRange(2)
            % Left shoulder
            leftCnt = leftCnt + 1;
        elseif dt > shRange(1) && dt <= shRange(2)
            % Right shoulder
            rightCnt = rightCnt + 1;
        end
        
        jb = jb + 1;
    end
end

baseCnt = [leftCnt, rightCnt];
baseBins = (shRange(2) - shRange(1)) / binS + 1;

end