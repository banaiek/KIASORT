function isRefr = refractoryMatrix(spikeIdx, labels, fs, unitIDs)
%--------------------------------------------------------------------------
% isRefr = refractoryMatrix(spikeIdx, labels, fs, unitIDs)
%
% INPUTS
%   spikeIdx   : column vector of sample indices (int64 or double)
%   labels     : column vector of same length giving the unit id of each spike
%   fs         : sampling frequency in Hz
%   unitIDs    : vector of unique unit labels you wish to test
%
% OUTPUT
%   isRefr     : logical N×N matrix.  isRefr(i,i) is true when unit i’s ACG
%                is refractory,  isRefr(i,j) (i≠j) is true when the CCG of
%                units i and j is refractory.
%
% PERFORMANCE
%   * No explicit histograms – a two‑pointer sweep is used for each pair.
%   * Memory footprint is O(∑ni) – fits millions of spikes easily.
%
%   Kia Supreme 27‑IV‑2025
%--------------------------------------------------------------------------
% constants (thresholds are from Kilosort4 paper, but can be tuned)
BIN_MS        = 1;                              % bin width (ms)
WIN_MS        = 500;                            % half‑width of correlogram (ms)
REF_MS        = 2;                              % central refractory half‑width (ms)
ACG_RATIO_THR = 0.10;                           % nk / E[nk]   threshold  (ACG)
CCG_RATIO_THR = 0.25;                           % nk / E[nk]   threshold  (CCG)
ACG_P_THR     = 0.20;                           % Poisson tail threshold   (ACG)
CCG_P_THR     = 0.05;                           % Poisson tail threshold   (CCG)

% pre‑compute constants in samples
binSamp  = round(BIN_MS * 1e-3 * fs);
winSamp  = round(WIN_MS * 1e-3 * fs);
refSamp  = round(REF_MS * 1e-3 * fs);
shoulderRange = [round(10e-3*fs) round(50e-3*fs)];  % 10‑50 ms baseline window

N  = numel(unitIDs);
isRefr = false(N,N);

% build cell array of spike vectors per unit, sorted (needed for sweep)
spkByUnit = cell(N,1);
for iu = 1:N
    spkByUnit{iu} = sort(spikeIdx(labels == unitIDs(iu)));
end

% main loop over unique pairs, including self‑pairs
for ii = 1:N
    s1 = spkByUnit{ii};
    for jj = ii:N
        s2 = spkByUnit{jj};
        
        % sweep: count coincidences inside ±winSamp at multiples of binSamp
        [centCnt, baseCnt, baseBins] = corrCoincidences(s1, s2, ...
                            winSamp, binSamp, refSamp, shoulderRange);
        if isempty(baseBins), continue; end          % not enough data
        
        % coincidence *rates*
        R = max(baseCnt ./ baseBins);                % conservative
        if R == 0,   continue; end                   % avoid div‑by‑0
        
        kList = 0:refSamp/binSamp;                  % central bins to test
        nk    = centCnt;                             % counts already per‑bin
        E_nk  = (2*kList + 1) * R;                   % expected coincidences
        
        ratio = min(nk ./ E_nk);                    % R12 in the paper
        % exact Poisson tail for each k
        pVals = poisscdf(nk, E_nk);
        qMin  = min(pVals);
        
        if ii == jj
            % ACG test
            pass = (ratio < ACG_RATIO_THR) && (qMin < ACG_P_THR);
        else
            % CCG test
            pass = (ratio < CCG_RATIO_THR) && (qMin < CCG_P_THR);
        end
        
        isRefr(ii,jj) = pass;
        isRefr(jj,ii) = pass;   % symmetry
    end
end
end
%--------------------------------------------------------------------------

function [centCnt, baseCnt, baseBins] = corrCoincidences(a,b,winS,binS,refS,shRange)
% Two‑pointer O(na+nb) sweep that returns:
%   centCnt  : vector of counts per central bin (width binS) from -refS .. +refS
%   baseCnt  : 2‑element vector [leftShoulder rightShoulder] total counts
%   baseBins : corresponding number of bins in each shoulder
%
if isempty(a) || isempty(b)
    centCnt  = 0; baseCnt = []; baseBins = [];
    return
end

% ensure 'a' is the shorter for efficiency
if numel(a) > numel(b), [a,b] = deal(b,a); end

centCnt = zeros(1,refS/binS + 1);
leftCnt = 0; rightCnt = 0;

ja = 1;
nb = numel(b);
for ia = 1:numel(a)
    tA = a(ia);
    % advance lower index in b
    while ja<=nb && b(ja) < tA - winS, ja = ja + 1; end
    jb = ja;
    while jb<=nb && b(jb) <= tA + winS
        dt = b(jb) - tA;
        if abs(dt) <= refS
            k = floor(abs(dt)/binS) + 1;
            centCnt(k) = centCnt(k) + 1;
        elseif dt < -shRange(1) && dt >= -shRange(2)
            leftCnt = leftCnt + 1;
        elseif dt >  shRange(1) && dt <=  shRange(2)
            rightCnt = rightCnt + 1;
        end
        jb = jb + 1;
    end
end
baseCnt  = [leftCnt rightCnt];
baseBins = diff(shRange)/binS + 1;  % number of bins in each shoulder
end