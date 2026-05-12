function [contRate, R, Q, baseRate] = estimateContamination(spk1, spk2, fs)
% Poisson-based refractory period contamination estimation.
%
% Computes how well a spike train (or pair of trains) respects the neural
% refractory period. Inspired by Kilosort 4 (Pachitariu et al., 2024).
%
% For a single train (ACG): estimates the contamination rate — the
% fraction of spikes likely originating from other neurons. A clean,
% well-isolated unit produces contRate near 0.
%
% For two trains (CCG): tests whether the pair exhibits a refractory dip,
% indicating they are oversplit fragments of the same neuron.
%
% Usage:
%   [contRate, R, Q, baseRate] = estimateContamination(spk, [], fs)   % ACG
%   [contRate, R, Q, baseRate] = estimateContamination(spk1, spk2, fs) % CCG
%
% Inputs:
%   spk1 - spike indices (samples) for unit 1
%   spk2 - spike indices for unit 2 (empty → autocorrelogram)
%   fs   - sampling frequency (Hz)
%
% Outputs:
%   contRate - contamination rate: min(observed/expected) across refractory
%              bins. 0 = perfect refractory, 1 = no refractory period.
%   R        - ratio vector per refractory bin (observed / expected)
%   Q        - Poisson CDF probability per refractory bin
%   baseRate - estimated baseline coincidence rate (spikes per bin)
%
% The refractory window is 0–2 ms, baseline is estimated from 10–50 ms
% shoulders, and a 1 ms bin width is used.

BIN_MS  = 1;     % bin width (ms)
WIN_MS  = 50;    % half-window for baseline (ms) — no need for full 500 ms
REF_MS  = 2;     % refractory half-width (ms)

binSamp  = round(BIN_MS  * 1e-3 * fs);
winSamp  = round(WIN_MS  * 1e-3 * fs);
refSamp  = round(REF_MS  * 1e-3 * fs);
shoulderLo = round(10e-3 * fs);   % 10 ms
shoulderHi = round(50e-3 * fs);   % 50 ms

isAuto = isempty(spk2);
if isAuto
    spk2 = spk1;
end

spk1 = sort(spk1(:));
spk2 = sort(spk2(:));

if numel(spk1) < 10 || numel(spk2) < 10
    contRate = 0;
    R = 0;
    Q = 0;
    baseRate = 0;
    return;
end

% ---- Count coincidences using efficient two-pointer sweep ----------------
nCentralBins = floor(refSamp / binSamp) + 1;
centCnt  = zeros(1, nCentralBins);
leftCnt  = 0;
rightCnt = 0;

% Always sweep shorter array against longer for efficiency
if numel(spk1) > numel(spk2)
    [spk1, spk2] = deal(spk2, spk1);
end

ja = 1;
nb = numel(spk2);

for ia = 1:numel(spk1)
    tA = spk1(ia);

    % Advance left pointer
    while ja <= nb && spk2(ja) < tA - winSamp
        ja = ja + 1;
    end

    jb = ja;
    while jb <= nb && spk2(jb) <= tA + winSamp
        dt = spk2(jb) - tA;
        absDt = abs(dt);

        % Skip zero-lag for autocorrelogram (same spike)
        if isAuto && dt == 0
            jb = jb + 1;
            continue;
        end

        if absDt <= refSamp
            k = floor(absDt / binSamp) + 1;
            if k <= nCentralBins
                centCnt(k) = centCnt(k) + 1;
            end
        elseif absDt >= shoulderLo && absDt <= shoulderHi
            if dt < 0
                leftCnt = leftCnt + 1;
            else
                rightCnt = rightCnt + 1;
            end
        end

        jb = jb + 1;
    end
end

% ---- Estimate baseline rate from shoulders -------------------------------
shoulderBins = (shoulderHi - shoulderLo) / binSamp + 1;
if shoulderBins <= 0
    contRate = 0; R = 0; Q = 0; baseRate = 0;
    return;
end

% Use the max of left/right shoulders (conservative — higher baseline
% means a harder test to pass as "clean")
baseRate = max(leftCnt, rightCnt) / shoulderBins;

if baseRate == 0
    % No baseline activity → cannot estimate contamination
    contRate = 0;
    R = zeros(1, nCentralBins);
    Q = zeros(1, nCentralBins);
    return;
end

% ---- Poisson model per refractory bin ------------------------------------
kList = 0:(nCentralBins - 1);
E_nk  = (2 * kList + 1) * baseRate;       % expected count per bin
E_nk  = max(E_nk, eps);

% R: observed / expected ratio (contamination rate per bin)
R = centCnt(1:nCentralBins) ./ E_nk;

% Q: Poisson CDF — probability of seeing <= nk by chance
Q = poisscdf(centCnt(1:nCentralBins), E_nk);

% Overall contamination rate: minimum ratio across bins
% (the bin with the deepest dip defines how refractory the unit is)
contRate = min(R);
contRate = max(0, min(contRate, 1));  % clamp to [0, 1]

end
