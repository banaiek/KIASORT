function Xd = remove_shared_noise(X, cfg)

% Decorrelates channels via noise-floor-equalised C^{-1/2}, estimated
% from low-amplitude ("quiet") time-points for robustness to spikes.
%
% Unlike standard ZCA which maps all eigenvector dimensions to unit
% variance (amplifying quiet dimensions and indirectly attenuating
% small spikes via RMS rescaling), this version equalises to the
% median eigenvalue (the noise floor).  Shared-noise dimensions are
% suppressed, independent-noise dimensions pass through at scale ≈ 1,
% and nothing is amplified — so spikes of all sizes are preserved.
%
%  X   – [nChannels x nTimepoints]
%  cfg – struct, optional fields:
%         .useGPU        – push to GPU                          (false)
%         .regFactor     – regularisation, fraction of noise
%                          floor eigenvalue                     (0.1)
%         .nPCs          – PCs to keep, 0 = all                 (0)
%         .covBatch      – max samples for cov estimate         (1e5)
%         .scaleOutput   – rescale to input RMS                 (true)
%         .quietPrctile  – %-ile threshold for quiet selection  (75)
%
%  Xd  – [nChannels x nTimepoints], same scale/mean as input

if nargin < 2, cfg = struct(); end

% ---- GPU handling (matches original interface) ----------------------
useGPU = false;
if isfield(cfg,'useGPU') && cfg.useGPU ...
        && gpuDeviceCount > 0 && ~isa(X,'gpuArray')
    useGPU = true;
    X      = gpuArray(X);
end

[nChan, nSamp] = size(X);

% ---- Config with defaults -------------------------------------------
regFactor    = getOr(cfg, 'regFactor',    0.1);
nPCs         = getOr(cfg, 'nPCs',         0);
covBatch     = getOr(cfg, 'covBatch',     min(nSamp, 1e5));
scaleOutput  = getOr(cfg, 'scaleOutput',  true);
quietPrctile = getOr(cfg, 'quietPrctile', 5);

if nPCs <= 0 || nPCs > nChan, nPCs = nChan; end

% ---- 1. Subsample for covariance estimation -------------------------
if covBatch < nSamp
    idx  = round(linspace(1, nSamp, covBatch));
    Xsub = X(:, idx);
else
    Xsub = X;
end

mu   = mean(Xsub, 2);
Xsub = Xsub - mu;                        % centre

% ---- 2. Keep only quiet time-points (robust to spikes) --------------
chanPow  = sum(Xsub.^2, 1);              % instantaneous power
th       = prctile(chanPow, quietPrctile);
quietIdx = chanPow <= th;
if sum(quietIdx) > 2*nChan                % need enough samples
    Xsub = Xsub(:, quietIdx);
end

% ---- 3. Noise covariance  [nChan x nChan] --------------------------
nS = size(Xsub, 2);
C  = (Xsub * Xsub') / nS;
C  = (C + C') * 0.5;                     % enforce perfect symmetry

% ---- 4. Eigendecomposition ------------------------------------------
[V, D] = eig(C, 'vector');
[D, ix] = sort(D, 'descend');
V = V(:, ix);

V = V(:, 1:nPCs);
D = max(D(1:nPCs), 0);                   % clamp negatives

% ---- 5. Noise-floor-equalised whitening weights ---------------------
%  Standard ZCA:   invSqrt = 1 ./ sqrt(D + reg)
%    → maps every dimension to unit variance
%    → amplifies low-variance (quiet) dimensions
%    → inflates total RMS → rescaling shrinks everything → kills small spikes
%
%  Noise-floor equalisation:   invSqrt = sqrt(noiseFloor ./ (D + reg))
%    → shared-noise dims (D >> floor): suppressed           ✓
%    → independent-noise dims (D ≈ floor): scale ≈ 1        ✓ untouched
%    → quiet dims (D < floor): scale ≈ 1, capped by reg     ✓ not amplified
%    → spikes ride through at original amplitude

noiseFloor = median(D);
reg        = regFactor * noiseFloor;
invSqrt    = sqrt(noiseFloor ./ (D + reg));

% ---- 6. Apply  Xd = V diag(invSqrt) V' (X - mu) -------------------
mu_full = mean(X, 2);
Xc      = X - mu_full;

%  two thin matrix multiplies – never forms the full W matrix
Xd = V * (invSqrt .* (V' * Xc));         % O(nChan·nPCs·nSamp)

% ---- 7. Rescale to original RMS (drop-in compatible) ----------------
if scaleOutput
    rms_in  = sqrt(mean(Xc.^2, 2)) + eps;
    rms_out = sqrt(mean(Xd.^2, 2)) + eps;
    Xd      = Xd .* (rms_in ./ rms_out);
end

Xd = Xd + mu_full;                       % restore DC offset

if useGPU, Xd = gather(Xd); end
end

% ---------------------------------------------------------------------
function v = getOr(s, f, d)
    if isfield(s, f), v = s.(f); else, v = d; end
end