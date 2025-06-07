function Xd = removeSharedNoise_fast2(X, ...
                                      noise_pct, ... % 75
                                      corr_th,   ... % 0.40
                                      amp_th,    ... % 0.05   (5 % of RMS)
                                      nIter,     ... % 2
                                      useGPU)        % false / true
tic
% Two‑pass percentile template with amplitude gate.
% No spike mask, no PC step — maximum speed, minimum artefacts.
%
% X : N×M double OR single.  If useGPU==true, X can be gpuArray.

% ---------- defaults ----------
if nargin<2||isempty(noise_pct), noise_pct = 75;  end
if nargin<3||isempty(corr_th),   corr_th   = 0.40;end
if nargin<4||isempty(amp_th),    amp_th    = 0.05;end
if nargin<5||isempty(nIter),     nIter     = 1;   end
if nargin<6||isempty(useGPU),    useGPU    = false;end

% ---------- move to GPU if requested ----------
if useGPU && ~isa(X,'gpuArray'), X = gpuArray(X); end

Xd = X;                                       % working copy
for it = 1:nIter
    % ----- shared‑noise template -----
    g  = prctile(Xd, noise_pct, 1);           % 1×M
    gc = g - mean(g);                         % centre once
    ng = sqrt(sum(gc.^2)) + eps;

    % ----- per‑channel centred copy -----
    Xc = Xd - mean(Xd,2);

    % ----- correlation -----
    r = (Xc*gc')./(sqrt(sum(Xc.^2,2))*ng);    % N×1

    % ----- regression gain -----
    beta = (Xd*g') / (g*g');                  % N×1

    % ----- amplitude gate -----
    sig_g = std(g);
    sig_x = std(Xd,0,2) + eps;                % N×1
    amp   = abs(beta).*sig_g ./ sig_x;        % N×1

    % ----- subtract only where both tests succeed -----
    keep          = (r > corr_th) & (amp > amp_th);
    Xd(keep,:)    = Xd(keep,:) - beta(keep).*g;
end

% ---------- back to CPU if needed ----------
if useGPU, Xd = gather(Xd); end
toc
end