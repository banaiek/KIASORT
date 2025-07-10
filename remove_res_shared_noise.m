function Xd = remove_res_shared_noise(X, cfg)
% shared-noise removal with signal preservation



p.useGPU = isfield(cfg,'useGPU') && cfg.useGPU && gpuDeviceCount>0;
p.maskThr = 3; % Base MAD multiplier
p.nPC = 5; % # PCs to evaluate
p.ds = 1; % Acceleration factor
p.varTol = 0.2; % Skip if low shared power
p.corrThr = 0.2; % Min correlation for shared noise
p.freqCutoff = 2000; % Hz, preserve above this
p.adaptiveMask = false; % Adaptive masking
p.iterative = true; % Iterative refinement
p.maxIter = 1; % Max iterations
% --------------------------------

if p.useGPU && ~isa(X,'gpuArray'), X = gpuArray(X); end


[n, M] = size(X);
Fs = cfg.samplingFrequency; % Sampling rate (Hz)

Xd = X;
mu = mean(X,2);
Xc = X - mu;

% frequency-based signal protection
if p.freqCutoff > 0 && Fs > 0
    % low-pass filter 
    [b, a] = butter(4, p.freqCutoff/(Fs/2), 'low');
    Xhigh = filtfilt(b, a, Xc')';
    Xlow = Xc - Xhigh;
else
    Xlow = Xc;
    Xhigh = zeros(size(Xc));
end

for iter = 1:p.maxIter
    
    mad = median(abs(Xlow - median(Xlow,2)),2) + eps;
    
    channelNoise = std(Xlow,0,2);
    channelSpike = std(Xhigh,0,2);
    snr = channelSpike ./ (channelNoise + eps);
    
    if p.adaptiveMask
        adaptThr = p.maskThr * (1 + 0.5 * tanh(2 - snr));
        mask = abs(Xlow) <= adaptThr .* mad;
    else
        mask = abs(Xlow) <= p.maskThr * mad;
    end
    
    if p.ds > 1
        ds_idx = 1:p.ds:M;
        mask_ds = mask(:,ds_idx);
        Xm = Xlow(:,ds_idx) .* mask_ds;
        Xlow_ds = Xlow(:,ds_idx);
    else
        mask_ds = mask;
        Xm = Xlow .* mask;
        Xlow_ds = Xlow;
    end
    
    W = sum(mask_ds,2);
    C = (Xm*Xm') ./ max(W,1);
    
    %  PC selection 
    try
        [U,S] = eigs(C,min(p.nPC,n-1),'largestreal','Tolerance',1e-3);
        S = diag(S); 
    catch
        [U,S,~] = svd(C);
        U = U(:,1:min(p.nPC,size(U,2)));
        S = diag(S);
        S = S(1:min(p.nPC,length(S)));
    end
    
    nPC_actual = min(size(U,2), length(S));
    U = U(:,1:nPC_actual);
    S = S(1:nPC_actual);
    
    projections = U' * Xlow_ds;
    validPC = false(nPC_actual,1);
    
    for pc = 1:nPC_actual
        pcWeights = abs(U(:,pc));
        
        if p.ds > 1
            % Upsample PC projection for correlation calculation
            pcProj = interp1(1:size(projections,2), projections(pc,:), ...
                            linspace(1,size(projections,2),M), 'linear');
        else
            pcProj = projections(pc,:);
        end
        
        % Count channels significantly correlated with this PC
        channelCorrs = zeros(n,1);
        for ch = 1:n
            if pcWeights(ch) > 0.1 % Only check channels with significant weight
                channelCorrs(ch) = abs(corr(Xlow(ch,:)', pcProj'));
            end
        end
        
        % PC is valid if correlated with multiple channels
        validPC(pc) = sum(channelCorrs > p.corrThr) > max(2, 0.3*n);
    end
    
    %  selective projection removal
    if any(validPC)
        Uvalid = U(:,validPC);
        
        weights = 1 ./ (1 + snr); 
        weights = weights / mean(weights);
        
        if p.ds > 1
            % Project and upsample
            proj_ds = Uvalid * (Uvalid' * Xlow_ds);
            proj = zeros(n, M);
            for ch = 1:n
                proj(ch,:) = interp1(1:size(proj_ds,2), proj_ds(ch,:), ...
                                   linspace(1,size(proj_ds,2),M), 'linear');
            end
        else
            proj = Uvalid * (Uvalid' * Xlow);
        end
        
        % Apply channel-specific weighting
        for ch = 1:n
            proj(ch,:) = proj(ch,:) * weights(ch);
        end
        
        Xlow = Xlow - proj;
    end
    
    % Check if we should continue iterating
    if ~p.iterative || iter >= p.maxIter
        break;
    end
    
    % Check convergence
    if any(validPC)
        sharedVar = sum(S(validPC)) / sum(var(Xlow,0,2));
    else
        sharedVar = 0;
    end
    if sharedVar < p.varTol
        break;
    end
end

Xd = Xlow + Xhigh + mu;


% Check for introduced artifacts
preVar = var(X,0,2);
postVar = var(Xd,0,2);
varRatio = postVar ./ preVar;


artifactChannels = varRatio > 1.5;
if any(artifactChannels)
    % Revert channels with artifacts
    Xd(artifactChannels,:) = X(artifactChannels,:);
    warning('Reverted %d channels due to artifact detection', sum(artifactChannels));
end

if p.useGPU, Xd = gather(Xd); end
end