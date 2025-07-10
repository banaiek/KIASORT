function Xd = remove_shared_noise(X, cfg)    

useGPU = false;
if isfield(cfg, 'useGPU') && cfg.useGPU && gpuDeviceCount > 0 && ~isa(X,'gpuArray')
    useGPU = true;
    X = gpuArray(X);
end

if isfield(cfg, 'noisePrctile')
    noise_pct = cfg.noisePrctile;
else
    noise_pct = 60;
end
if isfield(cfg, 'noiseCorr')
    corr_th = cfg.noiseCorr;
else
    corr_th = 0.4;
end

amp_th    = 0.05;
nIter     = 1;

Xd = X;           
for it = 1:nIter  
    g  = prctile(Xd, noise_pct, 1);           
    gc = g - mean(g);                         
    ng = sqrt(sum(gc.^2)) + eps;   
    Xc = Xd - mean(Xd,2);    
    r = (Xc*gc')./(sqrt(sum(Xc.^2,2))*ng);   
    beta = (Xd*g') / (g*g');                      
    sig_g = std(g);
    sig_x = std(Xd,0,2) + eps;                
    amp   = abs(beta).*sig_g ./ sig_x;            
    keep          = (r > corr_th) & (amp > amp_th);
    Xd(keep,:)    = Xd(keep,:) - beta(keep).*g;
end

if useGPU, Xd = gather(Xd); end

end