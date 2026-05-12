function [predLabels, matchScores, validKeep, bestVariantIdx] = kiaSort_template_matching(data, cfg)
% Match spikes to templates. bestVariantIdx is the matched variant index
% per spike, reused by performSecondRound to avoid re-searching variants.

if isfield(cfg, 'templateMatchMethod')
    method = cfg.templateMatchMethod;
else
    method = 'weighted_euclidean';
end

templates      = data.templateInfo.templates;
classLabels    = data.templateInfo.classLabels(:);
polarity       = data.templateInfo.polarity(:);
low_thr        = data.templateInfo.lowAmpThr(:);
high_thr       = data.templateInfo.highAmpThr(:);
ampVal         = data.amplitude(:);

[nClasses, nTemplatesPerClass, C, T] = size(templates);

validClassMask = classLabels ~= -1;
validClassIdx = find(validClassMask);
nValidClasses = length(validClassIdx);

if nValidClasses == 0
    error('No valid classes found (all labels are -1).');
end

fs           = cfg.samplingFrequency;
midPoint     = floor(cfg.spikeDuration * fs/(2*1000)) + 1;
spike_length = floor(cfg.clusteringSpikeDuration * fs/(2*1000));

waveform_full = data.waveform;
waveform_full(isnan(waveform_full)) = 0;
spikes = waveform_full(:, :, midPoint-spike_length:midPoint+spike_length);

[Nspikes, ~, T_spk] = size(spikes);

S = reshape(spikes, Nspikes, C*T_spk);

nTotalTemplates = nValidClasses * nTemplatesPerClass;
T_all = zeros(nTotalTemplates, C*T_spk);
templateClassLabels = zeros(nTotalTemplates, 1);
templateClassIdx = zeros(nTotalTemplates, 1);
templateVariantIdx = zeros(nTotalTemplates, 1);

for i = 1:nValidClasses
    origIdx = validClassIdx(i);
    for t = 1:nTemplatesPerClass
        flatIdx = (i-1)*nTemplatesPerClass + t;
        T_all(flatIdx, :) = reshape(squeeze(templates(origIdx, t, :, :)), 1, []);
        templateClassLabels(flatIdx) = classLabels(origIdx);
        templateClassIdx(flatIdx) = origIdx;
        templateVariantIdx(flatIdx) = t;
    end
end

nonZeroTemplates = any(T_all ~= 0, 2);
T_all = T_all(nonZeroTemplates, :);
templateClassLabels = templateClassLabels(nonZeroTemplates);
templateClassIdx = templateClassIdx(nonZeroTemplates);
templateVariantIdx = templateVariantIdx(nonZeroTemplates);
nTotalTemplates = size(T_all, 1);

if nTotalTemplates == 0
    predLabels = -ones(Nspikes, 1);
    matchScores = zeros(Nspikes, 1);
    validKeep = false(Nspikes, 1);
    bestVariantIdx = ones(Nspikes, 1);
    return;
end

bestTemplateIdx = [];

% GPU only for the matrix-multiply branches once the FLOP count clears
% the host<->device transfer overhead.
flopMM   = double(Nspikes) * double(nTotalTemplates) * double(C * T_spk);
useGPU_TM = cfg.useGPU && (gpuDeviceCount > 0) && ...
            ismember(lower(method), {'euclidean','weighted_euclidean','correlation','cosine'}) && ...
            flopMM >= 5e8;

switch lower(method)
    case 'euclidean'
        if useGPU_TM
            S_g = gpuArray(S); T_g = gpuArray(T_all);
        else
            S_g = S; T_g = T_all;
        end
        T_norm = sum(T_g.^2, 2)';
        S_norm = sum(S_g.^2, 2);
        distances = S_norm + T_norm - 2 * (S_g * T_g');
        distances = sqrt(max(distances, 0));

        [minDist, bestTemplateIdx] = min(distances, [], 2);
        if useGPU_TM
            minDist = gather(minDist);
            bestTemplateIdx = gather(bestTemplateIdx);
        end
        predLabels = templateClassLabels(bestTemplateIdx);
        matchScores = minDist;
        validKeep = true(Nspikes, 1);

    case 'weighted_euclidean'
        weights = computeMiddleChannelWeights(C, T_spk, cfg);
        if useGPU_TM
            S_w = gpuArray(S) .* gpuArray(weights);
            T_w = gpuArray(T_all) .* gpuArray(weights);
        else
            S_w = S .* weights;
            T_w = T_all .* weights;
        end

        T_norm = sum(T_w.^2, 2)';
        S_norm = sum(S_w.^2, 2);
        distances = S_norm + T_norm - 2 * (S_w * T_w');
        distances = sqrt(max(distances, 0));

        [minDist, bestTemplateIdx] = min(distances, [], 2);
        if useGPU_TM
            minDist = gather(minDist);
            bestTemplateIdx = gather(bestTemplateIdx);
        end
        predLabels = templateClassLabels(bestTemplateIdx);
        matchScores = minDist;
        validKeep = true(Nspikes, 1);
        
    case 'polynomial'
        [predLabels, matchScores, validKeep, bestTemplateIdx] = polynomialKernelMatching(...
            S, spikes, T_all, templates, templateClassLabels, templateClassIdx, ...
            validClassIdx, classLabels, C, T_spk, cfg);
        
    case 'poly_weighted'
        [predLabels, matchScores, validKeep, bestTemplateIdx] = polyWeightedMatching(...
            S, spikes, T_all, templates, templateClassLabels, templateClassIdx, ...
            validClassIdx, classLabels, C, T_spk, cfg);
        
    case 'mahalanobis'
        [predLabels, matchScores, validKeep, bestTemplateIdx] = mahalanobisMatching(S, T_all, templateClassLabels, cfg);
        
    case 'correlation'
        if useGPU_TM
            T_g = gpuArray(T_all); S_g = gpuArray(S);
        else
            T_g = T_all; S_g = S;
        end
        T_centered = T_g - mean(T_g, 2);
        S_centered = S_g - mean(S_g, 2);
        T_std = sqrt(sum(T_centered.^2, 2))';
        S_std = sqrt(sum(S_centered.^2, 2));
        T_std(T_std == 0) = eps;
        S_std(S_std == 0) = eps;

        corrMatrix = (S_centered * T_centered') ./ (S_std .* T_std);
        [maxCorr, bestTemplateIdx] = max(corrMatrix, [], 2);
        if useGPU_TM
            maxCorr = gather(maxCorr);
            bestTemplateIdx = gather(bestTemplateIdx);
        end
        predLabels = templateClassLabels(bestTemplateIdx);
        matchScores = maxCorr;
        if isfield(cfg, 'minCorrelation')
            minCorr = cfg.minCorrelation;
        else
            minCorr = 0.5;
        end
        validKeep = maxCorr >= minCorr;

    case 'cosine'
        if useGPU_TM
            T_g = gpuArray(T_all); S_g = gpuArray(S);
        else
            T_g = T_all; S_g = S;
        end
        T_norm = sqrt(sum(T_g.^2, 2))';
        S_norm = sqrt(sum(S_g.^2, 2));
        T_norm(T_norm == 0) = eps;
        S_norm(S_norm == 0) = eps;

        simMatrix = (S_g * T_g') ./ (S_norm .* T_norm);
        [maxSim, bestTemplateIdx] = max(simMatrix, [], 2);
        if useGPU_TM
            maxSim = gather(maxSim);
            bestTemplateIdx = gather(bestTemplateIdx);
        end
        predLabels = templateClassLabels(bestTemplateIdx);
        matchScores = maxSim;
        validKeep = maxSim >= 0.5;
        
    otherwise
        error('Unknown method: %s.', method);
end

bestVariantIdx = templateVariantIdx(bestTemplateIdx);

% Amplitude-based validation
uniquePred = unique(predLabels(validKeep));
for i = 1:length(uniquePred)
    if uniquePred(i) == -1
        continue;
    end
    
    origClassIdx = find(classLabels == uniquePred(i), 1);
    
    if isempty(origClassIdx)
        continue;
    end
    
    pred_idx = find(predLabels == uniquePred(i));
    
    if length(pred_idx) < 5
        continue;
    end
    
    if polarity(origClassIdx) == 1
        outliers = (ampVal(pred_idx) < low_thr(origClassIdx) | ampVal(pred_idx) > high_thr(origClassIdx));
    else
        outliers = (ampVal(pred_idx) > low_thr(origClassIdx) | ampVal(pred_idx) < high_thr(origClassIdx));
    end
    
    validKeep(pred_idx(outliers)) = false;
end

end


function weights = computeMiddleChannelWeights(C, T, cfg)
% Per-channel x per-time weighting biased toward the middle channel/peak.
middle_ch = ceil(C / 2);

if isfield(cfg, 'middleChannelWeight')
    midWeight = cfg.middleChannelWeight;
else
    midWeight = 1;
end

ch_weights = ones(1, C);
ch_weights(middle_ch) = midWeight;
for c = 1:C
    if c ~= middle_ch
        dist = abs(c - middle_ch);
        ch_weights(c) = 1 + (midWeight - 1) * exp(-dist^2 / 2);
    end
end

middle_t = ceil(T / 2);
t_sigma = max(T / 4, 3);
t_weights = 0.5 + 0.5 * exp(-((1:T) - middle_t).^2 / (2 * t_sigma^2));

[T_grid, C_grid] = meshgrid(t_weights, ch_weights);
weights = C_grid .* T_grid;
weights = weights / mean(weights(:));
weights = reshape(weights, 1, []);
end


function [predLabels, matchScores, validKeep, bestTemplateIdx] = polynomialKernelMatching(...
    S, spikes, T_all, templates, templateClassLabels, templateClassIdx, ...
    validClassIdx, classLabels, C, T, cfg)
% Polynomial kernel matching K(x,y) = (gamma * x.y + coef0)^degree.

[Nspikes, nFeatures] = size(S);
nTemplates = size(T_all, 1);
middle_ch = ceil(C / 2);

if isfield(cfg, 'polyDegree'), degree = cfg.polyDegree; else, degree = 3; end
if isfield(cfg, 'polyGamma'),  gamma  = cfg.polyGamma;  else, gamma  = 1 / nFeatures; end
if isfield(cfg, 'polyCoef0'),  coef0  = cfg.polyCoef0;  else, coef0  = 1; end

S_norm = sqrt(sum(S.^2, 2));
T_norm = sqrt(sum(T_all.^2, 2));
S_norm(S_norm == 0) = eps;
T_norm(T_norm == 0) = eps;
S_normalized = S ./ S_norm;
T_normalized = T_all ./ T_norm;

dotProducts = S_normalized * T_normalized';
kernelSim = (gamma * dotProducts + coef0).^degree;

[maxSim, bestTemplateIdx] = max(kernelSim, [], 2);
predLabels = templateClassLabels(bestTemplateIdx);
matchScores = maxSim;

if isfield(cfg, 'polyThreshold')
    threshold = cfg.polyThreshold;
else
    threshold = 0.5 * (gamma + coef0)^degree;
end

validKeep = maxSim >= threshold;

end


function [predLabels, matchScores, validKeep, bestTemplateIdx] = polyWeightedMatching(...
    S, spikes, T_all, templates, templateClassLabels, templateClassIdx, ...
    validClassIdx, classLabels, C, T, cfg)
% Polynomial kernel on the middle channel + weighted-cosine on the full
% waveform + a small derivative-similarity bonus. Per-class threshold is
% set from intra-class self-similarity; spikes also need an amplitude
% ratio in [0.2, 5] vs the matched template.

[Nspikes, nFeatures] = size(S);
nTemplates = size(T_all, 1);
middle_ch = ceil(C / 2);
middle_t = ceil(T / 2);

if isfield(cfg, 'polyDegree'),          degree    = cfg.polyDegree;          else, degree    = 3;   end
if isfield(cfg, 'middleChannelWeight'), midWeight = cfg.middleChannelWeight; else, midWeight = 1.5; end

S_3d = reshape(S, Nspikes, C, T);
T_3d = reshape(T_all, nTemplates, C, T);
S_mid = reshape(S_3d(:, middle_ch, :), Nspikes, T);
T_mid = reshape(T_3d(:, middle_ch, :), nTemplates, T);

% Normalise middle channel by peak (fall back to max-abs when peak ~ 0).
S_peak = S_mid(:, middle_t);
T_peak = T_mid(:, middle_t);
for i = 1:Nspikes
    if abs(S_peak(i)) < eps
        [~, idx] = max(abs(S_mid(i,:)));
        S_peak(i) = S_mid(i, idx);
    end
end
for i = 1:nTemplates
    if abs(T_peak(i)) < eps
        [~, idx] = max(abs(T_mid(i,:)));
        T_peak(i) = T_mid(i, idx);
    end
end
S_peak(abs(S_peak) < eps) = eps;
T_peak(abs(T_peak) < eps) = eps;
S_mid_norm = S_mid ./ abs(S_peak);
T_mid_norm = T_mid ./ abs(T_peak);

gamma_mid = 1 / T;
coef0 = 1;
dotProd_mid = S_mid_norm * T_mid_norm';
kernelSim_mid = (gamma_mid * dotProd_mid + coef0).^degree;

weights = computeMiddleChannelWeights(C, T, cfg);
S_w = S .* weights;
T_w = T_all .* weights;
S_w_norm = sqrt(sum(S_w.^2, 2));
T_w_norm = sqrt(sum(T_w.^2, 2));
S_w_norm(S_w_norm == 0) = eps;
T_w_norm(T_w_norm == 0) = eps;
S_w_normalized = S_w ./ S_w_norm;
T_w_normalized = T_w ./ T_w_norm;
cosineSim = S_w_normalized * T_w_normalized';

if isfield(cfg, 'midChannelKernelWeight')
    alpha = cfg.midChannelKernelWeight;
else
    alpha = 0.6;
end
combinedScore = alpha * kernelSim_mid + (1 - alpha) * cosineSim;

% Derivative-similarity bonus.
S_mid_deriv = diff(S_mid_norm, 1, 2);
T_mid_deriv = diff(T_mid_norm, 1, 2);
S_deriv_norm = sqrt(sum(S_mid_deriv.^2, 2));
T_deriv_norm = sqrt(sum(T_mid_deriv.^2, 2));
S_deriv_norm(S_deriv_norm == 0) = eps;
T_deriv_norm(T_deriv_norm == 0) = eps;
S_mid_deriv_n = S_mid_deriv ./ S_deriv_norm;
T_mid_deriv_n = T_mid_deriv ./ T_deriv_norm;
derivSim = S_mid_deriv_n * T_mid_deriv_n';
combinedScore = combinedScore + 0.1 * derivSim;

[maxScore, bestTemplateIdx] = max(combinedScore, [], 2);
predLabels = templateClassLabels(bestTemplateIdx);
matchScores = maxScore;

validKeep = true(Nspikes, 1);

% Per-class threshold = 0.1 * min intra-class self-similarity.
uniqueClasses = unique(templateClassLabels);
classThresholds = zeros(length(uniqueClasses), 1);
for i = 1:length(uniqueClasses)
    classLabel = uniqueClasses(i);
    templateMask = find(templateClassLabels == classLabel);
    if length(templateMask) > 1
        classMidNorm = T_mid_norm(templateMask, :);
        selfKernel = (gamma_mid * (classMidNorm * classMidNorm') + coef0).^degree;
        selfKernel(logical(eye(size(selfKernel)))) = inf;
        classThresholds(i) = 0.1 * min(selfKernel(:));
    else
        classThresholds(i) = 0.1;
    end
end

for i = 1:Nspikes
    predClass = predLabels(i);
    classIdx = find(uniqueClasses == predClass, 1);
    if ~isempty(classIdx)
        midKernelScore = kernelSim_mid(i, bestTemplateIdx(i));
        if midKernelScore < classThresholds(classIdx)
            validKeep(i) = false;
            continue;
        end
        ampRatio = abs(S_peak(i)) / abs(T_peak(bestTemplateIdx(i)));
        if ampRatio < 0.2 || ampRatio > 5.0
            validKeep(i) = false;
        end
    end
end

end


function [predLabels, matchScores, validKeep, bestTemplateIdx] = mahalanobisMatching(S, T_all, templateClassLabels, cfg)

[Nspikes, nFeatures] = size(S);
nTemplates = size(T_all, 1);

if isfield(cfg, 'mahalRegularization')
    lambda = cfg.mahalRegularization;
else
    lambda = 0.01;
end

if Nspikes > nFeatures
    covS = cov(S);
else
    covS = (S' * S) / max(Nspikes - 1, 1);
end

covS = covS + lambda * trace(covS) / nFeatures * eye(nFeatures);

try
    L = chol(covS, 'lower');
    S_white = (L \ S')';
    T_white = (L \ T_all')';
catch
    [U, D] = eig(covS);
    d = diag(D);
    d(d < eps) = eps;
    sqrtInvCov = U * diag(1 ./ sqrt(d)) * U';
    S_white = S * sqrtInvCov;
    T_white = T_all * sqrtInvCov;
end

T_norm = sum(T_white.^2, 2)';
S_norm = sum(S_white.^2, 2);
distances = S_norm + T_norm - 2 * (S_white * T_white');
distances = sqrt(max(distances, 0));

[minDist, bestTemplateIdx] = min(distances, [], 2);
predLabels = templateClassLabels(bestTemplateIdx);
matchScores = minDist;

if isfield(cfg, 'mahalThreshold')
    threshold = cfg.mahalThreshold;
else
    threshold = chi2inv(0.99, min(nFeatures, 50));
end
validKeep = minDist.^2 < threshold;

end