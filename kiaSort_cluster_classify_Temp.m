function [out, out_sampleFeatures] = kiaSort_cluster_classify_Temp(data, cfg, hp)

modelType               = cfg.modelType;
method                  = cfg.method;
sampleSpikeDuration     = cfg.spikeDuration;
fs                      = cfg.samplingFrequency;
clusteringSpikeDuration = cfg.clusteringSpikeDuration;
nComp                   = cfg.umapNComp;
testFraction            = cfg.testFraction;
usePCA                  = cfg.usePCA;
nPCAcomp                = cfg.nPCAcomp;
sample_dur = cfg.numSampleChunks * cfg.sampleChunkDuration;
sample_points = sample_dur * fs;
minClusterPoints = round(sample_dur * cfg.minRate);
maxClusterPoints = round(sample_dur * cfg.maxClusteringRate);

if isfield(cfg, 'modelType')
    if strcmp(cfg.modelType,'template')
        useTemplate = true;
    end
else
    useTemplate = false;
end

if isfield(cfg, 'numTemplatesPerCluster')
    numTemplatesPerCluster = cfg.numTemplatesPerCluster;
else
    numTemplatesPerCluster = 15;
end

midPoint = floor(sampleSpikeDuration * fs/(2*1000))+1;
spike_length = floor(clusteringSpikeDuration * fs/(2*1000));
spk_idx_full    = data.spk_idx_full;
spk_ID_full     = data.spk_ID_full;
rel_length = numel(spk_idx_full)/nPCAcomp;

while rel_length <= 1
    data.waveform_bp_full = [data.waveform_bp_full; data.waveform_bp_full];
    spk_idx_full = [spk_idx_full; spk_idx_full + 1000];
    spk_ID_full = [spk_ID_full; spk_ID_full];
    rel_length = numel(spk_idx_full)/nPCAcomp;
end

waveform_bp_full = data.waveform_bp_full;
waveform_bp_full(isnan(waveform_bp_full)) = 0;
waveform = waveform_bp_full(:,:,midPoint-spike_length:midPoint+spike_length);
[N, C, T] = size(waveform);
midChannel = ceil(C/2);
ampVals = waveform(:,midChannel,spike_length+1);

snr_vals = abs(ampVals./data.channel_thresholds_pos(midChannel));

channel_Var = var(waveform,[],[1,3]);
informative_Chans = (channel_Var./max(channel_Var) > 0.05);

umapWaveforms = waveform(:,find(informative_Chans),:);
[N2, C2, T2] = size(umapWaveforms);

flat_umapWaveforms = reshape(umapWaveforms, N2, C2*T2);



[~, C_full, T_full] = size(waveform_bp_full);
flattend_waveform_full = reshape(waveform_bp_full, N, C_full*T_full);
flattend_waveform = reshape(waveform, N, C*T);
hamWin = hamming(size(flattend_waveform,2)).^2;
hamWin = reshape(hamWin,C,T);
hamWin = hamWin(find(informative_Chans),:);
hamWin = reshape(hamWin, 1, C2*T2);
flat_umapWaveforms = flat_umapWaveforms.*hamWin;

if isfield(data.side_waveforms, 'waveform')
    if size(data.side_waveforms.waveform,3)>50
        side_waveforms = data.side_waveforms.waveform(:,:,midPoint-spike_length:midPoint+spike_length);
        side_waveforms(isnan(side_waveforms)) = 0;
        N_side = size(side_waveforms,1);
        flattend_side_waveforms = reshape(side_waveforms, N_side, C*T);
        [~,side_PCscores,~] = pca(flattend_side_waveforms,'Algorithm','eig','NumComponents',min(nPCAcomp,10));

        [epsilon, numPt] = estimate_dbscan_par(side_PCscores);
        numPt = max(min([ numPt, maxClusterPoints]),minClusterPoints);
        side_labels = dbscan(side_PCscores, epsilon, numPt,'Distance','minkowski','P',2);

        unique_side_labels = unique(side_labels);
        unique_side_labels(unique_side_labels==-1) = [];
        mean_side_waveforms = zeros(length(unique_side_labels), C, T);
        for i = 1:length(unique_side_labels)
            mean_side_waveforms(i,:,:) = mean(side_waveforms(side_labels == unique_side_labels(i), :, :),1 ,'omitmissing');
        end
    else
        mean_side_waveforms = [];
    end
else
    mean_side_waveforms = [];
end

PCA_waveform = flattend_waveform_full(:,sum(flattend_waveform_full,1)~=0);

warning('off', 'stats:pca:ColRankDefX');
warning('off', 'MATLAB:class:DynPropDuplicatesMethod');

umap_out = pythonUMAP(flat_umapWaveforms,nComp);
[PCA.coeff,PCA.score,PCA.latent,PCA.tsquared,PCA.explained,PCA.mu] = pca(PCA_waveform,'Algorithm','svd','NumComponents',nPCAcomp);
warning(warning);

PCA_score = PCA.score;

[umapNorm, ~] = mapminmax(umap_out', 0, 1);
[PCA_score_clustering, ~] = mapminmax([PCA_score(:,1:3)]', 0, 1);

umapNorm = umapNorm';
dataAll = [umapNorm, PCA_score_clustering'];

[epsilon, numPt] = estimate_dbscan_par(dataAll);
numPt = max(min([ numPt, maxClusterPoints]),minClusterPoints);
[nPoints, nDims] = size(dataAll);

minPctLimit    = round(0.001 * nPoints) + 5;
minPctDefault  = round(0.005 * nPoints) + 5;
minPtsStart    = max( min( minClusterPoints, minPctDefault ), minPctLimit);

maxPctLimit    = round(0.005 * nPoints) + 5;
maxPctDefault  = round(0.0075 * nPoints) + 5;
minPtsEnd      = min( max( maxClusterPoints, maxPctDefault ), maxPctLimit);

pt = estMinPts(dataAll);
numPt = min(max([min([ numPt, minPtsEnd]),minPtsStart, nDims+1]),pt);


initialLabels = noise_exclude_dbscan(dataAll, ampVals, data.channel_thresholds_pos(midChannel), epsilon, numPt);

labels = -ones(size(initialLabels));
globalCluster = 0;
uniqueLabels = unique(initialLabels);
initNumClasses = max(uniqueLabels);
for i = 1:length(uniqueLabels)
    if uniqueLabels(i)==-1, continue; end
    idx = find(initialLabels==uniqueLabels(i));
    [labels, globalCluster] = processCluster(idx, labels, fs, globalCluster, dataAll, spk_idx_full, cfg, 0, minPtsEnd, minPtsStart, sample_dur, numPt,initNumClasses, snr_vals);
end

if ~any(labels>=1)
    labels(:) = 1;
end

uniqueLabels = unique(labels);
numUniqueClusters = length(uniqueLabels);
class_polarity = zeros(numUniqueClusters,1);

trainingLbl = labels;

for i = 1:numUniqueClusters
    clusterIdx = find(labels == uniqueLabels(i));
    class_polarity(i,1) = mode(spk_ID_full(clusterIdx));

    misMatched_ID = spk_ID_full(clusterIdx) ~= class_polarity(i);
    trainingLbl(clusterIdx(misMatched_ID)) = -1;

    if size(PCA_score, 2) >= 3
        noise_idx = identifyOutliers(PCA_score(clusterIdx, 1:nPCAcomp));
        trainingLbl(clusterIdx(noise_idx)) = -1;
    end
end

clusterSampleCounts = zeros(numUniqueClusters,1);

meanClusterWaveform = zeros([numUniqueClusters, size(waveform,[2,3])]);

for i=1:numUniqueClusters
    idx = labels==uniqueLabels(i);
    clusterSampleCounts(i) = sum(idx);
    meanClusterWaveform(i,:,:)  = mean(waveform(idx,:,:), 1, 'omitmissing');
end

templateWaveforms = zeros([numUniqueClusters, numTemplatesPerCluster, C, T]);
templateWeights = zeros(numUniqueClusters, numTemplatesPerCluster);

for i = 1:numUniqueClusters
    if uniqueLabels(i) == -1
        templateWaveforms(i,:,:,:) = 0;
        templateWeights(i,:) = 0;
        continue;
    end

    idx = find(labels == uniqueLabels(i));
    clusterWaveforms = waveform(idx,:,:);
    clusterPCAscores = PCA_score(idx,:);

    if length(idx) >= numTemplatesPerCluster 
        [templates, weights] = generateMultipleTemplates(clusterWaveforms, numTemplatesPerCluster, clusterPCAscores);
        templateWaveforms(i,:,:,:) = templates;
        templateWeights(i,:) = weights;
    else
        meanWF = squeeze(meanClusterWaveform(i,:,:));
        for t = 1:numTemplatesPerCluster
            templateWaveforms(i,t,:,:) = meanWF;
        end
        templateWeights(i,:) = 1/numTemplatesPerCluster;
    end
end

[clusterSpikeDensity, ~, ~] = cluster_spike_density(spk_idx_full, labels, cfg);

[clusterRelabeling]  = kiaSort_process_clusters(meanClusterWaveform, clusterSpikeDensity, ...
    uniqueLabels, labels, PCA, spk_idx_full, mean_side_waveforms, cfg);

clusterRelabeling.originalLabels = uniqueLabels;
clusterRelabeling.mean_side_waveforms = mean_side_waveforms;
[clusterRelabeling] = realign_merge_Waveforms(meanClusterWaveform, clusterSampleCounts, clusterRelabeling);
updatedLabels = updateLabels(labels, uniqueLabels, clusterRelabeling.newLabels);


clusterData.channel_thresholds_pos = data.channel_thresholds_pos;
clusterData.channel_thresholds_neg = data.channel_thresholds_neg;
clusterData.waveform               = clusterRelabeling.newMeanWaveforms;
clusterData.unmerged_sampleCounts      = clusterSampleCounts;
clusterData.wavformChanelIdx = data.wavformChanelIdx;

[uniqueNewLabels, uniqLblID] = unique(clusterRelabeling.newLabels);
perRatio = clusterRelabeling.perRatio(uniqLblID);

[clusterSelection] = kiaSort_best_channel_detection(clusterData, 100, cfg, 1);
clusterSelection.classLabels = uniqueNewLabels;

clusterStatus = zeros(length(uniqueNewLabels), 1);
contaminationRate = zeros(length(uniqueNewLabels), 1);

ACG_R_CLEAN = 0.10;   % ratio threshold for "clean"
ACG_Q_CLEAN = 0.20;   % Poisson p-value threshold for "clean"
CONT_SEVERE = 0.50;   % contamination rate above which unit is suspect
ISI_BACKUP_THR1 = 0.1;   % ISI violation % at 1 ms (backup check)
ISI_BACKUP_THR2 = 0.2;   % ISI violation % at 2 ms (backup check)

for iMerged = 1:length(uniqueNewLabels)
    mergedLbl = uniqueNewLabels(iMerged);
    if mergedLbl == -1
        continue;
    end
    spk_in_merged = spk_idx_full(updatedLabels == mergedLbl);
    if numel(spk_in_merged) < 2
        continue;
    end

    [contRate, R, Q] = estimateContamination(spk_in_merged, [], fs);
    contaminationRate(iMerged) = contRate;

    if min(R) < ACG_R_CLEAN && min(Q) < ACG_Q_CLEAN
        
        clusterStatus(iMerged) = 0;
    elseif contRate >= CONT_SEVERE

        [~, ~, isv1] = getISIViolations(spk_in_merged, fs, 1);
        [~, ~, isv2] = getISIViolations(spk_in_merged, fs, 2);

        if isv1 > ISI_BACKUP_THR1 || isv2 > ISI_BACKUP_THR2
            clusterStatus(iMerged) = -4;
            clusterSelection.keep(iMerged) = 0;
        end
    end

end

if isfield(clusterSelection, 'lowAmpNotKept') && any(clusterSelection.lowAmpNotKept)
    lowAmpIdx = find(clusterSelection.lowAmpNotKept);
    for iLow = 1:length(lowAmpIdx)
        idx = lowAmpIdx(iLow);
        if uniqueNewLabels(idx) == -1
            continue;
        end
        if clusterStatus(idx) ~= 0
            continue;
        end
        if clusterSelection.mainNegativePolarity(idx)
            clusterStatus(idx) = -3;
        else
            clusterStatus(idx) = -2;
        end
    end
end

clusterSelection.clusterStatus = clusterStatus;
clusterSelection.contaminationRate = contaminationRate;
clusterRelabeling.clusterStatus = clusterStatus;

for i = 1:length(uniqueNewLabels)
    if (perRatio(i) <.5 && clusterSelection.rank(i) > ceil(C/4) ) || (perRatio(i) <.25 && clusterSelection.rank(i) > 5)
        trainingLbl(updatedLabels==uniqueNewLabels(i))=-1;
    end
end

low_thr = nan(length(uniqueLabels),1);
high_thr = nan(length(uniqueLabels),1);

stablePoints = clusterRelabeling.stablePoints;

for i = 1:length(uniqueLabels)
    trainingLbl(labels==uniqueLabels(i) & spk_idx_full < stablePoints(i,1) & spk_idx_full > stablePoints(i,2))=-1;
    class_idx = labels==uniqueLabels(i);
    lowPrc = prctile(ampVals(class_idx),5);
    highPrc = prctile(ampVals(class_idx),95);
    trainingLbl(labels==uniqueLabels(i) & (ampVals < lowPrc | ampVals > highPrc)) = -1;

    if class_polarity(i) == 1
        low_thr(i) = 0.25 * lowPrc;
        high_thr(i) = 2 * highPrc;
    else
        low_thr(i) = 0.25 * highPrc;
        high_thr(i) = 2 * lowPrc;
    end
end

net = [];
mdl = [];
classifierAccuracy = [];

if ~useTemplate
    if sum(trainingLbl~=-1)<50
        trainingLbl = labels;
    end

    numTest    = round(testFraction * N)+1;
    numTrain   = N - numTest;

    classificationLabels = uniqueLabels(uniqueLabels > 0);

    shuffledIdx = randperm(N);
    trainIdx = shuffledIdx(1:numTrain);
    testIdx = shuffledIdx(numTrain+1:end);
    trainIdx(trainingLbl(trainIdx) == -1) = [];
    testIdx(trainingLbl(testIdx) == -1) = [];

    count = 0 ;
    while any(~ismember(uniqueLabels(uniqueLabels~=-1),trainingLbl(trainIdx))) && count < 1
        shuffledIdx = randperm(N);
        trainIdx = shuffledIdx(1:numTrain);
        testIdx = shuffledIdx(numTrain+1:end);
        trainIdx(trainingLbl(trainIdx) == -1) = [];
        testIdx(trainingLbl(testIdx) == -1) = [];
        count = count + 1;
    end

    if ~isfield(cfg, 'method') || ~isfield(cfg, 'modelType')
        error('Configuration (cfg) must contain ''method'' and ''modelType'' fields.');
    end

    if usePCA
        Xinput = [PCA_score];
    else
        Xinput = PCA_waveform;
    end

    switch lower(method)
        case 'direct'
            switch lower(modelType)
                case 'cnn'
                    XTrain = waveform(trainIdx,:,:);
                    YTrain = labels(trainIdx,:);
                    XTestCNN  = waveform(testIdx,:,:);
                    XTestCNN = permute(XTestCNN, [2, 3, 1]);
                    XTest = reshape(XTestCNN, [C, T, 1, size(XTestCNN,3)]);
                    YTest = labels(testIdx,:);
                otherwise
                    XTrain = Xinput(trainIdx,:);
                    YTrain = labels(trainIdx,:);
                    XTest  = Xinput(testIdx,:);
                    YTest = labels(testIdx,:);
            end
        case 'indirect'
            XTrain1 = Xinput;
            YTrain1 = umapNorm;
    end

    switch lower(method)
        case 'direct'
            switch lower(modelType)
                case 'mlp'
                    net = mlpWaveformClassifier(XTrain, YTrain, hp);
                case 'cnn'
                    net = cnnWaveformClassifier(XTrain, YTrain, hp);
                case 'svm'
                    mdl = trainWaveformClassifier(XTrain, YTrain, modelType, hp);
                case 'gbmadaboost'
                    mdl = trainWaveformClassifier(XTrain, YTrain, modelType, hp);
                case 'gbmrusboost'
                    mdl = trainWaveformClassifier(XTrain, YTrain, modelType, hp);
                otherwise
                    error('Unsupported modelType ''%s'' for direct method.', cfg.modelType);
            end
        case 'indirect'
            net = mlpDimReduction(XTrain1, YTrain1, hp);
            umapPredNorm = predict(net, XTrain1);
            XTrain2 = umapPredNorm(trainIdx,:);
            YTrain2 = labels(trainIdx,:);
            XTest   = umapPredNorm(testIdx,:);
            YTest   = labels(testIdx,:);
            switch lower(modelType)
                case 'svm'
                    mdl = trainWaveformClassifier(XTrain2, YTrain2, modelType, hp);
                case 'gbmadaboost'
                    mdl = trainWaveformClassifier(XTrain2, YTrain2, modelType, hp);
                case 'gbmrusboost'
                    mdl = trainWaveformClassifier(XTrain2, YTrain2, modelType, hp);
                otherwise
                    error('Unsupported modelType ''%s'' for indirect method.', cfg.modelType);
            end
        otherwise
            error('Unsupported method ''%s''.', cfg.method);
    end

    switch lower(modelType)
        case 'mlp'
            predLabels = predict(net,XTest);
            predLabels = onehotdecode(predLabels,double(classificationLabels),2);
            classifierAccuracy = confusionmat(YTest, double(predLabels));
        case 'cnn'
            predLabels = predict(net, XTest);
            predLabels = onehotdecode(predLabels,double(classificationLabels),2);
            classifierAccuracy = confusionmat(YTest, double(predLabels));
        otherwise
            predLabels = predict(mdl,XTest);
            classifierAccuracy = confusionmat(YTest, predLabels);
    end
end

out_sampleFeatures.umapNorm        = umapNorm;
out_sampleFeatures.labels          = labels;
out_sampleFeatures.updatedLabels   = updatedLabels;
out_sampleFeatures.spk_idx         = spk_idx_full;
out_sampleFeatures.PCA_scores      = PCA_score;
out_sampleFeatures.classLabels     = uniqueLabels;
out_sampleFeatures.meanWaveform    = meanClusterWaveform;

out.clusteringInfo.epsilon          = epsilon;
out.clusteringInfo.PCA              = PCA;
out.clusteringInfo.numPt            = numPt;
out.clusteringInfo.clusterRelabeling = clusterRelabeling;
out.clusteringInfo.clusterSelection = clusterSelection;
out.clusteringInfo.classLabels      = uniqueLabels;

out.classifierInfo.valAccuracy      = classifierAccuracy;
out.classifierInfo.classLabels      = uniqueLabels;
out.classifierInfo.numClasses       = numUniqueClusters;
out.classifierInfo.hypPar           = hp;
out.classifierInfo.method           = cfg.method;
out.classifierInfo.modelType        = cfg.modelType;
out.classifierInfo.trainedNet       = net;
out.classifierInfo.trainedMdl       = mdl;
out.classifierInfo.class_polarity   = class_polarity;
out.classifierInfo.lowAmpThr        = low_thr;
out.classifierInfo.highAmpThr       = high_thr;

out.waveformInfo.size               = size(waveform);
out.waveformInfo.meanWaveform       = meanClusterWaveform;
out.waveformInfo.templateWaveforms  = templateWaveforms;
out.waveformInfo.templateWeights    = templateWeights;
out.waveformInfo.clusterSize        = clusterSampleCounts;
out.waveformInfo.informative_Chan   = informative_Chans;

out.cfg                             = cfg;
end

function [labels, globalCluster] = processCluster(idx, labels, fs, globalCluster, dataAll, spk_idx_full, cfg, depth, maxClusterPoints, minClusterPoints, sample_dur, numPt, initNumClasses, snr_vals)
[~,~,isi_viol] = getISIViolations(spk_idx_full(idx), fs, 2);
factor = length(idx)/size(dataAll,1);
if ((isi_viol<= 0.1 && depth >= 0)  && (initNumClasses>10 || depth > 0)) || depth >= 3
    globalCluster = globalCluster + 1;
    labels(idx) = globalCluster;
elseif isi_viol> 1 && length(idx) < minClusterPoints/2 && depth > 1
    labels(idx) = -1;
else
    [epsilon, ~] = estimate_dbscan_par(dataAll(idx,:));
    numPt = max(min([ max(factor*numPt,size(dataAll,2)) , factor*maxClusterPoints+5]),factor*minClusterPoints+5);
    numPt = max(numPt,size(dataAll,2));
    subLabels = dbscan(dataAll(idx,:), epsilon, round(numPt),'Distance','minkowski','P',1);

    uniqueSub = unique(subLabels);

    for j = 1:length(uniqueSub)
        if uniqueSub(j)==-1
            labels(idx(subLabels==uniqueSub(j))) = -1;
        else
            subIdx = idx(subLabels==uniqueSub(j));
            [labels, globalCluster] = processCluster(subIdx, labels, fs, globalCluster, dataAll, spk_idx_full, cfg, depth+1, maxClusterPoints, minClusterPoints, sample_dur, numPt, initNumClasses, snr_vals);
        end
    end
end
end

function noise_idx = identifyOutliers(data)
try
    mu = mean(data, 1);
    sigma = cov(data);
    if rcond(sigma) < 1e-10
        sigma = sigma + eye(size(sigma)) * 1e-6;
    end
    mahalDist = mahal(data, data);
    threshold = chi2inv(0.975, size(data, 2));
    noise_idx = mahalDist > threshold;
catch
    distFromCenter = sqrt(sum((data - mean(data)).^2, 2));
    threshold = prctile(distFromCenter, 97.5);
    noise_idx = distFromCenter > threshold;
end
end

function [templates, weights] = generateMultipleTemplates(clusterWaveforms, numTemplates, pcaScores)
[N, C, T] = size(clusterWaveforms);

if N <= numTemplates
    templates = zeros(numTemplates, C, T);
    for i = 1:N
        templates(i,:,:) = clusterWaveforms(i,:,:);
    end
    for i = N+1:numTemplates
        templates(i,:,:) = templates(N,:,:);
    end
    weights = ones(1, numTemplates) / numTemplates;
    return;
end

nPCAuse = min(50, size(pcaScores, 2));
pcaData = pcaScores(:, 1:nPCAuse);

try
    [subLabels, ~] = kmeans(pcaData, numTemplates, 'MaxIter', 200, 'Replicates', 3, 'Options', statset('UseParallel', false));
catch
    midC = ceil(C/2);
    midT = ceil(T/2);
    amps = clusterWaveforms(:, midC, midT);
    [~, sortIdx] = sort(amps);
    subLabels = zeros(N, 1);
    for i = 1:N
        subLabels(sortIdx(i)) = ceil(i / (N/numTemplates));
    end
    subLabels = min(subLabels, numTemplates);
    subLabels = max(subLabels, 1);
end

templates = zeros(numTemplates, C, T);
weights = zeros(1, numTemplates);

for k = 1:numTemplates
    subIdx = (subLabels == k);
    if sum(subIdx) > 0
        templates(k,:,:) = mean(clusterWaveforms(subIdx,:,:), 1, 'omitmissing');
        weights(k) = sum(subIdx);
    else
        templates(k,:,:) = mean(clusterWaveforms, 1, 'omitmissing');
        weights(k) = 1;
    end
end

weights = weights / sum(weights);
end

