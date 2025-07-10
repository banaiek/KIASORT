function [out, out_templateFeatures] = kiaSort_cluster_classify(data, cfg, hp)

%% Import parameters from config struct
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

%% Flatten the spike waveform matrix
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


ham_win = hamming(C)';
ham_win = repmat(ham_win,1,T);
flattend_waveform = reshape(waveform, N, C*T);
hamWin = hamming(size(flattend_waveform,2));
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

%% UMAP dimensionality reduction

PCA_waveform = flattend_waveform(:,sum(flattend_waveform,1)~=0);

warning('off', 'stats:pca:ColRankDefX');
warning('off', 'MATLAB:class:DynPropDuplicatesMethod');

% calculate umap and pca
umap_out = pythonUMAP(flat_umapWaveforms,nComp);
[PCA.coeff,PCA.score,PCA.latent,PCA.tsquared,PCA.explained,PCA.mu] = pca(PCA_waveform,'Algorithm','svd','NumComponents',nPCAcomp);
warning(warning);

% use the first PCA scores, construct matrix of UMAP, PCA, and main
% Amplitude
PCA_score = PCA.score;

[umapNorm, ~] = mapminmax(umap_out', 0, 1);
[PCA_score_clustering, ~] = mapminmax([PCA_score(:,1:3)]', 0, 1);

umapNorm = umapNorm';
dataAll = [umapNorm, PCA_score_clustering'];

% Estimate the epsilon and min points for dbscan clustering

[epsilon, numPt] = estimate_dbscan_par(dataAll);
numPt = max(min([ numPt, maxClusterPoints]),minClusterPoints);
[nPoints, nDims] = size(dataAll);

% Lower bound: the greater of 0.1% of nPoints and minClusterPoints vs. 0.5% of nPoints
minPctLimit    = round(0.001 * nPoints) + 5;
minPctDefault  = round(0.005 * nPoints) + 5;
minPtsStart    = max( min( minClusterPoints, minPctDefault ), minPctLimit);

% Upper bound: the lesser of 2% of nPoints and maxClusterPoints vs. 1% of nPoints
maxPctLimit    = round(0.005 * nPoints) + 5;
maxPctDefault  = round(0.0075 * nPoints) + 5;
minPtsEnd      = min( max( maxClusterPoints, maxPctDefault ), maxPctLimit);

numPt = max([min([ numPt, minPtsEnd]),minPtsStart, nDims+1]);


% cluster data using dbscan clustering  
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
    Id_idx = find(labels == uniqueLabels(i));
    class_polarity(i,1) = mode(spk_ID_full(Id_idx));
    misMatched_ID = spk_ID_full(Id_idx) ~= class_polarity(i,1);
    trainingLbl(Id_idx(misMatched_ID)) = -1;
    noise_idx = farthest_points(PCA_score(Id_idx,1:3));
    trainingLbl(Id_idx(noise_idx)) = -1;
end

clusterSampleCounts = zeros(numUniqueClusters,1);

% mean spike waveform for the full-length waveform
meanClusterWaveform = zeros([numUniqueClusters, size(waveform,[2,3])]);

for i=1:numUniqueClusters
    idx = labels==uniqueLabels(i);
    clusterSampleCounts(i) = sum(idx);
    meanClusterWaveform(i,:,:)  = mean(waveform(idx,:,:), 1, 'omitmissing');    
end

% calculating spike density for each cluster
[clusterSpikeDensity, ~, ~] = cluster_spike_density(spk_idx_full, labels, cfg);


% process clusters and find clusters to be merged.
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

% find the reference channel for each cluster
[clusterSelection] = kiaSort_best_channel_detection(clusterData, 100, cfg, 1);
clusterSelection.classLabels = uniqueNewLabels;


% exclude noisy clusters from the training process
for i = 1:length(uniqueNewLabels)
    if (perRatio(i) <.5 && clusterSelection.rank(i) > ceil(C/4) ) || (perRatio(i) <.25 && clusterSelection.rank(i) > 5) 
        trainingLbl(updatedLabels==uniqueNewLabels(i))=-1;
    end
end
% 

low_thr = nan(length(uniqueLabels),1);
high_thr = nan(length(uniqueLabels),1);

stablePoints = clusterRelabeling.stablePoints;

for i = 1:length(uniqueLabels)
        trainingLbl(labels==uniqueLabels(i) & spk_idx_full < stablePoints(i,1) & spk_idx_full > stablePoints(i,2))=-1;
        class_idx = labels==uniqueLabels(i);
        lowPrc = prctile(ampVals(class_idx),2.5);
        highPrc = prctile(ampVals(class_idx),97.5);
        trainingLbl(labels==uniqueLabels(i) & (ampVals < lowPrc | ampVals > highPrc)) = -1;
        
        if class_polarity(i) == 1
            low_thr(i) = 0.5 * lowPrc;
            high_thr(i) = 1.5 * highPrc;
        else
            low_thr(i) = 0.5 * highPrc;
            high_thr(i) = 1.5 * lowPrc;
        end

end


%% Method and Model selection for classification
if sum(trainingLbl~=-1)<50 
    trainingLbl = labels;
end

% Splitting the data into train and test sets
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

net = [];
mdl = [];

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

% method type
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



% test the model
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


out_templateFeatures.umapNorm        = umapNorm;
out_templateFeatures.labels          = labels;
out_templateFeatures.updatedLabels   = updatedLabels;
out_templateFeatures.spk_idx         = spk_idx_full;
out_templateFeatures.PCA_scores      = PCA_score;
out_templateFeatures.classLabels     = uniqueLabels;
out_templateFeatures.meanWaveform    = meanClusterWaveform;

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
out.waveformInfo.clusterSize        = clusterSampleCounts;
out.waveformInfo.informative_Chan   = informative_Chans;

out.cfg                             = cfg;
end

function [labels, globalCluster] = processCluster(idx, labels, fs, globalCluster, dataAll, spk_idx_full, cfg, depth, maxClusterPoints, minClusterPoints, sample_dur, numPt, initNumClasses, snr_vals)
[~,~,isi_viol] = getISIViolations(spk_idx_full(idx), fs, 2);
factor = length(idx)/size(dataAll,1);
if ((isi_viol<= 0.1 && depth >= 0) && (factor<0.25 || depth > 0) && (initNumClasses>10 || depth > 0)) || depth >= 3
    globalCluster = globalCluster + 1;
    labels(idx) = globalCluster;
elseif isi_viol> 2 && length(idx) < minClusterPoints/2 && depth > 1
    labels(idx) = -1;
else
    [epsilon, ~] = estimate_dbscan_par(dataAll(idx,:));
    numPt = max(min([ max(factor*numPt,size(dataAll,2)) , factor*maxClusterPoints+5]),factor*minClusterPoints+5);
    numPt = max(numPt,size(dataAll,2));
    subLabels = dbscan(dataAll(idx,:), epsilon, round(numPt),'Distance','minkowski','P',1);

    lowSNR = prctile(snr_vals(idx),2.5);
    highSNR = prctile(snr_vals(idx),97.5);
    medSNR = median(snr_vals(idx));

    if depth < 1 && factor > .49 && (length(idx)/sample_dur) > .5 && initNumClasses < 5 && medSNR<1.5
        try
            % tempLabels = kiaSort_comp_clustering(subLabels, dataAll(idx,:));
             tempLabels = kiaSort_graph_clustering(dataAll(idx,:), 'existingLabels', subLabels, 'qualityThreshold', .8);
            if (sum(tempLabels == -1) / length(idx)) <= .15
                subLabels = tempLabels;
            end
        catch
        end
    elseif depth < 1 && factor > .25 && (length(idx)/sample_dur) > .25 && initNumClasses < 5 && medSNR<1.5
        try
            % tempLabels = kiaSort_comp_clustering(subLabels, dataAll(idx,:));
             tempLabels = kiaSort_graph_clustering(dataAll(idx,:), 'existingLabels', subLabels, 'qualityThreshold', .8);
            if (sum(tempLabels == -1) / length(idx)) <= .15
                subLabels = tempLabels;
            end
        catch
        end

        %this part is added 
    elseif depth < 1  && (length(idx)/sample_dur) > .25 && ((medSNR-lowSNR) < (highSNR-medSNR)/1.25 || (medSNR-lowSNR) > 1.25 * (highSNR-medSNR)) 
        try
            tempLabels = kiaSort_graph_clustering(dataAll(idx,:), 'existingLabels', subLabels, 'qualityThreshold', .8);
            if (sum(tempLabels == -1) / length(idx)) < .1
                subLabels = tempLabels;
            end
        catch
        end
    end
    
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



