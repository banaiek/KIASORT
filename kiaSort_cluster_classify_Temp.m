function [out, out_sampleFeatures] = kiaSort_cluster_classify_Temp(data, cfg, hp)

modelType               = cfg.modelType;
method                  = cfg.method;
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
% Centre from the array itself so the crop tracks the extraction window
% regardless of how the two durations round.
midPoint = floor(size(waveform_bp_full,3) / 2) + 1;
waveform = waveform_bp_full(:,:,midPoint-spike_length:midPoint+spike_length);
[N, C, T] = size(waveform);
midChannel = ceil(C/2);
ampVals = waveform(:,midChannel,spike_length+1);

snr_vals = abs(ampVals./data.channel_thresholds_pos(midChannel));

sinkSnr    = cfgGet(cfg, 'noiseSinkSnr', 1.25);
enableSink = cfgGet(cfg, 'enableNoiseSink', true);
thrPosMid  = abs(data.channel_thresholds_pos(midChannel));
if isempty(data.channel_thresholds_neg) || numel(data.channel_thresholds_neg) < midChannel
    thrNegMid = thrPosMid;
else
    thrNegMid = abs(data.channel_thresholds_neg(midChannel));
end
sinkClassLabels = [];
% Margin around each class's own 1-99 percentile amplitude range.
ampBandLo = cfgGet(cfg, 'ampBandLowFactor',  0.5);
ampBandHi = cfgGet(cfg, 'ampBandHighFactor', Inf);

channel_Var = var(waveform,[],[1,3]);
informative_Chans = (channel_Var./max(channel_Var) > 0.05);

umapWaveforms = waveform(:,find(informative_Chans),:);
[N2, C2, T2] = size(umapWaveforms);

flat_umapWaveforms = reshape(umapWaveforms, N2, C2*T2);



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

% Clustering features come from the clustering crop only, so widening
% spikeDuration changes what is saved without moving any cluster.
PCA_waveform = flattend_waveform(:,sum(flattend_waveform,1)~=0);
nPCAcompEff  = max(1, min(nPCAcomp, size(PCA_waveform,2)));

warning('off', 'stats:pca:ColRankDefX');
warning('off', 'MATLAB:class:DynPropDuplicatesMethod');

umap_out = pythonUMAP(flat_umapWaveforms,nComp);
[PCA.coeff,PCA.score,PCA.latent,PCA.tsquared,PCA.explained,PCA.mu] = pca(PCA_waveform,'Algorithm','svd','NumComponents',nPCAcompEff);
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
splitLog = [];
for i = 1:length(uniqueLabels)
    if uniqueLabels(i)==-1, continue; end
    idx = find(initialLabels==uniqueLabels(i));
    [labels, globalCluster, splitLog] = processCluster(idx, labels, fs, globalCluster, dataAll, spk_idx_full, cfg, 0, minPtsEnd, minPtsStart, sample_dur, numPt,initNumClasses, snr_vals, splitLog);
end

if ~any(labels>=1)
    labels(:) = 1;
end

% Resolve the provisional ISI rejects. A cluster only becomes noise if it is
% above the keep gate AND its polarity still has another class; otherwise it
% is retained so low-amplitude noise has somewhere to go at sort time. The
% existing 1.25 amplitude gate keeps it out of the final units.
rejTags = unique(labels(labels <= -2));
if ~isempty(rejTags)
    jitGap = ceil(cfg.spikeDistance * fs / 1000);
    cIdx   = spike_length + 1;
    wStart = max(1, cIdx - jitGap);
    wEnd   = min(T, cIdx + jitGap);

    for r = 1:numel(rejTags)
        m   = labels == rejTags(r);
        pol = mode(spk_ID_full(m));
        if pol < 0, thr_p = thrNegMid; else, thr_p = thrPosMid; end

        wf      = reshape(mean(waveform(m, midChannel, :), 1, 'omitmissing'), 1, []);
        ampHigh = max(abs(wf(wStart:wEnd))) > sinkSnr * thr_p;

        nPol = 0;
        accepted = unique(labels(labels >= 1));
        for k = 1:numel(accepted)
            if mode(spk_ID_full(labels == accepted(k))) == pol
                nPol = nPol + 1;
            end
        end

        if enableSink && (~ampHigh || nPol < 2)
            globalCluster = globalCluster + 1;
            labels(m) = globalCluster;
            if ~ampHigh
                sinkClassLabels(end+1,1) = globalCluster; %#ok<AGROW>
            end
        else
            labels(m) = -1;
        end
    end
end

% Amplitude bimodality: a class holding a noise band plus a genuinely
% larger unit shows two amplitude modes that coexist in time. Drift also
% splits the amplitude histogram, but its modes are consecutive, so the
% same temporal-overlap veto separates the two cases.
if cfgGet(cfg, 'enableBimodalSplit', false)
    ampAcc = unique(labels(labels >= 1));
    nAmpSplit = 0;
    for a = 1:numel(ampAcc)
        if nAmpSplit >= cfgGet(cfg, 'ampSplitMaxPerChannel', 1), break; end
        m   = find(labels == ampAcc(a));
        if numel(m) < 4*cfgGet(cfg,'peelMinLowCount',20), continue; end
        v   = abs(double(ampVals(m)));
        [thrA, ~, ~, etaA] = kiaSort_otsu_split1d(v, cfgGet(cfg,'bimodalNumBins',64));
        if isnan(thrA) || etaA < cfgGet(cfg,'ampSplitSeparability',0.75)
            continue;
        end
        % The Otsu cut only locates the boundary; refine it so spikes near
        % the edge are assigned by proximity to the two modes rather than by
        % a bin edge. Modes here have unequal spread (tight noise band vs a
        % broader unit), which a raw histogram cut splits poorly.
        thrA = kiaSort_refine_1d_2means(v, thrA);
        lo = m(v <= thrA);  hi = m(v > thrA);
        minA = max(cfgGet(cfg,'peelMinLowCount',20), ...
                   ceil(cfgGet(cfg,'ampSplitMinChildFrac',0.15) * numel(m)));
        if numel(lo) < minA || numel(hi) < minA, continue; end

        tl = double(spk_idx_full(lo)); th = double(spk_idx_full(hi));
        ql = prctile(tl,[5 95]);       qh = prctile(th,[5 95]);
        spanA = min(ql(2)-ql(1), qh(2)-qh(1));
        if spanA <= 0, continue; end
        if (min(ql(2),qh(2)) - max(ql(1),qh(1)))/spanA < cfgGet(cfg,'bimodalMinTimeOverlap',0.5)
            continue;
        end

        globalCluster = globalCluster + 1;
        labels(lo)    = globalCluster;
        nAmpSplit     = nAmpSplit + 1;
    end
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
        noise_idx = identifyOutliers(PCA_score(clusterIdx, 1:min(nPCAcompEff, size(PCA_score,2))));
        trainingLbl(clusterIdx(noise_idx)) = -1;
    end
end

clusterSampleCounts = zeros(numUniqueClusters,1);

meanClusterWaveform = zeros([numUniqueClusters, size(waveform,[2,3])]);
% Full-length twin of meanClusterWaveform. Sorting decisions use the
% clustering-length one; only what gets saved uses this.
meanClusterWaveformFull = zeros([numUniqueClusters, size(waveform_bp_full,[2,3])]);

for i=1:numUniqueClusters
    idx = labels==uniqueLabels(i);
    clusterSampleCounts(i) = sum(idx);
    meanClusterWaveform(i,:,:)  = mean(waveform(idx,:,:), 1, 'omitmissing');
    meanClusterWaveformFull(i,:,:) = mean(waveform_bp_full(idx,:,:), 1, 'omitmissing');
end

% Drop small clusters whose per-spike maxima are concentrated on one
% channel but scattered in time -- likely formed from overlapping spikes
% of other (larger) classes.
sizes_valid = clusterSampleCounts(uniqueLabels >= 0 & ~ismember(uniqueLabels, sinkClassLabels));
if ~isempty(sizes_valid)
    sizeThr       = 0.10 * max(sizes_valid);
    chanFracThr   = 0.80;
    timeJitterThr = max(3, round(T / 5));

    for i = 1:numUniqueClusters
        L = uniqueLabels(i);
        if L < 0 || clusterSampleCounts(i) >= sizeThr || ismember(L, sinkClassLabels), continue; end
        idx = find(labels == L);
        if numel(idx) < 5, continue; end
        [maxPerCh, maxTimePerCh] = max(abs(waveform(idx, :, :)), [], 3);
        [~, mainCh] = max(maxPerCh, [], 2);
        modeCh = mode(mainCh);
        % With one channel the concentration test is vacuously true, which
        % would reduce this rule to a bare time-jitter test.
        if C < 2 || mean(mainCh == modeCh) < chanFracThr, continue; end
        times_on_dom = double(maxTimePerCh(mainCh == modeCh, modeCh));
        if numel(times_on_dom) < 5, continue; end
        if std(times_on_dom) <= timeJitterThr, continue; end
        labels(idx)      = -1;
        trainingLbl(idx) = -1;
    end

    survivingLbls = unique(labels(labels >= 0));
    if ~isempty(survivingLbls)
        newLbls     = -ones(size(labels));
        newTraining = -ones(size(trainingLbl));
        for k = 1:numel(survivingLbls)
            newLbls(labels == survivingLbls(k))           = k;
            newTraining(trainingLbl == survivingLbls(k))  = k;
        end
        labels      = newLbls;
        trainingLbl = newTraining;
    if ~isempty(sinkClassLabels)
        [~, sinkLoc] = ismember(sinkClassLabels, survivingLbls);
        sinkClassLabels = sinkLoc(sinkLoc > 0);
    end

        uniqueLabels        = unique(labels);
        numUniqueClusters   = length(uniqueLabels);
        class_polarity      = zeros(numUniqueClusters, 1);
        clusterSampleCounts = zeros(numUniqueClusters, 1);
        meanClusterWaveform = zeros([numUniqueClusters, size(waveform, [2, 3])]);
        meanClusterWaveformFull = zeros([numUniqueClusters, size(waveform_bp_full, [2, 3])]);
        for ii = 1:numUniqueClusters
            cMask = labels == uniqueLabels(ii);
            clusterSampleCounts(ii)     = sum(cMask);
            meanClusterWaveform(ii,:,:) = mean(waveform(cMask, :, :), 1, 'omitmissing');
            meanClusterWaveformFull(ii,:,:) = mean(waveform_bp_full(cMask, :, :), 1, 'omitmissing');
            if uniqueLabels(ii) >= 0
                class_polarity(ii, 1) = mode(spk_ID_full(cMask));
            end
        end
    end
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

    % A sink spans heterogeneous noise; the full variant count would give it
    % enough reach to win matches against real units.
    if ismember(uniqueLabels(i), sinkClassLabels)
        nT = max(1, min(numTemplatesPerCluster, cfgGet(cfg,'noiseSinkTemplates',3)));
        capPts = cfgGet(cfg,'noiseSinkMaxTemplatePts',5000);
        if numel(idx) > capPts
            idx = idx(round(linspace(1, numel(idx), capPts)));
        end
    else
        nT = numTemplatesPerCluster;
    end

    clusterWaveforms = waveform(idx,:,:);
    clusterPCAscores = PCA_score(idx,:);

    if length(idx) >= nT
        [templates, weights] = generateMultipleTemplates(clusterWaveforms, nT, clusterPCAscores);
        templateWaveforms(i,1:nT,:,:) = templates;
        templateWeights(i,1:nT) = weights;
    else
        meanWF = reshape(meanClusterWaveform(i,:,:), size(meanClusterWaveform,2), []);
        for t = 1:nT
            templateWaveforms(i,t,:,:) = meanWF;
        end
        templateWeights(i,1:nT) = 1/nT;
    end
end

[clusterSpikeDensity, ~, ~] = cluster_spike_density(spk_idx_full, labels, cfg);

sinkMask = ismember(uniqueLabels, sinkClassLabels);
[clusterRelabeling]  = kiaSort_process_clusters(meanClusterWaveform, clusterSpikeDensity, ...
    uniqueLabels, labels, PCA, spk_idx_full, mean_side_waveforms, cfg, sinkMask);

clusterRelabeling.originalLabels = uniqueLabels;
clusterRelabeling.mean_side_waveforms = mean_side_waveforms;
[clusterRelabeling] = realign_merge_Waveforms(meanClusterWaveform, meanClusterWaveformFull, clusterSampleCounts, clusterRelabeling);
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
isNoiseSink = false(length(uniqueNewLabels), 1);
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
            % Same rule as the iterative reject: only discard when the class
            % clears the keep gate and its polarity has another class left.
            ampLow = isfield(clusterSelection, 'lowAmpNotKept') && ...
                     numel(clusterSelection.lowAmpNotKept) >= iMerged && ...
                     clusterSelection.lowAmpNotKept(iMerged);
            polSelf = clusterSelection.mainNegativePolarity(iMerged);
            nPol    = sum(clusterSelection.mainNegativePolarity(:) == polSelf & ...
                          uniqueNewLabels(:) ~= -1);
            if ~ampLow && nPol >= 2
                clusterStatus(iMerged) = -4;
                clusterSelection.keep(iMerged) = 0;
            end
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

% A sink exists to absorb noise, never to be reported as a unit.
if ~isempty(sinkClassLabels)
    sinkRows = find(ismember(uniqueLabels, sinkClassLabels));
    for sR = 1:numel(sinkRows)
        idxS = find(uniqueNewLabels == clusterRelabeling.newLabels(sinkRows(sR)), 1);
        if isempty(idxS), continue; end
        isNoiseSink(idxS) = true;
        clusterSelection.keep(idxS) = 0;
        if clusterStatus(idxS) == 0
            if class_polarity(sinkRows(sR)) == -1
                clusterStatus(idxS) = -3;
            else
                clusterStatus(idxS) = -2;
            end
        end
    end
end
clusterSelection.clusterStatus = clusterStatus;
clusterSelection.isNoiseSink = isNoiseSink;
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
    lowPrc  = prctile(ampVals(class_idx), 1);
    highPrc = prctile(ampVals(class_idx), 99);
    trainingLbl(labels==uniqueLabels(i) & (ampVals < lowPrc | ampVals > highPrc)) = -1;

    % The sink owns the band below sinkSnr; the generic widening would let
    % it claim spikes well above its own distribution.
    if ismember(uniqueLabels(i), sinkClassLabels)
        low_thr(i) = 0;
        if class_polarity(i) == 1
            high_thr(i) =  sinkSnr * thrPosMid;
        else
            high_thr(i) = -sinkSnr * thrNegMid;
        end
    elseif class_polarity(i) == 1
        low_thr(i) = max(ampBandLo * lowPrc, 0);
        if isfinite(ampBandHi)
            high_thr(i) = max(ampBandHi * highPrc, 0);
        else
            high_thr(i) = Inf;      % no upper cap for a positive class
        end
    else
        low_thr(i) = min(ampBandLo * highPrc, 0);
        if isfinite(ampBandHi)
            high_thr(i) = min(ampBandHi * lowPrc, 0);
        else
            high_thr(i) = -Inf;     % no upper cap for a negative class
        end
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
out.clusteringInfo.bimodalSplitLog  = splitLog;
out.clusteringInfo.noiseSinkLabels  = sinkClassLabels;

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

function [labels, globalCluster, splitLog] = processCluster(idx, labels, fs, globalCluster, dataAll, spk_idx_full, cfg, depth, maxClusterPoints, minClusterPoints, sample_dur, numPt, initNumClasses, snr_vals, splitLog)
[~,~,isi_viol] = getISIViolations(spk_idx_full(idx), fs, 2);
factor = length(idx)/size(dataAll,1);
if ((isi_viol<= 0.1 && depth >= 0)  && (initNumClasses>10 || depth > 0)) || depth >= 3
    [parts, splitLog] = splitIfBimodal(idx, dataAll, spk_idx_full, fs, isi_viol, sample_dur, cfg, depth, splitLog);
    for p = 1:numel(parts)
        globalCluster = globalCluster + 1;
        labels(parts{p}) = globalCluster;
    end
elseif isi_viol> 1 && length(idx) < minClusterPoints/2 && depth > 1
    % Provisional reject. Whether it really becomes noise depends on its
    % amplitude and on how many classes its polarity has, neither of which
    % is known until the recursion finishes, so tag it and resolve later.
    labels(idx) = min([-1; labels(labels <= -2)]) - 1;
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
            [labels, globalCluster, splitLog] = processCluster(subIdx, labels, fs, globalCluster, dataAll, spk_idx_full, cfg, depth+1, maxClusterPoints, minClusterPoints, sample_dur, numPt, initNumClasses, snr_vals, splitLog);
        end
    end
end
end

function noise_idx = identifyOutliers(data)
% Restored on exit so the suppression does not leak into the session.
wsNear = warning('off', 'MATLAB:nearlySingularMatrix');
wsSing = warning('off', 'MATLAB:singularMatrix');
restoreWarn = onCleanup(@() warning([wsNear, wsSing])); %#ok<NASGU>
try
    mu = mean(data, 1);
    sigma = cov(data);
    % Ridge scaled to the data, then used directly: mahal() recomputes the
    % covariance itself, so calling it here would discard the regularised
    % sigma and invert the singular one. A cluster with fewer spikes than
    % components is rank-deficient by construction.
    ridge = 1e-6 * mean(diag(sigma));
    if ~isfinite(ridge) || ridge <= 0
        ridge = 1e-6;
    end
    if rcond(sigma) < 1e-10
        sigma = sigma + eye(size(sigma)) * ridge;
    end
    d = data - mu;
    mahalDist = sum((d / sigma) .* d, 2);
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



function v = cfgGet(cfg, name, dflt)
if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = dflt;
end
end


function [parts, splitLog] = splitIfBimodal(idx, dataAll, spk_idx_full, fs, isiParent, sample_dur, cfg, depth, splitLog)
% A cluster with clean ISI can still hold two neurons: a merged train of two
% cells that rarely co-fire looks refractory. Split only on a clear density
% valley, and only when the two sides overlap in time -- a drifting single
% unit also looks bimodal but its halves are consecutive, not interleaved.
parts = {idx};
if ~cfgGet(cfg, 'enableBimodalSplit', false), return; end
if ~isempty(splitLog) && sum([splitLog.split]) >= cfgGet(cfg, 'bimodalMaxSplits', 1), return; end

n        = numel(idx);
minUnit  = max(20, round(sample_dur * cfgGet(cfg, 'minRate', 0.15)));
minChild = max(minUnit, ceil(cfgGet(cfg, 'bimodalMinChildFrac', 0.20) * n));
if n < 2*minChild, return; end

% UMAP block only; the trailing 3 PCA columns are amplitude-dominated and
% carry the drift axis, which would cut a drifting unit in half.
dU = size(dataAll,2) - 3;
if dU < 2, return; end

maxTest = cfgGet(cfg, 'bimodalMaxTestPoints', 5000);
if n > maxTest
    sub = round(linspace(1, n, maxTest))';
else
    sub = (1:n)';
end
U  = dataAll(idx, 1:dU);
mu = mean(U(sub,:), 1);
Xc = U(sub,:) - mu;

[V, D] = eig(Xc'*Xc);
[~, k] = max(diag(D));
w = V(:,k);
t = Xc * w;

[thr, sep, valley, eta] = kiaSort_otsu_split1d(t, cfgGet(cfg, 'bimodalNumBins', 64));
rec = struct('n', n, 'depth', depth, 'sep', sep, 'valley', valley, 'eta', eta, ...
             'timeOverlap', NaN, 'medShift', NaN, 'isiParent', isiParent, ...
             'isiChild', [NaN NaN], 'split', false);

if isnan(thr) || eta < cfgGet(cfg,'bimodalSeparability',0.75)
    splitLog = appendSplitLog(splitLog, rec); return;
end

side = ((U - mu) * w) <= thr;
i1 = idx(side);  i2 = idx(~side);
if numel(i1) < minChild || numel(i2) < minChild
    splitLog = appendSplitLog(splitLog, rec); return;
end

t1 = double(spk_idx_full(i1));  t2 = double(spk_idx_full(i2));
q1 = prctile(t1, [5 95]);       q2 = prctile(t2, [5 95]);
span = min(q1(2)-q1(1), q2(2)-q2(1));
if span <= 0
    splitLog = appendSplitLog(splitLog, rec); return;
end
rec.timeOverlap = (min(q1(2),q2(2)) - max(q1(1),q2(1))) / span;
rec.medShift    = abs(median(t1) - median(t2)) / max(1, sample_dur*fs);
if rec.timeOverlap < cfgGet(cfg,'bimodalMinTimeOverlap',0.5) || ...
   rec.medShift    > cfgGet(cfg,'bimodalMaxMedianShift',0.3)
    splitLog = appendSplitLog(splitLog, rec); return;
end

isiCap = min(max(0.1, isiParent), 1);
[~,~,v1] = getISIViolations(spk_idx_full(i1), fs, 2);
[~,~,v2] = getISIViolations(spk_idx_full(i2), fs, 2);
rec.isiChild = [v1 v2];
if v1 > isiCap || v2 > isiCap
    splitLog = appendSplitLog(splitLog, rec); return;
end

rec.split = true;
splitLog  = appendSplitLog(splitLog, rec);
parts = {i1, i2};
end


function splitLog = appendSplitLog(splitLog, rec)
if isempty(splitLog)
    splitLog = rec;
else
    splitLog(end+1) = rec;
end
end
