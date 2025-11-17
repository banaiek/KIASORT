function [out] = kiaSort_process_clusters(meanWaveform, spikeDensity, ...
    initLabels, labels, pca, spike_idx, mean_side_waveforms, cfg)


sample_points = cfg.numSampleChunks * cfg.sampleChunkDuration * cfg.samplingFrequency;

spike_distance = round(cfg.spikeDistance .* cfg.samplingFrequency/1000);

S = size(meanWaveform, 1);
N = size(meanWaveform, 2);
M = size(meanWaveform, 3);
midChannel = round(N/2);

counts = zeros(S,1);
for i = 1:S
    if initLabels(i) == -1
        counts(i) = -Inf;
    else
        counts(i) = sum(labels == initLabels(i));
    end
end

merge_overlap = zeros(S,S);

% detecting overlapping waveforms

if ~isempty(mean_side_waveforms)
    merge_overlap = overlapping_waveform(meanWaveform, .4, counts, .05, 1, spike_distance+1, 4, mean_side_waveforms);
end

counts_merge = counts< 0.05*counts' | counts'< 0.05*counts;

corrThresholdHigh       = cfg.corrThresholdHigh;
corrThresholdMisAlign   = cfg.corrThresholdMisAlign;
refrac_threshold        = cfg.refrac_threshold;

thresholdISI = 2-cfg.spikeDistance;


fullLength = N*M;
validLength = numel(find(mean(meanWaveform,'omitmissing')~=0));

adj_corrThresholdHigh     = adjust_correlation_threshold(fullLength, validLength, corrThresholdHigh);
adj_corrThresholdMisAlign = adjust_correlation_threshold(fullLength, validLength, corrThresholdMisAlign);

ampVarThreshold         = cfg.ampVarThreshold;
ampXcorrVarThreshold    = cfg.ampXcorrVarThreshold;
corrThresholdSpkDensity = cfg.corrThresholdSpkDensity;

ampVar         = nan(S,S);
ampMedVar  = zeros(S, S);
ampVar2         = nan(S,S);
ampMedVar2  = zeros(S, S);
ampVarDrift = zeros(S, S);
midChannelVar = zeros(S, S);

maxXcorrVal = zeros(S, S);
maxXcorrLag = zeros(S, S);
isi_violation = zeros(S, S);
maxXcorrValNorm = zeros(S, S);
maxXcorrLagNorm = zeros(S, S);

[~, maxMidPoint] = max(squeeze(abs(meanWaveform(:,:,round(M/2)))),[],2);
similarKept = (maxMidPoint==round(M/2))==(maxMidPoint==round(M/2))';

% calculate amp. variation
maxWave = zeros(size(meanWaveform));
maxWave(:,:,round(M/2)-spike_distance-2:round(M/2)+spike_distance+2) = meanWaveform(:,:,round(M/2)-spike_distance-2:round(M/2)+spike_distance+2);
nonPeaks = ~(maxWave>0 & maxWave-circshift(maxWave,1,3)>0 & maxWave-circshift(maxWave,-1,3)>0) & ~(maxWave<0 & maxWave-circshift(maxWave,1,3)<0 & maxWave-circshift(maxWave,-1,3)<0);
maxWave(nonPeaks)=0;

maxAmp = max(abs(maxWave),[],3);
maxAmpP = abs(max((maxWave),[],3));
maxAmpN = abs(min((maxWave),[],3));
maxAmpP2 = squeeze(mean(abs(maxk((maxWave),2,3)),3,'omitmissing'));
maxAmpN2 = squeeze(mean(abs(mink((maxWave),2,3)),3,'omitmissing'));


midV  = maxAmp(:, round(N/2));
ranks = sum(maxAmp <= midV, 2);
rank_merge = ~(ranks < N-3 & ranks' < N-3);

% calculate bi and tri-phasic amplitude variation metric
for i = 1:S
    for j = i+1:S
        maxAmpPair = (max([maxAmp(i,:) , maxAmp(j,:)],[],'all'));

        ampDiffP = (maxAmpP(i,:) - maxAmpP(j,:))./maxAmpPair;
        ampDiffN = (maxAmpN(i,:) - maxAmpN(j,:))./maxAmpPair;
        ampDiff = abs(ampDiffP + ampDiffN)/2;

        ampDiffP2 = (maxAmpP2(i,:) - maxAmpP2(j,:))./maxAmpPair;
        ampDiffN2 = (maxAmpN2(i,:) - maxAmpN2(j,:))./maxAmpPair;
        ampDiff2 = abs(ampDiffP2 + ampDiffN2)/2;



        % validChan = (maxAmpN(i,:)>.1*maxAmpPair | maxAmpN(j,:)>.1*maxAmpPair);
        % 
        % % if sum(validChan) > 1 
        %     validChan = ones(size(validChan));
        % % end

        signedAmp = ampDiffP + ampDiffN;

        ampVarDrift(i,j) = mean(signedAmp(signedAmp~=0),'omitmissing');
        ampVarDrift(j,i) = ampVarDrift(i,j);

        ampVar(i,j)  = mean(ampDiff(ampDiff~=0 ),'omitmissing');
        ampVar(j,i)  = ampVar(i,j);
        ampMedVar(i,j)  = median(ampDiff(ampDiff~=0 ),'omitmissing');
        ampMedVar(j,i)  = ampMedVar(i,j);

        ampVar2(i,j)  = mean(ampDiff2(ampDiff2~=0 ),'omitmissing');
        ampVar2(j,i)  = ampVar2(i,j);
        ampMedVar2(i,j)  = median(ampDiff2(ampDiff2~=0 ),'omitmissing');
        ampMedVar2(j,i)  = ampMedVar2(i,j);

        if ampDiff2(midChannel) <.1 &&  ampDiff(midChannel) <.1
            midChannelVar(i, j) = 1;
            midChannelVar(j, i) = 1; 
        end


    end
end


% Flatten waveforms and compute pairwise max XCorr
flatWaves = zeros(S, N*M);
for sIdx = 1:S
    wave2d = squeeze(meanWaveform(sIdx, :, :));
    flatWaves(sIdx,:) = wave2d(:)';
end

% Compute the pairwise maximum cross-corr & best lag
for i = 1:S
    for j = i+1:S
        [bestCorr, bestLag] = max_half_corr(flatWaves(i,:), flatWaves(j,:), N, M, floor(M/3), 0);
        maxXcorrVal(i,j) = bestCorr;
        maxXcorrVal(j,i) = bestCorr;
        maxXcorrLag(i,j) = bestLag;
        maxXcorrLag(j,i) = -bestLag;

        [bestCorr, bestLag] = max_half_corr(flatWaves(i,:), flatWaves(j,:), N, M, floor(M/3), 1);
        maxXcorrValNorm(i,j) = bestCorr;
        maxXcorrValNorm(j,i) = bestCorr;
        maxXcorrLagNorm(i,j) = bestLag;
        maxXcorrLagNorm(j,i) = -bestLag;

    end
end

density_length = size(spikeDensity,2);
spikeDensityCorrVal = zeros(S,S);
startIdx = ones(S,1);
endIdx = density_length * ones(S,1);
perRatio = zeros(S,1);
adj_density = spikeDensity./max(spikeDensity,[],2);
for i = 1:S
    [startIdx(i), endIdx(i)] = findStableInterval(spikeDensity(i,:));
    adj_density(i,1:startIdx(i)) = 0;
    adj_density(i,endIdx(i):end) = 0;
end

% get isi violation for stable duration

for i = 1:S
    for j = i:S
        if i == j
            idx_i = (labels == initLabels(i)) & (spike_idx <= (endIdx(i)/density_length * sample_points)) & (spike_idx >= (startIdx(i)/density_length * sample_points));
            spike_idx_pair = spike_idx(idx_i);
            [~, ~, isi_viol] = getISIViolations(spike_idx_pair, cfg.samplingFrequency, thresholdISI);
            isi_violation (i, j) = isi_viol;

        else
            idx_i = (labels == initLabels(i)) & (spike_idx <= (endIdx(i)/density_length * sample_points)) & (spike_idx >= (startIdx(i)/density_length * sample_points));
            idx_j = (labels == initLabels(j)) & (spike_idx <= (endIdx(j)/density_length * sample_points)) & (spike_idx >= (startIdx(j)/density_length * sample_points));
            spike_idx_i = spike_idx(idx_i);
            spike_idx_j = spike_idx(idx_j);
            spike_idx_pair = [spike_idx_i; spike_idx_j];
            [~, ~, isi_viol] = getISIViolations(spike_idx_pair, cfg.samplingFrequency, thresholdISI);
            isi_violation (i, j) = isi_viol;
            isi_violation (j, i) = isi_violation (i, j);
        end

    end
end

% get spike density correlations for drift based merging
for i = 1:S
    for j = i+1:S
        c1 = corr(spikeDensity(i,:)', spikeDensity(j,:)', 'Type','Spearman');
        den_i = adj_density(i,:)';
        den_j = adj_density(j,:)';
        c2 = corr(den_i, den_j, 'Type','Spearman');
        if any(den_i~=0 | den_j~=0)
            c3 = corr(den_i(den_i~=0 | den_j~=0), den_j(den_i~=0 | den_j~=0), 'Type','Spearman');

        else
            c3 = 0;
        end
        c = min([c1,c2,c3]);
        spikeDensityCorrVal(i,j) = c;
        spikeDensityCorrVal(j,i) = c;
    end
end

pr_class = (endIdx-startIdx)/size(spikeDensity,2);
full_classes = pr_class./pr_class'>.25 | pr_class'./pr_class>.25;
merge_feat = kiaSort_merge_features(labels, pca, cfg);
merge_overlap = (merge_overlap+merge_overlap'>0);
corrMerge = (maxXcorrVal > adj_corrThresholdHigh);
XcorrMerge = (maxXcorrVal > adj_corrThresholdMisAlign & abs(maxXcorrLag) > 2 & abs(maxXcorrLag) <= 1.5 * spike_distance);
drift_based_Xcorr = maxXcorrVal > adj_corrThresholdMisAlign & spikeDensityCorrVal < corrThresholdSpkDensity;
ampMerge = ampVar < ampVarThreshold & midChannelVar & ampVar2 < ampVarThreshold;
ampMergeXcorr = ampVar < ampXcorrVarThreshold & ampVar2 < ampXcorrVarThreshold;
ampMergeFeat = ampMedVar < ampXcorrVarThreshold/2 & ampMedVar2 < ampXcorrVarThreshold/2;
ameMergeDrift = ampVar < ampXcorrVarThreshold & ampVar2 < ampXcorrVarThreshold & ampVarDrift < ampXcorrVarThreshold & ampVarDrift < .75 * ampVar;
d = diag(isi_violation);
isi_merge = (isi_violation - max(d, d') < refrac_threshold);
isRefr = refractoryMatrix(spike_idx, labels, cfg.samplingFrequency, initLabels);

merge_matrix = multiCond_Merge(meanWaveform);
merge1 =  (~isRefr | isi_merge) & (merge_matrix |(XcorrMerge & ampMergeXcorr));
merge2 = rank_merge & (~isRefr | isi_merge)  & ((merge_overlap  & similarKept & midChannelVar) | (merge_feat & ampMergeFeat & counts_merge  & similarKept & full_classes & midChannelVar) | (drift_based_Xcorr & ameMergeDrift));

[newLabels, changeType, timeLagChanged] = kiaSort_merge_cluster(merge1, merge2, maxXcorrLag, initLabels, counts, 1-maxXcorrVal);

for i = 1:S
    searchIdx = newLabels == newLabels(i);
    perRatio(i) = (max(endIdx(searchIdx))-min(startIdx(searchIdx)))/size(spikeDensity,2);
end

out.perRatio            = perRatio;
out.newLabels           = newLabels;
out.timeLagChanged      = timeLagChanged;
out.maxXcorrVal         = maxXcorrVal;
out.maxXcorrLag         = maxXcorrLag;
out.changeType          = changeType;
out.clusterSpikeDensity = spikeDensity;
out.stablePoints        = [startIdx, endIdx] .* sample_points/density_length;
end

function [startIdx, endIdx] = findStableInterval(d)
dropRate = 10;
n = numel(d);
if n < 3
    startIdx = 1; endIdx = n;
    return
end
cp = findchangepts(d, 'Statistic', 'rms','MinThreshold',10);

startIdx = 1; endIdx = n;
if isscalar(cp)
    if median(d(1:cp)) < median(d(cp+1:end))/dropRate
        startIdx = cp+1; endIdx = n;
    elseif median(d(1:cp))/dropRate > median(d(cp+1:end))
        startIdx = 1; endIdx = cp;
    end
elseif numel(cp)>=2
    for i = 1:numel(cp)
        if mean(d(1:cp(i))) < mean(d(cp(i)+1:end))/dropRate && mean(d(1:cp(i))) < mean(d(cp(i)+1:min(2*cp(i),n)))/dropRate
            startIdx = cp(i)+1;
        elseif mean(d(startIdx:cp(i)))/dropRate> mean(d(cp(i)+1:end)) && mean(d(startIdx:cp(i)))/dropRate > mean(d(cp(i)+1:min(2*cp(i),n)))
            endIdx = cp(i);
            break
        end
    end
end
end

