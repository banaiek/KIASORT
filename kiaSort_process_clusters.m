function [out] = kiaSort_process_clusters(meanWaveform, spikeDensity, ...
    initLabels, labels, pca, spike_idx, mean_side_waveforms, cfg)

sample_points = cfg.numSampleChunks * cfg.sampleChunkDuration * cfg.samplingFrequency;
S = size(meanWaveform, 1);
N = size(meanWaveform, 2);
M = size(meanWaveform, 3);
midChannel = round(N/2);
spike_distance = min(round(M/2)-3, round(cfg.spikeDistance .* cfg.samplingFrequency/1000));

% ---- Cluster counts ----
counts = zeros(S,1);
for i = 1:S
    if initLabels(i) == -1
        counts(i) = -Inf;
    else
        counts(i) = sum(labels == initLabels(i));
    end
end

% ---- Overlapping waveform detection ----
merge_overlap = zeros(S,S);
if ~isempty(mean_side_waveforms)
    merge_overlap = overlapping_waveform(meanWaveform, .2, counts, .05, 1, spike_distance+1, 4, mean_side_waveforms);
end

counts_merge = counts < 0.05*counts' | counts' < 0.05*counts;

% ---- Thresholds ----
corrThresholdMisAlign   = cfg.corrThresholdMisAlign;
refrac_threshold        = cfg.refrac_threshold;

thresholdISI   = 2 - cfg.spikeDistance;
isiViolation1  = .2;
thresholdISI_ms1 = 1;
isiViolation2  = .5;
thresholdISI_ms2 = 2;

fullLength  = N * M;
validLength = nnz(mean(meanWaveform, 'omitmissing') ~= 0);

adj_corrThresholdMisAlign = adjust_correlation_threshold(fullLength, validLength, corrThresholdMisAlign);

ampVarThreshold         = cfg.ampVarThreshold;
ampXcorrVarThreshold    = cfg.ampXcorrVarThreshold;
corrThresholdSpkDensity = cfg.corrThresholdSpkDensity;

% ---- Peak-based amplitude features (vectorized setup) ----
maxWave = zeros(size(meanWaveform));
halfBand = spike_distance + 2;
lo = round(M/2) - halfBand;
hi = round(M/2) + halfBand;
maxWave(:,:,lo:hi) = meanWaveform(:,:,lo:hi);

nonPeaks = ~(maxWave>0 & maxWave-circshift(maxWave,1,3)>0 & maxWave-circshift(maxWave,-1,3)>0) ...
         & ~(maxWave<0 & maxWave-circshift(maxWave,1,3)<0 & maxWave-circshift(maxWave,-1,3)<0);
maxWave(nonPeaks) = 0;

maxAmp  = max(abs(maxWave), [], 3);
maxAmpP = abs(max(maxWave, [], 3));
maxAmpN = abs(min(maxWave, [], 3));
maxAmpP2 = squeeze(mean(abs(maxk(maxWave,2,3)), 3, 'omitmissing'));
maxAmpN2 = squeeze(mean(abs(mink(maxWave,2,3)), 3, 'omitmissing'));

midV  = maxAmp(:, round(N/2));
ranks = sum(maxAmp <= midV, 2);
rank_merge = ~(ranks < N-3 & ranks' < N-3);

[~, maxMidPoint] = max(squeeze(abs(meanWaveform(:,:,round(M/2)))), [], 2);
similarKept = (maxMidPoint == round(M/2)) == (maxMidPoint == round(M/2))';

% ---- Amplitude variance (pairwise) ----
ampVar      = nan(S,S);
ampMedVar   = zeros(S,S);
ampVar2     = nan(S,S);
ampMedVar2  = zeros(S,S);
ampVarDrift = zeros(S,S);
midChannelVar = zeros(S,S);

for i = 1:S
    for j = i+1:S
        maxAmpPair = max([maxAmp(i,:), maxAmp(j,:)], [], 'all');

        ampDiffP = (maxAmpP(i,:) - maxAmpP(j,:)) ./ maxAmpPair;
        ampDiffN = (maxAmpN(i,:) - maxAmpN(j,:)) ./ maxAmpPair;
        ampDiff  = abs(ampDiffP + ampDiffN) / 2;

        ampDiffP2 = (maxAmpP2(i,:) - maxAmpP2(j,:)) ./ maxAmpPair;
        ampDiffN2 = (maxAmpN2(i,:) - maxAmpN2(j,:)) ./ maxAmpPair;
        ampDiff2  = abs(ampDiffP2 + ampDiffN2) / 2;

        signedAmp = ampDiffP + ampDiffN;
        ampVarDrift(i,j) = mean(signedAmp(signedAmp~=0), 'omitmissing');
        ampVarDrift(j,i) = ampVarDrift(i,j);

        nz = ampDiff ~= 0;
        ampVar(i,j)    = mean(ampDiff(nz), 'omitmissing');
        ampVar(j,i)    = ampVar(i,j);
        ampMedVar(i,j) = median(ampDiff(nz), 'omitmissing');
        ampMedVar(j,i) = ampMedVar(i,j);

        nz2 = ampDiff2 ~= 0;
        ampVar2(i,j)    = mean(ampDiff2(nz2), 'omitmissing');
        ampVar2(j,i)    = ampVar2(i,j);
        ampMedVar2(i,j) = median(ampDiff2(nz2), 'omitmissing');
        ampMedVar2(j,i) = ampMedVar2(i,j);

        if ampDiff2(midChannel) < .1 && ampDiff(midChannel) < .1
            midChannelVar(i,j) = 1;
            midChannelVar(j,i) = 1;
        end
    end
end

% ---- Cross-correlation (single pass — unnormalized only) ----
flatWaves = zeros(S, N*M);
for sIdx = 1:S
    flatWaves(sIdx,:) = reshape(squeeze(meanWaveform(sIdx,:,:)), 1, []);
end

maxXcorrVal = zeros(S,S);
maxXcorrLag = zeros(S,S);

for i = 1:S
    for j = i+1:S
        [bestCorr, bestLag] = max_half_corr(flatWaves(i,:), flatWaves(j,:), N, M, floor(M/3), 0);
        maxXcorrVal(i,j) = bestCorr;
        maxXcorrVal(j,i) = bestCorr;
        maxXcorrLag(i,j) = bestLag;
        maxXcorrLag(j,i) = -bestLag;
    end
end

% ---- Spike density analysis ----
density_length = size(spikeDensity, 2);
spikeDensityCorrVal = zeros(S,S);
startIdx = ones(S,1);
endIdx   = density_length * ones(S,1);
perRatio = zeros(S,1);
adj_density = spikeDensity ./ max(spikeDensity, [], 2);

for i = 1:S
    [startIdx(i), endIdx(i)] = findStableInterval(spikeDensity(i,:));
    adj_density(i, 1:startIdx(i)) = 0;
    adj_density(i, endIdx(i):end) = 0;
end

% ---- ISI violations & density correlation (combined pair loop) ----
isi_violation    = zeros(S,S);
isi_violation_1ms = zeros(S,S);
isi_violation_2ms = zeros(S,S);

for i = 1:S
    idx_i = (labels == initLabels(i)) ...
          & (spike_idx <= endIdx(i)/density_length * sample_points) ...
          & (spike_idx >= startIdx(i)/density_length * sample_points);
    si = spike_idx(idx_i);

    for j = i:S
        idx_j = (labels == initLabels(j)) ...
              & (spike_idx <= endIdx(j)/density_length * sample_points) ...
              & (spike_idx >= startIdx(j)/density_length * sample_points);

        if i == j
            spike_idx_pair = si;
        else
            spike_idx_pair = [si; spike_idx(idx_j)];
        end

        [~,~,viol] = getISIViolations(spike_idx_pair, cfg.samplingFrequency, thresholdISI);
        isi_violation(i,j) = viol;
        isi_violation(j,i) = viol;

        if i ~= j
            [~,~,v1] = getISIViolations(spike_idx_pair, cfg.samplingFrequency, thresholdISI_ms1);
            isi_violation_1ms(i,j) = v1;  isi_violation_1ms(j,i) = v1;

            [~,~,v2] = getISIViolations(spike_idx_pair, cfg.samplingFrequency, thresholdISI_ms2);
            isi_violation_2ms(i,j) = v2;  isi_violation_2ms(j,i) = v2;

            % --- Density correlation (merged into same loop) ---
            c1 = corr(spikeDensity(i,:)', spikeDensity(j,:)', 'Type','Spearman');
            den_i = adj_density(i,:)';
            den_j = adj_density(j,:)';
            c2 = corr(den_i, den_j, 'Type','Spearman');
            mask = den_i~=0 | den_j~=0;
            if any(mask)
                c3 = corr(den_i(mask), den_j(mask), 'Type','Spearman');
            else
                c3 = 0;
            end
            c = min([c1,c2,c3]);
            spikeDensityCorrVal(i,j) = c;
            spikeDensityCorrVal(j,i) = c;
        end
    end
end

isi_pass = (isi_violation_1ms < isiViolation1) & (isi_violation_2ms < isiViolation2);

% ---- Merge criteria ----
pr_class     = (endIdx - startIdx) / size(spikeDensity,2);
full_classes  = pr_class./pr_class' > .25 | pr_class'./pr_class > .25;
merge_feat    = kiaSort_merge_features(labels, pca, cfg);
merge_overlap = (merge_overlap + merge_overlap' > 0);

XcorrMerge       = maxXcorrVal > adj_corrThresholdMisAlign & abs(maxXcorrLag) > 2 & abs(maxXcorrLag) <= 1.5*spike_distance;
drift_based_Xcorr = maxXcorrVal > adj_corrThresholdMisAlign & spikeDensityCorrVal < corrThresholdSpkDensity;
ampMergeXcorr    = ampVar < ampXcorrVarThreshold & ampVar2 < ampXcorrVarThreshold;
ampMergeFeat     = ampMedVar < ampXcorrVarThreshold/2 & ampMedVar2 < ampXcorrVarThreshold/2;
ameMergeDrift    = ampVar < ampXcorrVarThreshold & ampVar2 < ampXcorrVarThreshold ...
                 & ampVarDrift < ampXcorrVarThreshold & ampVarDrift < .75*ampVar;

d = diag(isi_violation);
isi_merge = (isi_violation - max(d, d') < refrac_threshold);
isRefr    = refractoryMatrix(spike_idx, labels, cfg.samplingFrequency, initLabels);

merge_matrix = multiCond_Merge(meanWaveform);

merge1 = isi_pass & (~isRefr | isi_merge) & (merge_matrix | (XcorrMerge & ampMergeXcorr));
merge2 = isi_pass & rank_merge & (~isRefr | isi_merge) ...
       & ((merge_overlap & similarKept & midChannelVar) ...
        | (merge_feat & ampMergeFeat & counts_merge & similarKept & full_classes & midChannelVar) ...
        | (drift_based_Xcorr & ameMergeDrift));

[newLabels, changeType, timeLagChanged] = kiaSort_merge_cluster(merge1, merge2, maxXcorrLag, initLabels, counts, 1-maxXcorrVal);

for i = 1:S
    searchIdx = newLabels == newLabels(i);
    perRatio(i) = (max(endIdx(searchIdx)) - min(startIdx(searchIdx))) / size(spikeDensity,2);
end

% ---- Output ----
out.perRatio            = perRatio;
out.newLabels           = newLabels;
out.timeLagChanged      = timeLagChanged;
out.maxXcorrVal         = maxXcorrVal;
out.maxXcorrLag         = maxXcorrLag;
out.changeType          = changeType;
out.clusterSpikeDensity = spikeDensity;
out.stablePoints        = [startIdx, endIdx] .* sample_points/density_length;
end

%% ========================================================================
function [startIdx, endIdx] = findStableInterval(d)
dropRate = 10;
n = numel(d);
if n < 3
    startIdx = 1; endIdx = n;
    return
end
cp = findchangepts(d, 'Statistic','rms', 'MinThreshold',10);

startIdx = 1; endIdx = n;
if isscalar(cp)
    if median(d(1:cp)) < median(d(cp+1:end))/dropRate
        startIdx = cp+1;
    elseif median(d(1:cp))/dropRate > median(d(cp+1:end))
        endIdx = cp;
    end
elseif numel(cp) >= 2
    for i = 1:numel(cp)
        if mean(d(1:cp(i))) < mean(d(cp(i)+1:end))/dropRate ...
                && mean(d(1:cp(i))) < mean(d(cp(i)+1:min(2*cp(i),n)))/dropRate
            startIdx = cp(i)+1;
        elseif mean(d(startIdx:cp(i)))/dropRate > mean(d(cp(i)+1:end)) ...
                && mean(d(startIdx:cp(i)))/dropRate > mean(d(cp(i)+1:min(2*cp(i),n)))
            endIdx = cp(i);
            break
        end
    end
end
end