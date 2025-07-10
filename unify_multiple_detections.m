function [sampleFeatures, sortedSamples] = unify_multiple_detections(cfg, sampleFeatures, sortedSamples)
middleChannelIdx = cfg.num_channel_extract + 1;
searchChannel = round(cfg.num_channel_extract);
spikeDistance = 1.5 * cfg.spikeDistance * cfg.samplingFrequency/1000;
overlap_thr = cfg.overlap_thr;
numChannels = length(sortedSamples);
clusters = [];
for i = 1:numChannels
    if ~isempty(sampleFeatures{i,1})
        spk_idx = sampleFeatures{i,1}.spk_idx;
        spk_lbl = sampleFeatures{i,1}.updatedLabels;
        clust_sel = sortedSamples{i,1}.clusteringInfo.clusterSelection;
        keep_idx = find(clust_sel.keep);
        for k = 1:length(keep_idx)
            idx = keep_idx(k);
            lbl = clust_sel.classLabels(idx);
            spk = spk_idx(spk_lbl == lbl);
            clusters(end+1).channel = i;
            clusters(end).local_idx = idx;
            clusters(end).spikes = spk;
        end
    end
end
nClusters = length(clusters);
parent = 1:nClusters;
function r = find_parent(x)
    while parent(x) ~= x
        parent(x) = parent(parent(x));
        x = parent(x);
    end
    r = x;
end
function union_set(a, b)
    ra = find_parent(a);
    rb = find_parent(b);
    if ra ~= rb
        parent(rb) = ra;
    end
end
for p = 1:nClusters-1
    for q = p+1:nClusters
        if (clusters(q).channel - clusters(p).channel) <= searchChannel
            spk_p = clusters(p).spikes;
            spk_q = clusters(q).spikes;
            if isempty(spk_p) || isempty(spk_q)
                continue;
            end
            [d_pq, ~] = nearest_index(spk_p, spk_q);
            ratio_pq = sum(d_pq < spikeDistance) / numel(d_pq);
            [d_qp, ~] = nearest_index(spk_q, spk_p);
            ratio_qp = sum(d_qp < spikeDistance) / numel(d_qp);
            if (ratio_pq > overlap_thr && std(d_pq(d_pq < spikeDistance)) < 1.6) || (ratio_qp > overlap_thr && std(d_qp(d_qp < spikeDistance)) < 1.6)
                union_set(p, q);
            end
        end
    end
end
groups = containers.Map('KeyType','int32','ValueType','any');
for i = 1:nClusters
    r = find_parent(i);
    if isKey(groups, r)
        groups(r) = [groups(r), i];
    else
        groups(r) = i;
    end
end
group_keys = keys(groups);
for k = 1:length(group_keys)
    idxs = groups(group_keys{k});
    if numel(idxs) < 2
        continue;
    end
    best_idx = idxs(1);
    best_amp = sortedSamples{clusters(best_idx).channel,1}.clusteringInfo.clusterSelection.max_val(clusters(best_idx).local_idx);
    for j = idxs
        cur_amp = sortedSamples{clusters(j).channel,1}.clusteringInfo.clusterSelection.max_val(clusters(j).local_idx);
        if cur_amp > best_amp || (cur_amp == best_amp && clusters(j).channel < clusters(best_idx).channel)
            best_idx = j;
            best_amp = cur_amp;
        end
    end
    primary_channel = clusters(best_idx).channel;
    for j = idxs
        if j ~= best_idx
            ch = clusters(j).channel;
            if primary_channel > ch
                new_max = middleChannelIdx + (primary_channel - ch);
            elseif primary_channel < ch
                new_max = middleChannelIdx - (ch - primary_channel);
            else
                new_max = middleChannelIdx;
            end
            clust_sel = sortedSamples{ch,1}.clusteringInfo.clusterSelection;
            clust_sel.keep(clusters(j).local_idx) = 0;
            clust_sel.max_channelIdx(clusters(j).local_idx) = primary_channel;
            clust_sel.max_channel(clusters(j).local_idx) = new_max;
            sortedSamples{ch,1}.clusteringInfo.clusterSelection = clust_sel;
        end
    end
end
end