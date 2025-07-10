function [sampleFeatures, sortedSamples] = unify_cluster_relabling(cfg, sampleFeatures, sortedSamples)

numChannels = length(sortedSamples);


% run relabling three times to make sure it all channels are relabled correctly
for k = 1:3

    for i=1:numChannels

        if ~isempty(sampleFeatures{i, 1})

            spk_idx_i = sampleFeatures{i, 1}.spk_idx;
            spk_lbl_i = sampleFeatures{i, 1}.updatedLabels;
            clust_select_i = sortedSamples{i, 1}.clusteringInfo.clusterSelection;
            unique_lbl_i = sortedSamples{i, 1}.clusteringInfo.clusterSelection.classLabels;
            max_channelIdx_i = sortedSamples{i, 1}.clusteringInfo.clusterSelection.max_channelIdx;
            notKeep_i = find(~clust_select_i.keep);

            if any(notKeep_i)
                for iKeep = 1:length(notKeep_i)

                    lbl_i = clust_select_i.classLabels(notKeep_i(iKeep));
                    if lbl_i == -1
                        continue
                    end
                    clust_spk_i = spk_idx_i(spk_lbl_i == lbl_i);
                    relabeld_ch_i = max_channelIdx_i(notKeep_i(iKeep));


                    clust_spk_j = sampleFeatures{relabeld_ch_i, 1}.spk_idx;
                    spk_lbl_j = sampleFeatures{relabeld_ch_i, 1}.updatedLabels;
                    unique_lbl_j = sortedSamples{relabeld_ch_i, 1}.clusteringInfo.clusterSelection.classLabels;
                    max_channelIdx_j = sortedSamples{relabeld_ch_i, 1}.clusteringInfo.clusterSelection.max_channelIdx;

                    [d, idxij] = nearest_index(clust_spk_i, clust_spk_j);
                    valid_idx = d < 1.5*cfg.spikeDistance*cfg.samplingFrequency/1000;

                    if sum(valid_idx)/length(valid_idx) > .25

                        transfered_lbl = spk_lbl_j(idxij(valid_idx));
                        clusterSampleCounts = histc(transfered_lbl,unique_lbl_j);
                        [bestMatch_val,bestMatch_lbl] = max(clusterSampleCounts./length(clust_spk_i));

                        if max_channelIdx_j(bestMatch_lbl) == i
                            clust_spk_j_i = clust_spk_j(spk_lbl_j == unique_lbl_j(bestMatch_lbl));
                            [d_i, idxij_i] = nearest_index(clust_spk_j_i, spk_idx_i);
                            valid_idx_j_i = d_i < 1.5 * cfg.spikeDistance*cfg.samplingFrequency/1000;
                            if sum(valid_idx_j_i)/length(valid_idx_j_i) > .25
                                transfered_lbl_j_i = spk_lbl_i(idxij_i(valid_idx_j_i));
                                clusterSampleCounts_j_i = histc(transfered_lbl_j_i,unique_lbl_i);
                                [bestMatch_val_ji,bestMatch_lbl_ji] = max(clusterSampleCounts_j_i./length(clust_spk_j_i));
                                if bestMatch_lbl_ji == notKeep_i(iKeep)
                                    clust_select_i.keep(notKeep_i(iKeep)) = 1;
                                    relabeled_channel_max = clust_select_i.max_channel(notKeep_i(iKeep)) + max_channelIdx_j(bestMatch_lbl)-relabeld_ch_i;
                                    clust_select_i.max_channel(notKeep_i(iKeep)) = relabeled_channel_max;
                                    clust_select_i.max_channelIdx(notKeep_i(iKeep)) = max_channelIdx_j(bestMatch_lbl);
                                    sortedSamples{i, 1}.clusteringInfo.clusterSelection = clust_select_i;                                   
                                end                                
                            end
                            continue
                        else 
                            clust_select_i.max_channelIdx(notKeep_i(iKeep)) = max_channelIdx_j(bestMatch_lbl);
                            relabeled_channel_max = clust_select_i.max_channel(notKeep_i(iKeep)) + max_channelIdx_j(bestMatch_lbl)-relabeld_ch_i;
                            clust_select_i.max_channel(notKeep_i(iKeep)) = relabeled_channel_max;
                            sortedSamples{i, 1}.clusteringInfo.clusterSelection = clust_select_i;
                        end
                    end

                end
            end
        end
    end

end
end
