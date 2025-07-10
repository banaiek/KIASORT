function out = kiaSort_evaluate_labels(predLabels, clusterSelection, altMaxChannel, validKeep)

% Check which clusters to be kept, and find the reference channel for
% those not kept

classLabels = clusterSelection.classLabels;
keep = clusterSelection.keep;
keep(classLabels==-1) = 0;
numClasses = length(classLabels);
kept_idx = true(size(predLabels));
max_channel_idx = clusterSelection.max_channelIdx;


count = 0;
included_channels = [];
lookUpChannel_temp = [];
not_kept_idx_temp = [];
out = [];
for iClass = 1 : numClasses
    if keep(iClass) || ~any(predLabels == classLabels(iClass))
        continue;
    else

        not_kept_idx = find(predLabels == classLabels(iClass));
        kept_idx(not_kept_idx) = false;
        if  classLabels(iClass) ~=-1
            count = count + 1;
            not_kept_idx_temp{count, 1} = not_kept_idx;
            lookUpChannel_temp(count,1) = max_channel_idx(iClass);
        end
        included_channels = [included_channels ;max_channel_idx(iClass)];
    end
end

unique_channels = unique(lookUpChannel_temp);
count = 0;
for i = 1:length(unique_channels)
     count = count + 1;
     ch_id = find(lookUpChannel_temp == unique_channels(i));
     out.not_kept_idx{count, 1} = cell2mat(not_kept_idx_temp(ch_id,1));
     out.lookUpChannel(count,1) = unique_channels(i);
end

out.keep = kept_idx & validKeep;

for i = 1:count
    not_kept_idx_alt = find(predLabels > 0 & ~out.keep & altMaxChannel == out.lookUpChannel(i,1));
    if any(not_kept_idx_alt)
        all_not_kept = [out.not_kept_idx{i, 1}; not_kept_idx_alt];
        out.not_kept_idx{i, 1} = unique(all_not_kept);
    end
end

not_assigned = (predLabels > 0 & ~out.keep & ~ismember(altMaxChannel,included_channels));

if any(not_assigned)
    uniqe_Alt = unique(altMaxChannel(not_assigned));
    for i = 1:length(uniqe_Alt)
        not_kept_idx = find(predLabels > 0 & ~out.keep & altMaxChannel==uniqe_Alt(i));
        count = count + 1;
        out.not_kept_idx{count, 1} = not_kept_idx;
        out.lookUpChannel(count,1) = uniqe_Alt(i);
    end
end

out.nRelabling  = count;
end