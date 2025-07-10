function out = evaluate_labels(predLabels, clusterSelection, validKeep)

classLabels = clusterSelection.classLabels;
keep = clusterSelection.keep;
keep(classLabels==-1) = 0;
numClasses = length(classLabels);
kept_idx = true(size(predLabels));
max_channel_idx = clusterSelection.max_channelIdx;  


count = 0;

for iClass = 1 : numClasses
    if keep(iClass) || ~any(predLabels == classLabels(iClass))
        continue;
    else
        
        not_kept_idx = find(predLabels == classLabels(iClass));
        kept_idx(not_kept_idx) = false;
        if  classLabels(iClass) ~=-1
            count = count + 1;
            out.not_kept_idx{count, 1} = not_kept_idx;
            out.lookUpChannel(count,1) = max_channel_idx(iClass);
        end
    end
end
out.keep    = kept_idx & validKeep;
out.nRelabling  = count;
end