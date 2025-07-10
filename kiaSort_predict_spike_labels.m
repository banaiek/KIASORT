function [predLabels, lowDComponents, validKeep]  = kiaSort_predict_spike_labels(data, cfg)

% Import Clustering info
net          = data.classifierInfo.trainedNet;
mdl          = data.classifierInfo.trainedMdl;
uniqueLabels = data.classifierInfo.classLabels;
polarity     = data.classifierInfo.class_polarity;
low_thr      = data.classifierInfo.lowAmpThr;
high_thr     = data.classifierInfo.highAmpThr;
ampVal       = data.amplitude;

if cfg.usePCA
    nComp = cfg.nPCAcomp;
    centered_new_data = data.waveformNorm - data.PCA.mu;
    lowDComponents    = centered_new_data * data.PCA.coeff;
    Xinput = lowDComponents(:,1:nComp);
else
    centered_new_data = data.waveformNorm - data.PCA.mu;
    lowDComponents    = centered_new_data * data.PCA.coeff;
    Xinput = data.waveformNorm;
end

validKeep = true(size(Xinput,1),1);

switch lower(cfg.method)
    case 'direct'
        switch lower(cfg.modelType)
            case 'mlp'
                predLabels = predict(net,Xinput);
                predLabels = double(onehotdecode(predLabels,double(uniqueLabels),2));
            case 'cnn'
                predLabels = predict(net, Xinput);
                predLabels = double(onehotdecode(predLabels,double(uniqueLabels),2));
            otherwise
                [predLabels, score] = predict(mdl, Xinput);
                reliabilityScore    = max(score, [], 2);
                validKeep(predLabels < 0) = 0;
                validKeep(exp(reliabilityScore) < .75 & max(double(exp(score)),[],2)./median(maxk(exp(score),3,2),2) < 1.2 ) = 0;
                validKeep(exp(reliabilityScore) < .7) = 0;
        end


    case 'indirect'
        lowDComponents = predict(net, Xinput);
        predLabels = predict(mdl,lowDComponents);
    otherwise
        error('Unsupported method ''%s''.', cfg.method);
end

unique_pred = unique(predLabels(validKeep));
for i = 1:length(unique_pred)
    class_idx = find(uniqueLabels == unique_pred(i));
    pred_idx = find(predLabels == unique_pred(i));

    med_amp = median(ampVal(pred_idx));
    dev_amp = 3 * max(std(med_amp),mad(med_amp));
    
    mad_outlier = ampVal(pred_idx) < med_amp + dev_amp & ampVal(pred_idx) > med_amp - dev_amp;
    
    if polarity(class_idx) == 1
        outliers = (ampVal(pred_idx) < low_thr(class_idx) | ampVal(pred_idx) > high_thr(class_idx));
    else
        outliers = (ampVal(pred_idx) > low_thr(class_idx) | ampVal(pred_idx) < high_thr(class_idx));
    end

    validKeep(pred_idx(outliers)) = 0;
    validKeep(pred_idx(mad_outlier)) = 0;
end

lowDComponents = lowDComponents(:,1:3);

end
