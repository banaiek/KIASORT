function out = kiaSort_preprocess_waveforms(data, cfg)

fs              = cfg.samplingFrequency;
midPoint        = floor(cfg.spikeDuration * fs/(2*1000))+1;
spike_length    = floor(cfg.clusteringSpikeDuration * fs/(2*1000));

informative_Chan = data.waveformInfo.informative_Chan; 

waveform_full = data.waveform;
waveform_full(isnan(waveform_full)) = 0;
 waveform = waveform_full(:,:,midPoint-spike_length:midPoint+spike_length);
% waveform = waveform_full(:,find(informative_Chan),midPoint-spike_length:midPoint+spike_length);
[N, C, T] = size(waveform);
% flattened_waveform = reshape( permute(waveform, [1 3 2]),  N, C*T);
flattened_waveform = reshape(waveform,  N, C*T);
% flattened_waveform = flattened_waveform.*hamming(size(flattened_waveform,2))';
out.waveformNorm  = flattened_waveform(:,sum(flattened_waveform,1)~=0);

end