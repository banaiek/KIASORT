function [out] = kiaSort_signal_quality_check(inputSignal, cfg)

fs = cfg.samplingFrequency;

% Check for GPU usage.
useGPU = cfg.useGPU && (gpuDeviceCount > 0);
if useGPU
    inputSignal = gpuArray(inputSignal);
end

if ~useGPU
    setupParallel(cfg);
end

cfg.bandpass = [500 6000];
freq_pairs = equal_power_bands(cfg);
freq_peaks = [unique(freq_pairs(:)); round(fs/2.1)];
band_pairs = [freq_pairs; [freq_peaks(end-1), freq_peaks(end)]];
Ns_seq = round(0.5 * fs ./ movmean(freq_peaks, 2));
freq_Seq = freq_peaks(2:end);
Steps = length(freq_Seq);

decomposed_bands = zeros([Steps, size(inputSignal)], 'like', inputSignal);

if useGPU
    for i = 1:Steps
        decomposed_bands(i, :, :) = bandpass_filter(inputSignal', [freq_peaks(i) freq_peaks(i+1)], fs)';
    end
elseif cfg.parallelProcessing
    parfor i = 1:Steps
        decomposed_bands(i, :, :) = bandpass_filter(inputSignal', [freq_peaks(i) freq_peaks(i+1)], fs)';
    end
else
    for i = 1:Steps
        decomposed_bands(i, :, :) = bandpass_filter(inputSignal', [freq_peaks(i) freq_peaks(i+1)], fs)';
    end
end

rms_bands = std(decomposed_bands, [], 3, 'omitnan');

bandpass_signal = bandpass_filter(inputSignal', cfg.bandpass, fs);
thresh_MAD = median(abs(bandpass_signal-median(bandpass_signal)));

if useGPU
    if ~isempty(decomposed_bands)
        rms_bands = gather(rms_bands);
        thresh_MAD = gather(thresh_MAD);  
    end
end

power_factors = power_adjustment_factors(band_pairs);
rms_bands = power_factors .* rms_bands;

rel_rms = max(rms_bands([1, end],:)) ./ max(rms_bands([2, end-2],:));
scale_factor = rel_rms .* exp(0.5*(rel_rms - 1));

scale_factor = scale_factor/median(scale_factor);
scale_factor(thresh_MAD>mean(thresh_MAD)+ 5* std(thresh_MAD(scale_factor<3))) = 5;
scale_factor(thresh_MAD<mean(thresh_MAD)- 5* std(thresh_MAD(scale_factor<3))) = 5;
out.scale_factor = scale_factor;
out.rms_bands = rms_bands;
out.band_pairs = band_pairs;
out.Ns_seq = Ns_seq;
out.thresh_MAD = thresh_MAD(:);

end

function bp_data = bandpass_filter(signal, passband, fs)

[b, a] = butter(2, passband/(fs/2), 'bandpass');
if isa(signal, 'gpuArray')
    b = gpuArray(b);
    a = gpuArray(a);
end
bp_data = filtfilt(b, a, signal);
end
