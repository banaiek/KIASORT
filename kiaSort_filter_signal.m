function [out] = kiaSort_filter_signal(inputSignal, cfg)


passband = cfg.bandpass;
fs = cfg.samplingFrequency;


% Check for GPU usage.
useGPU = cfg.useGPU && (gpuDeviceCount > 0);
if useGPU
    inputSignal = gpuArray(inputSignal);
end

if ~useGPU
    setupParallel(cfg);
end

bandpass_signal = bandpass_filter(inputSignal, passband, fs)';

if useGPU
    bandpass_signal = gather(bandpass_signal);
end

out.bandpass_signal = bandpass_signal;
end

function bp_data = bandpass_filter(signal, passband, fs)

[b, a] = butter(4, passband/(fs/2), 'bandpass');
if isa(signal, 'gpuArray')
    b = gpuArray(b);
    a = gpuArray(a);
end
bp_data = filtfilt(b, a, signal);
end
