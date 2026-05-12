function [out] = kiaSort_signal_quality_check(inputSignal, cfg)
    fs = cfg.samplingFrequency;
    useGPU = cfg.useGPU && (gpuDeviceCount > 0);

    if useGPU
        inputSignal = gpuArray(inputSignal);
    else
        setupParallel(cfg);
    end

    cfg.bandpass = [500 6000];
    freq_pairs = equal_power_bands(cfg);
    freq_peaks = [unique(freq_pairs(:)); round(fs / 2.1)];
    band_pairs = [freq_pairs; [freq_peaks(end-1), freq_peaks(end)]];
    Ns_seq = round(0.5 * fs ./ movmean(freq_peaks, 2));
    Steps = length(freq_peaks) - 1;

    decomposed_bands = zeros([Steps, size(inputSignal)], 'like', inputSignal);

    if useGPU
        for i = 1:Steps
            decomposed_bands(i,:,:) = bandpass_filter(inputSignal', [freq_peaks(i) freq_peaks(i+1)], fs)';
        end
    elseif cfg.parallelProcessing
        parfor i = 1:Steps
            decomposed_bands(i,:,:) = bandpass_filter(inputSignal', [freq_peaks(i) freq_peaks(i+1)], fs)';
        end
    else
        for i = 1:Steps
            decomposed_bands(i,:,:) = bandpass_filter(inputSignal', [freq_peaks(i) freq_peaks(i+1)], fs)';
        end
    end

    rms_bands = std(decomposed_bands, [], 3, 'omitnan');

    bandpass_signal = squeeze(sum(decomposed_bands, 1));
    thresh_MAD = median(abs(bandpass_signal - median(bandpass_signal, 2)), 2);

    if useGPU
        rms_bands  = gather(rms_bands);
        thresh_MAD = gather(thresh_MAD);
        decomposed_bands = gather(decomposed_bands);
        bandpass_signal   = gather(bandpass_signal);
    end

    power_factors = power_adjustment_factors(band_pairs);
    rms_bands = power_factors .* rms_bands;

    rel_rms = max(rms_bands([1 end], :)) ./ max(rms_bands([2 end-2], :));
    scale_factor = rel_rms .* exp(0.5 * (rel_rms - 1));
    scale_factor = scale_factor / median(scale_factor);

    good_mask = scale_factor < 3;
    mad_mean = mean(thresh_MAD(good_mask));
    mad_std  = std(thresh_MAD(good_mask));
    scale_factor(thresh_MAD > mad_mean + 5 * mad_std) = 5;
    scale_factor(thresh_MAD < mad_mean - 5 * mad_std) = 5;

    [opt_bp, band_scores] = optimal_spike_band( ...
        bandpass_signal, scale_factor, fs, band_pairs(1,1), band_pairs(end,2));

    out.scale_factor     = scale_factor;
    out.rms_bands        = rms_bands;
    out.band_pairs       = band_pairs;
    out.Ns_seq           = Ns_seq;
    out.thresh_MAD       = thresh_MAD(:);
    out.optimal_bandpass = opt_bp;
    out.band_scores      = band_scores;
end

function [opt_bp, band_scores] = optimal_spike_band( ...
        broadband, scale_factor, fs, fmin, fmax)

    [~, nSamp] = size(broadband);

    good_ch = scale_factor < 3;
    if sum(good_ch) < 2, good_ch = true(size(scale_factor)); end
    ch_idx = find(good_ch);
    nGood  = length(ch_idx);

    nfft = 2^nextpow2(round(0.004 * fs));
    if nfft >= nSamp, nfft = 2^floor(log2(nSamp/2)); end
    half_nfft = nfft / 2;
    win_taper = hann(nfft);

    freqs = (0:half_nfft)' * (fs / nfft);
    band_mask = freqs >= fmin & freqs <= fmax;
    f_band = freqs(band_mask);
    nFreq  = length(f_band);

    refrac  = round(0.001 * fs);
    rng_idx = (-half_nfft+1:half_nfft)';

    spike_psd = zeros(nFreq, nGood);
    noise_psd = zeros(nFreq, nGood);
    valid_ch  = false(nGood, 1);

    for ci = 1:nGood
        ch = ch_idx(ci);
        trace = broadband(ch, :);

        noise_est = 1.4826 * median(abs(trace));
        thr = -4 * noise_est;

        cr = find(trace(1:end-1) > thr & trace(2:end) <= thr);
        if isempty(cr), continue; end

        last_kept = -inf;
        keep = false(size(cr));
        for k = 1:length(cr)
            if cr(k) - last_kept >= refrac
                keep(k) = true;
                last_kept = cr(k);
            end
        end
        cr = cr(keep);

        valid = cr > half_nfft & cr <= nSamp - half_nfft;
        spike_times = cr(valid);
        if length(spike_times) < 5, continue; end

        spike_idx = spike_times(:)' + rng_idx;
        spike_segs = trace(spike_idx);

        is_zone = false(1, nSamp);
        zi = spike_idx(:);
        zi = zi(zi >= 1 & zi <= nSamp);
        is_zone(zi) = true;

        cand = (half_nfft+1):half_nfft:(nSamp - half_nfft);
        noise_centers = cand(~is_zone(cand));
        if isempty(noise_centers), continue; end

        max_noise = max(2 * length(spike_times), 500);
        if length(noise_centers) > max_noise
            sel = round(linspace(1, length(noise_centers), max_noise));
            noise_centers = noise_centers(sel);
        end

        noise_segs = trace(noise_centers(:)' + rng_idx);
        valid_ch(ci) = true;

        Px_s = abs(fft(spike_segs .* win_taper, nfft, 1)).^2;
        Px_n = abs(fft(noise_segs .* win_taper, nfft, 1)).^2;
        Px_s = Px_s(1:half_nfft+1, :);
        Px_n = Px_n(1:half_nfft+1, :);

        spike_psd(:, ci) = mean(Px_s(band_mask, :), 2);
        noise_psd(:, ci) = mean(Px_n(band_mask, :), 2);
    end

    if ~any(valid_ch)
        opt_bp = [fmin, fmax];
        band_scores = struct('freq', [], 'score', []);
        return;
    end

    sp = spike_psd(:, valid_ch);
    np = noise_psd(:, valid_ch);

    snr    = median(sp ./ max(np, eps), 2);
    excess = max(sp - np, 0);
    efrac  = median(excess ./ max(sum(excess, 1), eps), 2);

    score = snr .* efrac;
    smooth_n = max(3, round(nFreq / 64));
    score = movmean(score, smooth_n);

    band_scores.freq  = f_band;
    band_scores.score = score;
    band_scores.snr   = snr;
    band_scores.efrac = efrac;

    if all(score <= 0) || ~any(isfinite(score))
        opt_bp = [fmin, fmax];
        return;
    end

    [peak_score, peak_idx] = max(score);
    lo_idx = find(score(1:peak_idx) >= 0.25 * peak_score, 1, 'first');
    hi_idx = peak_idx - 1 + find(score(peak_idx:end) >= 0.15 * peak_score, 1, 'last');
    if isempty(lo_idx), lo_idx = 1; end
    if isempty(hi_idx), hi_idx = nFreq; end

    opt_bp = round([f_band(lo_idx), f_band(hi_idx)]);
end

%% Butterworth bandpass with GPU support
function bp_data = bandpass_filter(signal, passband, fs)
    [b, a] = butter(2, passband / (fs/2), 'bandpass');
    if isa(signal, 'gpuArray')
        b = gpuArray(b);
        a = gpuArray(a);
    end
    bp_data = filtfilt(b, a, signal);
end