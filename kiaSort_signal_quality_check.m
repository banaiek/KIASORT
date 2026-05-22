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
        raw_for_psd       = gather(inputSignal);
    else
        raw_for_psd       = inputSignal;
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

    nyq = fs / 2.05;
    [opt_bp, band_scores] = optimal_spike_band( ...
        bandpass_signal, raw_for_psd, scale_factor, fs, 50, nyq);

    out.scale_factor     = scale_factor;
    out.rms_bands        = rms_bands;
    out.band_pairs       = band_pairs;
    out.Ns_seq           = Ns_seq;
    out.thresh_MAD       = thresh_MAD(:);
    out.optimal_bandpass = opt_bp;
    out.band_scores      = band_scores;
end

function [opt_bp, band_scores] = optimal_spike_band( ...
        broadband, raw_signal, scale_factor, fs, fmin, fmax)

    % broadband:  bandpass-filtered signal used ONLY for spike timing
    %             (MAD threshold + refractory). Reliable detection.
    % raw_signal: unfiltered input. The spike-window and noise-window
    %             FFTs are extracted from this so the PSD reflects the
    %             actual recording spectrum (including 1/f below the
    %             pre-filter cutoff). Without this, every recording's
    %             PSD looks the same shape because they all share the
    %             same pre-filter rolloff.

    [~, nSamp] = size(broadband);

    good_ch = scale_factor < 3;
    if sum(good_ch) < 2, good_ch = true(size(scale_factor)); end
    ch_idx = find(good_ch);
    nGood  = length(ch_idx);

    nfft = 2^nextpow2(round(0.016 * fs));
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
        bp_trace  = broadband(ch, :);
        raw_trace = raw_signal(ch, :);

        noise_est = 1.4826 * median(abs(bp_trace));
        thr = -4 * noise_est;

        cr = find(bp_trace(1:end-1) > thr & bp_trace(2:end) <= thr);
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
        spike_segs = raw_trace(spike_idx);

        is_zone = false(1, nSamp);
        zi = spike_idx(:);
        zi = zi(zi >= 1 & zi <= nSamp);
        is_zone(zi) = true;

        cand = (half_nfft+1):half_nfft:(nSamp - half_nfft);
        noise_centers = cand(~is_zone(cand));
        if isempty(noise_centers), continue; end

        max_noise = max(5 * length(spike_times), 2000);
        if length(noise_centers) > max_noise
            sel = round(linspace(1, length(noise_centers), max_noise));
            noise_centers = noise_centers(sel);
        end

        spike_segs = double(spike_segs);
        noise_segs = double(raw_trace(noise_centers(:)' + rng_idx));
        valid_ch(ci) = true;

        % Spike-triggered TIME-DOMAIN average ("template") then FFT.
        % LFP, thermal noise, and 60 Hz lines that are uncorrelated with
        % spike timing average out in the time domain, so the resulting
        % PSD is the pure spike spectrum |S(f)|^2 -- not |S|^2 + N.
        % Comparing this to the noise PSD makes sp - np truly reflect
        % spike-specific power; in LFP-heavy recordings sp - np goes
        % negative at low frequencies, where it should be.
        spike_template = mean(spike_segs, 2);
        Px_s = abs(fft(spike_template .* win_taper, nfft, 1)).^2;
        Px_n = mean(abs(fft(noise_segs   .* win_taper, nfft, 1)).^2, 2);
        Px_s = Px_s(1:half_nfft+1);
        Px_n = Px_n(1:half_nfft+1);

        spike_psd(:, ci) = Px_s(band_mask);
        noise_psd(:, ci) = Px_n(band_mask);
    end

    if ~any(valid_ch)
        opt_bp = [fmin, fmax];
        band_scores = struct('freq', [], 'score', []);
        return;
    end

    sp = spike_psd(:, valid_ch);
    np = noise_psd(:, valid_ch);

    bin_hz   = fs / nfft;
    smooth_n = max(3, round(150 / bin_hz));
    sp_med = movmedian(median(sp, 2), smooth_n);
    np_med = movmedian(median(np, 2), smooth_n);
    log_f  = log(max(f_band, eps));

    % Smooth the SIGNED difference first, THEN clip at zero. Clipping
    % before smoothing introduces a positive bias at frequencies where
    % the true spike power is ~0 but sp_med - np_med fluctuates -- the
    % positive fluctuations survive while the negative ones are killed,
    % so smoothed excess is artificially elevated at low f.
    diff_med = movmedian(sp_med - np_med, smooth_n);
    excess   = max(diff_med, 0);

    band_scores.freq       = f_band;
    band_scores.sp_med     = sp_med;
    band_scores.np_med     = np_med;
    band_scores.excess     = excess;
    band_scores.smooth_n   = smooth_n;
    band_scores.bin_hz     = bin_hz;

    if ~any(isfinite(excess)) || all(excess <= 0)
        opt_bp = [fmin, fmax];
        return;
    end

    % Use the peak of excess (which is recording-dependent) to set a
    % data-driven threshold. The lower edge is the lowest frequency in
    % [100, 800] where excess >= edge_frac * peak_excess; the upper
    % edge is the highest frequency in [1000, fmax] satisfying the same
    % test. Because the threshold scales with the peak, recordings with
    % a tall narrow prominence peak get a tight band; recordings with a
    % low broad peak get a wide band. The threshold is anchored to the
    % spike-spectrum excess, not to a fixed sp/np cutoff, so it varies
    % across datasets in a meaningful way.
    peak_search = find(f_band >= 500 & f_band <= 3000);
    if isempty(peak_search), peak_search = (1:nFreq)'; end
    [peak_excess, peak_rel] = max(excess(peak_search));
    peak_idx = peak_search(peak_rel);

    if peak_excess <= 0
        opt_bp = [fmin, fmax];
        return;
    end

    edge_frac = 0.5;
    edge_thr  = edge_frac * peak_excess;

    band_scores.peak_freq    = f_band(peak_idx);
    band_scores.peak_excess  = peak_excess;
    band_scores.edge_thr     = edge_thr;

    lo_search = find(f_band >= 100 & f_band <= 800);
    if isempty(lo_search), lo_search = (1:nFreq)'; end
    lo_freq = first_cross(excess, log_f, lo_search, edge_thr);

    hi_search = find(f_band >= 1000 & f_band <= fmax);
    if isempty(hi_search), hi_search = (1:nFreq)'; end
    hi_freq = last_cross(excess, log_f, hi_search, edge_thr);

    if lo_freq >= hi_freq
        lo_freq = log_f(1);
        hi_freq = log_f(end);
    end

    opt_bp = round([exp(lo_freq), exp(hi_freq)]);
end


function x_cross = first_cross(y, x, search_idx, thr)
    y = y(:); x = x(:);
    above = y(search_idx) >= thr;
    if ~any(above)
        x_cross = x(search_idx(end));
        return;
    end
    if above(1)
        x_cross = x(search_idx(1));
        return;
    end
    rel = find(above, 1, 'first');
    i = search_idx(rel);
    y1 = y(i-1); y2 = y(i);
    x1 = x(i-1); x2 = x(i);
    if y2 == y1
        x_cross = x1;
    else
        frac = (thr - y1) / (y2 - y1);
        x_cross = x1 + frac * (x2 - x1);
    end
end


function x_cross = last_cross(y, x, search_idx, thr)
    y = y(:); x = x(:);
    above = y(search_idx) >= thr;
    if ~any(above)
        x_cross = x(search_idx(1));
        return;
    end
    if above(end)
        x_cross = x(search_idx(end));
        return;
    end
    rel = find(above, 1, 'last');
    i = search_idx(rel);
    if i >= length(y)
        x_cross = x(end);
        return;
    end
    y1 = y(i); y2 = y(i+1);
    x1 = x(i); x2 = x(i+1);
    if y2 == y1
        x_cross = x1;
    else
        frac = (thr - y1) / (y2 - y1);
        x_cross = x1 + frac * (x2 - x1);
    end
end


function x_cross = walk_out(y, x, peak_idx, thr, side)
    % Walk outward from peak_idx until y drops below thr; return the
    % interpolated x at the crossing. Returns x(1) or x(end) if no
    % crossing exists.
    y = y(:); x = x(:);
    n = length(y);
    if strcmp(side, 'left')
        for i = peak_idx-1:-1:1
            if y(i) < thr
                y1 = y(i); y2 = y(i+1);
                x1 = x(i); x2 = x(i+1);
                if y2 == y1, x_cross = x1; return; end
                frac = (thr - y1) / (y2 - y1);
                x_cross = x1 + frac * (x2 - x1);
                return;
            end
        end
        x_cross = x(1);
    else
        for i = peak_idx+1:n
            if y(i) < thr
                y1 = y(i-1); y2 = y(i);
                x1 = x(i-1); x2 = x(i);
                if y2 == y1, x_cross = x1; return; end
                frac = (thr - y1) / (y2 - y1);
                x_cross = x1 + frac * (x2 - x1);
                return;
            end
        end
        x_cross = x(end);
    end
end


function x_cross = continuous_cross(y, x, thr, mode)
    % First (or last) value of x where y crosses thr, linearly
    % interpolated in x between the two adjacent samples around the
    % crossing. Returns NaN if no crossing exists.
    y = y(:); x = x(:);
    above = y >= thr;
    if ~any(above) || all(above)
        x_cross = NaN;
        return;
    end
    if strcmp(mode, 'first')
        i = find(above, 1, 'first');
        if i == 1, x_cross = x(1); return; end
        y1 = y(i-1); y2 = y(i);
        x1 = x(i-1); x2 = x(i);
    else
        i = find(above, 1, 'last');
        if i == numel(y), x_cross = x(end); return; end
        y1 = y(i); y2 = y(i+1);
        x1 = x(i); x2 = x(i+1);
    end
    if y2 == y1
        x_cross = x1;
    else
        frac = (thr - y1) / (y2 - y1);
        x_cross = x1 + frac * (x2 - x1);
    end
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