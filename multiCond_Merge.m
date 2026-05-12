function merge_matrix = multiCond_Merge(A, options)
% Robust waveform comparison with temporal alignment.
% A: N×Nch×T matrix (N groups, Nch channels, T time points)
%
% OPTIMIZATIONS:
%   - Shape similarity is shift-invariant (circular shift preserves max/min
%     over time) — computed once per pair, not per shift.
%   - Correlation computed via FFT cross-correlation across all shifts at
%     once, replacing the inner shift loop entirely.

    if nargin < 2
        options.corr_weight  = 0.5;
        options.shape_weight = 0.5;
        options.threshold    = 0.9;
        options.max_shift    = 15;
    end

    N   = size(A, 1);
    Nch = size(A, 2);
    T   = size(A, 3);

    % ------------------------------------------------------------------
    % 1. Shape similarity: shift-invariant, precompute per cluster
    %    circshift preserves max/min over dim-3, so peak-to-peak is constant.
    % ------------------------------------------------------------------
    peaks = squeeze(max(A,[],3) - min(A,[],3));   % N × Nch
    if N == 1, peaks = peaks(:)'; end              % force row

    peakNorms = sqrt(sum(peaks.^2, 2));            % N × 1

    max_shape_sim = zeros(N, N);
    for i = 1:N
        for j = i:N
            denom = peakNorms(i) + peakNorms(j);
            if denom > 0
                max_shape_sim(i,j) = 1 - norm(peaks(i,:) - peaks(j,:)) / denom;
            else
                max_shape_sim(i,j) = 1;
            end
            max_shape_sim(j,i) = max_shape_sim(i,j);
        end
    end

    % ------------------------------------------------------------------
    % 2. Correlation via FFT — all lags at once
    %    For circular shift s:
    %      score(s) = Σ_ch Σ_t a_norm(ch,t) · b_norm(ch, (t-s) mod T)
    %    which is the sum-of-channels circular cross-correlation, computed
    %    efficiently with fft/ifft.
    % ------------------------------------------------------------------

    % Pre-normalise each cluster (global zero-mean, unit-std across all
    % channels and time points — same as original code).
    A_norm2D = zeros(N, Nch, T);
    A_fft    = zeros(N, Nch, T);   % will hold complex FFTs

    for i = 1:N
        flat  = reshape(A(i,:,:), 1, []);
        mu    = mean(flat);
        sigma = std(flat);
        if sigma < eps, sigma = 1; end
        normed = (A(i,:,:) - mu) / sigma;
        A_norm2D(i,:,:) = normed;
        for ch = 1:Nch
            A_fft(i, ch, :) = fft(squeeze(normed(1, ch, :)));
        end
    end

    % Allowed shift indices in the FFT output vector (length T).
    %   shift s >= 0  →  index s + 1
    %   shift s <  0  →  index T + s + 1
    shifts   = -options.max_shift : options.max_shift;
    shiftIdx = mod(shifts, T) + 1;

    max_corr_sim = zeros(N, N);
    for i = 1:N
        fi = squeeze(A_fft(i,:,:));          % Nch × T  (complex)
        for j = i:N
            fj = squeeze(A_fft(j,:,:));      % Nch × T  (complex)

            % Cross-correlation summed over channels
            xc = sum(real(ifft(fi .* conj(fj), [], 2)), 1);   % 1 × T
            xc = xc / (Nch * T);

            max_corr_sim(i,j) = max(xc(shiftIdx));
            max_corr_sim(j,i) = max_corr_sim(i,j);
        end
    end

    % ------------------------------------------------------------------
    % 3. Combine & threshold
    % ------------------------------------------------------------------
    combined_sim = options.corr_weight  * max_corr_sim ...
                 + options.shape_weight * max_shape_sim;

    merge_matrix = combined_sim > options.threshold;
    merge_matrix(logical(eye(N))) = 0;
end