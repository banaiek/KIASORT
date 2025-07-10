function xcorr_vals = binary_xcorr(X_idx, Y_idx, total_samples, Fs, Fs_d, max_lag, step, smoothN)
    
    X = zeros(total_samples,1);
    Y = zeros(total_samples,1);
    X(X_idx) = 1;
    Y(Y_idx) = 1;

    % Downsample factor
    downsample_factor = Fs / Fs_d;
    assert(mod(downsample_factor, 1) == 0, 'Fs must be a multiple of Fs_d');

    % binary downsampling
    X_down = binary_downsample(X, downsample_factor);
    Y_down = binary_downsample(Y, downsample_factor);

    structElem = ones(1, 2*smoothN + 1);

    X_down = conv(X_down, structElem, 'same');
    Y_down = conv(Y_down, structElem, 'same');
    
    L = length(X_down);
    assert(length(Y_down) == L, 'Signals must have equal length after downsampling');
    assert(2*max_lag < L, 'max_lag exceeds signal length capabilities');
    start_idx = max_lag + 1;
    end_idx = L - max_lag;
    X_central = logical(X_down(start_idx:end_idx));
    
    lags = -max_lag:step:max_lag;
    num_lags = length(lags);
    xcorr_vals = zeros(1, num_lags);
    
    Y_logical = logical(Y_down);
    
    % Compute correlations for each lag
    for k = 1:num_lags
        i = lags(k);
        y_start = start_idx + i;
        y_end = end_idx + i;
        Y_segment = Y_logical(y_start:y_end);
        xcorr_vals(k) = sum(X_central & Y_segment);
    end
end

function downsampled = binary_downsample(signal, factor)
    
    orig_len = length(signal);
    new_len = ceil(orig_len / factor) * factor;
    
    % Pad with zeros 
    padded = false(new_len, 1);
    padded(1:orig_len) = logical(signal);
    
    reshaped = reshape(padded, factor, []);
    downsampled = any(reshaped, 1)';
end