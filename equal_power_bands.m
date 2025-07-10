function freq_pairs = equal_power_bands(cfg)
% equal_power_bands divides a frequency range into N bands of equal power for a 1/f PSD signal.
%

% Inputs: par.
%   - low_freq: Lower frequency bound in Hz 
%   - high_freq: Upper frequency bound in Hz 
%   - N: Number of desired frequency bands 
%
% Output:
%   - freq_pairs: An N-by-2 matrix where each row contains the lower and upper frequency bounds of a band
%
    % Import parameters
    low_freq = cfg.bandpass(1);
    high_freq = cfg.bandpass(2);
    N = cfg.numBands;

    % Ensure inputs are valid
    if low_freq <= 0 || high_freq <= 0
        error('Frequencies must be positive numbers.');
    end
    if low_freq >= high_freq
        error('low_freq must be less than high_freq.');
    end
    if N < 1 || floor(N) ~= N
        error('N must be a positive integer.');
    end

    % Calculate Total Power 
    P_total = log(high_freq / low_freq);

    % Calculate Power per Band 
    P_band = P_total / N;

    % Calculate Frequency Boundaries
    freq_bounds = zeros(N+1, 1);
    for i = 0:N
        freq_bounds(i+1) = low_freq * exp(i * P_band);
    end

    % Generate Frequency Pairs 
    freq_pairs = [freq_bounds(1:end-1), freq_bounds(2:end)];
end