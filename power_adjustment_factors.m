function factors = power_adjustment_factors(bands)
% compute_adjustment_factors Computes adjustment factors to equalize power across frequency bands.
%
% This function calculates adjustment factors for each frequency band to equalize
% the total power in each band, assuming a power spectral density that follows a 1/f power law.
%
% Usage:
%   factors = compute_adjustment_factors(bands)
%
% Input:
%   bands - An Nx2 matrix where each row represents a frequency band [f1, f2]
%
% Output:
%   factors - An Nx1 vector containing the adjustment factors for each band
%
% Example:
%   bands = [1, 2; 2, 4; 4, 8];
%   factors = compute_adjustment_factors(bands);

    % Number of frequency bands
    numBands = size(bands, 1);

    % Initialize array to store original power in each band
    originalPower = zeros(numBands, 1);

    % Compute the total power in each band by integrating 1/f over the band
    for i = 1:numBands
        f1 = bands(i, 1);
        f2 = bands(i, 2);
        
        % Check for valid frequency values
        if f1 <= 0 || f2 <= 0 || f2 <= f1
            error('Invalid frequency band: [%f, %f]', f1, f2);
        end
        
        % Calculate the total power in the band
        originalPower(i) = log(f2) - log(f1);
    end

    % Desired total power per band (mean of the original powers)
    desiredPower = mean(originalPower);

    % Compute adjustment factors for each band
    factors = sqrt(desiredPower ./ originalPower);
end