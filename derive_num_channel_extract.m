function n = derive_num_channel_extract(channel_locations, radius, max_n)
%DERIVE_NUM_CHANNEL_EXTRACT  Pick a global half-window from a y-radius.
%
%   n = derive_num_channel_extract(channel_locations, radius, max_n)
%
%   For every channel m the function counts how many channels with index
%   < m and > m fall within `radius` micrometres in y. The required
%   half-window per channel is max(before, after); the global value
%   returned is the median of those per-channel half-windows, capped at
%   max_n and floored at 1.
%
%   When channel_locations is empty or has fewer than 2 columns, the
%   configured max_n is returned unchanged.

    if nargin < 3 || isempty(max_n) || max_n < 1
        max_n = 1;
    end
    max_n = floor(max_n);

    if isempty(channel_locations) || size(channel_locations, 2) < 2 ...
            || isempty(radius) || radius <= 0
        n = max(1, max_n);
        return;
    end

    ys = channel_locations(:, 2);
    nCh = numel(ys);
    if nCh < 2
        n = max(1, max_n);
        return;
    end

    halfWin = zeros(nCh, 1);
    for m = 1:nCh
        within = abs(ys - ys(m)) <= radius;
        within(m) = false;
        idx = find(within);
        if isempty(idx)
            halfWin(m) = 1;
        else
            before = max(0, m - min(idx));
            after  = max(0, max(idx) - m);
            halfWin(m) = max(1, max(before, after));
        end
    end

    n = round(median(halfWin));
    n = max(1, min(n, max_n));
end
