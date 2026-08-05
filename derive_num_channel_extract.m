function n = derive_num_channel_extract(channel_locations, radius, max_n)
%DERIVE_NUM_CHANNEL_EXTRACT  Pick a global half-window from a y-radius.
%
%   n = derive_num_channel_extract(channel_locations, radius, max_n)
%
%   For every channel m the function looks at the channel indices that
%   fall within `radius` micrometres in y. The required half-window
%   per channel is max(before-span, after-span); the global value
%   returned is the LARGEST per-channel half-window across the probe,
%   capped at max_n and floored at 1.
%
%   max_n == 0 is passed through untouched: it selects main-channel-only
%   mode (a 1-channel window), so the radius is not consulted at all.
%
%   We use the max across channels (not the median, as an earlier
%   version did) because spike windows are sized once for the whole
%   probe but every spike must include all of its within-radius
%   neighbours -- a median was leaving the derived value one short of
%   the cap on uniform-pitch probes whenever per-channel half-windows
%   landed exactly on an integer just below the cap, so increasing
%   the radius felt like it had no effect. Using max guarantees that
%   bumping the radius drives the derived value upward up to max_n.
%
%   When channel_locations is empty or has fewer than 2 columns, the
%   configured max_n is returned unchanged.

    if nargin < 3 || isempty(max_n) || max_n < 0
        max_n = 1;
    end
    max_n = floor(max_n);
    if max_n == 0
        n = 0;
        return;
    end

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

    n = max(halfWin);
    n = max(1, min(n, max_n));
end
