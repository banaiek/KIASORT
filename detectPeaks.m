function [idx] = detectPeaks(signalIN,threshold,num_passing_points,polarity,tolerance)
% detect peaks that stay above threshold for n samples

if nargin < 3 || isempty(num_passing_points)
    num_passing_points = 0;
end

if nargin < 4 || isempty(polarity)
    polarity = 1;
end

if nargin < 5 || isempty(tolerance)
    tolerance = 0;
end

n = num_passing_points;
window_size = 2 * n + 1;

if polarity == 1

    if n
        is_above = signalIN > threshold;
        sum_in_window = movsum(double(is_above), [n n], 'Endpoints', 'fill');
        max_in_window = movmax(sum_in_window, [n n], 'Endpoints', 'fill');
        idx = find(max_in_window >= window_size-tolerance & signalIN-circshift(signalIN,1)>0 & signalIN-circshift(signalIN,-1)>0);
    else
        idx = find(signalIN > threshold & signalIN-circshift(signalIN,1)>0 & signalIN-circshift(signalIN,-1)>0);
    end


else
    if n
        is_above = signalIN < threshold;
        sum_in_window = movsum(double(is_above), [n n], 'Endpoints', 'fill');
        max_in_window = movmax(sum_in_window, [n n], 'Endpoints', 'fill');
        idx = find(max_in_window >= window_size-tolerance & signalIN-circshift(signalIN,1)<0 & signalIN-circshift(signalIN,-1)<0);
    else
        idx = find(signalIN < threshold & signalIN-circshift(signalIN,1)<0 & signalIN-circshift(signalIN,-1)<0);
    end

end

end