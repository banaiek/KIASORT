function [kept_x, kept_y, kept_label] = keep_highest_value(x, y, min_distance, label, min_distance_adj)

x = x(:);
y = y(:);

if nargin < 4 || isempty(label)
    label = ones(size(x));
end
if nargin < 5 || isempty(min_distance_adj)
    min_distance_adj = min_distance;
end

if isempty(x)
    kept_x = [];
    kept_y = [];
    kept_label = [];
    return
end

[x_sorted, sort_idx] = sort(x);
y_sorted     = y(sort_idx);
label_sorted = label(sort_idx);

N             = numel(x_sorted);
kept_indices  = [];
last_kept_idx = 1;
kept_indices(end+1) = last_kept_idx;

for i = 2:N
    distance = x_sorted(i) - x_sorted(last_kept_idx);
    if label_sorted(i) == label_sorted(last_kept_idx)
        thresh = min_distance_adj;
    else
        thresh = min_distance;
    end
    if distance > thresh
        kept_indices(end+1) = i;
        last_kept_idx       = i;
    else
        if y_sorted(i) > y_sorted(last_kept_idx)
            kept_indices(end) = i;
            last_kept_idx     = i;
        end
    end
end

kept_x     = x_sorted(kept_indices);
kept_y     = y_sorted(kept_indices);
kept_label = label_sorted(kept_indices);
end