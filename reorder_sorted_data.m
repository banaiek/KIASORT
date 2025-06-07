function varargout = reorder_sorted_data(idx, varargin)
% Reorder and sort data based on idx for vertical inputs.

if nargin < 1
    error('At least one input argument (idx) is required.');
end

[sorted_idx, sortOrder] = sort(idx);
[sorted_idx, uniqueIdx] = unique(sorted_idx, 'first');
sortedRows = sortOrder(uniqueIdx);

varargout{1} = sorted_idx;

for i = 1:length(varargin)
    if ndims(varargin{i}) == 3
        varargout{i + 1} = varargin{i}(sortedRows, :, :);
    else
        varargout{i + 1} = varargin{i}(sortedRows, :);
    end
end