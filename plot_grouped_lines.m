function ax_handle = plot_grouped_lines(N, Ng, ax_handle,varargin)



if length(varargin) < 3
    value = .85;
else
    value = varargin{3};
end

if length(varargin) < 2
    saturation = .85;
else
    saturation = varargin{2};
end

if length(varargin) < 1
    hue = .85;
else
    hue = varargin{1};
end

if Ng == 0
 colors_rgb  = hsv2rgb([hue saturation value]);   
 color_order = repmat(colors_rgb,N,1);
else
% group line plot colors
num_groups = ceil(N / Ng);

% distinct hues for colors
hues = linspace(0, 1, num_groups + 1);
hues(end) = [];

% HSV to RGB
colors_hsv = [hues', repmat(saturation, num_groups, 1), repmat(value, num_groups, 1)];
colors_rgb = hsv2rgb(colors_hsv);
color_order = repelem(colors_rgb, Ng, 1);
color_order = color_order(1:N, :);
end
ax_handle.ColorOrder = color_order;
ax_handle.NextPlot = 'replacechildren';

end