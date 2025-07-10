function cmap = polar_colormap(hsvColor, N)

    % default color map for crimson to white to yale
    if nargin < 1        
        pos_color   = [220, 20, 60] / 255;
        mid_color   = [1, 1, 1];
        neg_color   = [15, 77, 146] / 255;
    else
        pos_color    =  hsv2rgb(hsvColor);
        mid_color    =  hsv2rgb([1, 1-hsvColor(2) 1-hsvColor(3)]);        
        neg_color    =  hsv2rgb([1-hsvColor(1), hsvColor(2) hsvColor(3)]);
    end

    if nargin < 2
        N = 256; 
    end


    % Number of colors for each gradient segment
    N_half = floor(N / 2);

    % Generate gradient from Yale Blue to White
    r1 = linspace(neg_color(1), mid_color(1), N_half)';
    g1 = linspace(neg_color(2), mid_color(2), N_half)';
    b1 = linspace(neg_color(3), mid_color(3), N_half)';

    % Generate gradient from White to Crimson
    r2 = linspace(mid_color(1), pos_color(1), N - N_half)';
    g2 = linspace(mid_color(2), pos_color(2), N - N_half)';
    b2 = linspace(mid_color(3), pos_color(3), N - N_half)';

    r = [r1; r2];
    g = [g1; g2];
    b = [b1; b2];

    cmap = [r, g, b];

end