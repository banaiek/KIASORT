function y = smooth_square(L, trans, sharp)

    N = (L - 1) / 2;
    center = N + 1;
    x = 1:L;
    
    y = zeros(1, L);
    leftIdx = (x <= center);
    y(leftIdx) = 1 ./ (1 + exp(-sharp * (x(leftIdx) - (center - trans/2))));

    rightIdx = (x >= center);
    y(rightIdx) = 1 ./ (1 + exp( sharp * (x(rightIdx) - (center + trans/2))));
end