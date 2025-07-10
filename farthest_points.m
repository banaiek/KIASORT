function idx = farthest_points(pts)
    centroid = mean(pts, 1);
    dists = vecnorm(pts - centroid, 2, 2);
    n = ceil(0.25 * size(pts, 1));
    [~, idx] = maxk(dists, n);
end