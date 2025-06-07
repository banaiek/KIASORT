function merge_matrix = kiaSort_merge_features(labels, pca, cfg)

if isfield(cfg, 'nPCAcomp')
    X = pca.score(:, 1);
end
unique_labels = unique(labels);
nClusters = numel(unique_labels);
merge_matrix = zeros(nClusters);
overlap_matrix = zeros(nClusters);

if ~isfield(cfg, 'T_merge'), cfg.T_merge = 0.05; end
if ~isfield(cfg, 'G'), cfg.G = 2; end
if ~isfield(cfg, 'gamma'), cfg.gamma = 3; end

options = statset('MaxIter',1000, 'Display','off');
regVal = 1e-5;

for i = 1:nClusters-1
    idx_i = (labels == unique_labels(i));
    data_i = X(idx_i, :);
    center_i = mean(data_i, 1);
    for j = i+1:nClusters
        idx_j = (labels == unique_labels(j));
        data_j = X(idx_j, :);
        center_j = mean(data_j, 1);
        diff_center = center_j - center_i;
        norm_diff = norm(diff_center);
        if norm_diff == 0
            merge_matrix(i,j) = 1; 
            merge_matrix(j,i) = 1;
            continue;
        end
        d = diff_center / norm_diff;
        proj_i = data_i * d';
        proj_j = data_j * d';
        
        % Skip if there are too few samples for reliable fitting.
        if numel(proj_i) < 2 || numel(proj_j) < 2
            overlap = 0;
            overlap_matrix(i,j) = overlap;
            overlap_matrix(j,i) = overlap;
            continue;
        end
        
        nComp_i = cfg.G; 
        if numel(proj_i) < cfg.G, nComp_i = 1; end
        nComp_j = cfg.G; 
        if numel(proj_j) < cfg.G, nComp_j = 1; end
        
        try
            gm_i = fitgmdist(proj_i, nComp_i, 'RegularizationValue', regVal, 'Options', options);
        catch
            overlap = 0;
            overlap_matrix(i,j) = overlap;
            overlap_matrix(j,i) = overlap;
            continue;
        end
        try
            gm_j = fitgmdist(proj_j, nComp_j, 'RegularizationValue', regVal, 'Options', options);
        catch
            overlap = 0;
            overlap_matrix(i,j) = overlap;
            overlap_matrix(j,i) = overlap;
            continue;
        end
        
        [mu_i_all, sigma_i_all] = computeGMMStats(gm_i);
        [mu_j_all, sigma_j_all] = computeGMMStats(gm_j);
        
        % Define integration range based on union of supports.
        x_start = min(mu_i_all - cfg.gamma * sigma_i_all, mu_j_all - cfg.gamma * sigma_j_all);
        x_end   = max(mu_i_all + cfg.gamma * sigma_i_all, mu_j_all + cfg.gamma * sigma_j_all);
        x_range = linspace(x_start, x_end, 200);
        pdf_i = pdf(gm_i, x_range');
        pdf_j = pdf(gm_j, x_range');
        overlap = trapz(x_range, pdf_i .* pdf_j);
        overlap_matrix(i,j) = overlap;
        overlap_matrix(j,i) = overlap;
        
        if overlap >= cfg.T_merge
            merge_matrix(i,j) = 1;
            merge_matrix(j,i) = 1;
        end
    end
end

if (any(unique_labels==-1))
    overlap_matrix (1,:) = nan;
    overlap_matrix (:,1) = nan;
end
merge_matrix = overlap_matrix > mean(overlap_matrix(:),"omitmissing") + 1.68*std(overlap_matrix(:),"omitmissing");
merge_matrix(1:nClusters+1:end) = 0;
merge_matrix(1:nClusters+1:end) = 0;
end

function [mu_all, sigma_all] = computeGMMStats(gm)
% Compute the weighted mean and overall standard deviation for a 1-D GMM.
w = gm.ComponentProportion(:);
mu_vec = gm.mu(:);
mu_all = sum(w .* mu_vec);
var_vec = squeeze(gm.Sigma);
if ~isscalar(var_vec)
    var_vec = var_vec(:);
end
variance_all = sum(w .* (var_vec + mu_vec.^2)) - mu_all.^2;
sigma_all = sqrt(variance_all);
end