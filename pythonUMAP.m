function umapOut = pythonUMAP(X,nComp,seed)
% seed: [] (default) leaves UMAP unseeded; an integer makes the embedding
% reproducible, which an A/B of anything upstream of UMAP needs.
if nargin < 3, seed = []; end

try    
    numpy = py.importlib.import_module('numpy');
    umap_umap_ = py.importlib.import_module('umap.umap_');
    % reducer = umap_umap_.UMAP(pyargs('n_components', int32(nComp),'min_dist', 0.01));
       % )); 
       kw = py.dict(pyargs('p', 2));

        % reducer = umap_umap_.UMAP(pyargs(...
        %     'n_components', int32(nComp), ...
        %     'n_neighbors', int32(25), ...
        %     'min_dist', .25, ...
        %     'metric', 'minkowski','metric_kwds',kw));

        if isempty(seed)
            reducer = umap_umap_.UMAP(pyargs(...
                'n_components', int32(nComp), ...
                'n_neighbors', int32(25), ...
                'min_dist', .05));
        else
            reducer = umap_umap_.UMAP(pyargs(...
                'n_components', int32(nComp), ...
                'n_neighbors', int32(25), ...
                'min_dist', .05, ...
                'random_state', int64(seed)));
        end
    data_py = numpy.array(X);
    embedding = reducer.fit_transform(data_py);
    umapOut = double(embedding);
catch ME
    % Report the underlying Python error: a masked message cost a debugging
    % round when the failure was an argument type, not a missing package.
    error('pythonUMAP:failed', 'UMAP call failed: %s', ME.message);
end

end