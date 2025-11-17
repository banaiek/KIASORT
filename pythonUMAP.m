function umapOut = pythonUMAP(X,nComp)

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

                reducer = umap_umap_.UMAP(pyargs(...
            'n_components', int32(nComp), ...
            'n_neighbors', int32(25), ...
            'min_dist', .15));
    data_py = numpy.array(X);
    embedding = reducer.fit_transform(data_py);
    umapOut = double(embedding);
catch
    error('UMAP not installed or not found. Check the intruction for and make sure it is installed.');
    umapOut = [];
end

end