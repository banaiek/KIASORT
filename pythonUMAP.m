function umapOut = pythonUMAP(X,nComp)

try    
    numpy = py.importlib.import_module('numpy');
    umap_umap_ = py.importlib.import_module('umap.umap_');
    reducer = umap_umap_.UMAP(pyargs('n_components', int32(nComp))); 
    data_py = numpy.array(X);
    embedding = reducer.fit_transform(data_py);
    umapOut = double(embedding);
catch
    warning('UMAP not installed or not found. Skipping UMAP compare.');
    umapOut = [];
end

end