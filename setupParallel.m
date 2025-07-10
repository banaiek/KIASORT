function setupParallel(cfg)
    if cfg.parallelProcessing && isempty(gcp('nocreate'))
        if ischar(cfg.numParallelWorker) && strcmpi(cfg.numParallelWorker, 'auto')
            parpool;
        else
            cluster = parcluster();
            workers = cfg.numParallelWorker;
            if workers > cluster.NumWorkers
                workers = cluster.NumWorkers;
            end
            parpool(workers);
        end
    end
end