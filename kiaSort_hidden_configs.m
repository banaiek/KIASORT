function cfg = kiaSort_hidden_configs(cfg)

extra_cfg = struct(...
    'numBands',                 4, ...               % Number of frequency bands for multiband extraction
    'borderMargin',             2, ...               % Minimum spike margin from start/end (ms)
    'edgeSampleExclusion',      2, ...               % Exclusion (ms) for edge samples
    'spikeSampleDistance',      1, ...               % Spike sample distance (ms)
    'method',                   "direct", ...        % Sorting method ("direct" or "indirect")
    'overlap_thr',              0.5 ...             % Maximum allowed template overlap
    );

% Append any missing fields from extra_cfg to cfg
extra_fields = fieldnames(extra_cfg);
for i = 1:length(extra_fields)
    field = extra_fields{i};
    if ~isfield(cfg, field)
        cfg.(field) = extra_cfg.(field);
    end
end

extra_varName = struct(...
    'numBands',                'Num. Frequency Bands', ...
    'borderMargin',            'Border margin (ms)', ...
    'edgeSampleExclusion',     'Edge Exclusion (ms)', ...
    'spikeSampleDistance',     'Sample spk Dist. (ms)', ...
    'method',                  'Sorting: Method', ...
    'overlap_thr',             'Cluster Overlap Thr' ...
);

if ~isfield(cfg, 'varName')
    cfg.varName = struct();
end

% Append varName fields to cfg.varName
extra_var_fields = fieldnames(extra_varName);
for i = 1:length(extra_var_fields)
    key = extra_var_fields{i};
    if ~isfield(cfg.varName, key)
        cfg.varName.(key) = extra_varName.(key);
    end
end

allFields = fieldnames(cfg);
fieldsWithoutVarName = setdiff(allFields, {'varName'}, 'stable');
cfg = orderfields(cfg, [fieldsWithoutVarName; {'varName'}]);

end