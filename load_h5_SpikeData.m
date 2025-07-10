function out = load_h5_SpikeData(outputPath)

logFile = fullfile(outputPath, 'KIASort_GUI_log.txt');

fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');

% fields = {'upadatedLabels','labels','unifiedLabels','spike_idx', 'channelNum','features','amplitude'};
fields = {'unifiedLabels','spike_idx', 'channelNum','features','amplitude'};
out = cell2struct(cell(size(fields)), fields, 2);


sortedSpikeFolder = fullfile(outputPath, 'RES_Sorted');
if ~exist(sortedSpikeFolder,"dir")
    sortedSpikeFolder = outputPath;
end

for i = 1:numel(fields)
    fld = fields{i};
    h5File = [fld '.h5'];
    read_file = fullfile(sortedSpikeFolder,h5File);    
    if exist(read_file,"file")
    out.(fields{i}) = h5read(read_file, ['/' fld]);
    else
    fprintf(fid, ['Failed to load ', h5File, ' from sorted spike results.\n']);
    end
end

end
