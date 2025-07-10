function currentSize = saveh5SpikeData_gui(outputFolder, sorted_out, prevSize)


if prevSize ==0
    files = dir(fullfile(outputFolder, '*curated*.h5'));
for k = 1:length(files)
    fileToDelete = fullfile(outputFolder, files(k).name);
    delete(fileToDelete);
    fprintf('Deleted file: %s\n', fileToDelete);
end
end


% subfield names
fields = fieldnames(sorted_out);
currentSize = zeros(size(fields));

%  write each subfield to HDF5
for i = 1:numel(fields)
    fld = fields{i};
    data = cell2mat(sorted_out.(fld));
    h5File = [fld '.h5'];
    output_file = fullfile(outputFolder,h5File);
    
    currentSize(i) = length(data);
    
    
       dataSize = size(data);
       dataSize(1)  = inf;
       dataCount    = ones(1, ndims(data));
       
       dataCount(1) = dataCount(1) + prevSize(i);
       ChunkSize = dataSize;       
       ChunkSize(1) = min(size(data,1),100);
       datasetName = ['/' fld];
    if ~exist(output_file,"file") && ~isempty(data)
    h5create(output_file, datasetName, dataSize, 'ChunkSize', ChunkSize);
        if dataCount(1)>prevSize(i)
    h5write(output_file, ['/' fld], data, dataCount, size(data));
        end
    elseif ~isempty(data)
        if dataCount(1)>prevSize(i)
    h5write(output_file, ['/' fld], data, dataCount, size(data));
        end
    end
end

end