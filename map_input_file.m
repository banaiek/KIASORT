 function m = map_input_file(filename, cfg)
 
 if ~isfield(cfg, 'outputFolder')
     cfg.outputFolder = cfg.inputFolder;
 end

 logFile = fullfile(cfg.outputFolder, 'KIASort_log.txt');
 fid = fopen(logFile, 'a');
 if fid < 0
     error('Could not open log file: %s', logFile);
 end
 fprintf(fid, '\n=============================\n');
 fprintf(fid, 'Data loaded and mapped at %s\n', datetime);

        try
            fileinfo = dir(filename);
            total_bytes = fileinfo.bytes;
            total_samples = total_bytes / 2;
            num_samples = total_samples / cfg.numChannels;
            if mod(num_samples, 1) ~= 0
                fprintf(fid, 'ERROR: The number of samples per channel is not an integer. Check file size and channel count.\n');
                errorFlag = true;
                error('The number of samples per channel is not an integer. Check file size and channel count.');
            end
            num_samples = double(uint64(num_samples));
            m = memmapfile(filename, 'Format', {cfg.dataType, [cfg.numChannels, num_samples], 'data'});
            fprintf(fid, 'Memory-mapped file %s successfully.\n', filename);
        catch ME
            fprintf(fid, 'ERROR: Failed to memory-map file %s: %s\n', filename, ME.message);
            if ~isempty(ME.stack)
                fprintf(fid, '       In %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
            end
            errorFlag = true;
            error('Failed to memory-map file %s. Check the log file for details.', filename);
        end

end
