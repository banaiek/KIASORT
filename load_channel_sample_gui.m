function data = load_channel_sample_gui(outputPath, ch)


logFile = fullfile(outputPath, 'KIASort_GUI_log.txt');

fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');


inputFolder = fullfile(outputPath, 'RES_Samples');
inputFilename = fullfile(inputFolder, sprintf('channel_%d_results.mat', ch));
try
    data = load(inputFilename,'out').out;
    fprintf(fid, 'Samples from channel %d loaded. %s\n',ch, datetime);
catch ME
    fprintf(fid, 'ERROR: Failed to load channel %d sample from %s: %s.\n', ch, inputFilename, ME.message);
    sprintf('ERROR: Failed to load channel %d sample from %s: %s.\n', ch, inputFilename, ME.message);
end


end