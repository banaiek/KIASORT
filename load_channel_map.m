function [channel_mapping, channel_locations, channel_inclusion] = load_channel_map(inputPath, cfg)

if ~isfield(cfg,'outputFolder') && isfield(cfg,'inputFolder')
    cfg.outputFolder = cfg.inputFolder;
end

logFile = fullfile(cfg.outputFolder, 'KIASort_log.txt');
fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');



data = load(inputPath);
if isfield(data,'channel_mapping')
    channel_mapping = data.channel_mapping;
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'Channel map loaded! %s\n', datetime);
    end
elseif isfield(data,'chanMap')
    channel_mapping = data.chanMap;
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'Channel map loaded! %s\n', datetime);
    end
elseif isfield(data,'chanMap0ind')
    channel_mapping = data.chanMap0ind;
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'Channel map loaded! %s\n', datetime);
    end
else
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'No channel map was found! %s\n', datetime);
    end
    channel_mapping = 1 : cfg.numChannels;
end

if min(channel_mapping) == 0
    channel_mapping = channel_mapping +1;
end

if isfield(data,'channel_inclusion')
    channel_inclusion = data.channel_inclusion;
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'Channel inclusion info loaded! %s\n', datetime);
    end
elseif isfield(data,'connected')
    channel_inclusion = data.connected;
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'Channel inclusion info loaded! %s\n', datetime);
    end
else
    channel_inclusion = true(cfg.numChannels,1);
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'No channel inclusion info was found! %s\n', datetime);
    end
end

if isfield(data,'channel_locations')
    channel_locations = data.channel_locations;
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'Channel coordinates loaded! %s\n', datetime);
    end
elseif isfield(data,'xcoords') && isfield(data,'ycoords')
    channel_locations = [data.xcoords, data.ycoords, zeros(size(data.xcoords))];
    if isfield(data,'shankInd')
        channel_locations(:,3) = data.shankInd;
    end
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'Channel coordinates loaded! %s\n', datetime);
    end
else
    if isfield(cfg,'outputFolder')
        fprintf(fid, 'No channel coordinates were found! %s\n', datetime);
    end
    channel_locations      = zeros(cfg.numChannels, 3);
    channel_locations(:,2) = 1 : cfg.numChannels;
end

end