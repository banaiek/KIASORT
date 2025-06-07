function plot_single_channel_samples(cfg, parentPanel, figColor)

sampleRes = load_sorted_samples_gui(cfg.outputFolder);
included_channels = sampleRes.channel_info.channel_inclusion;
channelList = 1:min(cfg.numChannels, length(included_channels));
channelList = channelList(included_channels);


samplePlotGrid = uigridlayout(parentPanel, [1 3], ...
    'ColumnWidth', {'2x', '20x', '1.5x'}, ...
    'Padding', [5 5 5 5], ...
    'ColumnSpacing', 1,...
    'RowSpacing',1);
applyColorScheme(samplePlotGrid, figColor);


timeGrid = uigridlayout(samplePlotGrid, [14 1], ...
    'RowHeight', {'fit', 'fit', 'fit', 'fit', 'fit', 'fit','fit','fit','fit', 'fit','2x','1x','fit','6x'}, ...
    'Padding', [5 5 5 5], ...
    'ColumnSpacing', 1,...
    'RowSpacing',1);
timeGrid.Layout.Row = 1;
timeGrid.Layout.Column = 1;
applyColorScheme(timeGrid, figColor);

samplesPlotPanel = uipanel(samplePlotGrid,'BackgroundColor',figColor,'BorderColor',figColor);
samplesPlotPanel.Layout.Row = 1;
samplesPlotPanel.Layout.Column = 2;
applyColorScheme(samplesPlotPanel, figColor);

checkboxPanel = uipanel(samplePlotGrid,'BackgroundColor',figColor,'BorderColor',figColor);
checkboxPanel.Layout.Row = 1;
checkboxPanel.Layout.Column = 3;
checkboxPanel.Scrollable = 'on';
applyColorScheme(checkboxPanel, figColor);

tileforSamples = tiledlayout(samplesPlotPanel, 3, 5, 'Padding', 'compact', 'TileSpacing', 'compact');
axUMAP_All = nexttile(tileforSamples, 1);
axPCA_All = nexttile(tileforSamples, 2);
axXCorr = nexttile(tileforSamples, 3);
axUMAP_Time = nexttile(tileforSamples, 6, [2, 1]);
axPCA_Time = nexttile(tileforSamples, 7, [2, 1]);
axSPKds = nexttile(tileforSamples, 8, [2, 1]);
axAvgWaveform = nexttile(tileforSamples, 4, [3, 1]);

axIndWaveforms = nexttile(tileforSamples, 5, [3, 1]);


% features on the left pannel
exportFormatDD = uidropdown(timeGrid,...
    'Items',{'Image','Vector'},...
    'Value','Image',...
    'ValueChangedFcn',@(dd,ev)updateExportFormat(dd.Value));
exportFormatDD.Layout.Row = 1;
exportFormatDD.Layout.Column = 1;
applyColorScheme(exportFormatDD, figColor);
exportType = 'Image';

btnExportFig = uibutton(timeGrid,'Text','Export Fig.',...
    'ButtonPushedFcn',@(btn,ev)onExportFig());
btnExportFig.Layout.Row = 2;
btnExportFig.Layout.Column = 1;
applyColorScheme(btnExportFig, figColor);


channelListCell = arrayfun(@num2str, channelList, 'UniformOutput',false);
ddSampleChannel = uidropdown(timeGrid,...
    'Items',channelListCell,...
    'Value',channelListCell (1), 'ValueChangedFcn',@(dd,ev)updateChannel(dd.Value));
ddSampleChannel.Layout.Row = 4;
ddSampleChannel.Layout.Column = 1;
applyColorScheme(ddSampleChannel, figColor);

lblTimeStepSize = uilabel(timeGrid,'Text','Channel #:',...
    'HorizontalAlignment','right');
lblTimeStepSize.Layout.Row = 3;
lblTimeStepSize.Layout.Column = 1;
applyColorScheme(lblTimeStepSize, figColor);

waveformItems = {'1', '10', '25', '50', '100', '1000', 'Max'};
waveformDropdown = uidropdown(timeGrid, 'Items', waveformItems, 'Value', '25', ...
    'FontSize', 14, 'ValueChangedFcn', @(dd,ev)updateWaveforms(dd.Value));
waveformDropdown.Layout.Row = 6;
waveformDropdown.Layout.Column = 1;
applyColorScheme(waveformDropdown, figColor);

lblTimeStepSize = uilabel(timeGrid,'Text','# Waveform:',...
    'HorizontalAlignment','right');
lblTimeStepSize.Layout.Row = 5;
lblTimeStepSize.Layout.Column = 1;
applyColorScheme(lblTimeStepSize, figColor);


timeStepDropdown = uidropdown(timeGrid,...
    'Items',{'1','10','100','1000','full'},...
    'Value','full');
timeStepDropdown.Layout.Row = 8;
timeStepDropdown.Layout.Column = 1;
applyColorScheme(timeStepDropdown, figColor);

lblTimeStepSize = uilabel(timeGrid,'Text','Time window (s):',...
    'HorizontalAlignment','right');
lblTimeStepSize.Layout.Row = 7;
lblTimeStepSize.Layout.Column = 1;
applyColorScheme(lblTimeStepSize, figColor);


imagePlotDropdown = uidropdown(timeGrid,...
    'Items',{'Max XCorr.','Confusion Mat.'},...
    'Value','Max XCorr.', 'ValueChangedFcn', @(dd,ev)updateImgPlot(dd.Value));
imagePlotDropdown.Layout.Row = 10;
imagePlotDropdown.Layout.Column = 1;
applyColorScheme(imagePlotDropdown, figColor);

lblTimeStepSize = uilabel(timeGrid,'Text','Image plot:',...
    'HorizontalAlignment','right');
lblTimeStepSize.Layout.Row = 9;
lblTimeStepSize.Layout.Column = 1;
applyColorScheme(lblTimeStepSize, figColor);

timeSlider = uislider(timeGrid,'Value',0,...
    'Orientation','vertical','MajorTicksMode','auto');
timeSlider.Layout.Row = 14;
timeSlider.Layout.Column = 1;
applyColorScheme(timeSlider, figColor);

lblTimeStepSize = uilabel(timeGrid,'Text','Time step (s):',...
    'HorizontalAlignment','left');
lblTimeStepSize.Layout.Row = 13;
lblTimeStepSize.Layout.Column = 1;
applyColorScheme(lblTimeStepSize, figColor);

timeSlider.ValueChangedFcn     = @(sld,ev)updateXWindow(sld.Value, timeStepDropdown.Value);
timeStepDropdown.ValueChangedFcn = @(dd,ev)updateXWindow(timeSlider.Value, dd.Value);





currentChannelPlot = channelList(1);

channelData = load_channel_sample_gui(cfg.outputFolder, currentChannelPlot);
sortedSamples  = sampleRes.sortedSamples{currentChannelPlot,1};
tempFeat    = sampleRes.sampleFeatures{currentChannelPlot,1};

numPlotChannels = 2 * sortedSamples.cfg.num_channel_extract + 1;
trialLength = sampleRes.channel_info.num_samples;
plotChannelNumber = channelData.wavformChanelIdx;
spk_wf = channelData.waveform_bp_full;

umapComp = tempFeat.umapNorm;
pcScore = tempFeat.PCA_scores(:,1:3);
id_all = tempFeat.labels;
spk_inds = tempFeat.spk_idx;
tmax = min(sortedSamples.cfg.maxSampleSpan, trialLength/sortedSamples.cfg.samplingFrequency);
idxmax_chunk = sortedSamples.cfg.numSampleChunks * sortedSamples.cfg.sampleChunkDuration * sortedSamples.cfg.samplingFrequency;
t_spk = (spk_inds / idxmax_chunk) * tmax;
figColor = [.2 .2 .2];
spikeDensity = sortedSamples.clusteringInfo.clusterRelabeling.clusterSpikeDensity;
spikeDensityTime = linspace(0,tmax,size(spikeDensity,2));
uniqueLabels = sortedSamples.clusteringInfo.clusterRelabeling.originalLabels;
XcorrVal = sortedSamples.clusteringInfo.clusterRelabeling.maxXcorrVal;
channelMax = sortedSamples.clusteringInfo.clusterSelection.max_channel;
channelMaxLabels = sortedSamples.clusteringInfo.clusterSelection.classLabels;
updatedClusterLabels = sortedSamples.clusteringInfo.clusterRelabeling.newLabels;
waveformXaxis = -sortedSamples.cfg.spikeDuration/2 : 1000/(sortedSamples.cfg.samplingFrequency) : sortedSamples.cfg.spikeDuration/2;
confusionMat = sortedSamples.classifierInfo.valAccuracy;
confusionMat = 100 * confusionMat ./ sum(confusionMat, 2);


for i=1:length(uniqueLabels)
    updatedChannelMax(i,1) = channelMax(channelMaxLabels == updatedClusterLabels(i));
end

xWindowStart = 0;
xWindowEnd   = tmax;
xStepVal     = tmax;
timeSlider.Limits = [0 tmax];
imagePlotType = 'Max XCorr.';

maxVal = max(abs(spk_wf), [], [2,3]);
spk_wf(maxVal > prctile(maxVal, 99.9), :, :) = 0;

unique_ids = unique(id_all);
invalid_ids = (unique_ids == -1 | unique_ids == 0);
unique_ids(invalid_ids) = [];
unique_ids = unique_ids';
N = numel(unique_ids);
rgb_colors = generate_colors(N);
id_map = containers.Map(num2cell(unique_ids), num2cell(1:N));
special_ids = -1;
selectionOrder = [];
merged_id = updatedClusterLabels(~invalid_ids)';
groupNames = arrayfun(@num2str, unique_ids, 'UniformOutput', false);
mergedNames = arrayfun(@num2str, merged_id, 'UniformOutput', false);

if any(ismember(special_ids, unique(id_all)))
    groupNames = [strcat(groupNames, '-', mergedNames), {'-1'}];
else
    groupNames = strcat(groupNames, '-', mergedNames);
end
numGroups = length(groupNames);

gl = uigridlayout(checkboxPanel, [numGroups + 2, 1], ...
    'RowHeight', repmat({'fit'}, 1, numGroups + 2), ...
    'ColumnWidth', {'1x'}, ...
    'Padding', [10 10 10 10], ...
    'Scrollable','on',...
    'RowSpacing', 10);
applyColorScheme(gl, figColor);


chkAll = uicheckbox(gl, 'Text', 'All', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Value', true, 'ValueChangedFcn', @(src, event) selectAllGroups(src, event, gl));
addGroupCheckboxes(gl, numGroups, groupNames, rgb_colors);
applyColorScheme(chkAll, figColor);

numWaveforms = 100;

scatter3DPlot(axUMAP_All, umapComp, id_all, rgb_colors, id_map, special_ids, 'UMAP Scatter', [], []);
set(axUMAP_All, 'View', [30, 30]);
applyColorScheme(axUMAP_All, figColor);

scatter3DPlot(axPCA_All, pcScore, id_all, rgb_colors, id_map, special_ids, 'PCA Scatter', [], []);
set(axPCA_All, 'View', [30, 30]);


scatter3DPlot(axUMAP_Time, umapComp, id_all, rgb_colors, id_map, special_ids, 'UMAP Time Resolved', t_spk, []);
set(axUMAP_Time, 'View', [30, 30]);
scatter3DPlot(axPCA_Time, pcScore, id_all, rgb_colors, id_map, special_ids, 'PCA Time Resolved', t_spk, []);
set(axPCA_Time, 'View', [30, 30]);
plotAverageWaveforms(axAvgWaveform, spk_wf, id_all, rgb_colors, id_map, special_ids, 'Average Waveforms');
plotIndividualWaveforms(axIndWaveforms, spk_wf, id_all, rgb_colors, id_map, special_ids, numWaveforms, 'Individual Waveforms');
plotSpikeDensity(axSPKds, spikeDensity, uniqueLabels, spikeDensityTime, rgb_colors, id_map, special_ids, 'Spike Density');
plotXCorr('Max Xcorr R')




    function updateAllPlots()

        channelData = load_channel_sample_gui(cfg.outputFolder, currentChannelPlot);
        sortedSamples  = sampleRes.sortedSamples{currentChannelPlot,1};
        tempFeat    = sampleRes.sampleFeatures{currentChannelPlot,1};

        trialLength = sampleRes.channel_info.num_samples;
        plotChannelNumber = channelData.wavformChanelIdx;
        spk_wf = channelData.waveform_bp_full;        
        umapComp = tempFeat.umapNorm;
        pcScore = tempFeat.PCA_scores(:,1:3);
        id_all = tempFeat.labels;
        spk_inds = tempFeat.spk_idx;
        tmax = min(cfg.maxSampleSpan, trialLength/cfg.samplingFrequency);
        idxmax_chunk = cfg.numSampleChunks * cfg.sampleChunkDuration * cfg.samplingFrequency;
        t_spk = (spk_inds / idxmax_chunk) * tmax;
        figColor = [.2 .2 .2];
        spikeDensity = sortedSamples.clusteringInfo.clusterRelabeling.clusterSpikeDensity;
        spikeDensityTime = linspace(0,tmax,size(spikeDensity,2));
        uniqueLabels = sortedSamples.clusteringInfo.clusterRelabeling.originalLabels;
        XcorrVal = sortedSamples.clusteringInfo.clusterRelabeling.maxXcorrVal;
        channelMax = sortedSamples.clusteringInfo.clusterSelection.max_channel;
        channelMaxLabels = sortedSamples.clusteringInfo.clusterSelection.classLabels;
        updatedClusterLabels = sortedSamples.clusteringInfo.clusterRelabeling.newLabels;
        confusionMat = sortedSamples.classifierInfo.valAccuracy;
        confusionMat = 100 * confusionMat ./ sum(confusionMat, 2);

        for i=1:length(uniqueLabels)
            updatedChannelMax(i,1) = channelMax(channelMaxLabels == updatedClusterLabels(i));
        end

        xWindowStart = 0;
        xWindowEnd   = tmax;
        xStepVal     = tmax;
        timeSlider.Limits = [0 tmax];

        maxVal = max(abs(spk_wf), [], [2,3]);
        spk_wf(maxVal > prctile(maxVal, 99.9), :, :) = 0;

        unique_ids = unique(id_all);
        invalid_ids = (unique_ids == -1 | unique_ids == 0);
        unique_ids(invalid_ids) = [];
        unique_ids = unique_ids';
        N = numel(unique_ids);
        rgb_colors = generate_colors(N);
        id_map = containers.Map(num2cell(unique_ids), num2cell(1:N));
        special_ids = -1;
        selectionOrder = [];
        merged_id = updatedClusterLabels(~invalid_ids)';
        groupNames = arrayfun(@num2str, unique_ids, 'UniformOutput', false);
        mergedNames = arrayfun(@num2str, merged_id, 'UniformOutput', false);

        if any(ismember(special_ids, unique(id_all)))
            groupNames = [strcat(groupNames, '-', mergedNames), {'-1'}];
        else
            groupNames = strcat(groupNames, '-', mergedNames);
        end
        numGroups = length(groupNames);

        delete(allchild(checkboxPanel));

        gl = uigridlayout(checkboxPanel, [numGroups + 2, 1], ...
            'RowHeight', repmat({'fit'}, 1, numGroups + 2), ...
            'ColumnWidth', {'1x'}, ...
            'Padding', [10 10 10 10], ...
            'Scrollable','on',...
            'RowSpacing', 10);
        applyColorScheme(gl, figColor);


        chkAll = uicheckbox(gl, 'Text', 'All', 'FontSize', 14, 'FontWeight', 'bold', ...
            'Value', true, 'ValueChangedFcn', @(src, event) selectAllGroups(src, event, gl));
        addGroupCheckboxes(gl, numGroups, groupNames, rgb_colors);
        applyColorScheme(chkAll, figColor);



        scatter3DPlot(axUMAP_All, umapComp, id_all, rgb_colors, id_map, special_ids, 'UMAP Scatter', [], []);
        set(axUMAP_All, 'View', [30, 30]);
        applyColorScheme(axUMAP_All, figColor);

        scatter3DPlot(axPCA_All, pcScore, id_all, rgb_colors, id_map, special_ids, 'PCA Scatter', [], []);
        set(axPCA_All, 'View', [30, 30]);


        scatter3DPlot(axUMAP_Time, umapComp, id_all, rgb_colors, id_map, special_ids, 'UMAP Time Resolved', t_spk, []);
        set(axUMAP_Time, 'View', [30, 30]);
        scatter3DPlot(axPCA_Time, pcScore, id_all, rgb_colors, id_map, special_ids, 'PCA Time Resolved', t_spk, []);
        set(axPCA_Time, 'View', [30, 30]);
        plotAverageWaveforms(axAvgWaveform, spk_wf, id_all, rgb_colors, id_map, special_ids, 'Average Waveforms');
        plotIndividualWaveforms(axIndWaveforms, spk_wf, id_all, rgb_colors, id_map, special_ids, numWaveforms, 'Individual Waveforms');
        plotSpikeDensity(axSPKds, spikeDensity, uniqueLabels, spikeDensityTime, rgb_colors, id_map, special_ids, 'Spike Density');
        plotXCorr('Max Xcorr R')
    end

    function selectAllGroups(src, event, gridLayout)

        chkGroups = findall(gridLayout, 'Type', 'uicheckbox', 'Tag', 'groupCheckbox');

        if src.Value
            for j = 1:length(chkGroups)
                if ~chkGroups(j).Value
                    chkGroups(j).Value = true;
                    gid = extractGroupID(chkGroups(j).Text);
                    if ~ismember(gid, selectionOrder)
                        selectionOrder(end+1) = gid(end);
                    end
                end
            end
            if any(strcmp(groupNames, '-1'))
                idx = find(strcmp(groupNames, '-1'));
                if ~chkGroups(idx).Value && any(ismember(special_ids, unique(id_all)))
                    chkGroups(idx).Value = true;
                    selectionOrder = [selectionOrder, special_ids];
                end
            end
        else
            for j = 1:length(chkGroups)
                if chkGroups(j).Value
                    chkGroups(j).Value = false;
                    gid = extractGroupID(chkGroups(j).Text);
                    selectionOrder(selectionOrder == gid) = [];
                end
            end
            if any(strcmp(groupNames, '-1'))
                idx = find(strcmp(groupNames, '-1'));
                if chkGroups(idx).Value
                    chkGroups(idx).Value = false;
                    selectionOrder(selectionOrder == -1 | selectionOrder == 0) = [];
                end
            end
        end

        selectionOrder = unique(selectionOrder, 'stable');
        updatePlots(src);
    end

    function gid = extractGroupID(groupName)
        if strcmp(groupName, '-1')
            gid = [-1, 0];
        else
            tokens = regexp(groupName, '\d+', 'match');
            if ~isempty(tokens)
                gid = str2double(tokens{1});
            else
                gid = [];
            end
        end
    end

    function updatePlots(src, ~)


        if isgraphics(src, 'uicheckbox')
            groupName = src.Text;
            if ~strcmp(groupName, '-1')
                gid = str2double(groupName(1:strfind(groupName,'-')-1));
                if src.Value
                    if ~ismember(gid, selectionOrder)
                        selectionOrder(end+1) = gid;
                    end
                else
                    selectionOrder(selectionOrder == gid) = [];
                end
            else
                if src.Value
                    if ~any(ismember(special_ids, selectionOrder))
                        selectionOrder = [selectionOrder, special_ids];
                    end
                else
                    selectionOrder(selectionOrder == -1 | selectionOrder == 0) = [];
                end
            end
        end

        chkAll = findall(checkboxPanel, 'Type', 'uicheckbox', 'Text', 'All');
        chkGroups = findall(checkboxPanel, 'Type', 'uicheckbox', 'Tag', 'groupCheckbox');
        if any(~[chkGroups.Value])
            chkAll.Value = false;
        else
            chkAll.Value = true;
        end


        scatter3DPlot(axUMAP_All, umapComp, id_all, rgb_colors, id_map, special_ids, 'UMAP Scatter', [],selectionOrder);
        scatter3DPlot(axPCA_All, pcScore, id_all, rgb_colors, id_map, special_ids, 'PCA Scatter', [],selectionOrder);

        if strcmp(imagePlotType,'Max XCorr.')
            plotXCorr('Max Xcorr R')
        else
            plotConfusion('Classifier Confusion %')
        end

        replotifUpdated();
    end

    function replotifUpdated(~, ~)

        timeIndices = find(t_spk >= xWindowStart & t_spk <= xWindowEnd);
        selectedGroups = unique(selectionOrder, 'stable');

        scatter3DPlot(axUMAP_Time, umapComp(timeIndices, :), id_all(timeIndices), rgb_colors, id_map, special_ids, ...
            'UMAP Time Resolved', t_spk(timeIndices), selectedGroups);
        scatter3DPlot(axPCA_Time, pcScore(timeIndices, :), id_all(timeIndices), rgb_colors, id_map, special_ids, ...
            'PCA Time Resolved', t_spk(timeIndices), selectedGroups);
        plotAverageWaveforms(axAvgWaveform, spk_wf(timeIndices, :, :), id_all(timeIndices), rgb_colors, id_map, special_ids, ...
            'Average Waveforms');
        plotIndividualWaveforms(axIndWaveforms, spk_wf(timeIndices, :, :), id_all(timeIndices), rgb_colors, id_map, special_ids, ...
            numWaveforms, 'Individual Waveforms');
        plotSpikeDensity(axSPKds, spikeDensity, uniqueLabels, spikeDensityTime, rgb_colors, id_map, special_ids, 'Spike Density');
    end

    function updateWaveforms(newVal)
        numWaveforms = str2double(newVal);
        replotifUpdated();
    end

    function updateImgPlot(newVal)
        imagePlotType = newVal;
        if strcmp(imagePlotType,'Max XCorr.')
            plotXCorr('Max Xcorr R')
        else
            plotConfusion('Classifier Confusion %')
        end
    end

    function updateChannel(newVal)
        updatedChannelPlot = str2double(newVal);
        if updatedChannelPlot ~= currentChannelPlot
            currentChannelPlot = updatedChannelPlot;
            updateAllPlots();
        end
    end

    function save_callback(~, ~)
        [file, path] = uiputfile({'*.svg'; '*.eps'}, 'Save Figure As');
        if isequal(file,0)
            return;
        end
        [~, ~, ext] = fileparts(file);
        fullFileName = fullfile(path, file);
        switch lower(ext)
            case '.svg'
                saveas(gcf, fullFileName, 'svg');
            case '.eps'
                saveas(gcf, fullFileName, 'epsc');
            otherwise
                uialert(parentPanel.Parent, 'Unsupported file format.', 'Error');
        end
    end

    function scatter3DPlot(ax, data, labels, colors, id_map, special_ids, titleStr, t_spk_data, groupsToPlot)
        cla(ax);
        hold(ax, 'on');
        if isempty(groupsToPlot) && chkAll.Value
            groupsToPlot = unique(labels);
        end
        for i = 1:length(groupsToPlot)
            gid = groupsToPlot(i);
            idx = labels == gid;
            if any(idx)
                if ismember(gid, special_ids)
                    color = [0.5, 0.5, 0.5];
                else
                    color = colors(id_map(gid), :);
                end
                if ~isempty(t_spk_data)
                    if size(data, 2) >= 2
                        scatter3(ax, data(idx,1), data(idx,2), t_spk_data(idx), 24, color, 'filled');
                    else
                        scatter(ax, data(idx,1), t_spk_data(idx), 24, color, 'filled');
                    end
                else
                    if size(data, 2) >= 3
                        scatter3(ax, data(idx,1), data(idx,2), data(idx,3), 24, color, 'filled');
                    elseif size(data,2) >=2
                        scatter(ax, data(idx,1), data(idx,2), 24, color, 'filled');
                    else
                        scatter(ax, data(idx,1), zeros(size(data(idx,1))), 24, color, 'filled');
                    end
                end
            end
        end
        hold(ax, 'off');
        title(ax, titleStr, 'Color', 1-figColor);
        xlabel(ax, 'Comp. 1');
        ylabel(ax, 'Comp. 2');
        if ~isempty(t_spk_data)
            zlabel(ax, 'Time');
        elseif size(data, 2) >= 3
            zlabel(ax, 'Comp. 3');
        end
        set(ax, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor, 'ZColor',1-figColor, 'TickDir','out','box','on')
    end

    function plotAverageWaveforms(ax, waveforms, labels, colors, id_map, special_ids, titleStr)
        num_channels = size(waveforms,2);

        b0 = 5 * (num_channels:-1:1);
        cla(ax);
        hold(ax, 'on');

        groupsToPlotLabels = [];
        if isempty(selectionOrder) && chkAll.Value
        groupsToPlotLabels = uniqueLabels;
        elseif ~isempty(selectionOrder)
        groupsToPlotLabels = selectionOrder(ismember(selectionOrder,uniqueLabels));
        end

        for i = 1:length(groupsToPlotLabels)
            gid = groupsToPlotLabels(i);
            idx = labels == gid;
            if any(idx)
                avgWaveform = squeeze(mean(waveforms(idx, :, :),1,"omitmissing"));
                avgWaveform = .2 * avgWaveform / max(mean(abs(waveforms), 'all'));
                if ismember(gid, special_ids)
                    color = [0.5, 0.5, 0.5];
                else
                    color = colors(id_map(gid), :);
                end
                plot(ax, waveformXaxis, (b0' + avgWaveform)', 'Color', color, 'LineWidth', 2);
                channelMax_idx = updatedChannelMax(uniqueLabels==gid);
                if channelMax_idx <= numPlotChannels && channelMax_idx >= 1
                    plot(ax, waveformXaxis, (b0(channelMax_idx) + avgWaveform(channelMax_idx,:))', 'Color', color, 'LineWidth', 4);
                end
            end
        end
        hold(ax, 'off');
        title(ax, titleStr,"Color",1-figColor);
        xlabel(ax, 'Time');
        ylabel(ax, 'Channel');
        set(ax, 'YTick', fliplr(b0), 'YTickLabel', flipud(plotChannelNumber));
        set(ax, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor, 'TickDir','out',...
            'ylim',[b0(end)-5 b0(1)+5])


    end

    function plotIndividualWaveforms(ax, waveforms, labels, colors, id_map, special_ids, numWaveforms, titleStr)
        num_channels = size(waveforms,2);
        b0 = 5 * (num_channels:-1:1);
        cla(ax);
        hold(ax, 'on');
        
        groupsToPlotLabels = [];
        if isempty(selectionOrder) && chkAll.Value
        groupsToPlotLabels = uniqueLabels;
        elseif ~isempty(selectionOrder)
        groupsToPlotLabels = selectionOrder(ismember(selectionOrder,uniqueLabels));
        end

        for i = 1:length(groupsToPlotLabels)
            gid = groupsToPlotLabels(i);
            idx = find(labels == gid);
            if ~isempty(idx)
                numAvailable = length(idx);
                if numWaveforms < numAvailable
                    selectedIdx = idx(randperm(numAvailable, numWaveforms));
                else
                    selectedIdx = idx;
                end
                if ismember(gid, special_ids)
                    color = [0.5, 0.5, 0.5];
                else
                    color = colors(id_map(gid), :);
                end
                for j = selectedIdx'
                    indvWaveform = squeeze(waveforms(j, :, :));
                    if ~any(sum(indvWaveform, "all"))
                        continue;
                    end
                    indvWaveform = .2 *indvWaveform / max(mean(abs(waveforms), 'all'));
                    plot(ax, waveformXaxis, (b0' + indvWaveform)', 'Color', color);
                end
            end
        end
        hold(ax, 'off');
        title(ax, titleStr,"Color",1-figColor);
        xlabel(ax, 'Time');
        ylabel(ax, 'Channel');
        set(ax, 'YTick', fliplr(b0), 'YTickLabel', flipud(plotChannelNumber));
        set(ax, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor, 'TickDir','out',...
            'ylim',[b0(end)-5 b0(1)+5])
    end

    function rgb_colors = generate_colors(N)
        hues = linspace(0, 1, N+1);
        hues = hues(1:N);
        saturation = 0.8;
        value = 0.8;
        hsv_colors = [hues', repmat(saturation, N, 1), repmat(value, N, 1)];
        rgb_colors = hsv2rgb(hsv_colors);
    end

    function addGroupCheckboxes(gridLayout, numGroups, groupNames, rgb_colors)
        for i = 1:numGroups
            if ~strcmp(groupNames{i}, '-1')
                bgColor = rgb_colors(i, :);
            else
                bgColor = [0.5, 0.5, 0.5];
            end
            uicheckbox(gridLayout, 'Text', groupNames{i}, 'FontSize', 18, 'FontWeight', 'bold', ...
                'Value', false, 'ValueChangedFcn', @updatePlots, ...
                'FontColor', bgColor ,'Tag', 'groupCheckbox');
        end
    end

    function updateXWindow(sVal, stepStr)
        if strcmp(stepStr,'full')
            xStepVal = tmax;
        else
            xStepVal = str2double(stepStr);
        end
        xWindowStart = round(sVal);
        xWindowEnd   = xWindowStart + xStepVal;
        if xWindowEnd > timeSlider.Limits(2)
            xWindowEnd   = timeSlider.Limits(2);
            xWindowStart = xWindowEnd - xStepVal;
            if xWindowStart < timeSlider.Limits(1)
                xWindowStart = timeSlider.Limits(1);
            end
            timeSlider.Value = xWindowStart;
        end
        disp(['X-axis window => Start=',num2str(xWindowStart),...
            ', End=',num2str(xWindowEnd),...
            ', Step=',num2str(xStepVal)]);
        replotifUpdated();
    end

    function plotSpikeDensity(ax, spikeDensity, labels, spikeDensityTime, colors, id_map, special_ids, titleStr)
        cla(ax);
        hold(ax, 'on');
        densityYaxis = spikeDensityTime<xWindowEnd & spikeDensityTime>xWindowStart;
        densityXaxis = spikeDensity(:,densityYaxis);
        densityXaxis = densityXaxis./max(densityXaxis,[],2);

        groupsToPlotLabels = [];
        if isempty(selectionOrder) && chkAll.Value
        groupsToPlotLabels = uniqueLabels;
        elseif ~isempty(selectionOrder)
        groupsToPlotLabels = selectionOrder(ismember(selectionOrder,uniqueLabels));
        end

        if ~isempty(groupsToPlotLabels)

            for i = 1:length(groupsToPlotLabels)
                gid = groupsToPlotLabels(i);
                idx = find(labels == gid);
                if ~isempty(idx)
                    if ismember(gid, special_ids)
                        color = [0.5, 0.5, 0.5];
                    else
                        color = colors(id_map(gid), :);
                    end
                    plot(ax, densityXaxis(idx,:), spikeDensityTime(densityYaxis), 'Color', color,'LineWidth',2);
                end
            end
        end
        hold(ax, 'off');
        title(ax, titleStr,"Color",1-figColor);
        xlabel(ax, 'Norm. Spike Density');
        ylabel(ax, 'Time (s)');
        set(ax, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor, 'TickDir','out')

    end

    function plotXCorr(titleStr)
        cla(axXCorr);
        hold(axXCorr, 'on');
        groupsToPlotLabels = [];
        if isempty(selectionOrder) && chkAll.Value
        groupsToPlotLabels = uniqueLabels;
        elseif ~isempty(selectionOrder)
        groupsToPlotLabels = selectionOrder(ismember(selectionOrder,uniqueLabels));
        end

        if ~isempty(groupsToPlotLabels)
            idx = ismember(uniqueLabels, groupsToPlotLabels);
            imagesc(axXCorr, XcorrVal(idx,idx));
            cmap = polar_colormap;
            colormap(axXCorr, cmap)
            colorbar(axXCorr,'northoutside','Color',1-figColor)
            hold(axXCorr, 'off');
            title(axXCorr, titleStr,"Color",1-figColor);
            xlabel(axXCorr, 'Group #');
            ylabel(axXCorr, 'Group #');
            set(axXCorr, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor, 'TickDir','out',...
                'xlim',[1-.5 length(groupsToPlotLabels)+.5], 'ylim',[1-.5 length(groupsToPlotLabels)+.5], 'clim',[-1 1])
            set(axXCorr, 'YTick', 1:length(groupsToPlotLabels), 'YTickLabel', sort(groupsToPlotLabels), 'XTick', 1:length(groupsToPlotLabels), 'XTickLabel', sort(groupsToPlotLabels));
        else
            set(axXCorr, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor)
        end

    end

    function plotConfusion(titleStr)
        cla(axXCorr);
        hold(axXCorr, 'on');
        groupsToPlotLabels = [];
        if isempty(selectionOrder) && chkAll.Value
        groupsToPlotLabels = uniqueLabels(uniqueLabels>0);
        elseif ~isempty(selectionOrder)
        groupsToPlotLabels = selectionOrder(ismember(selectionOrder,uniqueLabels) & selectionOrder>0);
        end
        if ~isempty(groupsToPlotLabels)
            idx = ismember(uniqueLabels(uniqueLabels>0), groupsToPlotLabels);            
            imagesc(axXCorr, confusionMat(idx,idx));
            cmap = polar_colormap;
            colormap(axXCorr, cmap)
            colorbar(axXCorr,'northoutside','Color',1-figColor, 'Limits',[0 100])
            hold(axXCorr, 'off');
            title(axXCorr, titleStr,"Color",1-figColor);
            xlabel(axXCorr, 'Group #');
            ylabel(axXCorr, 'Group #');
            set(axXCorr, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor, 'TickDir','out',...
                'xlim',[1-.5 length(groupsToPlotLabels)+.5], 'ylim',[1-.5 length(groupsToPlotLabels)+.5], 'clim',[-100 100])
            set(axXCorr, 'YTick', 1:length(groupsToPlotLabels), 'YTickLabel', sort(groupsToPlotLabels), 'XTick', 1:length(groupsToPlotLabels), 'XTickLabel', sort(groupsToPlotLabels));
        else
            set(axXCorr, 'Color', figColor-.05, 'XColor',1-figColor, 'YColor',1-figColor)
        end

    end

    function onExportFig()
        saveFigPath = fullfile(cfg.outputFolder, 'exported_figures');
        if ~exist("saveFigPath",'dir')
            mkdir(saveFigPath);
        end

        channelBaseName =  sprintf('exported_Fig_Channel_%d_', currentChannelPlot);
        channelName = fullfile(saveFigPath, channelBaseName);

        
        files = dir([channelName, '*.eps']);

        if isempty(files)
            nextNum = 1;
        else
            nums = cellfun(@(s) sscanf(s, [channelBaseName,'%d.eps']), {files.name});
            nextNum = max(nums) + 1;
        end

        filename = fullfile(saveFigPath, sprintf([channelBaseName, '%d.eps'], nextNum));
        exportgraphics(tileforSamples, filename,'BackgroundColor','white','ContentType', exportType);
    end

    function updateExportFormat(ddVal)
        exportType = ddVal;
    end
end