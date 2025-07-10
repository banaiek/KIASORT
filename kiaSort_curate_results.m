function kiaSort_curate_results(cfg, parentPanel, figColor, parentFig)

logFile = fullfile(cfg.outputFolder, 'KIASort_GUI_log.txt');
fid = fopen(logFile, 'a');
if fid < 0
    error('Could not open log file: %s', logFile);
end
fprintf(fid, '\n=============================\n');
fprintf(fid, 'Sorting explorer was called \n');

sampleRes = load_sorted_samples_gui(cfg.outputFolder);

mappedData = [];

if isfield(cfg, 'altResFolder')
    if ~isempty(cfg.altResFolder)
        sortedRes   = load_h5_SpikeData(cfg.altResFolder);

        if isempty(mappedData)
            if isfield(cfg, 'inputFolder') && isfield(cfg, 'outputFolder')
                mappedData =  map_input_file(cfg.fullFilePath, cfg);
                num_Samples = size(mappedData.Data.data,2);
                trialLength = num_Samples/cfg.samplingFrequency;
            else
                disp('--- Data loading failed ---');
            end
        end
        

    else
        sortedRes   = load_h5_SpikeData(cfg.outputFolder);
        num_Samples = sampleRes.channel_info.num_samples;
        trialLength = sampleRes.channel_info.num_samples/(cfg.samplingFrequency);
    end
else
    sortedRes   = load_h5_SpikeData(cfg.outputFolder);
    num_Samples = sampleRes.channel_info.num_samples;
    trialLength = sampleRes.channel_info.num_samples/(cfg.samplingFrequency);
end

groupList   = sampleRes.crossChannelStats.unified_labels.label;
channelList = sampleRes.crossChannelStats.unified_labels.channelID;
sampleWaveform = sampleRes.crossChannelStats.unified_labels.meanWaveforms;
detectblity = sampleRes.crossChannelStats.unified_labels.detectblity;
mainPolarity = sampleRes.crossChannelStats.unified_labels.mainNegativePolarity;
sidePolarity = sampleRes.crossChannelStats.unified_labels.sideNegativePolarity;
numChannelPlot = cfg.num_channel_extract;
channelPlot = channelList + [-numChannelPlot : numChannelPlot];
channel_mapping = sampleRes.channel_info.channel_mapping;
chan_wave_inclusion = sampleRes.channel_info.chan_wave_inclusion;
numGroups   = length(groupList);
halfSpikeWaveDur = cfg.spikeDuration/2;
numSpikeSamples = round(halfSpikeWaveDur * cfg.samplingFrequency /1000);
spike_Xaxis = -numSpikeSamples:numSpikeSamples;
waveformXaxis = 1000 *spike_Xaxis/cfg.samplingFrequency;
xlocs = sampleRes.channel_info.channel_locations(:,1);
ylocs = sampleRes.channel_info.channel_locations(:,2);
mergedFlag = false(numGroups,1);
xNorm_adj = [];
yNorm_adj = [];
plotType = [];
selected_order = [];
selected_order_last = [];

% Initialize sparse matrices/cells for large data structures to improve memory usage
xcorr_vals = cell(numGroups,numGroups);
isiCounts = cell(numGroups,numGroups);
isiCenters = cell(numGroups,numGroups);
isiViolations = zeros(numGroups,numGroups);
unitIsolation = repmat({'NA'},[numGroups,1]);
DenCenters = cell(numGroups,numGroups);
DenCounts = cell(numGroups,numGroups);
presence_ratio = zeros(numGroups,numGroups);
numChannels = min(length(chan_wave_inclusion),length(channel_mapping));
selectUnits = groupList;
chanLabel = cell(numChannels, 1);
stable_length = [ones(numGroups,1), trialLength * ones(numGroups,1)];
maxSignal = [];

chkBoxSelect = false(numGroups,1);
chkBoxVis = false(numGroups,1);

% Pre-generate channel labels
chanLabel = arrayfun(@(x) sprintf('Ch %d',x), 1:numChannels, 'UniformOutput', false)';

% Initialize pagination variables for group list display
PAGE_SIZE = 50; % Number of groups to display per page
currentPage = 1;
totalPages = ceil(numGroups/PAGE_SIZE);

ccgLag = 100;
smoothN = 0;
thresholdISI = 1;

minDistanceY = [];
minDistanceX = [];
normalizeTimeAmp();

originalGroupList = groupList;
originalChannelList = channelList;
originalSampleWaveform = sampleWaveform;
originalUnifiedLabels = sortedRes.unifiedLabels;
originalChannelNum = sortedRes.channelNum;
originalSpikeIdx = sortedRes.spike_idx;

colorMapAll = generate_cluster_colors(numGroups);

disimlarityScore = sparse(numGroups, numGroups);
ampSimilarity = sparse(numGroups, numGroups);
PCA = cell(numChannels, 1);
preprocessed.firingRate = zeros(numGroups, 1);
preprocessed.logFiringRate = zeros(numGroups, 1);
preprocessed.isiViolation = zeros(numGroups, 1);
preprocessSortedData();
multiplePlotOrder = [];
multiPlotPanels = gobjects(7,1);
lastMultiPlotOrder = multiplePlotOrder;

parentFig.WindowKeyPressFcn = @keyPressHandler;

% Cache for frequently accessed data
persistent cachedFilteredData cachedFilteredRange
cachedFilteredData = [];
cachedFilteredRange = [];

% Find all components that support KeyPressFcn and assign the same callback
% allComponents = findall(parentFig, '-property', 'KeyPressFcn');
% for iComp = 1:numel(allComponents)
%     allComponents(iComp).KeyPressFcn = @keyPressHandler;
% end

% left panels
leftPanel = uipanel(parentPanel, ...
    'Title','Spike Group Controller', ...
    'BackgroundColor',figColor, ...
    'Scrollable','on');
leftPanel.Layout.Column = [1 5];
leftPanel.Layout.Row    = [2 3];
applyColorScheme(leftPanel, figColor);

leftGrid = uigridlayout(leftPanel, [4 1], ...
    'RowHeight',{'fit','fit','fit','1x'}, ...
    'ColumnWidth',{'1x'}, ...
    'Padding',10, ...
    'RowSpacing',10);
applyColorScheme(leftGrid, figColor);

% top left panel and its buttons
topPanel = uipanel(leftGrid, 'BorderType','none');
topPanel.Layout.Row = 1;
topPanel.Layout.Column = 1;
applyColorScheme(topPanel, figColor);

topButtonLayout = uigridlayout(topPanel,[2 9], ...
    'ColumnWidth',{'fit','fit','fit','fit','fit','fit','1.75x','1.5x','2x'}, ...
    'ColumnSpacing',7.5, ...
    'Padding',[0 0 0 0]);
applyColorScheme(topButtonLayout, figColor);

btnRemove = uibutton(topButtonLayout,'Text','Remove',...
    'ButtonPushedFcn',@(btn,ev)onRemove(), 'Tooltip', sprintf('Delete units: \n Ctrl+Shift+D or Cmd+Shift+D'));
btnRemove.Layout.Row = [1 2];
btnRemove.Layout.Column = 1;
applyColorScheme(btnRemove, figColor);

btnMerge = uibutton(topButtonLayout,'Text','Merge',...
    'ButtonPushedFcn',@(btn,ev)onMerge(), 'Tooltip', sprintf('Merge units: \n Ctrl+Shift+M or Cmd+Shift+M'));
btnMerge.Layout.Row = [1 2];
btnMerge.Layout.Column = 2;
applyColorScheme(btnMerge, figColor);

btnMUA = uibutton(topButtonLayout,'Text','Isolation',...
    'ButtonPushedFcn',@(btn,ev)onMUA(), 'Tooltip', sprintf('Isolation type: \n Ctrl+Shift+I or Cmd+Shift+I'));
btnMUA.Layout.Row = [1 2];
btnMUA.Layout.Column = 3;
applyColorScheme(btnMUA, figColor);

btnUndo = uibutton(topButtonLayout,'Text','Undo',...
    'ButtonPushedFcn',@(btn,ev)onUndo(), 'Tooltip', sprintf('Undo selected units: \n Ctrl+Shift+U or Cmd+Shift+U'));
btnUndo.Layout.Row = [1 2];
btnUndo.Layout.Column = 4;
applyColorScheme(btnUndo, figColor);

btnReset = uibutton(topButtonLayout,'Text','Reset',...
    'ButtonPushedFcn',@(btn,ev)onReset(), 'Tooltip', sprintf('Reset all units: \n Ctrl+Shift+R or Cmd+Shift+R'));
btnReset.Layout.Row = [1 2];
btnReset.Layout.Column = 5;
applyColorScheme(btnReset, figColor);

btnAutoCut = uibutton(topButtonLayout,'Text','Limit',...
    'ButtonPushedFcn',@(btn,ev)onAutoCut(), 'Tooltip', sprintf('Auto limit and crop the dropped rate: \n Ctrl+Shift+L or Cmd+Shift+L'));
btnAutoCut.Layout.Row = [1 2];
btnAutoCut.Layout.Column = 6;
applyColorScheme(btnAutoCut, figColor);

lblLimit = uilabel(topButtonLayout,'Text',sprintf('Drop Rate:'),...
    'FontWeight','bold','HorizontalAlignment','Right');
lblLimit.Layout.Row = 1;
lblLimit.Layout.Column = 7;
applyColorScheme(lblLimit, figColor);

uiLimit = uieditfield(topButtonLayout,'numeric','Value',5, ...
    'Tooltip',sprintf('Min drop of firing rate \n for auto limit.'),...
    'ValueChangedFcn',@(ui,ev)setDropVal(ui.Value));
uiLimit.Layout.Row = 2;
uiLimit.Layout.Column = 7;
applyColorScheme(uiLimit, figColor);
dropRate = 5;

lblSliderReAlign = uilabel(topButtonLayout,'Text',sprintf('Realign \n (ms):'),...
    'FontWeight','bold','HorizontalAlignment','Right');
lblSliderReAlign.Layout.Row = [1 2];
lblSliderReAlign.Layout.Column = 8;
applyColorScheme(lblSliderReAlign, figColor);

sliderReAlign = uislider(topButtonLayout,'Value', 0,...
    'Tooltip',sprintf('Realign the spike waveform \n for any selected unit.'),...
    'Limits',[-halfSpikeWaveDur halfSpikeWaveDur]);
sliderReAlign.MajorTicks = [-halfSpikeWaveDur 0 halfSpikeWaveDur];
sliderReAlign.Layout.Row = [1 2];
sliderReAlign.Layout.Column = 9;
applyColorScheme(sliderReAlign, figColor);
sliderReAlign.ValueChangedFcn = @(sld,ev)updateRealignSpikes(sld.Value);

btnSave = uibutton(parentPanel,'Text','Save',...
    'ButtonPushedFcn',@(btn,ev)onSave(), 'Tooltip', 'Ctrl+Shift+S or Cmd+Shift+S');
btnSave.Layout.Column = 4;
btnSave.Layout.Row = 1;
applyColorScheme(btnSave, figColor);

lblSave = uilabel(parentPanel,'Text',' ');
lblSave.Layout.Column = 5;
lblSave.Layout.Row = 1;
applyColorScheme(lblSave, figColor);

processingPanel = uipanel(leftGrid, 'title','Similarity Processing');
processingPanel.Layout.Row = 2;
processingPanel.Layout.Column = 1;
applyColorScheme(processingPanel, figColor);

% similarity Assessing panel
processingLayout = uigridlayout(processingPanel,[2 8], ...
    'ColumnWidth',{'fit','1x','1x','1x','1x','.5x','.75x','.75x'}, ...
    'ColumnSpacing',10, ...
    'Padding',[5 5 5 5]);
applyColorScheme(processingLayout, figColor);

btnPreProcess = uibutton(processingLayout,'Text','Process',...
    'Tooltip',sprintf('After you set the distance and the metric \n run similarity detection.'),...
    'ButtonPushedFcn',@(btn,ev)onSimilarityEstimation());
btnPreProcess.Layout.Row = 2;
btnPreProcess.Layout.Column = 1;
applyColorScheme(btnPreProcess, figColor);

lblProgressPros = uilabel(processingLayout,'Text','Not Processed.',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblProgressPros.Layout.Row=1;
lblProgressPros.Layout.Column = 1;
applyColorScheme(lblProgressPros , figColor);

lblEstType = uilabel(processingLayout,'Text','Metric:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblEstType.Layout.Row=1;
lblEstType.Layout.Column = 2;
applyColorScheme(lblEstType , figColor);

distanceDD = uidropdown(processingLayout,...
    'Items',{'XCorr','KL-div','Bhattacharyya'},...
    'Tooltip',sprintf('Select the metric for \n similarity detection.'),...
    'Value','XCorr',...
    'ValueChangedFcn',@(dd,ev)updateDistance(dd.Value));
distanceDD.Layout.Row = 2;
distanceDD.Layout.Column = 2;
applyColorScheme(distanceDD, figColor);
distanceEstType = 'XCorr';
distanceEstMeasureType = 'mean';

lblDist = uilabel(processingLayout,'Text','Similarity:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblDist.Layout.Row=1;
lblDist.Layout.Column = 4;
applyColorScheme(lblDist, figColor);

uiDist = uieditfield(processingLayout,'numeric','Value',0,"Limits",[0 1], ...
    'Tooltip',sprintf('Set min similarity \n between units.'),...
    "LowerLimitInclusive","on", ...
    "UpperLimitInclusive","on", ...
    'ValueChangedFcn',@(ui,ev)setDistVal(ui.Value,0));
uiDist.Layout.Row=2;
uiDist.Layout.Column=4;
applyColorScheme(uiDist, figColor);

lblAmpDist = uilabel(processingLayout,'Text','Amp. Var.:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAmpDist.Layout.Row=1;
lblAmpDist.Layout.Column = 5;
applyColorScheme(lblAmpDist, figColor);

uiAmpDist = uieditfield(processingLayout,'numeric','Value',0,"Limits",[0 1], ...
    'Tooltip',sprintf('Set mean Amp. Var. ratio.'),...
    "LowerLimitInclusive","on", ...
    "UpperLimitInclusive","on", ...
    'ValueChangedFcn',@(ui,ev)setDistVal(ui.Value,1));
uiAmpDist.Layout.Row=2;
uiAmpDist.Layout.Column=5;
applyColorScheme(uiAmpDist, figColor);

lblYDist = uilabel(processingLayout,'Text','Dist (µm):',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblYDist.Layout.Row=1;
lblYDist.Layout.Column = 3;
applyColorScheme(lblYDist, figColor);

uiYDist = uieditfield(processingLayout,'numeric','Value',0,"Limits",[0 range(ylocs)], ...
    'Tooltip',sprintf('Set max radius \n to search for similarity'),...
    "LowerLimitInclusive","on", ...
    "UpperLimitInclusive","on", ...
    'ValueChangedFcn',@(ui,ev)setDistYVal(ui.Value));
uiYDist.Layout.Row=2;
uiYDist.Layout.Column=3;
applyColorScheme(uiYDist, figColor);

ampThr = 0;
distThr = 0;
rowGroup = [];
colGroup = [];
lastGroup = 0;
yDistThr = 0;


if exist('Next.png','file')
    btnNextPair = uibutton(processingLayout, 'Icon','Next.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(1), 'Tooltip', sprintf('Next pair: \n Ctrl+Shift+Down or Cmd+Shift+Down'));
else
    btnNextPair = uibutton(processingLayout, 'Text','Next',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(1), 'Tooltip', sprintf('Next pair: \n Ctrl+Shift+Down or Cmd+Shift+Down'));
end
btnNextPair.Layout.Row = 2;
btnNextPair.Layout.Column = 6;
applyColorScheme(btnNextPair, figColor);

if exist('Prev.png','file')
    btnPrevPair = uibutton(processingLayout, 'Icon','Prev.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(-1), 'Tooltip', sprintf('Prev. pair: \n Ctrl+Shift+Up or Cmd+Shift+Up'));
else
    btnPrevPair = uibutton(processingLayout, 'Text','Previous',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(-1), 'Tooltip', sprintf('Prev. pair: \n Ctrl+Shift+Up or Cmd+Shift+Up'));
end
btnPrevPair.Layout.Row = 1;
btnPrevPair.Layout.Column = 6;
applyColorScheme(btnPrevPair, figColor);

pairNumber = uilabel(processingLayout,'Text','Pair #:',...
    'FontWeight','bold','HorizontalAlignment','left');
pairNumber.Layout.Row=1;
pairNumber.Layout.Column = 7;
applyColorScheme(pairNumber, figColor);

uifPair = uieditfield(processingLayout,'Value','0','Editable','on',...
    'ValueChangedFcn',@(ui,ev)changeuifPar(ui.Value,0));
uifPair.Layout.Row=2;
uifPair.Layout.Column=7;
applyColorScheme(uifPair, figColor);

uifPair2 = uieditfield(processingLayout,'Value','/0','Editable','off',...
    'ValueChangedFcn',@(ui,ev)changeuifParAll(ui.Value,0));
uifPair2.Layout.Row=2;
uifPair2.Layout.Column=8;
applyColorScheme(uifPair2, figColor);

% criteria selection panel
middlePanel = uipanel(leftGrid, 'title','Group Selection');
middlePanel.Layout.Row = 3;
middlePanel.Layout.Column = 1;
applyColorScheme(middlePanel, figColor);

middleButtonLayout = uigridlayout(middlePanel,[3 9], ...   
    'ColumnWidth',{'.75x','.75x','.75x','1x','.75x','1.25x','1x','.5x','.75x','.75x'}, ...
    'RowHeight',{'fit','fit','fit'},...
    'ColumnSpacing',10, ...
    'Padding',[5 5 5 5]);
applyColorScheme(middleButtonLayout, figColor);

lblRateInc = uilabel(middleButtonLayout,'Text','Rate:',...
    'FontWeight','bold','HorizontalAlignment','left');
lblRateInc.Layout.Row=1;
lblRateInc.Layout.Column = 1;
applyColorScheme(lblRateInc, figColor);

uifRate = uieditfield(middleButtonLayout,'numeric','Value',0,...
    'Tooltip','Set min firing rate',...
    'ValueChangedFcn',@(ui,ev)setRateIncVal(ui.Value));
uifRate.Layout.Row=2;
uifRate.Layout.Column=1;
applyColorScheme(uifRate, figColor);
rateInc = true(numGroups, 1);

lblISIInc = uilabel(middleButtonLayout,'Text','ISI:',...
    'FontWeight','bold','HorizontalAlignment','left');
lblISIInc.Layout.Row=1;
lblISIInc.Layout.Column = 2;
applyColorScheme(lblISIInc, figColor);

uifISI = uieditfield(middleButtonLayout,'numeric','Value',ceil(max(preprocessed.isiViolation)),...
    'Tooltip','Set max %ISI violation',...
    'ValueChangedFcn',@(ui,ev)setISIIncVal(ui.Value));
uifISI.Layout.Row=2;
uifISI.Layout.Column=2;
applyColorScheme(uifISI, figColor);
isiInc = true(numGroups, 1);

lblSNRInc = uilabel(middleButtonLayout,'Text','SNR:',...
    'FontWeight','bold','HorizontalAlignment','left');
lblSNRInc.Layout.Row=1;
lblSNRInc.Layout.Column = 3;
applyColorScheme(lblSNRInc, figColor);

uifSNR = uieditfield(middleButtonLayout,'numeric','Value',1,...
    'Tooltip','Set min SNR (Amp. rel. to Thr.)',...
    'ValueChangedFcn',@(ui,ev)setSNRIncVal(ui.Value));
uifSNR.Layout.Row=2;
uifSNR.Layout.Column=3;
applyColorScheme(uifSNR, figColor);
snrInc = true(numGroups, 1);

lblPlarityInc = uilabel(middleButtonLayout,'Text','Polarity:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblPlarityInc.Layout.Row=1;
lblPlarityInc.Layout.Column = 4;
applyColorScheme(lblPlarityInc, figColor);

polarityDD = uidropdown(middleButtonLayout,...
    'Items',{'All','Main Neg.', 'Main Pos.', 'All Pos.', '1 Chan Neg.'},...
    'Value','All',...
    'Tooltip','Select spike polarity',...
    'ValueChangedFcn',@(dd,ev)updatePolarity(dd.Value));
polarityDD.Layout.Row = 2;
polarityDD.Layout.Column = 4;
applyColorScheme(polarityDD, figColor);
polarityInc = true(numGroups, 1);

lblTglVis = uilabel(middleButtonLayout,'Text','Vis:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblTglVis.Layout.Row=1;
lblTglVis.Layout.Column = 5;
applyColorScheme(lblTglVis, figColor);

chkTogVis = uicheckbox(middleButtonLayout,'Text','',...
    'FontWeight','bold',...
    'Tooltip','Toggle Vis of selected units',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onToggleVis(cb.Value));
chkTogVis.Layout.Row = 2;
chkTogVis.Layout.Column = 5;
applyColorScheme(chkTogVis, figColor);

chkInclusion = uicheckbox(middleButtonLayout,'Text','Include',...
    'Tooltip',sprintf('Include units that\n pass all criteria.'),...
    'FontWeight','bold',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onINCEXC(cb.Value,1));
chkInclusion.Layout.Row = 1;
chkInclusion.Layout.Column = 6;
applyColorScheme(chkInclusion, figColor);
incChan =[];
exChan = [];

chkExclusion = uicheckbox(middleButtonLayout,'Text','Exclude',...
    'Tooltip',sprintf('Exclude units that\n do not pass all criteria.'),...
    'FontWeight','bold',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onINCEXC(cb.Value,0));
chkExclusion.Layout.Row = 2;
chkExclusion.Layout.Column = 6;
applyColorScheme(chkExclusion, figColor);


btnDeselectUnit = uibutton(middleButtonLayout, 'Text','Deselect',...
    'Tooltip',sprintf('Deselect all selected units.'),...
    'ButtonPushedFcn',@(btn,ev)onDeselect, 'Tooltip', sprintf('Deselect all selected units: \n Ctrl+Shift+H or Cmd+Shift+H'));
btnDeselectUnit.Layout.Row=2;
btnDeselectUnit.Layout.Column = 7;
applyColorScheme(btnDeselectUnit, figColor);

if exist('Next.png','file')
    btnNextUnit = uibutton(middleButtonLayout, 'Icon','Next.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(1), 'Tooltip', sprintf('Next unit: \n Ctrl+Shift+N or Cmd+Shift+N'));
else
    btnNextUnit = uibutton(middleButtonLayout, 'Text','Next',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(1), 'Tooltip', sprintf('Next unit: \n Ctrl+Shift+N or Cmd+Shift+N'));
end
btnNextUnit.Layout.Row = 2;
btnNextUnit.Layout.Column = 8;
applyColorScheme(btnNextUnit, figColor);

if exist('Prev.png','file')
    btnPrevUnit = uibutton(middleButtonLayout, 'Icon','Prev.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(-1), 'Tooltip', sprintf('Prev. unit: \n Ctrl+Shift+P or Cmd+Shift+P'));
else
    btnPrevUnit = uibutton(middleButtonLayout, 'Text','Previous',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(-1), 'Tooltip', sprintf('Prev. unit: \n Ctrl+Shift+P or Cmd+Shift+P'));
end
btnPrevUnit.Layout.Row = 1;
btnPrevUnit.Layout.Column = 8;
applyColorScheme(btnPrevUnit, figColor);

unitNumber = uilabel(middleButtonLayout,'Text','Unit #:',...
    'FontWeight','bold','HorizontalAlignment','left');
unitNumber.Layout.Row=1;
unitNumber.Layout.Column = 9;
applyColorScheme(unitNumber, figColor);

uifUnit = uieditfield(middleButtonLayout,'Value','0','Editable','on',...
    'ValueChangedFcn',@(ui,ev)changeuifUnit(ui.Value,0));
uifUnit.Layout.Row=2;
uifUnit.Layout.Column=9;
applyColorScheme(uifUnit, figColor);

uifUnit2 = uieditfield(middleButtonLayout,'Value',['/' num2str(numGroups)],'Editable','off',...
    'ValueChangedFcn',@(ui,ev)changeuifUnitAll(ui.Value,0));
uifUnit2.Layout.Row=2;
uifUnit2.Layout.Column=10;
applyColorScheme(uifUnit2, figColor);


% Add pagination controls
if totalPages > 1
    lblPage = uilabel(middleButtonLayout,'Text',sprintf('Page: %d/%d', currentPage, totalPages),...
        'FontWeight','bold','HorizontalAlignment','center');
    lblPage.Layout.Row = 3;
    lblPage.Layout.Column = 6;
    applyColorScheme(lblPage, figColor);

    lblPageSize = uilabel(middleButtonLayout,'Text','Size:',...
        'FontWeight','bold','HorizontalAlignment','center');
    lblPageSize.Layout.Row = 3;
    lblPageSize.Layout.Column = 4;
    applyColorScheme(lblPageSize, figColor);

    uifPageSize = uieditfield(middleButtonLayout,'Value','50',...
        'ValueChangedFcn',@(ui,ev)changeuifPageSize(ui.Value));
    uifPageSize.Layout.Row=3;
    uifPageSize.Layout.Column=5;
    applyColorScheme(uifPageSize, figColor);

    btnPrevPage = uibutton(middleButtonLayout,'Text','◀ Prev Page',...
        'ButtonPushedFcn',@(btn,ev)onChangePage(-1),'Tooltip', sprintf('Prev. unit: \n Ctrl+Shift+Left or Cmd+Shift+Left'));
    btnPrevPage.Layout.Row = 3;
    btnPrevPage.Layout.Column = [1 3];
    applyColorScheme(btnPrevPage, figColor);

    btnNextPage = uibutton(middleButtonLayout,'Text','Next Page ▶',...
        'ButtonPushedFcn',@(btn,ev)onChangePage(1),'Tooltip', sprintf('Prev. unit: \n Ctrl+Shift+Right or Cmd+Shift+Right'));
    btnNextPage.Layout.Row = 3;
    btnNextPage.Layout.Column = [7 9];
    applyColorScheme(btnNextPage, figColor);
end

lastUnit = 0;

%% Groups check boxes and radio panels
checkPanel = uipanel(leftGrid, ...
    'Title','Groups:', ...
    'BackgroundColor',figColor,...
    'Scrollable','on');
checkPanel.Layout.Row = 4;
checkPanel.Layout.Column = 1;
applyColorScheme(checkPanel, figColor);

% Calculate page index ranges for pagination
startIdx = (currentPage-1) * PAGE_SIZE + 1;
endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
visibleGroups = startIdx:endIdx;
visibleRows = length(visibleGroups) + 1;  % +1 for the "All" row

mainCheckGrid = uigridlayout(checkPanel, [visibleRows, 9], ...
    'RowHeight', repmat({'fit'},1,visibleRows), ...
    'ColumnWidth', repmat({'fit'},1,9), ...
    'RowSpacing',2, ...
    'Scrollable','on',...
    'ColumnSpacing',5, ...
    'Padding',[5 5 5 5]);
applyColorScheme(mainCheckGrid, figColor);

% Arrays for UI handles - only initialize for visible groups to save memory
lblGroupHandles       = gobjects(numGroups,1);
lblChannelHandles     = gobjects(numGroups,1);
groupVisCheckboxes    = gobjects(numGroups,1);
groupSelectCheckboxes = gobjects(numGroups,1);
groupRadioButtons     = gobjects(numGroups,1);
lblIsolation          = gobjects(numGroups,1);
lblRate               = gobjects(numGroups,1);
lblDetectablity       = gobjects(numGroups,1);
lblISIViolation       = gobjects(numGroups,1);

% Object for multiple plot
plotItemsNames  = {'Amplitude','Waveform', 'CCG', 'ISI', 'Density', 'PCs','Trace', 'Features'};
multipleSpinnerObj = gobjects(length(plotItemsNames),1);
multiplePlotObj = gobjects(length(plotItemsNames),1);
plotScaleFactor = ones(1, length(plotItemsNames));
reScaledFlag = false;

radioGroupBG = uibuttongroup(mainCheckGrid, ...
    'Title','', ...
    'Scrollable','on');
radioGroupBG.Layout.Row = [2 visibleRows];
radioGroupBG.Layout.Column = 9;
applyColorScheme(radioGroupBG, figColor);
rowPixelHeight = 25;
rowOffsetTop   = 5;

lblMerge = uilabel(mainCheckGrid,'Text','Merge To:',...
    'FontWeight','bold');
lblMerge.Layout.Row = 1;
lblMerge.Layout.Column = 9;
applyColorScheme(lblMerge, figColor);

% all selection row
lblAll = uilabel(mainCheckGrid,'Text','G #:',...
    'FontWeight','bold');
lblAll.Layout.Row = 1;
lblAll.Layout.Column = 1;
applyColorScheme(lblAll, figColor);

lblChannel = uilabel(mainCheckGrid,'Text','Ch#:','FontWeight','bold');
lblChannel.Layout.Row = 1;
lblChannel.Layout.Column = 2;
applyColorScheme(lblChannel, figColor);

lblunitIsolation = uilabel(mainCheckGrid,'Text','Iso:','FontWeight','bold');
lblunitIsolation.Layout.Row = 1;
lblunitIsolation.Layout.Column = 3;
applyColorScheme(lblunitIsolation, figColor);

lblunitRate = uilabel(mainCheckGrid,'Text','Rate:','FontWeight','bold');
lblunitRate.Layout.Row = 1;
lblunitRate.Layout.Column = 4;
applyColorScheme(lblunitRate, figColor);

lblSNR = uilabel(mainCheckGrid,'Text','SNR:','FontWeight','bold');
lblSNR.Layout.Row = 1;
lblSNR.Layout.Column = 5;
applyColorScheme(lblSNR, figColor);

lblISIV = uilabel(mainCheckGrid,'Text','ISIV %:');
lblISIV.Layout.Row = 1;
lblISIV.Layout.Column = 6;
applyColorScheme(lblISIV, figColor);

chkAllVis = uicheckbox(mainCheckGrid,'Text','Vis. All',...
    'FontWeight','bold',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onAllVisChecked(cb.Value));
chkAllVis.Layout.Row = 1;
chkAllVis.Layout.Column = 7;
applyColorScheme(chkAllVis, figColor);

chkAllSel = uicheckbox(mainCheckGrid,'Text','Sel. All',...
    'FontWeight','bold',...
    'ValueChangedFcn',@(cb,ev)onAllSelectChecked(cb.Value));
chkAllSel.Layout.Row = 1;
chkAllSel.Layout.Column = 8;
applyColorScheme(chkAllSel, figColor);

yRow1 = (visibleRows-1)*rowPixelHeight + rowOffsetTop;  % visibleRows-1 from the bottom
rdoNone = uiradiobutton(radioGroupBG, 'Text','none',...
    'FontWeight','bold','fontColor',1-figColor,...
    'Position',[10 yRow1 80 22]);
rdoNone.Value = true; % default selection

%% Create UI elements only for visible groups
createGroupUiElements();

%% Right plot and control panel
topRightPanel = uipanel(parentPanel,...
    'Title','Visualization Control',...
    'BackgroundColor',figColor);
topRightPanel.Layout.Column = 6;
topRightPanel.Layout.Row    = [1 2];
applyColorScheme(topRightPanel, figColor);

bottomRightPanel = uipanel(parentPanel,...
    'Title','Plots',...
    'BackgroundColor',figColor);
bottomRightPanel.Layout.Column = 6;
bottomRightPanel.Layout.Row    = 3;
applyColorScheme(bottomRightPanel, figColor);

topRightGrid = uigridlayout(topRightPanel,...
    'RowHeight',{'fit'},...
    'ColumnWidth',{'1x'},...
    'Padding',10);
applyColorScheme(topRightGrid, figColor);

bottomRightGrid = uigridlayout(bottomRightPanel,[5 3],...
    'RowHeight',{'1x','1x','12x','1x' ,'1x'},...
    'ColumnWidth',{'1x','12x','1x'},...
    'Padding',10);
applyColorScheme(bottomRightGrid, figColor);

% Row1: Buttons for Plot raw data, Plot filtered data, and dropdown
plotButtonPanel = uipanel(topRightGrid,'BorderType','none','BackgroundColor',figColor);
plotButtonPanel.Layout.Row = 1;
plotButtonPanel.Layout.Column = 1;
applyColorScheme(plotButtonPanel, figColor);

plotButtonGrid = uigridlayout(plotButtonPanel,[3 12],...
    'ColumnWidth',{'1x','1x','1x','1x','1x','1x','1x','1x','1x','1x','1x', '1x','1x','fit'},...
    'RowHeight',{'fit','fit','fit'},...
    'Padding',[0 0 0 0],...
    'RowSpacing', 1,...
    'ColumnSpacing',8);
applyColorScheme(plotButtonGrid, figColor);

plotTypeDD = uidropdown(plotButtonGrid,...
    'Items',{'Amplitude','Waveform', 'CCG', 'ISI', 'Density', 'PCs','Trace', 'Features','Multiple'},...
    'Tooltip',sprintf('Select plot type, \n select multiple for multiple types'),...
    'Value','Trace',...
    'ValueChangedFcn',@(dd,ev)updatePlotType(dd.Value));
plotTypeDD.Layout.Row = 2;
plotTypeDD.Layout.Column = 1;
applyColorScheme(plotTypeDD, figColor);

plotType = 'Trace';
lastPlotted = plotType;

plotFilterTypeDD = uidropdown(plotButtonGrid,...
    'Items',{'filtered','raw'},...
    'Tooltip',sprintf('Select Filter types'),...
    'Value','filtered',...
    'ValueChangedFcn',@(dd,ev)updatePlotFilterType(dd.Value));
plotFilterTypeDD.Layout.Row = 2;
plotFilterTypeDD.Layout.Column = 2;
applyColorScheme(plotFilterTypeDD, figColor);

plotFilterType = 'filtered';

chTypeDropdown = uidropdown(plotButtonGrid,...
    'Items',{'Config' ,'0', '1', '2', '3', '4', '5', '6', '7','8','9','10'},...
    'Tooltip',sprintf('Number of channels to show \n before/after the main channel.'),...
    'Value','Config');
chTypeDropdown.Layout.Row = 2;
chTypeDropdown.Layout.Column = 3;
applyColorScheme(chTypeDropdown, figColor);
chTypeDropdown.ValueChangedFcn      = @(sld,ev)updateChType(sld.Value);

numWaveformDropdown = uidropdown(plotButtonGrid,...
    'Items',{'1','5','10','25','50','100','1000'},...
    'Tooltip',sprintf('Number of waveforms to show.'),...
    'Value','25');
numWaveformDropdown.Layout.Row = 2;
numWaveformDropdown.Layout.Column = 4;
applyColorScheme(numWaveformDropdown, figColor);
numWaveformDropdown.ValueChangedFcn      = @(sld,ev)updateNumWaveforms(sld.Value);

numWaveforms = 25;

spikeWaveDurDD = uidropdown(plotButtonGrid,...
    'Items',{'Config','0.25','0.5','1','1.5','2'},...
    'Tooltip',sprintf('Select the waveform duration.'),...
    'Value','Config',...
    'ValueChangedFcn',@(dd,ev)updateWaveDur(dd.Value));
spikeWaveDurDD.Layout.Row = 2;
spikeWaveDurDD.Layout.Column = 5;
applyColorScheme(spikeWaveDurDD, figColor);

meanWaveformDD = uidropdown(plotButtonGrid,...
    'Items',{'Individual','Mean'},...
    'Tooltip',sprintf('Select to plot the mean or \n individual waveforms.'),...
    'Value','Individual',...
    'ValueChangedFcn',@(dd,ev)updateWavetype(dd.Value));
meanWaveformDD.Layout.Row = 2;
meanWaveformDD.Layout.Column = 6;
applyColorScheme(meanWaveformDD, figColor);
plotMean = 0;

ccgLagDD = uidropdown(plotButtonGrid,...
    'Items',{'25','50','75','100','200'},...
    'Tooltip',sprintf('Select CCG max time lag (ms).'),...
    'Value','100',...
    'ValueChangedFcn',@(dd,ev)updateCCGLag(dd.Value));
ccgLagDD.Layout.Row = 2;
ccgLagDD.Layout.Column = 7;
applyColorScheme(ccgLagDD, figColor);

ccgSmoothDD = uidropdown(plotButtonGrid,...
    'Items',{'0','1','2', '5', '10'},...
    'Tooltip',sprintf('Smoothing factor for \n CCG and density plots.'),...
    'Value','0',...
    'ValueChangedFcn',@(dd,ev)updateCCGSmoothN(dd.Value));
ccgSmoothDD.Layout.Row = 2;
ccgSmoothDD.Layout.Column = 8;
applyColorScheme(ccgSmoothDD, figColor);

isiThrDD = uidropdown(plotButtonGrid,...
    'Items',{'1','2','3','4','5'},...
    'Tooltip',sprintf('Select the thereshold to \n measure ISI violation.'),...
    'Value','1',...
    'ValueChangedFcn',@(dd,ev)updateISIThr(dd.Value));
isiThrDD.Layout.Row = 2;
isiThrDD.Layout.Column = 9;
applyColorScheme(isiThrDD, figColor);

ampSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 1,'MajorTicksMode','auto', 'Limits',[0 10]);
ampSlider.MajorTicks = [0 10];
ampSlider.Layout.Row = 2;
ampSlider.Layout.Column = 10;
applyColorScheme(ampSlider, figColor);
ampSlider.ValueChangedFcn      = @(sld,ev)updateScale(sld.Value);

lineWidthSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 1,'MajorTicksMode','auto', 'Limits',[0 5]);
lineWidthSlider.MajorTicks = [0 5];
lineWidthSlider.Layout.Row = 2;
lineWidthSlider.Layout.Column = 11;
applyColorScheme(lineWidthSlider, figColor);
lineWidthSlider.ValueChangedFcn      = @(sld,ev)updateLineWidth(sld.Value);
lineWidth = 1;

alphaSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 1,'MajorTicksMode','auto', 'Limits',[0 1]);
alphaSlider.MajorTicks = [0 1];
alphaSlider.Layout.Row = 2;
alphaSlider.Layout.Column = 12;
applyColorScheme(alphaSlider, figColor);
alphaSlider.ValueChangedFcn      = @(sld,ev)updateAlphaLevel(sld.Value);
alphaLevel = 1;

btnExportFig = uibutton(plotButtonGrid,'Text','Export Fig.',...
    'ButtonPushedFcn',@(btn,ev)onExportFig());
btnExportFig.Layout.Row = 2;
btnExportFig.Layout.Column = 13;
applyColorScheme(btnExportFig, figColor);

% labels for buttons
lblplotType = uilabel(plotButtonGrid,'Text','Plot Type:',...
    'HorizontalAlignment','left');
lblplotType.Layout.Row = 1;
lblplotType.Layout.Column = 1;
applyColorScheme(lblplotType, figColor);

lblFilterType = uilabel(plotButtonGrid,'Text','Data Type:',...
    'HorizontalAlignment','left');
lblFilterType.Layout.Row = 1;
lblFilterType.Layout.Column = 2;
applyColorScheme(lblFilterType, figColor);

lblchType = uilabel(plotButtonGrid,'Text','# Channel:',...
    'HorizontalAlignment','left');
lblchType.Layout.Row = 1;
lblchType.Layout.Column = 3;
applyColorScheme(lblchType, figColor);

lblNumSpikes = uilabel(plotButtonGrid,'Text','# Spikes:',...
    'HorizontalAlignment','left');
lblNumSpikes.Layout.Row = 1;
lblNumSpikes.Layout.Column = 4;
applyColorScheme(lblNumSpikes, figColor);

lblSpikeDur = uilabel(plotButtonGrid,'Text','Dur. (ms):',...
    'HorizontalAlignment','left');
lblSpikeDur.Layout.Row = 1;
lblSpikeDur.Layout.Column = 5;
applyColorScheme(lblSpikeDur, figColor);

lblMeanWaveform = uilabel(plotButtonGrid,'Text','Waveform:',...
    'HorizontalAlignment','left');
lblMeanWaveform.Layout.Row = 1;
lblMeanWaveform.Layout.Column = 6;
applyColorScheme(lblMeanWaveform, figColor);

lblCCGLag = uilabel(plotButtonGrid,'Text','Bins:',...
    'HorizontalAlignment','left');
lblCCGLag.Layout.Row = 1;
lblCCGLag.Layout.Column = 7;
applyColorScheme(lblCCGLag, figColor);

lblSmoothness = uilabel(plotButtonGrid,'Text','Smoothness:',...
    'HorizontalAlignment','left');
lblSmoothness.Layout.Row = 1;
lblSmoothness.Layout.Column = 8;
applyColorScheme(lblSmoothness, figColor);

lblisiThr = uilabel(plotButtonGrid,'Text','ISI Thr(ms):',...
    'HorizontalAlignment','left');
lblisiThr.Layout.Row = 1;
lblisiThr.Layout.Column = 9;
applyColorScheme(lblisiThr, figColor);

lblXAmpScale = uilabel(plotButtonGrid,'Text','Amp. Scale:',...
    'HorizontalAlignment','left');
lblXAmpScale.Layout.Row = 1;
lblXAmpScale.Layout.Column = 10;
applyColorScheme(lblXAmpScale, figColor);

lblXLineWidth = uilabel(plotButtonGrid,'Text','Line Width:',...
    'HorizontalAlignment','left');
lblXLineWidth.Layout.Row = 1;
lblXLineWidth.Layout.Column = 11;
applyColorScheme(lblXLineWidth, figColor);

lblAlphaLvl = uilabel(plotButtonGrid,'Text','Alpha:',...
    'HorizontalAlignment','left');
lblAlphaLvl.Layout.Row = 1;
lblAlphaLvl.Layout.Column = 12;
applyColorScheme(lblAlphaLvl, figColor);

exportFormatDD = uidropdown(plotButtonGrid,...
    'Items',{'Image','Vector'},...
    'Value','Image',...
    'ValueChangedFcn',@(dd,ev)updateExportFormat(dd.Value));
exportFormatDD.Layout.Row = 1;
exportFormatDD.Layout.Column = 13;
applyColorScheme(exportFormatDD, figColor);

exportType = 'Image';

%% Lower right pannels, Scrollers
xSlider = uislider(bottomRightGrid,'Value',0,...
    'MajorTicksMode','auto');
xSlider.Layout.Row = [4 5];
xSlider.Layout.Column = 2;
applyColorScheme(xSlider, figColor);

lblXAmpScale = uilabel(bottomRightGrid,'Text','Time window (s):',...
    'HorizontalAlignment','right');
lblXAmpScale.Layout.Row = 4;
lblXAmpScale.Layout.Column = 3;
applyColorScheme(lblXAmpScale, figColor);

xStepDropdown = uidropdown(bottomRightGrid,...
    'Items',{'0.001','0.01','0.1','1','2','5','10','100','1000','full'},...
    'Tooltip',sprintf('Window size for plots.'),...
    'Value','0.1');
xStepDropdown.Layout.Row = 5;
xStepDropdown.Layout.Column = 3;
applyColorScheme(xStepDropdown, figColor);

xSlider.ValueChangedFcn     = @(sld,ev)updateXWindow(sld.Value, xStepDropdown.Value);
xStepDropdown.ValueChangedFcn = @(dd,ev)updateXWindow(xSlider.Value, dd.Value);

% Y slider
channelValues = 2.^(0:10);
channelValues(channelValues > numChannels) = [];
tickIntIdx = nearest(channelValues, numChannels/10);
ySlider = uislider(bottomRightGrid,'Orientation','vertical','Value', 1,'MajorTicksMode','auto', 'Limits',[1 numChannels]);
ySlider.MajorTicks = [0 : channelValues(tickIntIdx): numChannels];

ySlider.Layout.Row = [3 5];
ySlider.Layout.Column = 1;
applyColorScheme(ySlider, figColor);

channelValues = arrayfun(@num2str, channelValues, 'UniformOutput', false);
yStepDropdown = uidropdown(bottomRightGrid,...
    'Items',[{'full'}, channelValues],...
    'Tooltip',sprintf('Range for channel visualization.'),...
    'Value','full');
yStepDropdown.Layout.Row = 2;
yStepDropdown.Layout.Column = 1;
applyColorScheme(yStepDropdown, figColor);

lblXAmpScale = uilabel(bottomRightGrid,'Text','# Channels:',...
    'HorizontalAlignment','right');
lblXAmpScale.Layout.Row = 1;
lblXAmpScale.Layout.Column = 1;
applyColorScheme(lblXAmpScale, figColor);

ySlider.ValueChangedFcn      = @(sld,ev)updateYWindow(sld.Value, yStepDropdown.Value);
yStepDropdown.ValueChangedFcn= @(dd,ev)updateYWindow(ySlider.Value, dd.Value);

axChannels = uiaxes(bottomRightGrid);
axChannels.Layout.Row = [1 3];
axChannels.Layout.Column = [2 3];
axChannels.Color   = [0.2 0.2 0.2];
axChannels.XColor  = [1 1 1];
axChannels.YColor  = [1 1 1];
axChannels.GridColor = [0.7 0.7 0.7];
axChannels.TickDir   =  "out";
title(axChannels,'On-signal Spike Explorer','Color','white');
xlabel(axChannels,'Time (s)','Color','white');
ylabel(axChannels,'Channel #','Color','white');

% Initialize X/Y window
xWindowStart = 0;
xWindowEnd   = 0.1;
xStepVal     = 0.1;
yWindowStart = 1;
yWindowEnd   = numChannels;
yStepVal     = numChannels;
ampScale     = 1;




%% Default run
plotOption();
normalizeTimeAmp();
refreshLabels();

%% =============== CALLBACKS ===============
    function keyPressHandler(~, event)
        persistent lastTime lastKey lastModifiers
        if isempty(lastTime)
            lastTime = datetime('now');
            lastKey = '';
            lastModifiers = {};
        end
        currentTime = datetime('now');
        if strcmpi(lastKey, event.Key) && isequal(lastModifiers, event.Modifier) ...
                && seconds(currentTime - lastTime) < 0.3
            return
        end
        lastTime = currentTime;
        lastKey = event.Key;
        lastModifiers = event.Modifier;

        if (ismember('control', event.Modifier) || ismember('command', event.Modifier)) ...
                && ismember('shift', event.Modifier)
            switch lower(event.Key)
                case 'r', onReset()
                case 'uparrow', onSelectSimilarPairs(-1)
                case 'downarrow', onSelectSimilarPairs(1)
                case 'leftarrow', onChangePage(-1)
                case 'rightarrow', onChangePage(1)
                case 'p', onSelectUnit(-1)
                case 'n', onSelectUnit(1)
                case 's', onSave()
                case 'd', onRemove()
                case 'm', onMerge()
                case 'i', onMUA()
                case 'u', onUndo()
                case 'l', onAutoCut()
                case 'h', onDeselect()
            end
        end
    end

    function changeuifPageSize(val)
        PAGE_SIZE = str2double(val);
        currentPage = 1;
        totalPages = ceil(numGroups/PAGE_SIZE);
        onChangePage(0);
    end
% Function to change page in pagination
    function onChangePage(direction)
        newPage = currentPage + direction;
        if newPage >= 1 && newPage <= totalPages
            currentPage = newPage;
            if exist('lblPage', 'var')
                lblPage.Text = sprintf('Page: %d/%d', currentPage, totalPages);
            end

            % Clean up old UI elements
            delete(findall(mainCheckGrid, 'Type', 'UILabel'));
            delete(findall(mainCheckGrid, 'Type', 'UICheckbox'));
            delete(findall(radioGroupBG, 'Type', 'UIRadioButton'));

            % Rebuild grid with new page size
            delete(mainCheckGrid);

            % Calculate new visible groups
            startIdx = (currentPage-1) * PAGE_SIZE + 1;
            endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
            visibleGroups = startIdx:endIdx;
            visibleRows = length(visibleGroups) + 1;  % +1 for the "All" row

            % Create new grid
            mainCheckGrid = uigridlayout(checkPanel, [visibleRows, 9], ...
                'RowHeight', repmat({'fit'},1,visibleRows), ...
                'ColumnWidth', repmat({'fit'},1,9), ...
                'RowSpacing',2, ...
                'Scrollable','on',...
                'ColumnSpacing',5, ...
                'Padding',[5 5 5 5]);
            applyColorScheme(mainCheckGrid, figColor);

            % Recreate radio group
            radioGroupBG = uibuttongroup(mainCheckGrid, ...
                'Title','', ...
                'Scrollable','on');
            radioGroupBG.Layout.Row = [2 visibleRows];
            radioGroupBG.Layout.Column = 9;
            applyColorScheme(radioGroupBG, figColor);

            % Recreate header row
            recreateHeaderRow();

            % Create none radio button
            yRow1 = (visibleRows-1)*rowPixelHeight + rowOffsetTop;
            rdoNone = uiradiobutton(radioGroupBG, 'Text','none',...
                'FontWeight','bold','fontColor',1-figColor,...
                'Position',[10 yRow1 80 22]);
            rdoNone.Value = true;

            % Create UI elements for the current page's groups
            createGroupUiElements();

            % Refresh state of checkboxes based on current selection
            refreshVisibility();
            refreshLabels();
        end
    end

% Function to create header row in the group list
    function recreateHeaderRow()
        lblAll = uilabel(mainCheckGrid,'Text','G #:',...
            'FontWeight','bold');
        lblAll.Layout.Row = 1;
        lblAll.Layout.Column = 1;
        applyColorScheme(lblAll, figColor);

        lblChannel = uilabel(mainCheckGrid,'Text','Ch#:','FontWeight','bold');
        lblChannel.Layout.Row = 1;
        lblChannel.Layout.Column = 2;
        applyColorScheme(lblChannel, figColor);

        lblunitIsolation = uilabel(mainCheckGrid,'Text','Iso:','FontWeight','bold');
        lblunitIsolation.Layout.Row = 1;
        lblunitIsolation.Layout.Column = 3;
        applyColorScheme(lblunitIsolation, figColor);

        lblunitRate = uilabel(mainCheckGrid,'Text','Rate:','FontWeight','bold');
        lblunitRate.Layout.Row = 1;
        lblunitRate.Layout.Column = 4;
        applyColorScheme(lblunitRate, figColor);

        lblSNR = uilabel(mainCheckGrid,'Text','SNR:','FontWeight','bold');
        lblSNR.Layout.Row = 1;
        lblSNR.Layout.Column = 5;
        applyColorScheme(lblSNR, figColor);

        lblISIV = uilabel(mainCheckGrid,'Text','ISIV %:');
        lblISIV.Layout.Row = 1;
        lblISIV.Layout.Column = 6;
        applyColorScheme(lblISIV, figColor);

        chkAllVis = uicheckbox(mainCheckGrid,'Text','Vis. All',...
            'FontWeight','bold',...
            'Value',false,...
            'ValueChangedFcn',@(cb,ev)onAllVisChecked(cb.Value));
        chkAllVis.Layout.Row = 1;
        chkAllVis.Layout.Column = 7;
        applyColorScheme(chkAllVis, figColor);

        chkAllSel = uicheckbox(mainCheckGrid,'Text','Sel. All',...
            'FontWeight','bold',...
            'ValueChangedFcn',@(cb,ev)onAllSelectChecked(cb.Value));
        chkAllSel.Layout.Row = 1;
        chkAllSel.Layout.Column = 8;
        applyColorScheme(chkAllSel, figColor);

        lblMerge = uilabel(mainCheckGrid,'Text','Merge To:',...
            'FontWeight','bold');
        lblMerge.Layout.Row = 1;
        lblMerge.Layout.Column = 9;
        applyColorScheme(lblMerge, figColor);
    end

% Create UI elements for visible groups
    function createGroupUiElements()
        startIdx = (currentPage-1) * PAGE_SIZE + 1;
        endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
        visibleGroups = startIdx:endIdx;

        for i = 1:length(visibleGroups)
            iG = visibleGroups(i);
            rowIdx = i + 1;  % +1 because row 1 is the header

            % group list
            lblGroupHandles(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('G: %d', groupList(iG)),...
                'FontWeight','bold');
            lblGroupHandles(iG).Layout.Row = rowIdx;
            lblGroupHandles(iG).Layout.Column = 1;

            % channel list
            lblChannelHandles(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('Ch %d', channelList(iG)),...
                'FontWeight','bold');
            lblChannelHandles(iG).Layout.Row = rowIdx;
            lblChannelHandles(iG).Layout.Column = 2;

            lblIsolation(iG) = uilabel(mainCheckGrid, ...
                'Text', unitIsolation{iG},...
                'FontWeight','bold');
            lblIsolation(iG).Layout.Row = rowIdx;
            lblIsolation(iG).Layout.Column = 3;

            lblRate(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('%.1f',preprocessed.firingRate(iG)),...
                'FontWeight','bold','FontColor',1-figColor);
            lblRate(iG).Layout.Row = rowIdx;
            lblRate(iG).Layout.Column = 4;

            lblDetectablity(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('%.1f',1 + detectblity(iG)),...
                'FontWeight','bold','FontColor',1-figColor);
            lblDetectablity(iG).Layout.Row = rowIdx;
            lblDetectablity(iG).Layout.Column = 5;

            lblISIViolation(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('%.1f',preprocessed.isiViolation(iG)),...
                'FontWeight','bold','FontColor',1-figColor);
            lblISIViolation(iG).Layout.Row = rowIdx;
            lblISIViolation(iG).Layout.Column = 6;

            % visualization checkboxes
            groupVisCheckboxes(iG) = uicheckbox(mainCheckGrid,...
                'Value',chkBoxVis(iG), 'Text','Vis.',...
                'ValueChangedFcn',@(cb,ev)onVisCheckboxChanged(iG));
            groupVisCheckboxes(iG).Layout.Row = rowIdx;
            groupVisCheckboxes(iG).Layout.Column = 7;

            % selection checkboxes
            groupSelectCheckboxes(iG) = uicheckbox(mainCheckGrid,...
                'Text','Select','Value',chkBoxSelect(iG),...
                'ValueChangedFcn',@(cb,ev)onSelectCheckboxChanged(iG));
            groupSelectCheckboxes(iG).Layout.Row = rowIdx;
            groupSelectCheckboxes(iG).Layout.Column = 8;

            % radio buttons
            yPos = (length(visibleGroups) - i )*rowPixelHeight + rowOffsetTop;
            groupRadioButtons(iG) = uiradiobutton(radioGroupBG,...
                'Text', sprintf('G%d',iG),...
                'FontWeight','bold',...
                'Position',[10, yPos, 80, 22]);
        end
    end

% Update visibility status based on current selection
    function refreshVisibility()
        startIdx = (currentPage-1) * PAGE_SIZE + 1;
        endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
        visibleGroups = startIdx:endIdx;

        for i = 1:length(visibleGroups)
            iG = visibleGroups(i);
            if any(selected_order == iG)
                groupVisCheckboxes(iG).Value = true;
                c = getGroupColor(groupList(iG));
                groupVisCheckboxes(iG).FontColor = c;
                groupVisCheckboxes(iG).FontWeight = 'Bold';
            else
                groupVisCheckboxes(iG).Value = false;
                groupVisCheckboxes(iG).FontColor = 1-figColor;
                groupVisCheckboxes(iG).FontWeight = 'normal';
            end
        end
    end

    function onAllVisChecked(val)
        % Batch update visibility checkboxes
        if val
            chkBoxVis(:) = 1;
            selected_order = 1:numGroups;
        else
            chkBoxVis(:) = 0;
            selected_order = [];
        end

        startIdx = (currentPage-1) * PAGE_SIZE + 1;
        endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
        visibleGroups = startIdx:endIdx;

        % Batch update UI elements
        for idx = 1:length(visibleGroups)
            k = visibleGroups(idx);
            groupVisCheckboxes(k).Value = val;

            if val
                c = getGroupColor(groupList(k));
                groupVisCheckboxes(k).FontColor = c;
                groupVisCheckboxes(k).FontWeight ='Bold';
                lblChannelHandles(k).BackgroundColor = c;
                lblChannelHandles(k).FontColor = bestTextColorFor(c);
            else
                groupVisCheckboxes(k).FontColor = 1-figColor;
                groupVisCheckboxes(k).FontWeight ='normal';
                lblChannelHandles(k).BackgroundColor = 'none';
                lblChannelHandles(k).FontColor = 1-figColor;
            end
        end

        drawnow limitrate;  % Update UI in batch
        plotOption();
    end

    function onAllSelectChecked(val)
        % Batch update selection checkboxes
        if val
            chkBoxSelect(:) = 1;
        else
            chkBoxSelect(:) = 0;
        end
        
        startIdx = (currentPage-1) * PAGE_SIZE + 1;
        endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
        visibleGroups = startIdx:endIdx;

        % Batch update UI elements
        for idx = 1:length(visibleGroups)
            k = visibleGroups(idx);
            groupSelectCheckboxes(k).Value = val;

            if val
                c = getGroupColor(groupList(k));
                groupSelectCheckboxes(k).FontColor = c;
                groupSelectCheckboxes(k).FontWeight ='Bold';
                lblChannelHandles(k).BackgroundColor = c;
                lblChannelHandles(k).FontColor = bestTextColorFor(c);
            else
                groupSelectCheckboxes(k).FontColor = 1-figColor;
                groupSelectCheckboxes(k).FontWeight ='normal';
                lblChannelHandles(k).BackgroundColor = 'none';
                lblChannelHandles(k).FontColor = 1-figColor;
            end
        end
        drawnow limitrate;  % Update UI in batch
    end

    function onVisCheckboxChanged(idx)
        c = getGroupColor(groupList(idx));
        if groupVisCheckboxes(idx).Value
            chkBoxVis(idx) = 1;
            selected_order = [selected_order,idx];
            groupVisCheckboxes(idx).FontColor = c;
            groupVisCheckboxes(idx).FontWeight ='Bold';
        else
            chkBoxVis(idx) = 0;
            selected_order(selected_order==idx) = [];
            groupVisCheckboxes(idx).FontColor = 1-figColor;
            groupVisCheckboxes(idx).FontWeight ='normal';
        end

        if isprop(lblChannelHandles(idx), 'BackgroundColor') && (groupSelectCheckboxes(idx).Value || groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).BackgroundColor = c;
        elseif isprop(lblChannelHandles(idx), 'BackgroundColor') && (~groupSelectCheckboxes(idx).Value && ~groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).BackgroundColor = 'none';
        end

        if isprop(lblChannelHandles(idx), 'FontColor') && (groupSelectCheckboxes(idx).Value || groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).FontColor = bestTextColorFor(c);
        elseif isprop(lblChannelHandles(idx), 'FontColor') && (~groupSelectCheckboxes(idx).Value && ~groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).FontColor = 1-figColor;
        end
        drawnow limitrate;
        selected_order = unique(selected_order,'stable');
        plotOption();
    end

    function onSelectCheckboxChanged(idx)
        c = getGroupColor(groupList(idx));
        if groupSelectCheckboxes(idx).Value
            chkBoxSelect(idx) = 1;
            groupRadioButtons(idx).Value = true;
            groupSelectCheckboxes(idx).FontColor = c;
            groupSelectCheckboxes(idx).FontWeight ='Bold';
        else
            chkBoxSelect(idx) = 0;
            groupSelectCheckboxes(idx).FontColor = 1-figColor;
            groupSelectCheckboxes(idx).FontWeight ='normal';
        end

        if isprop(lblChannelHandles(idx), 'BackgroundColor') && (groupSelectCheckboxes(idx).Value || groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).BackgroundColor = c;
        elseif isprop(lblChannelHandles(idx), 'BackgroundColor') && (~groupSelectCheckboxes(idx).Value && ~groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).BackgroundColor = 'none';
        end

        if isprop(lblChannelHandles(idx), 'FontColor') && (groupSelectCheckboxes(idx).Value || groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).FontColor = bestTextColorFor(c);
        elseif isprop(lblChannelHandles(idx), 'FontColor') && (~groupSelectCheckboxes(idx).Value && ~groupVisCheckboxes(idx).Value)
            lblChannelHandles(idx).FontColor = 1-figColor;
        end
        drawnow limitrate;
        sliderReAlign.Value = 0;
        updateRealignSpikes(0);
        drawnow limitrate;
    end

    function onRemove()
        % For each group that is selected => set label = NaN => color => gray
        selectedIdx = find(chkBoxSelect);
        for k = selectedIdx'
            sortedRes.unifiedLabels(sortedRes.unifiedLabels==groupList(k)) = NaN;
            sortedRes.channelNum(sortedRes.unifiedLabels==groupList(k)) = NaN;
            groupList(k) = NaN;
            channelList(k) = NaN;
        end
        refreshLabels();
        plotOption();
    end

    function onMerge()
        % Find which radio is selected
        selIdx = [];
        for i = 1:numGroups
            if isvalid(groupRadioButtons(i)) && isprop(groupRadioButtons(i),'Value')
                if groupRadioButtons(i).Value
                    selIdx = i;
                    break;
                end
            end
        end

        if isempty(selIdx)
            uialert(parentFig,'No valid group selected for merging.','Merge Warning');
            return;
        end

        mergeLabel = groupList(selIdx);
        mergeChannel = channelList(selIdx);

        if isnan(mergeLabel)
            uialert(parentFig,'Cannot merge into a removed (NaN) group.','Merge Error');
            return;
        end

        if ~any([groupSelectCheckboxes(groupList == selIdx).Value] == 1) && ~any([groupSelectCheckboxes(originalGroupList == selIdx).Value] == 1)
            uialert(parentFig,'Inconsistent merging! Merge into one of the selected groups.','Merge Error');
            return;
        end

        % Perform the merge => all selected become the same label
        selectedIdx = find(chkBoxSelect);
        for k = selectedIdx'
            mergedFlag(k) = true;
            sortedRes.unifiedLabels(sortedRes.unifiedLabels==groupList(k)) = mergeLabel;
            sortedRes.channelNum(sortedRes.unifiedLabels==groupList(k))    = mergeChannel;
            groupList(k)   = mergeLabel;
            channelList(k) = mergeChannel;
            stable_length(k,:) = stable_length(selIdx,:);
        end
        updateChannelPlot();
        refreshLabels();
        plotOption();
    end

    function onReset()
        groupList = originalGroupList;
        channelList = originalChannelList;
        sampleWaveform = originalSampleWaveform;
        sortedRes.unifiedLabels = originalUnifiedLabels;
        sortedRes.channelNum = originalChannelNum;
        sortedRes.spike_idx = originalSpikeIdx;
        mergedFlag = false(numGroups,1);
        chkBoxSelect(:) = 0;
        chkBoxVis(:) = 0;
        lastGroup = 0;
        lastUnit = 0;
        stable_length = [ones(numGroups,1), trialLength * ones(numGroups,1)];
        selected_order = [];

        % Uncheck all selection and visibility for current page
        startIdx = (currentPage-1) * PAGE_SIZE + 1;
        endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
        for k = startIdx:endIdx
            if isvalid(groupSelectCheckboxes(k))
                groupSelectCheckboxes(k).Value = false;
                groupVisCheckboxes(k).Value = false;
            end
        end




        selected_order = [];
        unitIsolation = repmat({'NA'},[numGroups,1]);
        refreshDensity();
        refreshISI();
        refreshCCG();
        refreshLabels();
        chkAllSel.Value = false;
        chkAllVis.Value = false;
        rdoNone.Value = true;
        updateChannelPlot();
        plotOption();
    end

    function onMUA()
        visibleGroups = find(chkBoxSelect);
        if ~isempty(visibleGroups)
            currentType = unitIsolation{visibleGroups(1)};
            typeMap = {'NA', 'SUA+'; 'SUA+', 'SUA'; 'SUA', 'MUA+'; 'MUA+', 'MUA'; 'MUA', 'NA'};
            nextType = typeMap{strcmp(typeMap(:,1), currentType), 2};
            
            % Batch update
            for i = 1:length(visibleGroups)
                unitIsolation{visibleGroups(i)} = nextType;
            end
        end
        refreshLabels();
    end

    function onUndo()
        selectedIdx = find(chkBoxSelect);
        for k = selectedIdx'
            groupList(k) = originalGroupList(k);
            channelList(k) = originalChannelList(k);
            sampleWaveform(k,:,:) = originalSampleWaveform(k,:,:);
            sortedRes.unifiedLabels(originalUnifiedLabels==groupList(k)) = groupList(k);
            sortedRes.channelNum (originalUnifiedLabels==groupList(k)) = channelList(k);
            sortedRes.spike_idx(originalUnifiedLabels==groupList(k))   = originalSpikeIdx(originalUnifiedLabels==groupList(k));
            mergedFlag(k,1) = false;
            unitIsolation{k} = 'NA';
            stable_length(k,:) = [1, trialLength];
        end
        refreshDensity();
        refreshISI();
        refreshCCG();
        refreshLabels();
        plotOption();
    end

    function onAutoCut()
        selectedIdx = find(chkBoxSelect);
        for k = selectedIdx'
            stable_length(k,:) = [1, trialLength];
            if ~isempty(DenCounts{k,k})
                [startIdx, endIdx] = findStableInterval(DenCounts{k,k});
                x_start = DenCenters{k,k}(startIdx);
                stable_length(k,1) = trialLength * (x_start-DenCenters{k,k}(1)) / DenCenters{k,k}(end);
                if stable_length(k,1) < 1
                    stable_length(k,1) = 1;
                elseif stable_length(k,1) > trialLength
                    stable_length(k,1) = trialLength;
                end
                x_end = DenCenters{k,k}(endIdx);
                stable_length(k,2) = trialLength * (x_end) / DenCenters{k,k}(end);
                if stable_length(k,2) < 1
                    stable_length(k,2) = 1;
                elseif stable_length(k,2) > trialLength
                    stable_length(k,2) = trialLength;
                end
            end
        end
        plotOption();

        function [startIdx, endIdx] = findStableInterval(d)
            n = numel(d);
            if n < 3
                startIdx = 1; endIdx = n;
                return
            end
            cp = findchangepts(d, 'Statistic', 'rms','MinThreshold',5);

            startIdx = 1; endIdx = n;
            if isscalar(cp)
                if mean(d(1:cp)) < mean(d(cp+1:end))/dropRate
                    startIdx = cp+1; endIdx = n;
                elseif mean(d(1:cp))/dropRate > mean(d(cp+1:end))
                    startIdx = 1; endIdx = cp;
                end
            elseif numel(cp)>=2
                for i = 1:numel(cp)
                    if mean(d(1:cp(i))) < mean(d(cp(i)+1:end))/dropRate && mean(d(1:cp(i))) < mean(d(cp(i)+1:min(2*cp(i),n)))/dropRate
                        startIdx = cp(i)+1;
                    elseif mean(d(startIdx:cp(i)))/dropRate> mean(d(cp(i)+1:end)) && mean(d(startIdx:cp(i)))/dropRate > mean(d(cp(i)+1:min(2*cp(i),n)))
                        endIdx = cp(i);
                        break
                    end
                end
            end
        end
    end

    function onSave()
        lblSave.Text = "Saving in Progress...";
        drawnow;
        persistent validSpikes
        validSpikes = zeros(size(sortedRes.spike_idx));
        excludeSpikes();
        outputFolder = fullfile(cfg.outputFolder,'RES_Sorted');

        % Create the output directory if it doesn't exist
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end

        sorted_out.unifiedLabels_curated{1,1} = sortedRes.unifiedLabels;
        sorted_out.channelNum_curated{1,1} = sortedRes.channelNum;
        sorted_out.spike_idx_curated{1,1} = sortedRes.spike_idx;
        sorted_out.inclusion_curated{1,1} = validSpikes;
        sorted_out.waveforms_curated{1,1} = extractCuratedWaveforms();
        numFields = numel(fieldnames(sorted_out));
        saveh5SpikeData_gui(outputFolder, sorted_out, zeros(numFields,1));

        curatedSamples.unifiedLabels = unique(groupList);
        curatedSamples.channelNum    = channelList(unique(groupList));
        curatedSamples.waveform      = sampleWaveform(unique(groupList),:,:);
        curatedSamples.unitIsolation = unitIsolation(unique(groupList));
        curatedSamples.validInterval = stable_length(unique(groupList),:);
        curatedSamples.spikePolarity = mainPolarity(unique(groupList),:);
        curatedSamples.original.unifiedLabels = originalUnifiedLabels;
        curatedSamples.original.channelNum    = originalChannelNum;
        curatedSamples.original.waveform      = originalSampleWaveform;
        curatedSamples.original.unitIsolation = unitIsolation;

        filename = fullfile(outputFolder,'curated_sample.mat');
        save(filename,'curatedSamples', '-v7.3');

        uialert(parentFig,'Curated sorted data saved successfully!','Success','Icon','success')
        disp('Saved Curated Data');
        lblSave.Text = "Saved Curated Data!";
        drawnow;

        function excludeSpikes()
            % Vectorized exclusion
            for i = 1:numGroups
                idx = find(sortedRes.unifiedLabels == groupList(i));
                spk_idx = sortedRes.spike_idx(idx);
                valid_idx = spk_idx <= stable_length(i, 2) & spk_idx >= stable_length(i, 1);
                validSpikes(idx(valid_idx)) = 1;
            end
        end

        function waveforms = extractCuratedWaveforms()
            numSpikeSamples = round(halfSpikeWaveDur * cfg.samplingFrequency /1000);
            waveformLength = 2*numSpikeSamples + 1;

            if isempty(mappedData)
                mappedData = map_input_file(cfg.fullFilePath, cfg);
            end

            BATCH_SIZE = 5000;  % Increased batch size for better performance
            total_spikes = length(sortedRes.spike_idx);
            num_batches = ceil(total_spikes / BATCH_SIZE);
            waveforms = zeros(total_spikes, waveformLength);

            for batch = 1:num_batches
                start_idx = (batch-1)*BATCH_SIZE + 1;
                end_idx = min(start_idx + BATCH_SIZE - 1, total_spikes);
                batch_indices = start_idx:end_idx;

                % Vectorized operations
                spike_idx_batch = sortedRes.spike_idx(batch_indices);
                channel_batch = channel_mapping(sortedRes.channelNum(batch_indices));
                
                % Extract waveforms more efficiently
                for i = 1:length(batch_indices)
                    idx_range = spike_idx_batch(i) + (-numSpikeSamples:numSpikeSamples);
                    if all(idx_range > 0 & idx_range <= size(mappedData.Data.data, 2))
                        waveforms(batch_indices(i), :) = mappedData.Data.data(channel_batch(i), idx_range);
                    end
                end
            end
        end
    end

    function updateScale(sVal)
        ampScale = sVal;
        plotOption();
    end

    function updateLineWidth(sVal)
        if sVal == 0
            sVal = eps;
        end
        lineWidth = sVal;
        plotOption();
    end

    function updateAlphaLevel(sVal)
        if sVal == 0
            sVal = eps;
        end
        alphaLevel = sVal;
        plotOption();
    end

    function updateRealignSpikes(sVal)

        selectedGroups = find(chkBoxSelect);

        if isempty(selectedGroups)
            return;
        end

        shiftVal = round(sVal * cfg.samplingFrequency/1000);
        
        for k = selectedGroups'
            % Vectorized spike realignment
            shiftedSpikesIdx = sortedRes.unifiedLabels == groupList(k);
            sortedRes.spike_idx(shiftedSpikesIdx) = sortedRes.spike_idx(shiftedSpikesIdx) - shiftVal;
            
            % Update sample waveform
            if shiftVal ~= 0
                sampleWaveform(k,:,:) = circshift(sampleWaveform(k,:,:), shiftVal, 3);
            end
        end

        if sVal~= 0
            sliderReAlign.Value = 0;
            updateRealignSpikes(0);
            drawnow limitrate;
        end

        plotOption();
    end

    function updateExportFormat(ddVal)
        exportType = ddVal;
    end

    function setDropVal(val)
        dropRate = val;
    end

    function onExportFig()
        saveFigPath = fullfile(cfg.outputFolder, 'exported_figures');
        if ~exist(saveFigPath,'dir')
            mkdir(saveFigPath);
        end
        files = dir(fullfile(saveFigPath, ['exported_Fig_' plotType '*.eps']));
        if isempty(files)
            nextNum = 1;
        else
            nums = cellfun(@(s) sscanf(s, ['exported_Fig_' plotType '%d.eps']), {files.name});
            nextNum = max(nums) + 1;
        end
        filename = fullfile(saveFigPath, sprintf(['exported_Fig_' plotType '%d.eps'], nextNum));
        exportgraphics(axChannels, filename,'BackgroundColor','white','ContentType', exportType);
    end

    function setDistYVal(val)
        yDistThr = val;
    end


    function setDistVal(val,type)

disimlarityScore(15,14)
        if nargin>0
            if type == 0
                if strcmp(distanceEstType,'XCorr')
                    distThr = val;
                else
                    distThr = prctile(disimlarityScore,(1-val)*100,"all");
                end
            else
                ampThr = val;
            end
        end

        switch distanceEstType
            case 'KL-div'
                [rowGroup, colGroup] = find(disimlarityScore <= distThr & ampSimilarity <= ampThr);
            case 'Bhattacharyya'
                [rowGroup, colGroup] = find(disimlarityScore <= distThr & ampSimilarity <= ampThr);
            case 'XCorr'
                [rowGroup, colGroup] = find(disimlarityScore >= distThr & ampSimilarity <= ampThr);
        end

        if chkInclusion.Value
            keepPair = incChan(rowGroup) & incChan(colGroup);
            rowGroup = rowGroup(keepPair);
            colGroup = colGroup(keepPair);
        elseif chkExclusion.Value
            keepPair = exChan(rowGroup) & exChan(colGroup);
            rowGroup = rowGroup(keepPair);
            colGroup = colGroup(keepPair);
        end

        pairNum = ['/' num2str(length(colGroup))];
        changeuifPar('0',1);
        changeuifParAll(pairNum);
    end

    function changeuifPar(val,sv)
        if ~sv
            inGroup = str2double(val);
            inGroup = min(max(inGroup,0), length(colGroup));
            lastGroup = inGroup;
            onSelectSimilarPairs(0);
        else
            uifPair.Value = val;
        end
        drawnow limitrate;
    end

    function changeuifParAll(val)
        uifPair2.Value = val;
        drawnow limitrate;
    end

    function setRateIncVal(val)
        rateInc = preprocessed.firingRate > val;
        updateIncExcCheck();
    end

    function setISIIncVal(val)
        isiInc = preprocessed.isiViolation < val;
        updateIncExcCheck();
    end

    function setSNRIncVal(val)
        snrInc = 1+detectblity > val;
        updateIncExcCheck();
    end

    function onINCEXC(val, type)
        if ~type
            if val
                chkInclusion.Value = false;
            end
            drawnow limitrate;
        else
            if val
                chkExclusion.Value = false;
                drawnow limitrate;
            end
        end
        updateIncExcCheck();
        setDistVal();
    end

    function onToggleVis(val)
        updateIncExcCheck();
        if ~val
            onAllVisChecked(val);
        end
    end

    function updateIncExcCheck()
        if chkInclusion.Value
            lastUnit = 0;
            incChan = isiInc & rateInc & snrInc & polarityInc;
            selectUnits = find(incChan);

            startIdx = (currentPage-1) * PAGE_SIZE + 1;
            endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
            visibleGroups = startIdx:endIdx;
            chkBoxSelect(:) = 0;
            chkBoxSelect(selectUnits) = 1;
            if chkTogVis.Value
                chkBoxVis(:) = 0;
                chkBoxVis(selectUnits) = 1;
                selected_order = selectUnits';
            else
                selected_order = [];
            end
            
            for i = visibleGroups
                if isvalid(groupSelectCheckboxes(i))
                    groupSelectCheckboxes(i).Value = incChan(i);
                    if chkTogVis.Value
                        groupVisCheckboxes(i).Value = incChan(i);
                    end
                end
            end
            unitNum = sprintf('/%d',length(selectUnits));
            changeuifUnitAll(unitNum);
            changeuifUnit(num2str(lastUnit),1);
            
        elseif ~chkExclusion.Value
            selected_order = [];
            onAllSelectChecked(false);
            onAllVisChecked(false);
            selectUnits = groupList;
            lastUnit = 0;
            unitNum = sprintf('/%d',length(selectUnits));
            changeuifUnitAll(unitNum);
            changeuifUnit(num2str(lastUnit),1);

        end

        if chkExclusion.Value
            lastUnit = 0;
            exChan = ~isiInc | ~rateInc | ~snrInc & polarityInc;
            selectUnits = find(exChan);

            startIdx = (currentPage-1) * PAGE_SIZE + 1;
            endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);
            visibleGroups = startIdx:endIdx;
            chkBoxSelect(:) = 0;
            chkBoxSelect(selectUnits) = 1;
            if chkTogVis.Value
                chkBoxVis(:) = 0;
                chkBoxVis(selectUnits) = 1;
            end
            for i = visibleGroups
                if isvalid(groupSelectCheckboxes(i))
                    groupSelectCheckboxes(i).Value = exChan(i);
                    if chkTogVis.Value
                        groupVisCheckboxes(i).Value = exChan(i);
                        if groupVisCheckboxes(i).Value
                            selected_order = [selected_order, i];
                        else
                            selected_order(selected_order==i) = [];
                        end
                    end
                end
            end
            unitNum = sprintf('/%d',length(selectUnits));
            changeuifUnitAll(unitNum);
            changeuifUnit(num2str(lastUnit),1);
        elseif ~chkInclusion.Value
            selected_order = [];
            onAllSelectChecked(false);
            onAllVisChecked(false);
            selectUnits = groupList;
            lastUnit = 0;
            unitNum = sprintf('/%d',length(selectUnits));
            changeuifUnitAll(unitNum);
            changeuifUnit(num2str(lastUnit),1);
        end

        selected_order = unique(selected_order,'stable');
        if chkTogVis.Value
            plotOption();
        end

        if ~isempty(rowGroup)
            setDistVal();
        end
    end


    function changeuifUnit(val,sv)
        if ~sv
            inGroup = str2double(val);
            inGroup = min(max(inGroup,0), length(selectUnits));
            lastUnit = inGroup;
            onSelectUnit(0)
        else
            uifUnit.Value = val;
        end
        drawnow limitrate;
    end

    function changeuifUnitAll(val)
        uifUnit2.Value = val;
        drawnow limitrate;
    end


    function onDeselect
        onAllSelectChecked(false);
        onAllVisChecked(false);
        chkAllSel.Value = false;
        chkAllVis.Value = false;
        selected_order = [];
        plotOption();
        refreshLabels();
    end

    function onSelectUnit(val)
        selected_order = [];
        onAllSelectChecked(false);
        onAllVisChecked(false);
        thisUnit = lastUnit + val;
        if thisUnit >= 1 && thisUnit <= length(selectUnits)
            thisCol = selectUnits(thisUnit);

            % Calculate which page contains this unit
            targetPage = ceil(thisCol / PAGE_SIZE);
            if targetPage ~= currentPage
                currentPage = targetPage;
                onChangePage(0); % Refresh page without changing page number
            end

            % Make sure the unit is now visible
            if isvalid(groupVisCheckboxes(thisCol))
                groupVisCheckboxes(thisCol).Value = true;
                groupSelectCheckboxes(thisCol).Value = true;
                chkBoxSelect(thisCol) = 1;
                chkBoxVis(thisCol) = 1;
                selected_order = [thisCol];
                groupRadioButtons(thisCol).Value = true;
                c = getGroupColor(groupList(thisCol));
                groupVisCheckboxes(thisCol).FontColor = c;
                groupVisCheckboxes(thisCol).FontWeight ='Bold';
                groupSelectCheckboxes(thisCol).FontColor = c;
                groupSelectCheckboxes(thisCol).FontWeight ='Bold';
                lblChannelHandles(thisCol).BackgroundColor = c;
                lblChannelHandles(thisCol).FontColor = bestTextColorFor(c);
                lastUnit = thisUnit;
            else
                lastUnit = 0;
            end
        else
            lastUnit = 0;
        end
        if val~=0
        unitNum = sprintf('/%d',length(selectUnits));
        changeuifUnit(num2str(lastUnit),1);
        changeuifUnitAll(unitNum);
        end
        plotOption()
    end

    function onSelectSimilarPairs(val)
        selected_order = [];
        onAllSelectChecked(false);
        onAllVisChecked(false);
        thisGroup = lastGroup + val;
        if thisGroup >= 1 && thisGroup <= length(colGroup)
            thisCol = colGroup(thisGroup);
            thisRow = rowGroup(thisGroup);

            % Calculate which page contains these units
            targetPageCol = ceil(thisCol / PAGE_SIZE);
            targetPageRow = ceil(thisRow / PAGE_SIZE);

            % If both units are on different pages, navigate to the first one's page
            if targetPageCol ~= currentPage && targetPageRow ~= currentPage
                currentPage = max(targetPageCol,targetPageRow);
                onChangePage(0); % Refresh page without changing page number
            elseif targetPageCol ~= currentPage
                currentPage = max(targetPageCol,targetPageRow);
                onChangePage(0);
            elseif targetPageRow ~= currentPage
                currentPage = max(targetPageCol,targetPageRow);
                onChangePage(0);
            end

            chkBoxVis(thisCol) = 1;
            chkBoxVis(thisRow) = 1;
            chkBoxSelect(thisCol) = 1;
            chkBoxSelect(thisRow) = 1;
            selected_order = [thisCol,thisRow];

            if isvalid(groupVisCheckboxes(thisCol)) && isprop(groupVisCheckboxes(thisCol),'Value')

                groupVisCheckboxes(thisCol).Value    = true;
                groupSelectCheckboxes(thisCol).Value = true;

                c = getGroupColor(groupList(thisCol));
                groupVisCheckboxes(thisCol).FontColor    = c;
                groupVisCheckboxes(thisCol).FontWeight   = 'Bold';
                groupSelectCheckboxes(thisCol).FontColor  = c;
                groupSelectCheckboxes(thisCol).FontWeight = 'Bold';
                lblChannelHandles(thisCol).BackgroundColor = c;
                lblChannelHandles(thisCol).FontColor       = bestTextColorFor(c);
            end

            if isvalid(groupVisCheckboxes(thisRow)) && isprop(groupVisCheckboxes(thisRow),'Value')
                % Make checkboxes checked                 
                groupVisCheckboxes(thisRow).Value    = true;
                groupSelectCheckboxes(thisRow).Value = true;

                c = getGroupColor(groupList(thisRow));
                groupVisCheckboxes(thisRow).FontColor    = c;
                groupVisCheckboxes(thisRow).FontWeight   = 'Bold';
                groupSelectCheckboxes(thisRow).FontColor  = c;
                groupSelectCheckboxes(thisRow).FontWeight = 'Bold';
                lblChannelHandles(thisRow).BackgroundColor = c;
                lblChannelHandles(thisRow).FontColor       = bestTextColorFor(c);
                
            end

            selected_order = [thisCol,thisRow];
            if preprocessed.firingRate(thisCol) > preprocessed.firingRate(thisRow)
                if isvalid(groupRadioButtons(thisCol))  && isprop(groupRadioButtons(thisCol),'Value')
                    groupRadioButtons(thisCol).Value = true;
                end
            else
                if isvalid(groupRadioButtons(thisRow)) && isprop(groupRadioButtons(thisRow),'Value')
                    groupRadioButtons(thisRow).Value = true;
                end
            end
            lastGroup = thisGroup;

        else
            lastGroup = 0;
        end
        
        if val~=0
        pairNum = sprintf('/%d',length(colGroup));
        changeuifPar(num2str(lastGroup),1);
        changeuifParAll(pairNum);
        end
        plotOption()
    end

    function updatePolarity(ddVal)
        switch ddVal
            case 'All'
                polarityInc = true(numGroups,1);
            case 'Main Neg.'
                polarityInc = mainPolarity;
            case 'Main Pos.'
                polarityInc = ~mainPolarity;
            case 'All Pos.'
                polarityInc = ~sidePolarity & ~mainPolarity;
            case '1 Chan Neg.'
                polarityInc = sidePolarity;
        end
        updateIncExcCheck()
    end

    function updateDistance(ddVal)
        distanceEstType = ddVal;
    end

%% All supporting functions
    function refreshLabels()
        startIdx = (currentPage-1) * PAGE_SIZE + 1;
        endIdx = min(startIdx + PAGE_SIZE - 1, numGroups);

        for k = startIdx:endIdx
            if ~isvalid(lblGroupHandles(k))
                continue;
            end

            lblVal = groupList(k);
            if isnan(lblVal)
                txt = 'G: NaN';
            else
                txt = sprintf('G: %g', lblVal);
            end
            lblGroupHandles(k).Text = txt;

            lblChannelVal = channelList(k);
            if isnan(lblChannelVal)
                txt = 'Ch NaN';
            else
                txt = sprintf('Ch %g', lblChannelVal);
            end
            lblChannelHandles(k).Text = txt;

            if isprop(lblIsolation(k),'FontColor')
                switch unitIsolation{k}
                    case 'NA'
                        lblIsolation(k).FontColor = [.5 .5 .5];
                    case 'SUA+'
                        lblIsolation(k).FontColor = [.2 .9 .2];
                    case 'SUA'
                        lblIsolation(k).FontColor = [.6 .9 .2];
                    case 'MUA+'
                        lblIsolation(k).FontColor = [.9 .6 .2];
                    case 'MUA'
                        lblIsolation(k).FontColor = [.9 .2 .2];
                end
                lblIsolation(k).Text = unitIsolation{k};
            end

            txt = sprintf('%.1f%%', preprocessed.isiViolation(k));
            lblISIViolation(k).Text = txt;

            c = getGroupColor(lblVal);
            if isprop(lblGroupHandles(k), 'BackgroundColor')
                lblGroupHandles(k).BackgroundColor = c;
            end
            if isprop(lblGroupHandles(k), 'FontColor')
                lblGroupHandles(k).FontColor = bestTextColorFor(c);
            end

            % Channel label
            if isprop(lblChannelHandles(k), 'BackgroundColor') && ...
                    (isvalid(groupSelectCheckboxes(k)) && (groupSelectCheckboxes(k).Value || groupVisCheckboxes(k).Value))
                lblChannelHandles(k).BackgroundColor = c;
            else
                lblChannelHandles(k).BackgroundColor = figColor;
            end

            if isprop(lblChannelHandles(k), 'FontColor') && ...
                    (isvalid(groupSelectCheckboxes(k)) && (groupSelectCheckboxes(k).Value || groupVisCheckboxes(k).Value))
                lblChannelHandles(k).FontColor = bestTextColorFor(c);
            else
                lblChannelHandles(k).FontColor = 1-figColor;
            end

            if isprop(groupVisCheckboxes(k), 'FontColor')
                if groupVisCheckboxes(k).Value
                    groupVisCheckboxes(k).FontColor = c;
                else
                    groupVisCheckboxes(k).FontColor = 1-figColor;
                end
            end

            if isprop(groupSelectCheckboxes(k), 'FontColor')
                if groupSelectCheckboxes(k).Value
                    groupSelectCheckboxes(k).FontColor = c;
                else
                    groupSelectCheckboxes(k).FontColor = 1-figColor;
                end
            end

            if isprop(groupRadioButtons(k), 'FontColor')
                groupRadioButtons(k).FontColor = c;
            end
        end
    end

    function plotOption()
        switch plotType
            case 'Trace'
                multipleCheckbox('off');
                plot_on_Signal();
            case 'Amplitude'
                multipleCheckbox('off');
                plotAmplitude();
            case 'Features'
                multipleCheckbox('off');
                plotFeatures();
            case 'Waveform'
                multipleCheckbox('off');
                plotWaveforms();
            case 'CCG'
                multipleCheckbox('off');
                plotCCG();
            case 'ISI'
                multipleCheckbox('off');
                plotISI();
            case 'Density'
                multipleCheckbox('off');
                plotDensity();
            case 'PCs'
                multipleCheckbox('off');
                plotChannelPCA();
            case 'Multiple'
                multipleCheckbox('on');
                plotMultiple();
        end
    end

    function multipleCheckbox(stateVal)
        if strcmp(stateVal,'on') && ~isprop(multiplePlotObj(1),'Value')
            for i = 1:length(multiplePlotObj)
                multiplePlotObj(i) = uicheckbox(plotButtonGrid,...
                    'Value',false, 'Text',plotItemsNames{i},...
                    'ValueChangedFcn',@(cb,ev)onPlotCheckboxChanged(i));
                multiplePlotObj(i).Visible = 'on';
                multiplePlotObj(i).Layout.Row = 3;
                multiplePlotObj(i).Layout.Column = i;
                multiplePlotObj(i).FontColor = 1-figColor;

                multipleSpinnerObj(i) = uispinner(plotButtonGrid, ...
                    'Limits', [1 10], ...
                    'Value', 1, ...
                    'ValueChangedFcn', @(sp,ev) updateScale(i, sp.Value));
                multipleSpinnerObj(i).Visible = 'on';
                multipleSpinnerObj(i).Layout.Row = 4;  % Positioned below the checkbox
                multipleSpinnerObj(i).Layout.Column = i;
                applyColorScheme(multipleSpinnerObj(i),figColor);
            end
        elseif strcmp(stateVal,'off')
            delete(multiplePlotObj(isvalid(multiplePlotObj)));
            delete(multipleSpinnerObj(isvalid(multipleSpinnerObj)));
            multiplePlotOrder = [];
            multipleSpinnerObj = gobjects(length(plotItemsNames),1);
            multiplePlotObj = gobjects(length(plotItemsNames),1);
            plotScaleFactor = ones(length(plotItemsNames),1);
            reScaledFlag = false;
        end

        function updateScale(idx, Val)
            plotScaleFactor(idx) = Val;
            reScaledFlag = true;
            plotOption();
        end
    end

    function onPlotCheckboxChanged(idx)
        if multiplePlotObj(idx).Value
            multiplePlotOrder = [multiplePlotOrder,idx];
        else
            multiplePlotOrder(multiplePlotOrder==idx) = [];
        end
        multiplePlotOrder = unique(multiplePlotOrder,'stable');
        plotOption();
    end

    function plotMultiple()
        % Optimize multiple plot display
        disp('--- Plot Multiple called ---');

        if strcmp(lastPlotted,'Multiple')
            if ~isequal(lastMultiPlotOrder,multiplePlotOrder) || reScaledFlag
                delete(allchild(axChannels));
            end
        else
            delete(allchild(axChannels));
            delete(axChannels);
            axChannels = uipanel(bottomRightGrid,'BorderType','none',...
                'BackgroundColor',figColor);
            axChannels.Layout.Column = [2 3];
            axChannels.Layout.Row    = [1 3];
            applyColorScheme(axChannels, figColor);
        end
        lastPlotted = plotType;

        countPlot = length(multiplePlotOrder);
        if countPlot>0
            if ~isequal(lastMultiPlotOrder,multiplePlotOrder) || reScaledFlag
                if countPlot == 1
                    gridXY = [1, 1];
                elseif countPlot == 2
                    gridXY = [1, 2];
                elseif countPlot > 2 && countPlot < 5
                    gridXY = [2, 2];
                elseif countPlot >=5
                    gridXY = [2, 3];
                end

                scaleFactors = plotScaleFactor(multiplePlotOrder(1:countPlot));
                rows = zeros(countPlot,1);
                cols = zeros(countPlot,1);
                for k = 1:countPlot
                    [rows(k), cols(k)] = ind2sub([gridXY(1), gridXY(2)], k);
                end

                rowMax = ones(gridXY(1),1);
                for r = 1:gridXY(1)
                    mask = rows == r;
                    if any(mask)
                        rowMax(r) = max(scaleFactors(mask));
                    end
                end

                colMax = ones(1, gridXY(2));
                for c = 1:gridXY(2)
                    mask = cols == c;
                    if any(mask)
                        colMax(c) = max(scaleFactors(mask));
                    end
                end

                rowHeights = arrayfun(@(x)sprintf('%.fx',x),rowMax,'UniformOutput',false);
                colWidths = arrayfun(@(x)sprintf('%.fx',x),colMax,'UniformOutput',false);

                gridMultiple = uigridlayout(axChannels, gridXY,...
                    'RowHeight',rowHeights, 'ColumnWidth',colWidths,...
                    'Padding',0, 'RowSpacing',1,'ColumnSpacing',1);
                applyColorScheme(gridMultiple, figColor);

                multiPlotPanels = gobjects(countPlot,1);
                for k = 1:countPlot
                    multiPlotPanels(k) = uipanel(gridMultiple,"BorderColor",1-figColor);
                    [row,col] = ind2sub([gridXY(1), gridXY(2)],k);
                    multiPlotPanels(k).Layout.Row = row;
                    multiPlotPanels(k).Layout.Column = col;
                    applyColorScheme(multiPlotPanels(k), figColor);
                end
            end

            for k = 1:countPlot
                mltiPlotIdx = multiplePlotOrder(k);
                delete(allchild(multiPlotPanels(k)));
                switch mltiPlotIdx
                    case 1
                        multiPlotPanels(k).Title = 'Amplitude';
                        plotAmplitude(multiPlotPanels(k))
                    case 2
                        multiPlotPanels(k).Title = 'Waveform';
                        plotWaveforms(multiPlotPanels(k))
                    case 3
                        multiPlotPanels(k).Title = 'Correlogram:';
                        plotCCG(multiPlotPanels(k))
                    case 4
                        multiPlotPanels(k).Title = 'ISI Violation %:';
                        plotISI(multiPlotPanels(k))
                    case 5
                        multiPlotPanels(k).Title = 'Presence Ratio:';
                        plotDensity(multiPlotPanels(k))
                    case 6
                        multiPlotPanels(k).Title = 'Channel PC';
                        plotChannelPCA(multiPlotPanels(k))
                    case 7
                        multiPlotPanels(k).Title = 'Trace';
                        plot_on_Signal(multiPlotPanels(k))
                    case 8
                        multiPlotPanels(k).Title = 'Features';
                        plotFeatures(multiPlotPanels(k))
                end
            end
        end

        lastMultiPlotOrder = multiplePlotOrder;
        reScaledFlag = false;
    end

    function normalizeTimeAmp()
        [xNorm, ~] = mapminmax(xlocs, -1, 1);
        [yNorm, ~] = mapminmax(ylocs, 0, 1);

        sortedX = sort(xNorm(:));
        differences = diff(sortedX);
        nonZeroDiffs = differences(differences > 0);
        if isempty(nonZeroDiffs)
            minDistanceX = 1;
        else
            minDistanceX = min(nonZeroDiffs);
        end

        sortedY = sort(yNorm(:));
        differencesY = diff(sortedY);
        nonZeroDiffsY = differencesY(differencesY > 0);
        if isempty(nonZeroDiffsY)
            minDistanceY = 1;
        else
            minDistanceY = min(nonZeroDiffsY);
        end

        xNorm_adj = (xNorm * range(waveformXaxis) * 2) ./ minDistanceX;
        yNorm_adj = (yNorm ) ./ minDistanceY;
    end

    function c = getGroupColor(lbl)
        if isnan(lbl)
            c = [0.6 0.6 0.6]; % gray
        else
            idx = mod(round(lbl)-1, numGroups) + 1;
            c   = colorMapAll(idx,:);
        end
    end

    function updateXWindow(sVal, stepStr)
        if strcmp(stepStr,'full')
            xStepVal = trialLength;
        else
            xStepVal = str2double(stepStr);
        end
        xWindowStart = round(sVal);
        xWindowEnd   = xWindowStart + xStepVal;
        if xWindowEnd > xSlider.Limits(2)
            xWindowEnd   = xSlider.Limits(2);
            if xWindowStart < xSlider.Limits(1)
                xWindowStart = xSlider.Limits(1);
            end
            xSlider.Value = xWindowStart;
        end
        cfg.XWindowStart = xWindowStart;
        cfg.XWindowEnd   = xWindowEnd;
        disp(['X-axis window => Start=',num2str(xWindowStart),...
            ', End=',num2str(xWindowEnd),...
            ', Step=',num2str(xStepVal)]);
        plotOption();
    end

    function updateYWindow(sVal, stepStr)
        if strcmp(stepStr,'full')
            yStepVal = numChannels;
        else
            yStepVal = str2double(stepStr);
        end
        yWindowStart = round(sVal);
        yWindowEnd   = yWindowStart + yStepVal;
        if yWindowEnd > ySlider.Limits(2)
            yWindowEnd   = ySlider.Limits(2);
            yWindowStart = yWindowEnd - yStepVal;
            if yWindowStart < ySlider.Limits(1)
                yWindowStart = ySlider.Limits(1);
            end
            ySlider.Value = yWindowStart;
        end
        cfg.YWindowStart = yWindowStart;
        cfg.YWindowEnd   = yWindowEnd;
        disp(['Y-axis window => Start=',num2str(yWindowStart),...
            ', End=',num2str(yWindowEnd),...
            ', Step=',num2str(yStepVal)]);
        plotOption();
    end


    function plot_on_Signal(panelOption)
        disp('--- Plot filtered data function called ---');

        if isempty(mappedData)
            if isfield(cfg, 'inputFolder') && isfield(cfg, 'outputFolder')
                mappedData = map_input_file(cfg.fullFilePath, cfg);
                trialLength = size(mappedData.Data.data,2) / cfg.samplingFrequency;
                xSlider.Limits = [0 trialLength];
                roundTicks = 25:25:trialLength;
                tickInterval = roundTicks(nearest(roundTicks, floor(trialLength/10)));
                xSlider.MajorTicks = 0:tickInterval:trialLength;
                disp('--- Data loaded ---');
            else
                disp('--- Data loading failed ---');
                return;
            end
        end

        % Reuse or create axes
        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'Trace')
                cla(axChannels, 'reset');
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uiaxes(bottomRightGrid);
                axChannels.Layout.Row = [1 3];
                axChannels.Layout.Column = [2 3];
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
            maxLimit = 100;
        else
            plotAxis = uiaxes(panelOption, 'Units', 'normalized', 'Position', [0 0 1 1]);
            maxLimit = 0.1;
        end

        % Set up axes properties
        view(plotAxis, 2);
        set(plotAxis, 'tickdir', 'out', 'Color', figColor, 'XColor', 1 - figColor, 'YColor', 1 - figColor);
        box(plotAxis, 'off');

        disp(['Plot Type = ', plotType]);


        lastPlotted = plotType;

        [sy, sx] = size(mappedData.Data.data);

        xDataRange = round(xWindowStart * cfg.samplingFrequency)+1 : min([round(xWindowEnd * cfg.samplingFrequency), round((xWindowStart+maxLimit) * cfg.samplingFrequency), sx]);
        yDataRange = max(1, yWindowStart) : min(yWindowEnd, sy);
        yDataRange = yDataRange(chan_wave_inclusion(yDataRange));

        minYAx = 1E6;
        maxYAx = -1E6;
        for k = selected_order
            minYAx = min(minYAx,min(channelPlot(k,:)));
            maxYAx = max(maxYAx,max(channelPlot(k,:)));
        end

        if strcmp(plotType,'Multiple') && minYAx ~= 1E6 && maxYAx ~= -1E6
            yDataRange(yDataRange<minYAx | yDataRange>maxYAx ) = [];
        end


        % Find spikes in range
        inRange_mask = sortedRes.spike_idx > xDataRange(1) & sortedRes.spike_idx < xDataRange(end);
        inRange_idx = find(inRange_mask);
        unifiedLabels = sortedRes.unifiedLabels(inRange_idx);
        spike_idx = sortedRes.spike_idx(inRange_idx);

        % Cache data extraction and filtering
        persistent lastXDataRange lastYDataRange lastFilteredData
        if isempty(lastFilteredData) || ~isequal(lastXDataRange, xDataRange) || ~isequal(lastYDataRange, yDataRange)
            inputData = double(mappedData.Data.data(channel_mapping(yDataRange), xDataRange));
            
            if strcmp(plotFilterType,'filtered')
                inputData = bandpass_filter_GUI(inputData', cfg.bandpass, cfg.samplingFrequency)';
            end

            if cfg.denoising
                inputData = remove_shared_noise(inputData, cfg);
            end
            
            lastXDataRange = xDataRange;
            lastYDataRange = yDataRange;
            lastFilteredData = inputData;
        else
            inputData = lastFilteredData;
        end

        waveformMask = nan(size(inputData));


        normIn = inputData./max(abs(inputData),[],'all');
        plot(plotAxis, xDataRange/cfg.samplingFrequency, (yDataRange' + normIn * ampScale/2)', 'LineWidth', lineWidth/2);
        plotAxis = plot_grouped_lines(numel(yDataRange), 0, plotAxis, 0, 0, 0.5);
        hold(plotAxis,'on')
        

        % Trace spike waveforms on top - optimized
        if any(spike_idx) && ~isempty(selected_order)
            for k = selected_order
                if any(ismember(channelPlot(k,:),yDataRange))
                    if chkBoxVis(k)
                        plot_spike_idx = spike_idx(unifiedLabels == groupList(k));
                        if ~isempty(plot_spike_idx)
                            chPlot = channelPlot(k,:);
                            chPlot = (ismember(yDataRange,chPlot));
                            c = getGroupColor(groupList(k));
                            
                            % Vectorized spike waveform plotting
                            plottingIdx = unique(spike_Xaxis + plot_spike_idx(:));
                            validIdx = ismember(xDataRange,plottingIdx);
                            classWaveform = waveformMask;
                            classWaveform(chPlot,validIdx) = normIn(chPlot,validIdx);
                            plot(plotAxis, xDataRange/cfg.samplingFrequency, (yDataRange' + classWaveform * ampScale/2)', 'Color', c, 'LineWidth', lineWidth);
 
                        end
                    end
                end
            end
        end
        hold(plotAxis,'off')
        if ~strcmp(plotType,'Multiple')
            title(plotAxis,'Trace Spike Explorer','Color','white');
            xlabel(plotAxis, 'Time (s)');
            ylabel(plotAxis, 'Channel #');
        end
        set(plotAxis,'tickdir','out');
        plotAxis.Color   = [0.1 0.1 0.1];
        plotAxis.XColor  = [1 1 1];
        plotAxis.YColor  = [1 1 1];
        rangePlot = [xDataRange(1)/cfg.samplingFrequency xDataRange(end)/cfg.samplingFrequency];

        if strcmp(plotType,'Multiple')
            minYAx((minYAx)==1E6) = yWindowStart;
            maxYAx((maxYAx)==-1E6) = yWindowEnd;
            axis(plotAxis,[rangePlot(1) rangePlot(2) minYAx-.5 maxYAx+.5]);
            yTicksIntervals = round((maxYAx-minYAx)/10);
            yTicksIntervals = max(yTicksIntervals,1);
            set(plotAxis,'ytick',[minYAx:yTicksIntervals:maxYAx])
        else
            axis(plotAxis,[rangePlot(1) rangePlot(2) yWindowStart-.5 yWindowEnd+.5]);
            yTicksIntervals = round((yWindowEnd-yWindowStart)/10);
            yTicksIntervals = max(yTicksIntervals,1);
            set(plotAxis,'ytick',[yWindowStart:yTicksIntervals:yWindowEnd])

        end
            
        axis(plotAxis,'tight')
        box(plotAxis, 'off');


    end


    function plotFeatures(panelOption)

        disp('--- Plot Features called ---');

        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'Features')
                cla(axChannels,'reset');
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uiaxes(bottomRightGrid,'Color',figColor);
                axChannels.Layout.Row = [1 3];
                axChannels.Layout.Column = [2 3];
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = uiaxes(panelOption, 'Units', 'normalized', 'Position', [0 0 1 1],'Color',figColor);
        end

        view(plotAxis, 3);
        if ~strcmp(plotType,'Multiple')
            title(plotAxis,'Feature across all channels','Color','white');
        end
        xlabel(plotAxis, 'Time (s)');
        ylabel(plotAxis, 'Norm. PC1');
        zlabel(plotAxis, 'Norm. PC2 on Channel #');
        set(plotAxis,'tickdir','out');
        plotAxis.Color   = figColor;
        plotAxis.XColor  = 1-figColor;
        plotAxis.YColor  = 1-figColor;
        plotAxis.ZColor  = 1-figColor;


        if ~isempty(selected_order) 
            minSliderVal = min(channelPlot(selected_order(end),:))-1;
            maxSliderVal = max(channelPlot(selected_order(end),:))+1;
            if ~isempty(minSliderVal)
                if minSliderVal(end) < yWindowStart
                    yWindowStart = max(floor(minSliderVal(end)),1);
                    yWindowEnd = min(yWindowStart + str2double(yStepDropdown.Value), numChannels);
                    ySlider.Value = yWindowStart;
                    drawnow limitrate;
                elseif  maxSliderVal(end) > yWindowEnd
                    yWindowEnd = min(ceil(maxSliderVal(end)), numChannels);
                    yWindowStart = max(yWindowEnd-str2double(yStepDropdown.Value),1);
                    ySlider.Value = yWindowStart;
                    drawnow limitrate;
                end
            end
        end

        xDataRange = round(xWindowStart * cfg.samplingFrequency)+1 : min(round(xWindowEnd * cfg.samplingFrequency), num_Samples);
        length(xDataRange)
        yDataRange = max(1, yWindowStart) : min(yWindowEnd, numChannels);
        yDataRange = yDataRange(chan_wave_inclusion(yDataRange));

        % Optimized spike finding
        inRange_mask = sortedRes.spike_idx > xDataRange(1) & sortedRes.spike_idx < xDataRange(end);
        inRange_idx = find(inRange_mask);

        unifiedLabels = sortedRes.unifiedLabels(inRange_idx);
        spike_idx = sortedRes.spike_idx(inRange_idx);
        spike_features = sortedRes.features(inRange_idx,:);
        norm_spike_features = mapminmax(spike_features', -.5, .5)';

        hold(plotAxis,'on')

        if any(spike_idx) && ~isempty(selected_order)
            for k = selected_order
                scatterIdx = find(unifiedLabels == groupList(k));
                if length(scatterIdx) > 1E4
                    scatterIdx = scatterIdx(randperm(length(scatterIdx), 1E4));
                end
                plot_spike_idx = spike_idx(scatterIdx);
                if ~isempty(plot_spike_idx)
                    c = getGroupColor(groupList(k));
                    scatter3(plotAxis, plot_spike_idx/cfg.samplingFrequency,norm_spike_features(scatterIdx,1),channelList(k)+norm_spike_features(scatterIdx,2),...
                        25,'markerfacecolor',c,'markeredgecolor','none')
                end
            end
        end

        hold(plotAxis,'off')


        axis(plotAxis,'tight',[xWindowStart xWindowEnd -.5 .5 yWindowStart-.5 yWindowEnd+.5 ]);
        yTicksIntervals = round((yWindowEnd-yWindowStart)/10);
        yTicksIntervals = max(yTicksIntervals,1);
        set(plotAxis,'ztick',[yWindowStart:yTicksIntervals:yWindowEnd])
        set(plotAxis,'ytick',[-.5,0, .5],'ytickLabel',[-1,0, 1])


        selected_order_last = selected_order;
    end


    function plotWaveforms(panelOption)
tic
        disp('--- Plot waveforms called ---');

        if isempty(mappedData)
            if isfield(cfg,'inputFolder') && isfield(cfg,'outputFolder')
                mappedData  = map_input_file(cfg.fullFilePath, cfg);
                trialLength = size(mappedData.Data.data,2) / cfg.samplingFrequency;
                xSlider.Limits = [0 trialLength];
                roundTicks     = 25:25:trialLength;
                tickInterval   = roundTicks(nearest(roundTicks, floor(trialLength/10)));
                xSlider.MajorTicks = 0:tickInterval:trialLength;
                disp('--- Data loaded ---');
            else
                disp('--- Data loading failed ---');
                return;
            end
        end


        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'Waveform')
                cla(axChannels,'reset');
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uiaxes(bottomRightGrid);
                axChannels.Layout.Row    = [1 3];
                axChannels.Layout.Column = [2 3];
            end
            plotAxis   = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = uiaxes(panelOption,'Units','normalized','Position',[0 0 1 1]);
        end

        view(plotAxis,2);
        if ~strcmp(plotType,'Multiple')
            title(plotAxis,'Spike Waveforms','Color','white');
        end
        set(plotAxis,'TickDir','out','Color',figColor,'XColor',1-figColor,'YColor',1-figColor);
        box(plotAxis,'off');
        axis(plotAxis,'off');

        if ~isempty(selected_order) && ~isequal(selected_order, selected_order_last)
            lastSel = selected_order(end);
            chVals  = channelPlot(lastSel,:);
            minVal  = min(chVals) - 1;
            maxVal  = max(chVals) + 1;
            if ~isempty(minVal)
                if minVal(end) < yWindowStart
                    yWindowStart = max(floor(minVal(end)),1);
                    yWindowEnd   = min(yWindowStart + str2double(yStepDropdown.Value), numChannels);
                    ySlider.Value = yWindowStart;
                elseif maxVal(end) > yWindowEnd
                    yWindowEnd   = min(ceil(maxVal(end)), numChannels);
                    yWindowStart = max(yWindowEnd - str2double(yStepDropdown.Value), 1);
                    ySlider.Value = yWindowStart;
                end
                drawnow limitrate;
            end
        end

        xDataRange = round(xWindowStart*cfg.samplingFrequency)+1 : ...
            min(round(xWindowEnd*cfg.samplingFrequency), num_Samples);

        if numel(xDataRange) > 10*cfg.samplingFrequency
            totalBlocks     = floor(numel(xDataRange)/cfg.samplingFrequency);
            selectedBlocks  = sort(randsample(totalBlocks, min(10,totalBlocks)));
            subsampledRange = [];
            for i = 1:numel(selectedBlocks)
                idx = (selectedBlocks(i)-1)*cfg.samplingFrequency + (1:cfg.samplingFrequency);
                subsampledRange = [subsampledRange, xDataRange(idx(idx<=numel(xDataRange)))];
            end
            xDataRange = subsampledRange(subsampledRange>=xDataRange(1) & subsampledRange<=xDataRange(end));
        end


        candidateRange = max(1, yWindowStart) : min(yWindowEnd, numChannels);
        candidateRange = candidateRange(chan_wave_inclusion(candidateRange));


        if ~isempty(selected_order)
            plotChans = unique(channelPlot(selected_order, :));   % channels used by selected units
            plotChans = plotChans(plotChans>0);                  % strip any zeros
            yDataRange = intersect(candidateRange, plotChans, 'stable');
            if isempty(yDataRange)                               % fallback if nothing intersects
                yDataRange = candidateRange;
            end
        else
            yDataRange = candidateRange;
        end


        mask = false(num_Samples,1);
        mask(xDataRange) = true;
        convMask = conv(double(mask), ones(2*numSpikeSamples+1,1), 'same');
        inRange = sortedRes.spike_idx > numSpikeSamples & ...
            sortedRes.spike_idx <= num_Samples-numSpikeSamples & ...
            convMask(sortedRes.spike_idx) >= (2*numSpikeSamples+1);

        unifiedLabels = sortedRes.unifiedLabels(inRange);
        spike_idx     = sortedRes.spike_idx(inRange);
        [~, spike_idx] = ismember(spike_idx, xDataRange);

        persistent cachedInputData lastXDataRange lastYDataRange
        if isempty(cachedInputData) || ~isequal(lastXDataRange,xDataRange) || ~isequal(lastYDataRange,yDataRange)
            inputData       = double(mappedData.Data.data(channel_mapping(yDataRange), xDataRange));
            lastXDataRange  = xDataRange;
            lastYDataRange  = yDataRange;
            cachedInputData = inputData;
        else
            inputData = cachedInputData;
        end

        if strcmp(plotFilterType,'filtered')
            inputData = bandpass_filter_GUI(inputData', cfg.bandpass, cfg.samplingFrequency)';
        end

        if isempty(maxSignal)
            sampleSize = min(size(mappedData.Data.data,2), 10000);
            sampleIdx  = randsample(size(mappedData.Data.data,2), sampleSize);
            maxSignal  = max(abs(double(mappedData.Data.data(channel_mapping(chan_wave_inclusion), sampleIdx))), [], 'all');
        end
        normIn = inputData ./ maxSignal;

        chanPlotList = [];
        isplotIdx    = [];
        hold(plotAxis,'on');

        plotCount = 0; minAx = inf; maxAx = -inf;
        max_YAx = -inf; min_YAx = inf;
        uniqPlots = numel(unique(groupList(selected_order)));

        if any(spike_idx) && ~isempty(selected_order)
            for k = selected_order
                if ismember(groupList(k), isplotIdx)
                    continue;
                end
                isplotIdx = [isplotIdx; groupList(k)];

                plot_spike_idx = spike_idx(unifiedLabels == groupList(k));
                if isempty(plot_spike_idx)
                    chPlot = channelPlot(k, :);
                    chPlot = chPlot(ismember(chPlot, yDataRange));
                    min_t  = waveformXaxis + (uniqPlots/3) * min(xNorm_adj(chPlot));
                    min_t  = min_t + 1.1 * plotCount * range(min_t);
                    max_t  = waveformXaxis + (uniqPlots/3) * max(xNorm_adj(chPlot));
                    max_t  = max_t + 1.1 * plotCount * range(max_t);
                    plotCount = plotCount + 1;
                    continue;
                end

                chPlot = channelPlot(k, :);
                chPlot = chPlot(ismember(chPlot, yDataRange));
                if ~any(chPlot)
                    continue;
                end
                chanPlotList = [chanPlotList, chPlot];
                c = getGroupColor(groupList(k));

                spike_indices = plot_spike_idx(:);
                numSpikes     = numel(spike_indices);
                if ~plotMean && numSpikes > numWaveforms
                    spike_indices = spike_indices(randperm(numSpikes, numWaveforms));
                    numSpikes     = numel(spike_indices);
                end

                all_indices = spike_indices + spike_Xaxis;
                numChPlots  = numel(chPlot);

                if ~plotMean
                    Waveforms = zeros(numSpikes, numChPlots, numel(spike_Xaxis));
                else
                    Waveforms = nan(numSpikes, numChPlots, numel(spike_Xaxis));
                end

                for cIdx = 1:numChPlots
                    chan  = chPlot(cIdx);
                    rowIdx = find(yDataRange == chan, 1);
                    validIndices = all_indices > 0 & all_indices <= size(inputData,2);
                    for sIdx = 1:numSpikes
                        if all(validIndices(sIdx,:))
                            Waveforms(sIdx, cIdx, :) = ampScale * normIn(rowIdx, all_indices(sIdx,:));
                        end
                    end
                end

                if plotMean
                    Waveforms = mean(Waveforms, 1, 'omitnan');
                end

                for cIdx = 1:numChPlots
                    chan = chPlot(cIdx);
                    t = waveformXaxis + (uniqPlots) * xNorm_adj(chan);
                    t = t + 1.1 * plotCount * range(t);
                    minAx = min(minAx, min(t));
                    maxAx = max(maxAx, max(t));
                    wf = squeeze(Waveforms(:, cIdx, :)) + yNorm_adj(chan);

                    min_YAx = min(min_YAx, yNorm_adj(chan)-1);
                    max_YAx = max(max_YAx, yNorm_adj(chan)+1);
                    if chan == channelList(k)
                        plot(plotAxis, t, wf', 'Color', [c, alphaLevel], 'LineWidth', lineWidth*2);
                    else
                        plot(plotAxis, t, wf', 'Color', [c, alphaLevel], 'LineWidth', lineWidth);
                    end
                    plot(plotAxis, [t(round(end/2)) t(round(end/2))], ...
                        [yNorm_adj(chan)-1 yNorm_adj(chan)+1], ...
                        'Color', 1-figColor, 'LineWidth', lineWidth/10, 'LineStyle', '-');
                end
                plotCount = plotCount + 1;
            end
        end

        chanPlotList = unique(chanPlotList);

        if isinf(maxAx),   maxAx   = 1; end
        if isinf(minAx),   minAx   = 0; end
        if isinf(max_YAx), max_YAx = 1; end
        if isinf(min_YAx), min_YAx = 0; end

        for i = 1:numel(chanPlotList)
            chan = chanPlotList(i);
            xPos = ((uniqPlots) * xNorm_adj(chan) + waveformXaxis(1)) - 0.05*(maxAx-minAx);
            yPos = yNorm_adj(chan);
            text(plotAxis, xPos, yPos, chanLabel(chan), ...
                'Color', 1-figColor, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                'FontSize', 9, 'FontWeight', 'bold', 'BackgroundColor', [figColor 0.7]);
        end

        hold(plotAxis,'off');

        axis(plotAxis, 'tight', [minAx - .1*(maxAx-minAx), maxAx + .05*(maxAx-minAx), ...
            min_YAx - .1*(max_YAx-min_YAx), max_YAx + .05*(max_YAx-min_YAx)]);

        selected_order_last = selected_order;
    toc
    end

    function plotCCG(panelOption)
        disp('--- Plot CCG called ---');


        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'CCG')
                delete(allchild(axChannels));
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uipanel(bottomRightGrid,'BorderType','none',...
                    'BackgroundColor',figColor);
                axChannels.Layout.Column = [1 5];
                axChannels.Layout.Row    = [1 5];
                applyColorScheme(axChannels, figColor);
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = panelOption;
        end


        plotList = unique(groupList(selected_order));
        if strcmp(plotType,'Multiple')
        numCCG = min(numel(plotList),8);
        else
        numCCG = min(numel(plotList),10);
        end
        ccgStep = 1;
        if numCCG
            plotingGroups = plotList(numel(plotList)-numCCG+1:end);

            tileforCCG = tiledlayout(plotAxis, numCCG, numCCG, 'Padding', 'tight', 'TileSpacing', "tight");


            for iCCG = 1:numCCG
                for jCCG = 1:numCCG

                    groupID_i = (plotingGroups(iCCG));
                    groupID_j = (plotingGroups(jCCG));
                    if groupID_i>0 && groupID_j>0
                        c = getGroupColor(groupList(groupID_i));

                        if isempty(xcorr_vals{groupID_i,groupID_j}) || mergedFlag(groupID_i) || mergedFlag(groupID_j)
                            spike_idx_i = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_i);
                            spike_idx_j = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_j);
                            xcorr_vals{groupID_i,groupID_j} = binary_xcorr(spike_idx_i, spike_idx_j, num_Samples, cfg.samplingFrequency, 1000, ccgLag, ccgStep, smoothN);
                            if iCCG == jCCG
                                xcorr_vals{groupID_i,groupID_j}(ccgLag+1-2*smoothN : ccgLag+1+2*smoothN) = 0;
                            end
                        end

                        tiledAX{iCCG,jCCG} = nexttile(tileforCCG);
                        applyColorScheme(tiledAX{iCCG,jCCG}, figColor);
                        bar(tiledAX{iCCG,jCCG},-ccgLag:ccgStep:ccgLag,xcorr_vals{groupID_i,groupID_j},'facecolor',c,'edgecolor',c)
                        if groupID_i == groupID_j
                            title(tiledAX{iCCG,jCCG}, sprintf('G%g',...
                                groupID_i),'color',1-figColor);
                        else
                            title(tiledAX{iCCG,jCCG}, sprintf('G%g-G%g',...
                                groupID_i,groupID_j),'color',1-figColor);
                        end
                        axis(tiledAX{iCCG,jCCG},'tight','off')
                    end
                end
            end
        end


    end

    function plotISI(panelOption)
        disp('--- Plot ISI called ---');

        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'ISI')
                delete(allchild(axChannels));
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uipanel(bottomRightGrid,'BorderType','none',...
                    'BackgroundColor',figColor);
                axChannels.Layout.Column = [1 5];
                axChannels.Layout.Row    = [1 5];
                applyColorScheme(axChannels, figColor);
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = panelOption;
        end

        plotList = unique(groupList(selected_order));
        numISI = min(numel(plotList),10);

        if numISI
            plotingGroups = plotList(numel(plotList)-numISI+1:end);

            tileforISI = tiledlayout(plotAxis, numISI, numISI, 'Padding', 'tight', 'TileSpacing', "tight");

            for iISI = 1:numISI
                for jISI = 1:numISI

                    groupID_i = (plotingGroups(iISI));
                    groupID_j = (plotingGroups(jISI));
                    if groupID_i>0 && groupID_j>0
                        c = getGroupColor(groupList(groupID_i));
                        if isempty(isiCounts{groupID_i,groupID_j}) || mergedFlag(groupID_i) || mergedFlag(groupID_j)
                            spike_idx_i = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_i);
                            spike_idx_j = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_j);
                            if iISI==jISI
                                spike_idx = unique(spike_idx_i);
                            else
                                spike_idx = [unique(spike_idx_i); unique(spike_idx_j)];
                            end

                            [isiCounts{groupID_i,groupID_j}, isiCenters{groupID_i,groupID_j}, isiViolations(groupID_i,groupID_j)] = getISIViolations(spike_idx, cfg.samplingFrequency , thresholdISI);
                        end

                        tiledAX{iISI,jISI} = nexttile(tileforISI);
                        applyColorScheme(tiledAX{iISI,jISI}, figColor);
                        b = bar(tiledAX{iISI,jISI},isiCenters{groupID_i,groupID_j},isiCounts{groupID_i,groupID_j},...
                            'FaceColor','flat','edgecolor','flat');
                        b.CData(1,:) = [1 .1 .1];
                        b.CData(2:end,:) = repmat(c,[length(isiCenters{groupID_i,groupID_j})-1,1]);
                        b.LineWidth = 5 * thresholdISI/1000;

                        if thresholdISI > 1

                            if isiViolations(groupID_i,groupID_j)< 1
                                isiColor = 1-figColor;
                            elseif isiViolations(groupID_i,groupID_j)>= 1 && isiViolations(groupID_i,groupID_j) < 2.5
                                isiColor = [.8 .5 .15];
                            else
                                isiColor = [.8 .15 .15];
                            end

                        else
                            if isiViolations(groupID_i,groupID_j)< .25
                                isiColor = 1-figColor;
                            elseif isiViolations(groupID_i,groupID_j)>= .25 && isiViolations(groupID_i,groupID_j) < 0.5
                                isiColor = [.8 .5 .15];
                            else
                                isiColor = [.8 .15 .15];
                            end
                        end

                        if iISI == jISI
                            title(tiledAX{iISI,jISI}, sprintf('G%g:  %.2f%%',...
                                groupID_i,isiViolations(groupID_i,groupID_j)),'color',isiColor);
                        else
                            title(tiledAX{iISI,jISI}, sprintf('%.2f%%',...
                                isiViolations(groupID_i,groupID_j)),'color',isiColor);
                        end
                        axis(tiledAX{iISI,jISI},'tight','off')
                    end
                end
            end
        end

    end


    function plotDensity(panelOption)
        disp('--- Plot density called ---');

        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'Density')
                delete(allchild(axChannels));
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uipanel(bottomRightGrid, 'BorderType', 'none',...
                    'BackgroundColor', figColor);
                axChannels.Layout.Column = [1 5];
                axChannels.Layout.Row    = [1 5];
                applyColorScheme(axChannels, figColor);
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = panelOption;
        end

        plotList = unique(groupList(selected_order));
        numDen = min(numel(plotList), 10);

        if numDen
            plotingGroups = plotList(numel(plotList)-numDen+1:end);
            tileforDen = tiledlayout(plotAxis, numDen, numDen, 'Padding', 'tight', 'TileSpacing', "tight");

            for iDen = 1:numDen
                for jDen = 1:numDen

                    groupID_i = plotingGroups(iDen);
                    groupID_j = plotingGroups(jDen);
                    if groupID_i > 0 && groupID_j > 0
                        c = getGroupColor(groupList(groupID_i));
                        if isempty(DenCounts{groupID_i, groupID_j}) || mergedFlag(groupID_i) || mergedFlag(groupID_j)
                            spike_idx_i = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_i);
                            spike_idx_j = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_j);
                            if iDen == jDen
                                spike_idx = unique(spike_idx_i);
                            else
                                spike_idx = [unique(spike_idx_i); unique(spike_idx_j)];
                            end
                            [DenCenters{groupID_i, groupID_j}, DenCounts{groupID_i, groupID_j}, presence_ratio(groupID_i, groupID_j)] = ...
                                presenceRatio(spike_idx, cfg.samplingFrequency, trialLength, round(trialLength/(2*ccgLag)), smoothN);
                        end

                        tiledAX{iDen,jDen} = nexttile(tileforDen);
                        applyColorScheme(tiledAX{iDen,jDen}, figColor);
                        set(tiledAX{iDen,jDen}, 'Color', figColor);  % Ensure background color is figColor

                        b = bar(tiledAX{iDen,jDen}, DenCenters{groupID_i, groupID_j}, DenCounts{groupID_i, groupID_j},...
                            'FaceColor', c, 'EdgeColor', c);

                        if presence_ratio(groupID_i, groupID_j) > 0.9
                            DenColor = 1 - figColor;
                        elseif presence_ratio(groupID_i, groupID_j) >= 0.5 && presence_ratio(groupID_i, groupID_j) < 0.9
                            DenColor = [.8 .5 .15];
                        else
                            DenColor = [.8 .15 .15];
                        end

                        if iDen == jDen
                            title(tiledAX{iDen,jDen}, sprintf('G%g: %.2f', groupID_i, presence_ratio(groupID_i, groupID_j)), 'color', DenColor);
                        else
                            title(tiledAX{iDen,jDen}, sprintf('%.2f', presence_ratio(groupID_i, groupID_j)), 'color', DenColor);
                        end

                        % For diagonal plots, add the interactive red vertical line
                        if iDen == jDen
                            hold(tiledAX{iDen,jDen}, 'on');
                            yl = get(tiledAX{iDen,jDen}, 'YLim');
                            set(tiledAX{iDen,jDen}, 'XLimMode', 'manual', 'YLimMode', 'manual');
                            default_x1 = stable_length(groupID_i,1) * DenCenters{groupID_i, groupID_j}(end) / trialLength;
                            default_x2 = stable_length(groupID_i,2) * DenCenters{groupID_i, groupID_j}(end) / trialLength;
                            pos1 = [default_x1, yl(1); default_x1, yl(2)];
                            pos2 = [default_x2, yl(1); default_x2, yl(2)];
                            hGreenLine = imline(tiledAX{iDen,jDen}, pos1);
                            hRedLine = imline(tiledAX{iDen,jDen}, pos2);
                            setColor(hGreenLine, [.1 .9 .1]);
                            setColor(hRedLine, [.9 .1 .1]);
                            hLine = findobj(hRedLine, 'Type', 'line');
                            set(hLine, 'LineWidth', 2);
                            setPositionConstraintFcn(hGreenLine, @(pos) [pos(1,1) yl(1); pos(1,1) yl(2)]);
                            addNewPositionCallback(hGreenLine, @(p) greenLineCallback(p, groupID_i));
                            setPositionConstraintFcn(hRedLine, @(pos) [pos(1,1) yl(1); pos(1,1) yl(2)]);
                            addNewPositionCallback(hRedLine, @(p) redLineCallback(p, groupID_i));
                            hold(tiledAX{iDen,jDen}, 'off');
                            xlim(tiledAX{iDen,jDen},[min(DenCenters{groupID_i, groupID_j})-range(DenCenters{groupID_i, groupID_j})/20 max(DenCenters{groupID_i, groupID_j})+range(DenCenters{groupID_i, groupID_j})/20]);
                        else
                            xlim(tiledAX{iDen,jDen},[min(DenCenters{groupID_i, groupID_j})-range(DenCenters{groupID_i, groupID_j})/20 max(DenCenters{groupID_i, groupID_j})+range(DenCenters{groupID_i, groupID_j})/20]);
                            axis(tiledAX{iDen,jDen},'off')
                        end
                        set(tiledAX{iDen,jDen}, 'Color', figColor);

                    end
                end
            end
        end


        function greenLineCallback(pos, group_idx)
            new_x = pos(1,1);
            fprintf('New x position for group %d: %f\n', group_idx, new_x);
            if new_x > DenCenters{group_idx, group_idx}(end)
                new_x = DenCenters{group_idx, group_idx}(end);
            end
            stable_length(group_idx,1) = trialLength * (new_x-DenCenters{group_idx, group_idx}(1)) / DenCenters{group_idx, group_idx}(end);

            if stable_length(group_idx,1) > trialLength
                stable_length(group_idx,1) = trialLength;
            elseif stable_length(group_idx,1) < 1
                stable_length(group_idx,1) = 1;
            end
        end

        function redLineCallback(pos, group_idx)
            new_x = pos(1,1);
            fprintf('New x position for group %d: %f\n', group_idx, new_x);
            if new_x > DenCenters{group_idx, group_idx}(end)
                new_x = DenCenters{group_idx, group_idx}(end);
            end
            stable_length(group_idx,2) = trialLength * new_x / DenCenters{group_idx, group_idx}(end);
            if stable_length(group_idx,2) > trialLength
                stable_length(group_idx,2) = trialLength;
            elseif stable_length(group_idx,2) < 1
                stable_length(group_idx,2) = 1;
            end
        end

    end



    function plotChannelPCA(panelOption)
        disp('--- Plot waveforms called ---');

        if isempty(mappedData)
            if isfield(cfg, 'inputFolder') && isfield(cfg, 'outputFolder')
                mappedData =  map_input_file(cfg.fullFilePath, cfg);
                trialLength = size(mappedData.Data.data,2)/cfg.samplingFrequency;
                xSlider.Limits = [0 trialLength];
                roundTicks = 25:25:trialLength;
                nearestRoundTick = nearest(roundTicks,floor(trialLength/10));
                tickInterval = roundTicks(nearestRoundTick);
                xSlider.MajorTicks = 0:tickInterval:trialLength;
                disp('--- Data loaded ---');
            else
                disp('--- Data loading failed ---');
            end
        end

        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'PCs')
                cla(axChannels,'reset');
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uiaxes(bottomRightGrid);
                axChannels.Layout.Row = [1 3];
                axChannels.Layout.Column = [2 3];
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = uiaxes(panelOption, 'Units', 'normalized', 'Position', [0 0 1 1]);
        end

        if ~strcmp(lblProgressPros.Text,'Complete!')
            uialert(parentFig,'Similarity Processing needs to be performed first for this plot type!','Error');
            return,
        end

        view(plotAxis, 2);
        title(plotAxis,'PCs on Channels','Color','white');
        set(plotAxis,'tickdir','out');
        plotAxis.Color   = figColor;
        plotAxis.XColor  = 1-figColor;
        plotAxis.YColor  = 1-figColor;
        box(plotAxis,'off');
        axis(plotAxis,'off');
        if ~isempty(selected_order) & ~isequal(selected_order, selected_order_last)
            minSliderVal = min(channelPlot(selected_order(end),:))-1;
            maxSliderVal = max(channelPlot(selected_order(end),:))+1;
            if ~isempty(minSliderVal)
                if minSliderVal(end) < yWindowStart
                    yWindowStart = max(floor(minSliderVal(end)),1);
                    yWindowEnd = min(yWindowStart + str2double(yStepDropdown.Value), numChannels);
                    ySlider.Value = yWindowStart;
                    drawnow limitrate;
                elseif  maxSliderVal(end) > yWindowEnd
                    yWindowEnd = min(ceil(maxSliderVal(end)), numChannels);
                    yWindowStart = max(yWindowEnd-str2double(yStepDropdown.Value),1);
                    ySlider.Value = yWindowStart;
                    drawnow limitrate;
                end
            end
        end


        xDataRange = round(xWindowStart * cfg.samplingFrequency) + 1 : min(round(xWindowEnd * cfg.samplingFrequency), num_Samples);
        if length(xDataRange) > 10 * cfg.samplingFrequency
            totalBlocks = floor(length(xDataRange) / cfg.samplingFrequency);
            selectedBlocks = sort(randsample(totalBlocks, 10));
            subsampledRange = [];
            for i = 1:numel(selectedBlocks)
                idx = (selectedBlocks(i)-1) * cfg.samplingFrequency + (1:cfg.samplingFrequency);
                subsampledRange = [subsampledRange, xDataRange(idx)];
            end
            xDataRange = subsampledRange(subsampledRange<=xDataRange(end) & subsampledRange>=xDataRange(1));
        end

        yDataRange = max(1, yWindowStart) : min(yWindowEnd, numChannels);
        yDataRange = yDataRange(chan_wave_inclusion(yDataRange));

        mask = false(num_Samples, 1);
        mask(xDataRange) = true;
        convMask = conv(double(mask), ones(2*numSpikeSamples+1,1), 'same');
        inRange = sortedRes.spike_idx > numSpikeSamples & sortedRes.spike_idx <= num_Samples-numSpikeSamples & ...
            convMask(sortedRes.spike_idx) >= (2*numSpikeSamples+1);

        unifiedLabels = sortedRes.unifiedLabels(inRange);
        spike_idx = sortedRes.spike_idx(inRange);
        [~, spike_idx] = ismember(spike_idx, xDataRange);

        persistent cachedInputData lastXDataRange lastYDataRange
        if isempty(cachedInputData) || ~isequal(lastXDataRange, xDataRange) || ~isequal(lastYDataRange, yDataRange)
            inputData = double(mappedData.Data.data(channel_mapping(yDataRange), xDataRange));
            lastXDataRange = xDataRange;
            lastYDataRange = yDataRange;
            cachedInputData = inputData;
        else
            inputData = cachedInputData;
        end

        if strcmp(plotFilterType, 'filtered')
            inputData = bandpass_filter_GUI(inputData', cfg.bandpass, cfg.samplingFrequency)';
        end



        chanPlotList = [];
        minAx=inf;
        maxAx=-inf;

        max_YAx = -inf;
        min_YAx = inf;

        hold(plotAxis,'on')

        if any(spike_idx) && ~isempty(selected_order)
            for k = selected_order
                if chkBoxVis(k)
                    plot_spike_idx = spike_idx(unifiedLabels == groupList(k));
                    if ~isempty(plot_spike_idx)

                        chPlot = channelPlot(k,:);
                        chPlot = chPlot(ismember(chPlot,yDataRange));
                        chanPlotList = [chanPlotList, chPlot];
                        c = getGroupColor(groupList(k));
                        spike_indices = plot_spike_idx(:);
                        numSpikes = numel(spike_indices);
                        if ~plotMean
                            if numSpikes > 4*numWaveforms
                                selected_ids = randperm(numSpikes,4*numWaveforms);
                                spike_indices = spike_indices(selected_ids);
                                numSpikes = numel(spike_indices);
                            end
                            all_indices = spike_indices + spike_Xaxis ;
                            numChannelsPlots = numel(chPlot);
                            pcScores = zeros(numSpikes, numChannelsPlots, 2);

                            for cIdx = 1:numChannelsPlots
                                chan = chPlot(cIdx);
                                sampled_waveforms = reshape(inputData(yDataRange==chan, all_indices), [numSpikes, numel(spike_Xaxis)]);
                                waveform_scores = (sampled_waveforms - PCA{chan,1}.mu) * PCA{chan,1}.coeff;
                                pcScores(:,cIdx,:) = 3 * waveform_scores(:,1:2)./PCA{chan,1}.max;
                            end
                        else
                            all_indices = spike_indices + spike_Xaxis ;
                            numChannelsPlots = numel(chPlot);
                            pcScores = nan(numSpikes, numChannelsPlots, 2);
                            for cIdx = 1:numChannelsPlots
                                chan = chPlot(cIdx);
                                sampled_waveforms = reshape(inputData(yDataRange==chan, all_indices), [numSpikes, numel(spike_Xaxis)]);
                                waveform_scores = (sampled_waveforms - PCA{chan,1}.mu) * PCA{chan,1}.coeff;
                                pcScores(:,cIdx,:) = 3 *waveform_scores(:,1:2)./PCA{chan,1}.max;
                            end
                            pcScores = mean(pcScores,1,'omitmissing');
                        end

                        for cIdx = 1:numChannelsPlots
                            chan = chPlot(cIdx);
                            xOff = xNorm_adj(chan);
                            yOff = yNorm_adj(chan);
                            xScatter = squeeze(pcScores(:, cIdx, 1)) + xOff;
                            yScatter = squeeze(pcScores(:, cIdx, 2))/2 + yOff;

                            minAx = min(minAx,min(xScatter));
                            maxAx = max(maxAx,max(xScatter));

                            min_YAx = min(min_YAx,yNorm_adj(chan)-1);
                            max_YAx = max(max_YAx,yNorm_adj(chan)+1);

                            if chan == channelList(k)
                                scatter(plotAxis,xScatter,yScatter,lineWidth*10,'filled','MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerFaceAlpha',alphaLevel, 'MarkerEdgeAlpha', alphaLevel)
                            else
                                scatter(plotAxis,xScatter,yScatter,lineWidth*5,'filled','MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerFaceAlpha',alphaLevel, 'MarkerEdgeAlpha', alphaLevel)
                            end

                        end

                    end
                end
            end
        end

        chanPlotList = unique(chanPlotList);


        hold(plotAxis,'off')


        maxAx(isinf(maxAx)) =1;
        minAx(isinf(minAx)) =0;

        min_YAx(isinf(min_YAx)) =0;
        max_YAx(isinf(max_YAx)) =1;
        %
        for i = 1:length(chanPlotList)
            chan = chanPlotList(i);
            xPos = minAx + (xNorm_adj(chan)) - 0.05*(maxAx-minAx);
            yPos = yNorm_adj(chan);

            text(plotAxis, xPos, yPos, chanLabel(chan), 'color', 1 - figColor, ...
                'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
                'FontSize', 9, 'FontWeight', 'bold', 'BackgroundColor', [figColor 0.7]);
        end


        axis(plotAxis, 'tight',[minAx-.2*(maxAx-minAx), maxAx + .05*(maxAx-minAx), ...
            min_YAx-.1*(max_YAx-min_YAx), max_YAx + .05*(max_YAx-min_YAx)]);


        disp(['Time range: ', num2str(xWindowStart), ' to ', num2str(xWindowEnd)]);
        disp(['Channel range: ', num2str(yDataRange(1)), ' to ', num2str(yDataRange(end))]);
        selected_order_last = selected_order;

    end


    function plotAmplitude(panelOption)

        disp('--- Plot Features called ---');

        if ~strcmp(plotType,'Multiple')
            if strcmp(lastPlotted,'Amplitude')
                cla(axChannels,'reset');
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uiaxes(bottomRightGrid,'Color',figColor);
                axChannels.Layout.Row = [1 3];
                axChannels.Layout.Column = [2 3];
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = uiaxes(panelOption, 'Units', 'normalized', 'Position', [0 0 1 1],'Color',figColor);
        end

        if ~strcmp(plotType,'Multiple')
            title(plotAxis,'Spike Amplitude','Color','white');
        end
        xlabel(plotAxis, 'Time (s)');
        ylabel(plotAxis, 'Amp');
        set(plotAxis,'tickdir','out');
        plotAxis.Color   = figColor;
        plotAxis.XColor  = 1-figColor;
        plotAxis.YColor  = 1-figColor;



        xDataRange = round(1 : num_Samples);


        inRange_idx = find(sortedRes.spike_idx > xDataRange(1) & sortedRes.spike_idx < xDataRange(end));

        unifiedLabels = sortedRes.unifiedLabels(inRange_idx);
        spike_idx = sortedRes.spike_idx(inRange_idx);
        spike_amplitude = sortedRes.amplitude(inRange_idx);
        norm_spike_amplitude = sign(spike_amplitude).*sqrt(abs(spike_amplitude/max(abs(spike_amplitude),[],"all")));

        hold(plotAxis,'on')

        if any(spike_idx) && ~isempty(selected_order)
            for k = selected_order
                classPolarity = double(mainPolarity(k));
                classPolarity(classPolarity==1) = -1;
                classPolarity(classPolarity==0) = 1;
                scatterIdx = find(unifiedLabels == groupList(k));
                if length(scatterIdx) > 5E4
                    scatterIdx = scatterIdx(randperm(length(scatterIdx), 5E4));
                end
                plot_spike_idx = spike_idx(scatterIdx);
                if ~isempty(plot_spike_idx)
                    c = getGroupColor(groupList(k));
                    scatter(plotAxis, plot_spike_idx/cfg.samplingFrequency,classPolarity * abs(norm_spike_amplitude(scatterIdx)),...
                        lineWidth * 2.5,'markerfacecolor',c,'markeredgecolor','none')
                end
            end
        end

        hold(plotAxis,'off')


        axis(plotAxis,'tight',[0 num_Samples/cfg.samplingFrequency -1.1 1.1 ]);
        yTicksIntervals = 0.5;
        set(plotAxis,'ytick',[-1:yTicksIntervals:1]);
    end

    function updatePlotType(newVal)
        plotType = newVal;
        disp(['Plot type set to ', plotType]);
        plotOption();
    end

    function updatePlotFilterType(newVal)
        plotFilterType = newVal;
        disp(['Plot filter type set to ', plotFilterType]);
        plotOption();
    end

% CCG plot functions
    function updateISIThr(ddVal)
        thresholdISI = str2double(ddVal);
        preprocessSortedData();
        refreshLabels();
        disp(['Plot CCG Step set to ', ddVal]);
        refreshISI();
        setISIIncVal(uifISI.Value);
        plotOption();
    end

    function updateCCGLag(ddVal)
        ccgLag = str2double(ddVal);
        disp(['Plot CCG Lag set to ', ddVal]);
        refreshCCG();
        refreshDensity();
        plotOption();
    end

    function updateCCGSmoothN(ddVal)
        smoothN = str2double(ddVal);
        disp(['Plot CCG smoothness set to ', ddVal]);
        refreshCCG();
        refreshDensity();
        plotOption();
    end

    function refreshCCG()
        xcorr_vals = cell(numGroups,numGroups);
    end

    function refreshISI()
        isiCounts = cell(numGroups,numGroups);
        isiCenters = cell(numGroups,numGroups);
        isiViolations = zeros(numGroups,numGroups);
    end

    function refreshDensity()
        DenCenters = cell(numGroups,numGroups);
        DenCounts = cell(numGroups,numGroups);
        presence_ratio = zeros(numGroups,numGroups);
    end




    function updateChType(sVal)
        if strcmp(sVal,'Config')
            numChannelPlot = cfg.num_channel_extract;
        else
            numChannelPlot = str2double(sVal);
        end
        channelPlot = channelList + [-numChannelPlot:numChannelPlot];
        plotOption();
    end

    function updateChannelPlot()
        channelPlot = channelList + [-numChannelPlot:numChannelPlot];
        plotOption();
    end

    function updateNumWaveforms(sldVal)
        numWaveforms = str2double(sldVal);
        disp(['Number of wavforms set to ', num2str(numWaveforms)]);
        plotOption();
    end

    function updateWaveDur(sldVal)

        if strcmp(sldVal,'Config')
            halfSpikeWaveDur = cfg.clusteringSpikeDuration/2;
        else
            halfSpikeWaveDur = str2double(sldVal)/2;
            numSpikeSamples = round(halfSpikeWaveDur * cfg.samplingFrequency /1000);
            spike_Xaxis = -numSpikeSamples:numSpikeSamples;
            waveformXaxis = 1000 *spike_Xaxis/cfg.samplingFrequency;
            sliderReAlign.Limits = [-halfSpikeWaveDur halfSpikeWaveDur];
            sliderReAlign.MajorTicks = [-halfSpikeWaveDur 0 halfSpikeWaveDur];
        end

        disp(['Number of wavforms set to ', num2str(halfSpikeWaveDur)]);
        normalizeTimeAmp();
        plotOption();
    end

    function updateWavetype(ddVal)
        if strcmp(ddVal,'Mean')
            plotMean = 1;
        else
            plotMean = 0;
        end
        disp(['Plot ', ddVal,'Waveforms']);
        plotOption();
    end


    function preprocessSortedData()
        % Vectorized preprocessing
        for i = 1:numGroups
            groupID_i = groupList(i);
            spike_idx_i = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_i);
            [~, ~, preprocessed.isiViolation(i,1)] = getISIViolations(spike_idx_i, cfg.samplingFrequency, thresholdISI);
            preprocessed.firingRate(i,1) = numel(spike_idx_i)./trialLength;
        end
        preprocessed.logFiringRate = log(preprocessed.firingRate);
    end




    function onSimilarityEstimation()

        lblProgressPros.Text = 'Processing...';
        drawnow;
        groupSimScore = cell(numGroups, numGroups);
        ampVar = cell(numGroups, numGroups);
        if isempty(mappedData)
            if isfield(cfg, 'inputFolder') && isfield(cfg, 'outputFolder')
                mappedData =  map_input_file(cfg.fullFilePath, cfg);
                trialLength = size(mappedData.Data.data,2)/cfg.samplingFrequency;
                xSlider.Limits = [0 trialLength];
                roundTicks = 25:25:trialLength;
                nearestRoundTick = nearest(roundTicks,floor(trialLength/10));
                tickInterval = roundTicks(nearestRoundTick);
                xSlider.MajorTicks = 0:tickInterval:trialLength;
                disp('--- Data loaded ---');
            else
                disp('--- Data loading failed ---');
            end
        end


        group_idx = cell(numChannels, 1);
        PCA = cell(numChannels, 1);
        lastGroup = 0;

        for iChan = 1:numChannels

            if chan_wave_inclusion(iChan)

                channels = find(channelList<= iChan + numChannelPlot & channelList >= iChan - numChannelPlot);
                sampled_waveforms = [];
                mean_sample_waveforms = [];

                if ~isempty(channels)
                    for jChan = 1:length(channels)
                        spike_idx_i = sortedRes.spike_idx(sortedRes.unifiedLabels == groupList(channels(jChan)));
                        numSamples = min(length(spike_idx_i),100);
                        sampled_idx = datasample(spike_idx_i, numSamples,'Replace',false);
                        waveform_idx = sampled_idx(:) + [-numSpikeSamples:numSpikeSamples];
                        inputData = double(mappedData.Data.data(channel_mapping(iChan), waveform_idx'));
                        class_waveform = reshape(inputData, 2*numSpikeSamples +1, numSamples)';
                        sampled_waveforms = [sampled_waveforms; class_waveform];
                        mean_sample_waveforms = [mean_sample_waveforms; mean(class_waveform,1,'omitmissing')];
                        group_idx{iChan,1} = [group_idx{iChan,1}; repmat(channels(jChan),[numSamples,1])];
                    end
                end
                [PCA{iChan,1}.coeff,PCA{iChan,1}.score,~,~,~,PCA{iChan,1}.mu] = ...
                    pca(sampled_waveforms,'Algorithm','eig',  'Centered','on','NumComponents',2);
                PCA{iChan,1}.max = max(abs(PCA{iChan,1}.score));

                for g1 = 1:length(channels)
                    for g2 = 1:length(channels)
                        if g2 < g1 && abs(ylocs(channelList(channels(g2))) - ylocs(channelList(channels(g1)))) <= yDistThr

                            mw1 = mean_sample_waveforms(g1,:);
                            mw2 = mean_sample_waveforms(g2,:);
                            g1_idx = group_idx{iChan,1} == channels(g1);
                            g2_idx = group_idx{iChan,1} == channels(g2);
                            if sum(g1_idx)>1 && sum(g2_idx)>1
                                pc_g1 = PCA{iChan,1}.score(g1_idx,:);
                                pc_g2 = PCA{iChan,1}.score(g2_idx,:);
                                L1 = size(pc_g1,1);
                                L2 = size(pc_g2,1);
                                logicalGroups = [true(L1,1); false(L2,1)];
                                switch distanceEstType
                                    case 'KL-div'
                                        similarityScore = relativeEntropy([pc_g1; pc_g2],logicalGroups);
                                    case 'Bhattacharyya'
                                        similarityScore = bhattacharyyaDistance([pc_g1; pc_g2],logicalGroups);
                                    case 'XCorr'
                                        [similarityScore, ~] = max_half_corr(mw1, mw2, 1, 2*numSpikeSamples+1, round(numSpikeSamples/2),0);
                                end
                            else
                                similarityScore = nan;
                            end
                            maxAmp = max(abs([mw1; mw2]),[],2);
                            maxAmpP = abs(max(([mw1; mw2]),[],2));
                            maxAmpN = abs(min(([mw1; mw2]),[],2));
                            ampDiffP = (maxAmpP(1,:) - maxAmpP(2,:))./(max([maxAmp(1,:) , maxAmp(2,:)],[],'all'));
                            ampDiffN = (maxAmpN(1,:) - maxAmpN(2,:))./(max([maxAmp(1,:) , maxAmp(2,:)],[],'all'));
                            ampDiff = (abs(ampDiffP) + abs(ampDiffN))/2;
                            ampDiff(ampDiff==0) = nan;
                            ampVar{channels(g1), channels(g2)} = ...
                                [ampVar{channels(g1), channels(g2)}; ampDiff];
                            groupSimScore{channels(g1), channels(g2)} = ...
                                [groupSimScore{channels(g1), channels(g2)}; similarityScore ];

                        end

                    end                    
                end

            end
            
            lblProgressPros.Text = sprintf('Processing %3.0d',round(100*iChan/numChannels));
            drawnow;
        end


        if strcmp(distanceEstMeasureType,'median')
            disimlarityScore = cellfun(@(x) median(x,'all','omitmissing'),groupSimScore);
            ampSimilarity = cellfun(@(x) median(x,'all','omitmissing'),ampVar);
        else
            disimlarityScore = cellfun(@(x) mean(x,'all','omitmissing'),groupSimScore);
            ampSimilarity = cellfun(@(x) mean(x,'all','omitmissing'),ampVar);
        end

        lblProgressPros.Text = 'Complete!';
        drawnow;
        plotOption();
    end

end



function colorMap = generate_cluster_colors(N)

colorMap = zeros(N,3);
goldenRatio = 0.618;
h = 0;
for i = 1:N
    h = mod(h + goldenRatio, 1);
    s = 0.9;
    v = 0.9;
    rgb_color = hsv2rgb([h, s ,v]);
    brightness = sum(rgb_color .* [0.299, 0.587, 0.114]);
    while brightness < 0.4
        s = rand(1) ;
        rgb_color = hsv2rgb([h, s ,v]);
        brightness = sum(rgb_color .* [0.299, 0.587, 0.114]);
    end
    colorMap(i,:) = rgb_color;
end
end



function tc = bestTextColorFor(bg)

brightness = sum(bg .* [0.299, 0.587, 0.114]);
if brightness < 0.5
    tc = [1 1 1];
else
    tc = [0 0 0];
end
end

