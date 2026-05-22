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

% Cache of waveform line handles + the unscaled data behind them so that
% scale / alpha / line-width sliders can update the existing plot in
% place instead of triggering a full plotWaveforms redraw. The cache
% is populated in plotWaveforms (single-plot mode) and consumed by the
% updateScale / updateAlphaLevel / updateLineWidth fast paths.
wfPlotCache = struct('lines', {{}}, 'unscaledY', {{}}, ...
                     'yOff', [], 'isMain', logical([]), ...
                     'rgb', zeros(0,3), 'axis', []);

% Most recent uitable handle from plotTable. paintSelectedTableRows
% repaints per-row background colours on this handle whenever the
% selection changes (table-click toggle, walker step, etc.) so the
% user can see at a glance which row is currently selected.
tableHandle = [];

% Standalone "curation table" window. Lives in a separate uifigure
% so a per-unit selection doesn't clear / re-render it (which is
% what used to happen when Table was a plot type). Click "Open
% table" to launch; the window persists across selection changes,
% preserves its column sort, and can be refreshed on demand.
tableFigure       = [];   % uifigure handle, or empty when closed
tableContentPanel = [];   % uipanel inside the figure that hosts the uitable

% Initialize sparse matrices/cells for large data structures to improve memory usage
xcorr_vals = cell(numGroups,numGroups);
isiCounts = cell(numGroups,numGroups);
isiCenters = cell(numGroups,numGroups);
isiViolations = zeros(numGroups,numGroups);
unitIsolation = repmat({'NA'},[numGroups,1]);
% Free-text note per unit. Lives in memory while the GUI is open and
% gets serialised into curated_sample.mat / curated_metrics.csv on
% Save. A user can type "good FS interneuron, layer 5, drifted"
% etc.
unitNotes     = repmat({''}, [numGroups, 1]);

% Split-button options. The dropdown / edit field next to the Split
% button mirror these closure variables; onSplit reads them whenever
% the user triggers a split.
splitMethod      = 'K-means';
splitNumClusters = 2;
DenCenters = cell(numGroups,numGroups);
DenCounts = cell(numGroups,numGroups);
presence_ratio = zeros(numGroups,numGroups);
numChannels = min(length(chan_wave_inclusion),length(channel_mapping));
selectUnits = groupList;
% Unit# navigation walks the shown (nonzero-rate) groups by default; a
% channel inclusion/exclusion filter overrides it with its own set.
chanFilterActive = false;
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
displayedGroups = (1:numGroups)';
totalPages = max(1, ceil(numel(displayedGroups)/PAGE_SIZE));

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

% --- Resume the last saved curation, if any ----------------------
% The Save button writes RES_Sorted/curated_sample.mat with a
% .session sub-struct that captures every per-unit and per-spike
% array we need. On launch we check for that file and apply the
% session on top of the freshly-loaded raw sort data; the
% original* vars above are left as the true pre-curation state
% so Reset / Undo still walk back to the sorter output.
%
% The restore is best-effort: a missing, corrupt, or dimension-
% mismatched file is silently ignored (with a console hint) and
% the GUI starts on the raw data. A fresh dataset that hasn't
% been saved yet just hits the "no file" branch.
try
    curatedSampleFile = fullfile(cfg.outputFolder, 'RES_Sorted', 'curated_sample.mat');
    if exist(curatedSampleFile, 'file')
        loaded = load(curatedSampleFile, 'curatedSamples');
        if isfield(loaded, 'curatedSamples') && ...
                isfield(loaded.curatedSamples, 'session')
            sess = loaded.curatedSamples.session;
            nSpikeFile  = numel(sortedRes.unifiedLabels);
            nSpikeSaved = 0;
            if isfield(sess,'spikeLabels')
                nSpikeSaved = numel(sess.spikeLabels);
            end
            nUnitSaved = 0;
            if isfield(sess,'groupList')
                nUnitSaved = numel(sess.groupList);
            end

            % We can apply per-unit state when nUnitSaved >= numGroups.
            % nUnitSaved > numGroups happens when the saved session
            % had Split-created units; we extend the per-unit arrays
            % with the same NaN-original markers onSplit uses.
            if nUnitSaved >= numGroups
                if nUnitSaved > numGroups
                    extra = nUnitSaved - numGroups;
                    wfShape = size(originalSampleWaveform);
                    originalGroupList(end+1:nUnitSaved,1)   = NaN;
                    originalChannelList(end+1:nUnitSaved,1) = NaN;
                    originalSampleWaveform(end+1:nUnitSaved,:,:) = ...
                        nan(extra, wfShape(2), wfShape(3));
                    sampleWaveform(end+1:nUnitSaved,:,:) = ...
                        nan(extra, wfShape(2), wfShape(3));
                    detectblity(end+1:nUnitSaved,1)  = 0;
                    mainPolarity(end+1:nUnitSaved,1) = false;
                    sidePolarity(end+1:nUnitSaved,1) = false;
                    mergedFlag(end+1:nUnitSaved,1)   = false;
                    unitIsolation(end+1:nUnitSaved,1) = {'NA'};
                    unitNotes(end+1:nUnitSaved,1)    = {''};
                    stable_length(end+1:nUnitSaved,:) = ...
                        repmat([1, trialLength], extra, 1);
                    chkBoxSelect(end+1:nUnitSaved,1) = false;
                    chkBoxVis(end+1:nUnitSaved,1)    = false;
                    selectUnits = [selectUnits; nan(extra,1)];
                    numGroups   = nUnitSaved;
                    displayedGroups = (1:numGroups)';
                    totalPages  = max(1, ceil(numel(displayedGroups) / PAGE_SIZE));
                    % Per-pair caches grow on demand so just size them now.
                    xcorr_vals(numGroups,numGroups)   = {[]};
                    isiCounts(numGroups,numGroups)    = {[]};
                    isiCenters(numGroups,numGroups)   = {[]};
                    isiViolations(numGroups,numGroups)= 0;
                    DenCenters(numGroups,numGroups)   = {[]};
                    DenCounts(numGroups,numGroups)    = {[]};
                    presence_ratio(numGroups,numGroups)= 0;
                end

                % Apply per-unit arrays. Each isfield + length match
                % so a partial save (older format) doesn't crash.
                if numel(sess.groupList)     == numGroups, groupList     = sess.groupList(:);     end
                if isfield(sess,'channelList')   && numel(sess.channelList)   == numGroups, channelList   = sess.channelList(:);   end
                if isfield(sess,'unitIsolation') && numel(sess.unitIsolation) == numGroups, unitIsolation = sess.unitIsolation(:); end
                if isfield(sess,'stable_length') && size(sess.stable_length,1) == numGroups, stable_length = sess.stable_length;    end
                if isfield(sess,'mainPolarity')  && numel(sess.mainPolarity)  == numGroups, mainPolarity  = sess.mainPolarity(:);  end
                if isfield(sess,'sidePolarity')  && numel(sess.sidePolarity)  == numGroups, sidePolarity  = sess.sidePolarity(:);  end
                if isfield(sess,'mergedFlag')    && numel(sess.mergedFlag)    == numGroups, mergedFlag    = sess.mergedFlag(:);    end
                if isfield(sess,'unitNotes')     && numel(sess.unitNotes)     == numGroups, unitNotes     = sess.unitNotes(:);     end
                if isfield(sess,'detectblity')   && numel(sess.detectblity)   == numGroups, detectblity   = sess.detectblity(:);   end
                % channelPlot is derived from channelList -- recompute.
                channelPlot = channelList + [-numChannelPlot:numChannelPlot];
            else
                fprintf('Saved curation has %d units, current sort has %d -- skipping per-unit restore.\n', nUnitSaved, numGroups);
            end

            % Apply per-spike arrays only when the spike count
            % matches exactly. A re-sort would change spike count,
            % in which case we keep the freshly-loaded sortedRes.
            if nSpikeSaved == nSpikeFile && nSpikeSaved > 0
                sortedRes.unifiedLabels = sess.spikeLabels;
                if isfield(sess,'spikeIdx')      && numel(sess.spikeIdx)      == nSpikeFile, sortedRes.spike_idx  = sess.spikeIdx;      end
                if isfield(sess,'spikeChannels') && numel(sess.spikeChannels) == nSpikeFile, sortedRes.channelNum = sess.spikeChannels; end
            elseif nSpikeSaved > 0
                fprintf('Saved curation has %d spikes, current sort has %d -- skipping spike-level restore.\n', nSpikeSaved, nSpikeFile);
            end
        end
    end
catch ME
    fprintf('Could not restore saved curation: %s\n', ME.message);
end

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

% Interactive cell drag state: gridMultiple is the uigridlayout that
% hosts the multi-plot panels (assigned inside plotMultiple). The
% swap button click pattern is two-step: click cell A's grip, then
% click cell B's grip and they swap. The resize grip enters a
% mouse-tracking mode -- click once to start dragging the cell
% boundaries, click again (or release) to commit.
gridMultiple        = [];
multiSwapPendingK   = 0;
multiResizeDrag     = struct('active',false,'k',0,'startMouse',[0 0], ...
    'startRH',{{}},'startCW',{{}},'rowIdx',0,'colIdx',0, ...
    'cellW',0,'cellH',0,'origMotionFcn',[],'origDownFcn',[]);

parentFig.WindowKeyPressFcn = @keyPressHandler;

% Best-effort persistent user prefs. Loaded here (defaults applied
% later when the relevant widgets exist), saved on Save button
% press AND on figure close. Lives next to the curated outputs
% (cfg.outputFolder/RES_Sorted/kiaSort_curate_prefs.mat) so the
% tuning is per-dataset; the legacy home-directory file is still
% read as a fallback so older installs migrate cleanly the first
% time they hit Save.
projectPrefsDir  = fullfile(cfg.outputFolder, 'RES_Sorted');
projectPrefsFile = fullfile(projectPrefsDir, 'kiaSort_curate_prefs.mat');
legacyPrefsFile  = fullfile(getKiaPrefsDir(), 'kiaSort_curate_prefs.mat');
loadedPrefs = struct();
prefsFile = '';
if exist(projectPrefsFile, 'file')
    prefsFile = projectPrefsFile;
elseif exist(legacyPrefsFile, 'file')
    prefsFile = legacyPrefsFile;
end
if ~isempty(prefsFile)
    try
        loadedPrefs = load(prefsFile);
        if isfield(loadedPrefs, 'prefs')
            loadedPrefs = loadedPrefs.prefs;
        end
    catch
        loadedPrefs = struct();
    end
end

% Wire up an idempotent save-on-close hook. The previous CloseRequestFcn
% is preserved and chained, so we don't accidentally swallow another
% caller's logic.
try
    prevCloseFcn = parentFig.CloseRequestFcn;
    parentFig.CloseRequestFcn = @(src,evt) onParentClose(src, evt, prevCloseFcn);
catch
end

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

leftGrid = uigridlayout(leftPanel, [5 1], ...
    'RowHeight',{'fit','fit','fit','fit','1x'}, ...
    'ColumnWidth',{'1x'}, ...
    'Padding',10, ...
    'RowSpacing',10);
applyColorScheme(leftGrid, figColor);

% top left panel and its buttons
topPanel = uipanel(leftGrid, 'BorderType','none');
topPanel.Layout.Row = 1;
topPanel.Layout.Column = 1;
applyColorScheme(topPanel, figColor);

% topButtonLayout: 2 rows x 9 cols.
%   Cols 1-4 hold the eight action buttons in a 2x4 grid:
%       Row 1: Remove   Merge   Split   Limit
%       Row 2: Isolation Undo   Reset   Notes
%   Cols 5-6 are the compact Split options (Method + #Clust), sized
%   to their label widths.
%   Col 7  Drop Rate label / edit (slimmer than before).
%   Col 8  Realign label (left of slider), col 9 Realign slider --
%   the slider spans both rows so it has room to breathe.
topButtonLayout = uigridlayout(topPanel,[2 9], ...
    'ColumnWidth',{'fit','fit','fit','fit',60,55,'1.25x','fit','2x'}, ...
    'ColumnSpacing',7.5, ...
    'Padding',[0 0 0 0]);
applyColorScheme(topButtonLayout, figColor);

% Row 1 buttons.
btnRemove = uibutton(topButtonLayout,'Text','Remove',...
    'ButtonPushedFcn',@(btn,ev)onRemove(), ...
    'Tooltip','Remove the selected units. Ctrl+D');
btnRemove.Layout.Row = 1;
btnRemove.Layout.Column = 1;
applyColorScheme(btnRemove, figColor);

btnMerge = uibutton(topButtonLayout,'Text','Merge',...
    'ButtonPushedFcn',@(btn,ev)onMerge(), ...
    'Tooltip','Merge selected units into the radio-button choice. Ctrl+M');
btnMerge.Layout.Row = 1;
btnMerge.Layout.Column = 2;
applyColorScheme(btnMerge, figColor);

btnSplit = uibutton(topButtonLayout,'Text','Split',...
    'ButtonPushedFcn',@(btn,ev)onSplit(), ...
    'Tooltip','Split the selected unit(s). Ctrl+T');
btnSplit.Layout.Row = 1;
btnSplit.Layout.Column = 3;
applyColorScheme(btnSplit, figColor);

btnAutoCut = uibutton(topButtonLayout,'Text','Limit',...
    'ButtonPushedFcn',@(btn,ev)onAutoCut(), ...
    'Tooltip','Trim selected units to their stable interval. Ctrl+L');
btnAutoCut.Layout.Row = 1;
btnAutoCut.Layout.Column = 4;
applyColorScheme(btnAutoCut, figColor);

% Row 2 buttons.
btnMUA = uibutton(topButtonLayout,'Text','Isolation',...
    'ButtonPushedFcn',@(btn,ev)onMUA(), ...
    'Tooltip','Cycle isolation label (SUA+/SUA/MUA+/MUA). Ctrl+K');
btnMUA.Layout.Row = 2;
btnMUA.Layout.Column = 1;
applyColorScheme(btnMUA, figColor);

btnUndo = uibutton(topButtonLayout,'Text','Undo',...
    'ButtonPushedFcn',@(btn,ev)onUndo(), ...
    'Tooltip','Revert selected units to their pre-edit state. Ctrl+U');
btnUndo.Layout.Row = 2;
btnUndo.Layout.Column = 2;
applyColorScheme(btnUndo, figColor);

btnReset = uibutton(topButtonLayout,'Text','Reset',...
    'ButtonPushedFcn',@(btn,ev)onReset(), ...
    'Tooltip','Reset all edits. Ctrl+R');
btnReset.Layout.Row = 2;
btnReset.Layout.Column = 3;
applyColorScheme(btnReset, figColor);

btnNotes = uibutton(topButtonLayout,'Text','Notes',...
    'ButtonPushedFcn',@(btn,ev)onEditNotes(), ...
    'Tooltip','Edit the unit''s note. Ctrl+E');
btnNotes.Layout.Row = 2;
btnNotes.Layout.Column = 4;
applyColorScheme(btnNotes, figColor);

% Split options: clustering method + number of sub-clusters. Sit
% right next to the Split button. Sized via fixed pixel widths so
% they take only as much space as their labels need.
lblSplitMethod = uilabel(topButtonLayout,'Text','Method:',...
    'FontWeight','bold','HorizontalAlignment','Center');
lblSplitMethod.Layout.Row = 1;
lblSplitMethod.Layout.Column = 5;
applyColorScheme(lblSplitMethod, figColor);

ddSplitMethod = uidropdown(topButtonLayout, ...
    'Items',{'K-means','GMM','Graph','UMAP+DBSCAN'}, ...
    'Value', splitMethod, ...
    'Tooltip','Clustering method used by Split. UMAP+DBSCAN falls back to PCA+DBSCAN when UMAP is unavailable.', ...
    'ValueChangedFcn', @(dd,ev) setSplitMethod(dd.Value));
ddSplitMethod.Layout.Row = 2;
ddSplitMethod.Layout.Column = 5;
applyColorScheme(ddSplitMethod, figColor);

lblSplitN = uilabel(topButtonLayout,'Text','#Clust:',...
    'FontWeight','bold','HorizontalAlignment','Center');
lblSplitN.Layout.Row = 1;
lblSplitN.Layout.Column = 6;
applyColorScheme(lblSplitN, figColor);

uiSplitN = uieditfield(topButtonLayout,'numeric','Value', splitNumClusters, ...
    'Limits',[2 10], 'LowerLimitInclusive','on', 'UpperLimitInclusive','on', ...
    'RoundFractionalValues','on', ...
    'Tooltip','Number of sub-clusters for Split (2-10).', ...
    'ValueChangedFcn', @(ui,ev) setSplitNumClusters(ui.Value));
uiSplitN.Layout.Row = 2;
uiSplitN.Layout.Column = 6;
applyColorScheme(uiSplitN, figColor);

% Right-hand side: Drop Rate (col 7) + Realign label/slider stacked
% in col 8 (label on top, slider underneath).
lblLimit = uilabel(topButtonLayout,'Text',sprintf('Drop Rate:'),...
    'FontWeight','bold','HorizontalAlignment','Right');
lblLimit.Layout.Row = 1;
lblLimit.Layout.Column = 7;
applyColorScheme(lblLimit, figColor);

uiLimit = uieditfield(topButtonLayout,'numeric','Value',5, ...
    'Tooltip','Fold drop in firing rate Limit needs to trim. Higher = more conservative.',...
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
    'Tooltip','Shift selected units'' spike times by this many ms.',...
    'Limits',[-halfSpikeWaveDur halfSpikeWaveDur]);
sliderReAlign.MajorTicks = [-halfSpikeWaveDur 0 halfSpikeWaveDur];
sliderReAlign.Layout.Row = [1 2];
sliderReAlign.Layout.Column = 9;
applyColorScheme(sliderReAlign, figColor);
sliderReAlign.ValueChangedFcn = @(sld,ev)updateRealignSpikes(sld.Value);

btnSave = uibutton(parentPanel,'Text','Save',...
    'ButtonPushedFcn',@(btn,ev)onSave(), ...
    'Tooltip','Save the curation. Ctrl+S');
btnSave.Layout.Column = 4;
btnSave.Layout.Row = 1;
applyColorScheme(btnSave, figColor);

lblSave = uilabel(parentPanel,'Text',' ');
lblSave.Layout.Column = 5;
lblSave.Layout.Row = 1;
applyColorScheme(lblSave, figColor);

processingPanel = uipanel(leftGrid, 'title','Similarity Processing');
processingPanel.Layout.Row = 3;
processingPanel.Layout.Column = 1;
applyColorScheme(processingPanel, figColor);

% similarity Assessing panel
processingLayout = uigridlayout(processingPanel,[2 8], ...
    'ColumnWidth',{'fit','1x','1x','1x','1x','.5x','.75x','.75x'}, ...
    'ColumnSpacing',10, ...
    'Padding',[5 5 5 5]);
applyColorScheme(processingLayout, figColor);

btnPreProcess = uibutton(processingLayout,'Text','Process',...
    'Tooltip','Compute pairwise similarity.',...
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
    'Tooltip','Similarity metric.',...
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
    'Tooltip','Walker similarity threshold. 0 = no filter.',...
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
    'Tooltip','Max amplitude mismatch for the pair walker. 0 = no filter.',...
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
    'Tooltip','Max y-distance (µm) for Process / pair walker.',...
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

% Auto-curation defaults (configurable via the "Auto curation" panel).
% These are deliberately kept independent of the Similarity Processing
% values so that running Auto does not depend on whatever the user has
% typed in the manual fields. Defaults below are overridden by the
% persistent-prefs file (loadedPrefs, populated near the top of this
% function) when the user has saved a previous session's tuning.
autoDistVal     = 200;     % search radius (um) for similar-pair gating
autoSimVal      = 0.85;    % similarity threshold (XCorr default; KL/Bhatta uses <=0.10)
autoAmpVal      = 0.15;    % amplitude variance ratio gate
autoIsiVal      = 1.0;     % max ISI violation %
autoOverlapFrac = 0.10;    % coincidence overlap fraction (Phase: removal)
autoIsiBudget   = 1.5;     % merged-ISI budget (Phase: merge)
autoPCVal       = 0.30;    % PC-distance gate (Phase: merge)
autoCCGRatio    = 5;       % CCG peak / baseline-median ratio (Phase: cleaning)
% Timing-consistency gate (Phase: cleaning). After CCG peak-at-zero
% triggers, the contaminated side's coincident-spike lags (relative to
% the nearest clean-side spike) must cluster tightly around their
% median, and that median must itself sit near zero. Random coincidences
% spread the lags evenly inside the +-0.5 ms window; a real sorter-
% introduced duplicate produces a sharp sub-ms peak. Only the
% consistent-lag subset is stripped.
autoLagTightMs         = 0.2;
autoLagConsistencyFrac = 0.75;
% Cleaning acts only on similar-quality pairs (|SNR_k - SNR_j| <
% autoCCGSnrDiff); larger gaps belong to overlap removal. The lower
% firing-rate side is the one stripped, capped at autoCCGMaxStripFrac
% of its spikes per pair.
autoCCGSnrDiff         = 0.3;
autoCCGMaxStripFrac    = 0.5;

% Pull saved Auto-panel values from the prefs file so the user gets
% their tuning back on the next launch. Each guard verifies the
% loaded value is finite and inside the editfield's accepted range
% before overwriting the default; a corrupt prefs entry can never
% break the GUI's bring-up.
if isfield(loadedPrefs,'autoDistVal') && isfinite(loadedPrefs.autoDistVal) && ...
        loadedPrefs.autoDistVal >= 0
    autoDistVal = loadedPrefs.autoDistVal;
end
if isfield(loadedPrefs,'autoSimVal') && isfinite(loadedPrefs.autoSimVal) && ...
        loadedPrefs.autoSimVal >= 0 && loadedPrefs.autoSimVal <= 1
    autoSimVal = loadedPrefs.autoSimVal;
end
if isfield(loadedPrefs,'autoAmpVal') && isfinite(loadedPrefs.autoAmpVal) && ...
        loadedPrefs.autoAmpVal >= 0 && loadedPrefs.autoAmpVal <= 1
    autoAmpVal = loadedPrefs.autoAmpVal;
end
if isfield(loadedPrefs,'autoIsiVal') && isfinite(loadedPrefs.autoIsiVal) && ...
        loadedPrefs.autoIsiVal >= 0 && loadedPrefs.autoIsiVal <= 100
    autoIsiVal = loadedPrefs.autoIsiVal;
end
if isfield(loadedPrefs,'autoOverlapFrac') && isfinite(loadedPrefs.autoOverlapFrac) && ...
        loadedPrefs.autoOverlapFrac >= 0 && loadedPrefs.autoOverlapFrac <= 1
    autoOverlapFrac = loadedPrefs.autoOverlapFrac;
end
if isfield(loadedPrefs,'autoIsiBudget') && isfinite(loadedPrefs.autoIsiBudget) && ...
        loadedPrefs.autoIsiBudget >= 1 && loadedPrefs.autoIsiBudget <= 10
    autoIsiBudget = loadedPrefs.autoIsiBudget;
end
if isfield(loadedPrefs,'autoPCVal') && isfinite(loadedPrefs.autoPCVal) && ...
        loadedPrefs.autoPCVal >= 0 && loadedPrefs.autoPCVal <= 5
    autoPCVal = loadedPrefs.autoPCVal;
end
if isfield(loadedPrefs,'autoCCGRatio') && isfinite(loadedPrefs.autoCCGRatio) && ...
        loadedPrefs.autoCCGRatio >= 0.5 && loadedPrefs.autoCCGRatio <= 10
    autoCCGRatio = loadedPrefs.autoCCGRatio;
end
if isfield(loadedPrefs,'autoLagTightMs') && isfinite(loadedPrefs.autoLagTightMs) && ...
        loadedPrefs.autoLagTightMs > 0 && loadedPrefs.autoLagTightMs <= 2
    autoLagTightMs = loadedPrefs.autoLagTightMs;
end
if isfield(loadedPrefs,'autoLagConsistencyFrac') && isfinite(loadedPrefs.autoLagConsistencyFrac) && ...
        loadedPrefs.autoLagConsistencyFrac >= 0 && loadedPrefs.autoLagConsistencyFrac <= 1
    autoLagConsistencyFrac = loadedPrefs.autoLagConsistencyFrac;
end

% Audit trail of pairs/units touched by the most recent Auto run.
% Each pair of vectors is matched 1:1 so the pair walker on the Auto
% panel can step through them. autoPairType (closure variable,
% default 'Merged') selects which set the Next / Prev buttons walk.
autoMergedPrimary    = [];
autoMergedAbsorbed   = [];
autoOverlapDroppedK  = [];   % unit dropped by the >autoOverlapFrac test
autoOverlapTriggerJ  = [];   % the neighbour that caused the drop
autoCCGStripCont     = [];   % unit whose coincident spikes were stripped
autoCCGStripClean    = [];   % the cleaner partner of that pair
autoPairType         = 'Merged';
if isfield(loadedPrefs,'autoPairType') && (ischar(loadedPrefs.autoPairType) || isstring(loadedPrefs.autoPairType)) ...
        && any(strcmp(char(loadedPrefs.autoPairType), ...
        {'Merged','Overlap removed','SUA+','SUA','MUA+','MUA','NaN'}))
    autoPairType = char(loadedPrefs.autoPairType);
end
lastAutoPair         = 0;


if exist('Next.png','file')
    btnNextPair = uibutton(processingLayout, 'Icon','Next.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(1), 'Tooltip','Next similar pair. Ctrl+]');
else
    btnNextPair = uibutton(processingLayout, 'Text','Next',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(1), 'Tooltip','Next similar pair. Ctrl+]');
end
btnNextPair.Layout.Row = 2;
btnNextPair.Layout.Column = 6;
applyColorScheme(btnNextPair, figColor);

if exist('Prev.png','file')
    btnPrevPair = uibutton(processingLayout, 'Icon','Prev.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(-1), 'Tooltip','Previous similar pair. Ctrl+[');
else
    btnPrevPair = uibutton(processingLayout, 'Text','Previous',...
        'ButtonPushedFcn',@(btn,ev)onSelectSimilarPairs(-1), 'Tooltip','Previous similar pair. Ctrl+[');
end
btnPrevPair.Layout.Row = 1;
btnPrevPair.Layout.Column = 6;
applyColorScheme(btnPrevPair, figColor);

pairNumber = uilabel(processingLayout,'Text','Pair #:',...
    'FontWeight','bold','HorizontalAlignment','left',...
    'Tooltip','Similar-pair walker. Step with Ctrl+[ and Ctrl+].');
pairNumber.Layout.Row=1;
pairNumber.Layout.Column = 7;
applyColorScheme(pairNumber, figColor);

uifPair = uieditfield(processingLayout,'Value','0','Editable','on',...
    'Tooltip','Index of the displayed similar pair. Step with Ctrl+[ and Ctrl+].',...
    'ValueChangedFcn',@(ui,ev)changeuifPar(ui.Value,0));
uifPair.Layout.Row=2;
uifPair.Layout.Column=7;
applyColorScheme(uifPair, figColor);

uifPair2 = uieditfield(processingLayout,'Value','/0','Editable','off',...
    'Tooltip','Total similar pairs. Step with Ctrl+[ and Ctrl+].',...
    'ValueChangedFcn',@(ui,ev)changeuifParAll(ui.Value,0));
uifPair2.Layout.Row=2;
uifPair2.Layout.Column=8;
applyColorScheme(uifPair2, figColor);

% Auto curation panel (parameters + Auto button)
autoPanel = uipanel(leftGrid, 'title','Auto Curation');
autoPanel.Layout.Row = 2;
autoPanel.Layout.Column = 1;
applyColorScheme(autoPanel, figColor);

% Auto panel layout: 4 rows x 10 cols. Each parameter cell is
% label-on-top + edit-below (2 rows tall). The top half holds
% CCG, Sim, Amp.Var., PC Dist; the bottom half holds Overlap,
% ISI, Budget, Dist. CCG sits at col 1 directly above Overlap so
% the two coincidence-related gates are paired top-to-bottom in
% the leftmost column. Cols 5 and 9 are 0-width spacers that
% double the visual gap before the pair walker (col 6+) and
% before the Auto button (col 10) -- keeping the parameter
% block, walker, and trigger button as three clearly separated
% groups without depending on per-column-pair spacing (which
% uigridlayout doesn't support).
%
%       Col 1   Col 2  Col 3  Col 4    Col 5   Col 6  Col 7   Col 8  Col 9   Col 10 (Auto)
% R1:   Dist:   Sim:   Amp:   PC Dist: ----- gap PairType: -- (spans 6-8) -- gap   Auto
% R2:   [v]     [v]    [v]    [v]      ----- gap [dropdown] - (spans 6-8) -- gap   (Auto)
% R3:   Overlap:ISI %: Budget:CCG:     ----- gap Prev    Pair #:  (empty)   gap   (Auto)
% R4:   [v]     [v]    [v]    [v]      ----- gap Next    [v]      /0        gap   (Auto)
autoLayout = uigridlayout(autoPanel,[4 10], ...
    'ColumnWidth',{60, 75, 65, 65, 0, 28, 50, 30, 0, 75}, ...
    'RowHeight',{'fit','fit','fit','fit'}, ...
    'ColumnSpacing',10, ...
    'RowSpacing',4, ...
    'Padding',[5 5 5 5]);
applyColorScheme(autoLayout, figColor);

% --- CCG peak/median ratio (bottom-right slot) ---
% Lag-zero CCG peak must exceed autoCCGRatio times the median of
% the off-zero bins for cleanup to fire on a candidate pair.
lblAutoCCG = uilabel(autoLayout,'Text','CCG:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoCCG.Layout.Row = 3;
lblAutoCCG.Layout.Column = 4;
applyColorScheme(lblAutoCCG, figColor);

uiAutoCCG = uieditfield(autoLayout,'numeric','Value',autoCCGRatio, ...
    'Limits',[0.5 10], ...
    'Tooltip','CCG lag-0 peak / baseline-median ratio. Higher = stricter. Default 5.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setAutoCCGRatio(ui.Value));
uiAutoCCG.Layout.Row = 4;
uiAutoCCG.Layout.Column = 4;
applyColorScheme(uiAutoCCG, figColor);

lblAutoSim = uilabel(autoLayout,'Text','Similarity:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoSim.Layout.Row = 1;
lblAutoSim.Layout.Column = 2;
applyColorScheme(lblAutoSim, figColor);

uiAutoSim = uieditfield(autoLayout,'numeric','Value',autoSimVal, ...
    'Limits',[0 1], ...
    'Tooltip','Waveform similarity gate.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setAutoSimVal(ui.Value));
uiAutoSim.Layout.Row = 2;
uiAutoSim.Layout.Column = 2;
applyColorScheme(uiAutoSim, figColor);

lblAutoAmp = uilabel(autoLayout,'Text','Amp. Var.:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoAmp.Layout.Row = 1;
lblAutoAmp.Layout.Column = 3;
applyColorScheme(lblAutoAmp, figColor);

uiAutoAmp = uieditfield(autoLayout,'numeric','Value',autoAmpVal, ...
    'Limits',[0 1], ...
    'Tooltip','Max amplitude mismatch between two units. Lower = stricter.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setAutoAmpVal(ui.Value));
uiAutoAmp.Layout.Row = 2;
uiAutoAmp.Layout.Column = 3;
applyColorScheme(uiAutoAmp, figColor);

% --- PC Dist (col 4 top, above Dist) -----------------------------
lblAutoPC = uilabel(autoLayout,'Text','PC Dist:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoPC.Layout.Row = 1;
lblAutoPC.Layout.Column = 4;
applyColorScheme(lblAutoPC, figColor);

uiAutoPC = uieditfield(autoLayout,'numeric','Value',autoPCVal, ...
    'Limits',[0 5], ...
    'Tooltip','Normalised PC shape distance gate. Lower = stricter.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setAutoPCVal(ui.Value));
uiAutoPC.Layout.Row = 2;
uiAutoPC.Layout.Column = 4;
applyColorScheme(uiAutoPC, figColor);


% Walker selector for the Auto panel. Two flavours of target:
%   PAIR options (highlight TWO related units per click):
%     * 'Merged'          -> autoMergedPrimary / autoMergedAbsorbed
%     * 'Overlap removed' -> autoOverlapDroppedK / autoOverlapTriggerJ
%   UNIT options (highlight ONE unit per click, by isolation class):
%     * 'SUA+', 'SUA', 'MUA+', 'MUA' -> units of that isolation class.
%     * 'NaN'                        -> dropped/unclassified units
%                                       (groupList NaN or iso 'NA').
% (CCG-cleaned pairs do not get their own walker entry -- the cleanup
% phase only strips coincident spikes, it does not change unit
% identity, so there is nothing to step through visually.)
% Spans cols 5-7 on the top two rows so the dropdown sits exactly
% on top of the Prev/Next + Pair# walker below.
lblAutoPairType = uilabel(autoLayout,'Text','Pair type:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoPairType.Layout.Row = 1;
lblAutoPairType.Layout.Column = [6 8];
applyColorScheme(lblAutoPairType, figColor);

autoPairTypeItems = {'Merged','Overlap removed', ...
                     'SUA+','SUA','MUA+','MUA','NaN'};
% Make sure autoPairType matches one of the dropdown's items
% (a stale value from an earlier build would otherwise throw when
% uidropdown tries to apply it as Value).
if ~ismember(autoPairType, autoPairTypeItems)
    autoPairType = 'Merged';
end
ddAutoPairType = uidropdown(autoLayout, ...
    'Items', autoPairTypeItems, ...
    'Value', autoPairType, ...
    'Tooltip','Walker target: pair classes step through pairs; isolation classes step through single units.', ...
    'ValueChangedFcn', @(dd,ev) setAutoPairType(dd.Value));
ddAutoPairType.Layout.Row = 2;
ddAutoPairType.Layout.Column = [6 8];
applyColorScheme(ddAutoPairType, figColor);

% Auto button: rightmost column, spans all four rows so it reads as
% the pipeline trigger rather than just another parameter cell.
% Two-line label keeps the button compact while making the action
% explicit ('Auto curation' instead of just 'Auto').
btnAutoCurate = uibutton(autoLayout,'Text',sprintf('Auto\ncuration'),...
    'Tooltip','Run automatic curation.',...
    'ButtonPushedFcn',@(btn,ev)onAutoCurate());
btnAutoCurate.Layout.Row = [1 4];
btnAutoCurate.Layout.Column = 10;
applyColorScheme(btnAutoCurate, figColor);

% --- Bottom half (rows 3-4): drop / quality params + Dist ---------
% Overlap sits in col 1 directly under CCG so the two
% coincidence-related gates are paired top-to-bottom.
lblAutoOverlap = uilabel(autoLayout,'Text','Overlap:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoOverlap.Layout.Row = 3;
lblAutoOverlap.Layout.Column = 1;
applyColorScheme(lblAutoOverlap, figColor);

uiAutoOverlap = uieditfield(autoLayout,'numeric','Value',autoOverlapFrac, ...
    'Limits',[0 1], ...
    'Tooltip','Drop a unit if this fraction of its spikes overlap with a neighbour. Lower = stricter.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setOverlapVal(ui.Value));
uiAutoOverlap.Layout.Row = 4;
uiAutoOverlap.Layout.Column = 1;
applyColorScheme(uiAutoOverlap, figColor);

lblAutoIsi = uilabel(autoLayout,'Text','ISI %:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoIsi.Layout.Row = 3;
lblAutoIsi.Layout.Column = 2;
applyColorScheme(lblAutoIsi, figColor);

uiAutoIsi = uieditfield(autoLayout,'numeric','Value',autoIsiVal, ...
    'Limits',[0 100], ...
    'Tooltip','Max ISI violation %% for the merged unit.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setAutoIsiVal(ui.Value));
uiAutoIsi.Layout.Row = 4;
uiAutoIsi.Layout.Column = 2;
applyColorScheme(uiAutoIsi, figColor);

lblAutoBudget = uilabel(autoLayout,'Text','Budget:',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoBudget.Layout.Row = 3;
lblAutoBudget.Layout.Column = 3;
applyColorScheme(lblAutoBudget, figColor);

uiAutoBudget = uieditfield(autoLayout,'numeric','Value',autoIsiBudget, ...
    'Limits',[1 10], ...
    'Tooltip','Merged-ISI budget vs parents'' average. 1 = no slack.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setAutoBudgetVal(ui.Value));
uiAutoBudget.Layout.Row = 4;
uiAutoBudget.Layout.Column = 3;
applyColorScheme(uiAutoBudget, figColor);

% Dist sits in the top-left slot.
lblAutoDist = uilabel(autoLayout,'Text','Dist (µm):',...
    'FontWeight','bold','HorizontalAlignment','Left');
lblAutoDist.Layout.Row = 1;
lblAutoDist.Layout.Column = 1;
applyColorScheme(lblAutoDist, figColor);

uiAutoDist = uieditfield(autoLayout,'numeric','Value',autoDistVal, ...
    'Limits',[0 range(ylocs)], ...
    'Tooltip','Max y-distance (µm) between two units. Smaller = stricter.',...
    'LowerLimitInclusive','on','UpperLimitInclusive','on', ...
    'ValueChangedFcn',@(ui,ev)setAutoDistVal(ui.Value));
uiAutoDist.Layout.Row = 2;
uiAutoDist.Layout.Column = 1;
applyColorScheme(uiAutoDist, figColor);

% Pair walker (mirrors the Similarity Processing panel's Next/Prev
% buttons but iterates only over pairs the most recent Auto run
% produced for the selected pair type). Sits under the dropdown:
% Prev/Next stacked in col 6, Pair# label/edit stacked in col 7,
% /total parked in col 8 row 4.
if exist('Prev.png','file')
    btnPrevAutoPair = uibutton(autoLayout, 'Icon','Prev.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectMergedPairs(-1), ...
        'Tooltip','Previous Auto pair / unit. Ctrl+Up');
else
    btnPrevAutoPair = uibutton(autoLayout, 'Text','Prev',...
        'ButtonPushedFcn',@(btn,ev)onSelectMergedPairs(-1), ...
        'Tooltip','Previous Auto pair / unit. Ctrl+Up');
end
btnPrevAutoPair.Layout.Row = 3;
btnPrevAutoPair.Layout.Column = 6;
applyColorScheme(btnPrevAutoPair, figColor);

if exist('Next.png','file')
    btnNextAutoPair = uibutton(autoLayout, 'Icon','Next.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectMergedPairs(1), ...
        'Tooltip','Next Auto pair / unit. Ctrl+Down');
else
    btnNextAutoPair = uibutton(autoLayout, 'Text','Next',...
        'ButtonPushedFcn',@(btn,ev)onSelectMergedPairs(1), ...
        'Tooltip','Next Auto pair / unit. Ctrl+Down');
end
btnNextAutoPair.Layout.Row = 4;
btnNextAutoPair.Layout.Column = 6;
applyColorScheme(btnNextAutoPair, figColor);

% Counter title flips between '#Pair' and '#Unit' depending on
% whether the dropdown is sitting on a pair-mode or unit-mode item.
% Initialise it from autoPairType so a stale unit-class value
% (re-loaded session, hot-reload) renders the right label even
% before the user touches the dropdown.
if any(strcmp(autoPairType, {'SUA+','SUA','MUA+','MUA','NaN'}))
    lblAutoPairText = 'Unit #:';
else
    lblAutoPairText = 'Pair #:';
end
lblAutoPair = uilabel(autoLayout,'Text',lblAutoPairText,...
    'FontWeight','bold','HorizontalAlignment','left',...
    'Tooltip','Auto pair / unit walker. Step with Ctrl+Up and Ctrl+Down.');
lblAutoPair.Layout.Row = 3;
lblAutoPair.Layout.Column = 7;
applyColorScheme(lblAutoPair, figColor);

uifAutoPair = uieditfield(autoLayout,'Value','0','Editable','on',...
    'Tooltip','Index of the displayed Auto pair / unit. Step with Ctrl+Up and Ctrl+Down.',...
    'ValueChangedFcn',@(ui,ev)changeAutoPair(ui.Value));
uifAutoPair.Layout.Row = 4;
uifAutoPair.Layout.Column = 7;
applyColorScheme(uifAutoPair, figColor);

uifAutoPair2 = uieditfield(autoLayout,'Value','/0','Editable','off',...
    'Tooltip','Total Auto pairs / units. Step with Ctrl+Up and Ctrl+Down.');
uifAutoPair2.Layout.Row = 4;
uifAutoPair2.Layout.Column = 8;
applyColorScheme(uifAutoPair2, figColor);

% criteria selection panel
middlePanel = uipanel(leftGrid, 'title','Group Selection');
middlePanel.Layout.Row = 4;
middlePanel.Layout.Column = 1;
applyColorScheme(middlePanel, figColor);

middleButtonLayout = uigridlayout(middlePanel,[3 10], ...
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
    'Tooltip','Min firing rate (Hz) for the Include filter.',...
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
    'Tooltip','Max ISI violation %% for the Include filter.',...
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
    'Tooltip','Min SNR for the Include filter.',...
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
    'Tooltip','Spike-polarity filter for Include.',...
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
    'Tooltip','Toggle Vis on every selected unit.',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onToggleVis(cb.Value));
chkTogVis.Layout.Row = 2;
chkTogVis.Layout.Column = 5;
applyColorScheme(chkTogVis, figColor);

chkInclusion = uicheckbox(middleButtonLayout,'Text','Include',...
    'Tooltip','Tick units that pass all filters above.',...
    'FontWeight','bold',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onINCEXC(cb.Value,1));
chkInclusion.Layout.Row = 1;
chkInclusion.Layout.Column = 6;
applyColorScheme(chkInclusion, figColor);
incChan =[];
exChan = [];

chkExclusion = uicheckbox(middleButtonLayout,'Text','Exclude',...
    'Tooltip','Tick units that fail any filter above.',...
    'FontWeight','bold',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onINCEXC(cb.Value,0));
chkExclusion.Layout.Row = 2;
chkExclusion.Layout.Column = 6;
applyColorScheme(chkExclusion, figColor);


btnDeselectUnit = uibutton(middleButtonLayout, 'Text','Deselect',...
    'ButtonPushedFcn',@(btn,ev)onDeselect, ...
    'Tooltip','Untick all selected units. Ctrl+H');
btnDeselectUnit.Layout.Row=2;
btnDeselectUnit.Layout.Column = 7;
applyColorScheme(btnDeselectUnit, figColor);

if exist('Next.png','file')
    btnNextUnit = uibutton(middleButtonLayout, 'Icon','Next.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(1), 'Tooltip','Next unit. Ctrl+N');
else
    btnNextUnit = uibutton(middleButtonLayout, 'Text','Next',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(1), 'Tooltip','Next unit. Ctrl+N');
end
btnNextUnit.Layout.Row = 2;
btnNextUnit.Layout.Column = 8;
applyColorScheme(btnNextUnit, figColor);

if exist('Prev.png','file')
    btnPrevUnit = uibutton(middleButtonLayout, 'Icon','Prev.png', 'Text','',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(-1), 'Tooltip','Previous unit. Ctrl+P');
else
    btnPrevUnit = uibutton(middleButtonLayout, 'Text','Previous',...
        'ButtonPushedFcn',@(btn,ev)onSelectUnit(-1), 'Tooltip','Previous unit. Ctrl+P');
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

uifUnit2 = uieditfield(middleButtonLayout,'Value',['/' num2str(length(selectUnits))],'Editable','off',...
    'ValueChangedFcn',@(ui,ev)changeuifUnitAll(ui.Value,0));
uifUnit2.Layout.Row=2;
uifUnit2.Layout.Column=10;
applyColorScheme(uifUnit2, figColor);


% Pagination row 3 layout (10-col middleButtonLayout):
%   Col 1   First (⏮)
%   Cols 2-3 Prev (◀)
%   Col 4   "Size:" label
%   Col 5   page-size edit
%   Col 6   editable current-page number
%   Col 7   "/Y" total-pages label
%   Cols 8-9 Next (▶)
%   Col 10  Last (⏭)
if totalPages > 1
    btnFirstPage = uibutton(middleButtonLayout,'Text','⏮',...
        'ButtonPushedFcn',@(btn,ev)goToPage(1), ...
        'Tooltip','First page.');
    btnFirstPage.Layout.Row = 3;
    btnFirstPage.Layout.Column = 1;
    applyColorScheme(btnFirstPage, figColor);

    btnPrevPage = uibutton(middleButtonLayout,'Text','◀ Prev',...
        'ButtonPushedFcn',@(btn,ev)onChangePage(-1), ...
        'Tooltip','Previous page. Ctrl+Left');
    btnPrevPage.Layout.Row = 3;
    btnPrevPage.Layout.Column = [2 3];
    applyColorScheme(btnPrevPage, figColor);

    lblPageSize = uilabel(middleButtonLayout,'Text','Size:',...
        'FontWeight','bold','HorizontalAlignment','center');
    lblPageSize.Layout.Row = 3;
    lblPageSize.Layout.Column = 4;
    applyColorScheme(lblPageSize, figColor);

    uifPageSize = uieditfield(middleButtonLayout,'Value','50',...
        'Tooltip','Number of unit rows shown per page.', ...
        'ValueChangedFcn',@(ui,ev)changeuifPageSize(ui.Value));
    uifPageSize.Layout.Row=3;
    uifPageSize.Layout.Column=5;
    applyColorScheme(uifPageSize, figColor);

    % Editable current page. The total ("/Y") sits next to it as a
    % label so only the page number itself is user-typeable.
    uifPageNum = uieditfield(middleButtonLayout,'numeric', ...
        'Value', currentPage, ...
        'Limits',[1 max(1,totalPages)], ...
        'LowerLimitInclusive','on', 'UpperLimitInclusive','on', ...
        'RoundFractionalValues','on', ...
        'Tooltip','Jump to page.', ...
        'ValueChangedFcn', @(ui,ev) goToPage(ui.Value));
    uifPageNum.Layout.Row = 3;
    uifPageNum.Layout.Column = 6;
    applyColorScheme(uifPageNum, figColor);

    lblPage = uilabel(middleButtonLayout, ...
        'Text', sprintf('/%d', totalPages), ...
        'FontWeight','bold','HorizontalAlignment','left');
    lblPage.Layout.Row = 3;
    lblPage.Layout.Column = 7;
    applyColorScheme(lblPage, figColor);

    btnNextPage = uibutton(middleButtonLayout,'Text','Next ▶',...
        'ButtonPushedFcn',@(btn,ev)onChangePage(1), ...
        'Tooltip','Next page. Ctrl+Right');
    btnNextPage.Layout.Row = 3;
    btnNextPage.Layout.Column = [8 9];
    applyColorScheme(btnNextPage, figColor);

    btnLastPage = uibutton(middleButtonLayout,'Text','⏭',...
        'ButtonPushedFcn',@(btn,ev)goToPage(totalPages), ...
        'Tooltip','Last page.');
    btnLastPage.Layout.Row = 3;
    btnLastPage.Layout.Column = 10;
    applyColorScheme(btnLastPage, figColor);
end

lastUnit = 0;

%% Groups check boxes and radio panels
checkPanel = uipanel(leftGrid, ...
    'Title','Groups:', ...
    'BackgroundColor',figColor,...
    'Scrollable','on');
checkPanel.Layout.Row = 5;
checkPanel.Layout.Column = 1;
applyColorScheme(checkPanel, figColor);

% Calculate page index ranges for pagination
startIdx = (currentPage-1) * PAGE_SIZE + 1;
endIdx = min(startIdx + PAGE_SIZE - 1, numel(displayedGroups));
visibleGroups = displayedGroups(max(1,startIdx):endIdx);
visibleRows = length(visibleGroups) + 1;  % +1 for the "All" row

mainCheckGrid = uigridlayout(checkPanel, [visibleRows, 10], ...
    'RowHeight', repmat({'fit'},1,visibleRows), ...
    'ColumnWidth', repmat({'fit'},1,10), ...
    'RowSpacing',2, ...
    'Scrollable','on',...
    'ColumnSpacing',5, ...
    'Padding',[5 5 5 5]);
applyColorScheme(mainCheckGrid, figColor);

% Arrays for UI handles - only initialize for visible groups to save memory
lblGroupHandles       = gobjects(numGroups,1);
lblIdHandles          = gobjects(numGroups,1);
lblChannelHandles     = gobjects(numGroups,1);
groupVisCheckboxes    = gobjects(numGroups,1);
groupSelectCheckboxes = gobjects(numGroups,1);
groupRadioButtons     = gobjects(numGroups,1);
lblIsolation          = gobjects(numGroups,1);
lblRate               = gobjects(numGroups,1);
lblDetectablity       = gobjects(numGroups,1);
lblISIViolation       = gobjects(numGroups,1);

% Object for multiple plot
plotItemsNames  = {'Amplitude','Waveform', 'CCG', 'ISI', 'Density', 'Amp Dist', 'PCs', 'Trace', 'Features'};
multipleSpinnerObj = gobjects(length(plotItemsNames),1);
multiplePlotObj = gobjects(length(plotItemsNames),1);
plotScaleFactor = ones(1, length(plotItemsNames));
reScaledFlag = false;

radioGroupBG = uibuttongroup(mainCheckGrid, ...
    'Title','', ...
    'Scrollable','on');
radioGroupBG.Layout.Row = [2 visibleRows];
radioGroupBG.Layout.Column = 10;
applyColorScheme(radioGroupBG, figColor);
rowPixelHeight = 25;
rowOffsetTop   = 5;

lblMerge = uilabel(mainCheckGrid,'Text','Merge To:',...
    'FontWeight','bold');
lblMerge.Layout.Row = 1;
lblMerge.Layout.Column = 10;
applyColorScheme(lblMerge, figColor);

% all selection row
lblIdHdr = uilabel(mainCheckGrid,'Text','ID:','FontWeight','bold');
lblIdHdr.Layout.Row = 1;
lblIdHdr.Layout.Column = 1;
applyColorScheme(lblIdHdr, figColor);

lblAll = uilabel(mainCheckGrid,'Text','G #:',...
    'FontWeight','bold');
lblAll.Layout.Row = 1;
lblAll.Layout.Column = 2;
applyColorScheme(lblAll, figColor);

lblChannel = uilabel(mainCheckGrid,'Text','Ch#:','FontWeight','bold');
lblChannel.Layout.Row = 1;
lblChannel.Layout.Column = 3;
applyColorScheme(lblChannel, figColor);

lblunitIsolation = uilabel(mainCheckGrid,'Text','Iso:','FontWeight','bold');
lblunitIsolation.Layout.Row = 1;
lblunitIsolation.Layout.Column = 4;
applyColorScheme(lblunitIsolation, figColor);

lblunitRate = uilabel(mainCheckGrid,'Text','Rate:','FontWeight','bold');
lblunitRate.Layout.Row = 1;
lblunitRate.Layout.Column = 5;
applyColorScheme(lblunitRate, figColor);

lblSNR = uilabel(mainCheckGrid,'Text','SNR:','FontWeight','bold');
lblSNR.Layout.Row = 1;
lblSNR.Layout.Column = 6;
applyColorScheme(lblSNR, figColor);

lblISIV = uilabel(mainCheckGrid,'Text','ISIV %:');
lblISIV.Layout.Row = 1;
lblISIV.Layout.Column = 7;
applyColorScheme(lblISIV, figColor);

chkAllVis = uicheckbox(mainCheckGrid,'Text','Vis. All',...
    'FontWeight','bold',...
    'Value',false,...
    'ValueChangedFcn',@(cb,ev)onAllVisChecked(cb.Value));
chkAllVis.Layout.Row = 1;
chkAllVis.Layout.Column = 8;
applyColorScheme(chkAllVis, figColor);

chkAllSel = uicheckbox(mainCheckGrid,'Text','Sel. All',...
    'FontWeight','bold',...
    'ValueChangedFcn',@(cb,ev)onAllSelectChecked(cb.Value));
chkAllSel.Layout.Row = 1;
chkAllSel.Layout.Column = 9;
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

plotButtonGrid = uigridlayout(plotButtonPanel,[3 14],...
    'ColumnWidth',{'1x','1x','1x','1x','1x','1x','1x','1x','1x','1x','1x', '1x','1x','fit'},...
    'RowHeight',{'fit','fit','fit'},...
    'Padding',[0 0 0 0],...
    'RowSpacing', 1,...
    'ColumnSpacing',8);
applyColorScheme(plotButtonGrid, figColor);

plotTypeDD = uidropdown(plotButtonGrid,...
    'Items',{'Amplitude','Waveform', 'CCG', 'ISI', 'Density', 'Amp Dist', 'PCs','Trace','Features','Multiple'},...
    'Tooltip','Select a plot, or Multiple to show several at once.',...
    'Value','Trace',...
    'ValueChangedFcn',@(dd,ev)updatePlotType(dd.Value));
plotTypeDD.Layout.Row = 2;
plotTypeDD.Layout.Column = 1;
applyColorScheme(plotTypeDD, figColor);

plotType = 'Trace';
lastPlotted = plotType;

plotFilterTypeDD = uidropdown(plotButtonGrid,...
    'Items',{'filtered','raw'},...
    'Tooltip','Bandpass-filtered or raw signal.',...
    'Value','filtered',...
    'ValueChangedFcn',@(dd,ev)updatePlotFilterType(dd.Value));
plotFilterTypeDD.Layout.Row = 2;
plotFilterTypeDD.Layout.Column = 2;
applyColorScheme(plotFilterTypeDD, figColor);

plotFilterType = 'filtered';

chTypeDropdown = uidropdown(plotButtonGrid,...
    'Items',{'Config' ,'0', '1', '2', '3', '4', '5', '6', '7','8','9','10'},...
    'Tooltip','How many channels above/below the main one.',...
    'Value','Config');
chTypeDropdown.Layout.Row = 2;
chTypeDropdown.Layout.Column = 3;
applyColorScheme(chTypeDropdown, figColor);
chTypeDropdown.ValueChangedFcn      = @(sld,ev)updateChType(sld.Value);

numWaveformDropdown = uidropdown(plotButtonGrid,...
    'Items',{'1','5','10','25','50','100','1000'},...
    'Tooltip','Max waveforms drawn per unit.',...
    'Value','25');
numWaveformDropdown.Layout.Row = 2;
numWaveformDropdown.Layout.Column = 4;
applyColorScheme(numWaveformDropdown, figColor);
numWaveformDropdown.ValueChangedFcn      = @(sld,ev)updateNumWaveforms(sld.Value);

numWaveforms = 25;

spikeWaveDurDD = uidropdown(plotButtonGrid,...
    'Items',{'Config','0.25','0.5','1','1.5','2'},...
    'Tooltip','Half-window of waveform plot (ms).',...
    'Value','Config',...
    'ValueChangedFcn',@(dd,ev)updateWaveDur(dd.Value));
spikeWaveDurDD.Layout.Row = 2;
spikeWaveDurDD.Layout.Column = 5;
applyColorScheme(spikeWaveDurDD, figColor);

meanWaveformDD = uidropdown(plotButtonGrid,...
    'Items',{'Individual','Mean'},...
    'Tooltip','Show every spike or just the mean.',...
    'Value','Individual',...
    'ValueChangedFcn',@(dd,ev)updateWavetype(dd.Value));
meanWaveformDD.Layout.Row = 2;
meanWaveformDD.Layout.Column = 6;
applyColorScheme(meanWaveformDD, figColor);
plotMean = 0;

ccgLagDD = uidropdown(plotButtonGrid,...
    'Items',{'25','50','75','100','200'},...
    'Tooltip','CCG half-window (ms).',...
    'Value','100',...
    'ValueChangedFcn',@(dd,ev)updateCCGLag(dd.Value));
ccgLagDD.Layout.Row = 2;
ccgLagDD.Layout.Column = 7;
applyColorScheme(ccgLagDD, figColor);

ccgSmoothDD = uidropdown(plotButtonGrid,...
    'Items',{'0','1','2', '5', '10'},...
    'Tooltip','Boxcar smoothing for CCG / density. 0 = none.',...
    'Value','0',...
    'ValueChangedFcn',@(dd,ev)updateCCGSmoothN(dd.Value));
ccgSmoothDD.Layout.Row = 2;
ccgSmoothDD.Layout.Column = 8;
applyColorScheme(ccgSmoothDD, figColor);

isiThrDD = uidropdown(plotButtonGrid,...
    'Items',{'1','2','3','4','5'},...
    'Tooltip','ISI violation window (ms).',...
    'Value','1',...
    'ValueChangedFcn',@(dd,ev)updateISIThr(dd.Value));
isiThrDD.Layout.Row = 2;
isiThrDD.Layout.Column = 9;
applyColorScheme(isiThrDD, figColor);

ampSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 5,'MajorTicksMode','auto', 'Limits',[0 10], ...
    'Tooltip','Waveform vertical scale.');
ampSlider.MajorTicks = [0 10];
ampSlider.Layout.Row = 2;
ampSlider.Layout.Column = 10;
applyColorScheme(ampSlider, figColor);
ampSlider.ValueChangedFcn      = @(sld,ev)updateScale(sld.Value);

lineWidthSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 0.15,'MajorTicksMode','auto', 'Limits',[0 5], ...
    'Tooltip','Waveform line thickness.');
lineWidthSlider.MajorTicks = [0 5];
lineWidthSlider.Layout.Row = 2;
lineWidthSlider.Layout.Column = 11;
applyColorScheme(lineWidthSlider, figColor);
lineWidthSlider.ValueChangedFcn      = @(sld,ev)updateLineWidth(sld.Value);
lineWidth = 0.15;

alphaSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 0.5,'MajorTicksMode','auto', 'Limits',[0 1], ...
    'Tooltip','Trace opacity.');
alphaSlider.MajorTicks = [0 1];
alphaSlider.Layout.Row = 2;
alphaSlider.Layout.Column = 12;
applyColorScheme(alphaSlider, figColor);
alphaSlider.ValueChangedFcn      = @(sld,ev)updateAlphaLevel(sld.Value);
alphaLevel = 0.5;

chkOverlay = uicheckbox(plotButtonGrid, 'Text','', ...
    'Value', false, ...
    'Tooltip','Overlay all groups at the same location.', ...
    'ValueChangedFcn', @(cb,ev)updatePlotOverlay(cb.Value));
chkOverlay.Layout.Row = 2;
chkOverlay.Layout.Column = 13;
applyColorScheme(chkOverlay, figColor);
plotOverlay = false;

lblOverlay = uilabel(plotButtonGrid,'Text','Overlay:',...
    'HorizontalAlignment','left');
lblOverlay.Layout.Row = 1;
lblOverlay.Layout.Column = 13;
applyColorScheme(lblOverlay, figColor);

btnExportFig = uibutton(plotButtonGrid,'Text','Export Fig.',...
    'ButtonPushedFcn',@(btn,ev)onExportFig(), ...
    'Tooltip','Export the main plot as a vector file.');
btnExportFig.Layout.Row = 2;
btnExportFig.Layout.Column = 14;
applyColorScheme(btnExportFig, figColor);

% Open the curation table in a separate window. The table used to
% live as a plot type / multi-plot panel, which meant every
% selection click rebuilt it (and lost the user's column sort).
% In a standalone uifigure it survives across selection changes
% and can be refreshed on demand.
% Shortcut help button -- pops the same dialog as Ctrl+/.
% Sits underneath the Overlay checkbox (column 13) so it's
% discoverable without the user having to know the keyboard
% shortcut, and it slots into the empty cell that already
% existed in the third row of the plot toolbar.
btnShortcutHelp = uibutton(plotButtonGrid,'Text','?',...
    'ButtonPushedFcn',@(btn,ev) showShortcutsHelp(), ...
    'Tooltip','Show the list of keyboard shortcuts. Ctrl+/');
btnShortcutHelp.Layout.Row = 3;
btnShortcutHelp.Layout.Column = 13;
applyColorScheme(btnShortcutHelp, figColor);

btnTableWindow = uibutton(plotButtonGrid,'Text','Open table',...
    'ButtonPushedFcn',@(btn,ev) openTableWindow(), ...
    'Tooltip','Open the curation table in a standalone window.');
btnTableWindow.Layout.Row = 3;
btnTableWindow.Layout.Column = 14;
applyColorScheme(btnTableWindow, figColor);

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
    'Tooltip','Export format: raster image or vector.',...
    'ValueChangedFcn',@(dd,ev)updateExportFormat(dd.Value));
exportFormatDD.Layout.Row = 1;
exportFormatDD.Layout.Column = 14;
applyColorScheme(exportFormatDD, figColor);

exportType = 'Image';

% Apply persistent user prefs to the just-built widgets. Each block
% is guarded so a missing / corrupt field falls back to the default.
if ~isempty(fieldnames(loadedPrefs))
    try
        if isfield(loadedPrefs,'plotType') && any(strcmpi(loadedPrefs.plotType, plotTypeDD.Items))
            plotTypeDD.Value = loadedPrefs.plotType;
            plotType = loadedPrefs.plotType;
        end
        if isfield(loadedPrefs,'plotFilterType') && any(strcmpi(loadedPrefs.plotFilterType, plotFilterTypeDD.Items))
            plotFilterTypeDD.Value = loadedPrefs.plotFilterType;
            plotFilterType = loadedPrefs.plotFilterType;
        end
        if isfield(loadedPrefs,'numWaveforms') && ismember(num2str(loadedPrefs.numWaveforms), numWaveformDropdown.Items)
            numWaveformDropdown.Value = num2str(loadedPrefs.numWaveforms);
            numWaveforms = loadedPrefs.numWaveforms;
        end
        if isfield(loadedPrefs,'plotMean') && ismember(loadedPrefs.plotMean, [0 1])
            plotMean = loadedPrefs.plotMean;
            if plotMean == 1
                meanWaveformDD.Value = 'Mean';
            else
                meanWaveformDD.Value = 'Individual';
            end
        end
        if isfield(loadedPrefs,'plotOverlay')
            plotOverlay = logical(loadedPrefs.plotOverlay);
            chkOverlay.Value = plotOverlay;
        end
        if isfield(loadedPrefs,'ccgLag') && ismember(num2str(loadedPrefs.ccgLag), ccgLagDD.Items)
            ccgLagDD.Value = num2str(loadedPrefs.ccgLag);
            ccgLag = loadedPrefs.ccgLag;
        end
        if isfield(loadedPrefs,'smoothN') && ismember(num2str(loadedPrefs.smoothN), ccgSmoothDD.Items)
            ccgSmoothDD.Value = num2str(loadedPrefs.smoothN);
            smoothN = loadedPrefs.smoothN;
        end
        if isfield(loadedPrefs,'thresholdISI') && ismember(num2str(loadedPrefs.thresholdISI), isiThrDD.Items)
            isiThrDD.Value = num2str(loadedPrefs.thresholdISI);
            thresholdISI = loadedPrefs.thresholdISI;
        end
        if isfield(loadedPrefs,'ampScale') && isfinite(loadedPrefs.ampScale)
            ampScale = loadedPrefs.ampScale;
            ampSlider.Value = max(min(ampScale, ampSlider.Limits(2)), ampSlider.Limits(1));
        end
        if isfield(loadedPrefs,'lineWidth') && isfinite(loadedPrefs.lineWidth)
            lineWidth = loadedPrefs.lineWidth;
            lineWidthSlider.Value = max(min(lineWidth, lineWidthSlider.Limits(2)), lineWidthSlider.Limits(1));
        end
        if isfield(loadedPrefs,'alphaLevel') && isfinite(loadedPrefs.alphaLevel)
            alphaLevel = loadedPrefs.alphaLevel;
            alphaSlider.Value = max(min(alphaLevel, alphaSlider.Limits(2)), alphaSlider.Limits(1));
        end
        if isfield(loadedPrefs,'exportType') && any(strcmpi(loadedPrefs.exportType, exportFormatDD.Items))
            exportFormatDD.Value = loadedPrefs.exportType;
            exportType = loadedPrefs.exportType;
        end
    catch
        % Any pref restoration failure is non-fatal -- defaults stand.
    end
end

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
    'Tooltip','Time-window length (s). ''full'' = whole recording.',...
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
    'Tooltip','Channels in the vertical viewport.',...
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
ampScale     = 5;




%% Default run
plotOption();
normalizeTimeAmp();
refreshLabels();

%% =============== CALLBACKS ===============
    function keyPressHandler(src, event)
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

        % Unmodified Space inside the table window toggles the
        % focused unit's Vis/Select pair. We only fire when the
        % event came from the table figure -- pressing Space in
        % the main GUI shouldn't fight with text-entry widgets.
        if isempty(event.Modifier) && strcmpi(event.Key, 'space')
            if ~isempty(tableFigure) && isgraphics(tableFigure) && ...
                    isequal(src, tableFigure)
                toggleTableFocusedSelection();
            end
            return;
        end

        % Curation shortcuts fire on the Control key alone (no Cmd,
        % no Alt/Option, no Shift). Mac users press the dedicated
        % "control" / ⌃ key; Windows/Linux users press Ctrl. The
        % same key name on every platform.
        % Letters are mnemonic (D=Delete, M=Merge, T=spliT,
        % K=isolation tier, U=Undo, R=Reset, L=Limit, H=deselect,
        % E=Edit notes, S=Save, P/N=Prev/Next unit).
        if ismember('control', event.Modifier) ...
                && ~ismember('command', event.Modifier) ...
                && ~ismember('alt', event.Modifier) ...
                && ~ismember('shift', event.Modifier)
            % Track whether this keystroke actually mutated the
            % curation state. Only mutating actions need a full
            % table rebuild; navigation / help / save / bare-Ctrl
            % auto-repeat events do not -- and rebuilding for
            % those would silently throw away the user's column
            % sort. Bare-modifier events (event.Key = 'control')
            % don't match any case and leave dataChanged false,
            % which is exactly what stops them from re-rendering
            % the standalone table.
            dataChanged = false;
            switch lower(event.Key)
                case 'r', onReset();              dataChanged = true;
                case 'uparrow',   onSelectMergedPairs(-1)   % Auto walker prev
                case 'downarrow', onSelectMergedPairs(1)    % Auto walker next
                case {'leftbracket','open_bracket'}, onSelectSimilarPairs(-1)  % Similarity walker prev
                case {'rightbracket','close_bracket'}, onSelectSimilarPairs(1) % Similarity walker next
                case 'leftarrow', onChangePage(-1)
                case 'rightarrow', onChangePage(1)
                case 'p', onSelectUnit(-1)
                case 'n', onSelectUnit(1)
                case 's', onSave()
                case 'd', onRemove();             dataChanged = true;
                case 'm', onMerge();              dataChanged = true;
                case 'k', onMUA();                dataChanged = true;
                case 'u', onUndo();               dataChanged = true;
                case 'l', onAutoCut();            dataChanged = true;
                case 'h', onDeselect()
                case 't', onSplit();              dataChanged = true;   % T = spliT
                case 'e', onEditNotes();          dataChanged = true;   % E = Edit notes
                case 'slash', showShortcutsHelp()    % Ctrl+/ (also '?')
                case 'questionmark', showShortcutsHelp()
            end
            % Repaint selection colours on the standalone table so
            % shortcuts that change which units are ticked
            % (P/N walker steps, walker arrows, Deselect) update
            % the row tinting without rebuilding -- preserving
            % the user's column sort. Full rebuild is reserved
            % for actions that mutate the underlying data.
            if dataChanged && ~isempty(tableFigure) && isgraphics(tableFigure) && ...
                    isequal(src, tableFigure)
                try, refreshTableWindow(); catch, end
            elseif ~isempty(tableHandle) && isgraphics(tableHandle)
                try, paintSelectedTableRows(); catch, end
            end
        end
    end

    function toggleTableFocusedSelection()
        % Helper for the Space key. Reads the currently selected
        % cell in the standalone table, looks up the unit at that
        % row, and runs the Vis/Select toggle. Re-paints the row
        % colour and triggers a plot refresh so the rest of the
        % GUI tracks the change.
        if isempty(tableHandle) || ~isgraphics(tableHandle), return; end
        sel = tableHandle.Selection;
        if isempty(sel), return; end
        if size(sel, 2) >= 2
            row = sel(1, 1);
        elseif numel(sel) >= 1
            row = sel(1);
        else
            return;
        end
        if ~isfinite(row) || row < 1, return; end
        if isa(tableHandle.Data, 'table')
            if row > height(tableHandle.Data), return; end
        else
            if row > size(tableHandle.Data, 1), return; end
        end
        unitIdx = tableUnitIdx(tableHandle.Data{row, 1}, tableHandle.Data{row, 2});
        toggleUnitSelect(unitIdx);
        try, paintSelectedTableRows(); catch, end
        try, plotOption(); catch, end
    end

    function showShortcutsHelp()
        % Centralised reference for every keyboard shortcut wired into
        % the GUI. Triggered by Ctrl+? (or "/"). All shortcuts use
        % the Control key as the sole modifier (Mac users press
        % the dedicated "control" / ⌃ key, separate from Cmd) so
        % we don't collide with macOS Cmd+Shift+letter bindings.
        msg = sprintf([ ...
            'Keyboard shortcuts (Ctrl + key):\n\n', ...
            '  D    Remove selected units (Delete)\n', ...
            '  M    Merge selected units\n', ...
            '  T    Split selected unit\n', ...
            '  K    Cycle isolation tier of selected units\n', ...
            '  U    Undo selected units (also dissolves merge group)\n', ...
            '  R    Reset every edit (with confirmation)\n', ...
            '  L    Apply Limit (auto-trim) to selected units\n', ...
            '  E    Edit note for the radio-selected unit\n', ...
            '  H    Deselect every selected unit\n', ...
            '  N    Select next unit\n', ...
            '  P    Select previous unit\n', ...
            '  S    Save curated results\n', ...
            '  Up   / Down    Walk previous / next Auto pair (or unit)\n', ...
            '  [    / ]       Walk previous / next similar pair\n', ...
            '  Left / Right   Walk previous / next page\n', ...
            '  ?    Show this list\n\n', ...
            'Inside the standalone Curation Table window:\n', ...
            '  Space   Toggle Vis/Select for the focused row\n', ...
            '  Click   Set the focused row as the merge target\n', ...
            '  All Ctrl+key shortcuts above also work here.\n']);
        try
            uialert(parentFig, msg, 'Keyboard shortcuts', 'Icon', 'info');
        catch
            disp(msg);
        end
    end

    function changeuifPageSize(val)
        PAGE_SIZE = str2double(val);
        currentPage = 1;
        recomputeDisplayedGroups();
        onChangePage(0);
    end

    function recomputeDisplayedGroups()
        mask = preprocessed.firingRate(:) > 0 & ~isnan(groupList(:));
        displayedGroups = find(mask);
        totalPages = max(1, ceil(numel(displayedGroups) / PAGE_SIZE));
        if currentPage > totalPages, currentPage = totalPages; end
        % Default Unit# navigation walks exactly the shown groups, in the
        % same order as their displayed Group ID (1..#shown).
        if ~chanFilterActive
            selectUnits = displayedGroups;
        end
    end

    function vg = currentPageVisibleGroups()
        s = (currentPage-1) * PAGE_SIZE + 1;
        e = min(s + PAGE_SIZE - 1, numel(displayedGroups));
        if s > numel(displayedGroups)
            vg = displayedGroups([]);
        else
            vg = displayedGroups(s:e);
        end
    end

    function goToPage(target)
        % Jump straight to a specific page. Used by First / Last
        % buttons and the editable page number field. Clamps to
        % [1, totalPages] so a typo can't crash the navigation.
        if isempty(target) || ~isfinite(target), return; end
        target = round(target);
        if target < 1, target = 1; end
        if target > totalPages, target = totalPages; end
        currentPage = target;
        onChangePage(0);
    end

% Function to change page in pagination
    function onChangePage(direction)
        % Bail out if the main figure has been closed but a callback
        % (e.g. the table window's keyboard handler) still fires.
        % Without this, uigridlayout(checkPanel,...) below crashes
        % with "'Parent' cannot be set to a deleted object" when
        % checkPanel was destroyed alongside parentFig.
        if isempty(checkPanel) || ~isgraphics(checkPanel) || ~isvalid(checkPanel)
            return;
        end
        newPage = currentPage + direction;
        if newPage >= 1 && newPage <= totalPages
            currentPage = newPage;
            % Keep the editable page-number field and the "/Y" total
            % label in sync no matter how the page changed (button
            % click, keyboard shortcut, programmatic call).
            if exist('lblPage', 'var') && isvalid(lblPage)
                lblPage.Text = sprintf('/%d', totalPages);
            end
            if exist('uifPageNum','var') && isvalid(uifPageNum)
                if uifPageNum.Limits(2) ~= max(1,totalPages)
                    uifPageNum.Limits = [1 max(1,totalPages)];
                end
                uifPageNum.Value = currentPage;
            end

            % Clean up old UI elements
            delete(findall(mainCheckGrid, 'Type', 'UILabel'));
            delete(findall(mainCheckGrid, 'Type', 'UICheckbox'));
            delete(findall(radioGroupBG, 'Type', 'UIRadioButton'));

            % Rebuild grid with new page size
            delete(mainCheckGrid);

            % Calculate new visible groups
            visibleGroups = currentPageVisibleGroups();
            visibleRows = length(visibleGroups) + 1;  % +1 for the "All" row

            % Create new grid
            mainCheckGrid = uigridlayout(checkPanel, [visibleRows, 10], ...
                'RowHeight', repmat({'fit'},1,visibleRows), ...
                'ColumnWidth', repmat({'fit'},1,10), ...
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
            radioGroupBG.Layout.Column = 10;
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

            % If a previously-ticked group sits on this freshly-built
            % page, point the buttongroup's SelectedObject at its radio
            % so onMerge has a live target even when the merge partner
            % lives on another page.
            visGrpSel = currentPageVisibleGroups();
            for kSel = visGrpSel(:)'
                if chkBoxSelect(kSel) && isgraphics(groupRadioButtons(kSel))
                    radioGroupBG.SelectedObject = groupRadioButtons(kSel);
                    break;
                end
            end
            drawnow;

            % Refresh state of checkboxes based on current selection
            refreshVisibility();
            refreshLabels();
        end
    end

% Function to create header row in the group list
    function recreateHeaderRow()
        lblIdHdr = uilabel(mainCheckGrid,'Text','ID:','FontWeight','bold');
        lblIdHdr.Layout.Row = 1;
        lblIdHdr.Layout.Column = 1;
        applyColorScheme(lblIdHdr, figColor);

        lblAll = uilabel(mainCheckGrid,'Text','G #:',...
            'FontWeight','bold');
        lblAll.Layout.Row = 1;
        lblAll.Layout.Column = 2;
        applyColorScheme(lblAll, figColor);

        lblChannel = uilabel(mainCheckGrid,'Text','Ch#:','FontWeight','bold');
        lblChannel.Layout.Row = 1;
        lblChannel.Layout.Column = 3;
        applyColorScheme(lblChannel, figColor);

        lblunitIsolation = uilabel(mainCheckGrid,'Text','Iso:','FontWeight','bold');
        lblunitIsolation.Layout.Row = 1;
        lblunitIsolation.Layout.Column = 4;
        applyColorScheme(lblunitIsolation, figColor);

        lblunitRate = uilabel(mainCheckGrid,'Text','Rate:','FontWeight','bold');
        lblunitRate.Layout.Row = 1;
        lblunitRate.Layout.Column = 5;
        applyColorScheme(lblunitRate, figColor);

        lblSNR = uilabel(mainCheckGrid,'Text','SNR:','FontWeight','bold');
        lblSNR.Layout.Row = 1;
        lblSNR.Layout.Column = 6;
        applyColorScheme(lblSNR, figColor);

        lblISIV = uilabel(mainCheckGrid,'Text','ISIV %:');
        lblISIV.Layout.Row = 1;
        lblISIV.Layout.Column = 7;
        applyColorScheme(lblISIV, figColor);

        chkAllVis = uicheckbox(mainCheckGrid,'Text','Vis. All',...
            'FontWeight','bold',...
            'Value',false,...
            'ValueChangedFcn',@(cb,ev)onAllVisChecked(cb.Value));
        chkAllVis.Layout.Row = 1;
        chkAllVis.Layout.Column = 8;
        applyColorScheme(chkAllVis, figColor);

        chkAllSel = uicheckbox(mainCheckGrid,'Text','Sel. All',...
            'FontWeight','bold',...
            'ValueChangedFcn',@(cb,ev)onAllSelectChecked(cb.Value));
        chkAllSel.Layout.Row = 1;
        chkAllSel.Layout.Column = 9;
        applyColorScheme(chkAllSel, figColor);

        lblMerge = uilabel(mainCheckGrid,'Text','Merge To:',...
            'FontWeight','bold');
        lblMerge.Layout.Row = 1;
        lblMerge.Layout.Column = 10;
        applyColorScheme(lblMerge, figColor);
    end

% Create UI elements for visible groups
    function createGroupUiElements()
        visibleGroups = currentPageVisibleGroups();

        for i = 1:length(visibleGroups)
            iG = visibleGroups(i);
            rowIdx = i + 1;  % +1 because row 1 is the header

            % contiguous group id (rank among shown / nonzero groups)
            lblIdHandles(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('%d', (currentPage-1)*PAGE_SIZE + i),...
                'FontWeight','bold','FontColor',1-figColor);
            lblIdHandles(iG).Layout.Row = rowIdx;
            lblIdHandles(iG).Layout.Column = 1;

            % group list
            lblGroupHandles(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('G: %d', groupList(iG)),...
                'FontWeight','bold');
            lblGroupHandles(iG).Layout.Row = rowIdx;
            lblGroupHandles(iG).Layout.Column = 2;

            % channel list
            lblChannelHandles(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('Ch %d', channelList(iG)),...
                'FontWeight','bold');
            lblChannelHandles(iG).Layout.Row = rowIdx;
            lblChannelHandles(iG).Layout.Column = 3;

            lblIsolation(iG) = uilabel(mainCheckGrid, ...
                'Text', unitIsolation{iG},...
                'FontWeight','bold');
            lblIsolation(iG).Layout.Row = rowIdx;
            lblIsolation(iG).Layout.Column = 4;

            lblRate(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('%.1f',preprocessed.firingRate(iG)),...
                'FontWeight','bold','FontColor',1-figColor);
            lblRate(iG).Layout.Row = rowIdx;
            lblRate(iG).Layout.Column = 5;

            lblDetectablity(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('%.1f',1 + detectblity(iG)),...
                'FontWeight','bold','FontColor',1-figColor);
            lblDetectablity(iG).Layout.Row = rowIdx;
            lblDetectablity(iG).Layout.Column = 6;

            lblISIViolation(iG) = uilabel(mainCheckGrid, ...
                'Text', sprintf('%.1f',preprocessed.isiViolation(iG)),...
                'FontWeight','bold','FontColor',1-figColor);
            lblISIViolation(iG).Layout.Row = rowIdx;
            lblISIViolation(iG).Layout.Column = 7;

            % visualization checkboxes
            groupVisCheckboxes(iG) = uicheckbox(mainCheckGrid,...
                'Value',chkBoxVis(iG), 'Text','Vis.',...
                'ValueChangedFcn',@(cb,ev)onVisCheckboxChanged(iG));
            groupVisCheckboxes(iG).Layout.Row = rowIdx;
            groupVisCheckboxes(iG).Layout.Column = 8;

            % selection checkboxes
            groupSelectCheckboxes(iG) = uicheckbox(mainCheckGrid,...
                'Text','Select','Value',chkBoxSelect(iG),...
                'ValueChangedFcn',@(cb,ev)onSelectCheckboxChanged(iG));
            groupSelectCheckboxes(iG).Layout.Row = rowIdx;
            groupSelectCheckboxes(iG).Layout.Column = 9;

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
        visibleGroups = currentPageVisibleGroups();

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

        visibleGroups = currentPageVisibleGroups();

        for idx = 1:length(visibleGroups)
            k = visibleGroups(idx);
            if ~isgraphics(groupVisCheckboxes(k)), continue; end
            groupVisCheckboxes(k).Value = val;
            if val
                c = getGroupColor(groupList(k));
                groupVisCheckboxes(k).FontColor = c;
                groupVisCheckboxes(k).FontWeight ='Bold';
                if isgraphics(lblChannelHandles(k))
                    lblChannelHandles(k).BackgroundColor = c;
                    lblChannelHandles(k).FontColor = bestTextColorFor(c);
                end
            else
                groupVisCheckboxes(k).FontColor = 1-figColor;
                groupVisCheckboxes(k).FontWeight ='normal';
                if isgraphics(lblChannelHandles(k))
                    lblChannelHandles(k).BackgroundColor = 'none';
                    lblChannelHandles(k).FontColor = 1-figColor;
                end
            end
        end

        drawnow limitrate;
        plotOption();
    end

    function onAllSelectChecked(val)
        % Batch update selection checkboxes
        if val
            chkBoxSelect(:) = 1;
        else
            chkBoxSelect(:) = 0;
        end

        visibleGroups = currentPageVisibleGroups();

        for idx = 1:length(visibleGroups)
            k = visibleGroups(idx);
            if ~isgraphics(groupSelectCheckboxes(k)), continue; end
            groupSelectCheckboxes(k).Value = val;
            if val
                c = getGroupColor(groupList(k));
                groupSelectCheckboxes(k).FontColor = c;
                groupSelectCheckboxes(k).FontWeight ='Bold';
                if isgraphics(lblChannelHandles(k))
                    lblChannelHandles(k).BackgroundColor = c;
                    lblChannelHandles(k).FontColor = bestTextColorFor(c);
                end
            else
                groupSelectCheckboxes(k).FontColor = 1-figColor;
                groupSelectCheckboxes(k).FontWeight ='normal';
                if isgraphics(lblChannelHandles(k))
                    lblChannelHandles(k).BackgroundColor = 'none';
                    lblChannelHandles(k).FontColor = 1-figColor;
                end
            end
        end
        drawnow limitrate;
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
        % Mark each selected group as removed (groupList = NaN) but
        % keep its entries in sortedRes.unifiedLabels intact so the
        % unit's spikes can still be visualised via effLabel(k).
        selectedIdx = find(chkBoxSelect);
        for k = selectedIdx'
            groupList(k)   = NaN;
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

        % Validate via chkBoxSelect (closure state) rather than the
        % per-page UI handles. groupSelectCheckboxes(g) is a
        % GraphicsPlaceholder when g lives on a different page, so
        % indexing-and-reading-.Value through a multi-page mask
        % crashes with "Unrecognized property 'Value' for class
        % 'matlab.graphics.GraphicsPlaceholder'". chkBoxSelect
        % carries the same logical state without that hazard.
        sel_now  = chkBoxSelect(:) & (groupList(:) == selIdx);
        sel_orig = chkBoxSelect(:) & (originalGroupList(:) == selIdx);
        if ~any(sel_now) && ~any(sel_orig)
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

    function onReset(skipConfirm)
        % Confirm before wiping every curation edit -- Reset is the
        % only single click that can undo Auto sweeps + manual merges
        % + removals + isolation labels at once, so a misclick should
        % not be silent. The skipConfirm flag is used by onAutoCurate
        % (which calls onReset internally before each Auto run); the
        % public callbacks always confirm.
        if nargin < 1, skipConfirm = false; end
        if ~skipConfirm
            try
                answer = uiconfirm(parentFig, ...
                    'Reset will undo all your curation. Are you sure?', ...
                    'Confirm Reset', ...
                    'Options',{'Yes','No'}, ...
                    'DefaultOption','No', ...
                    'CancelOption','No', ...
                    'Icon','warning');
            catch
                % If the figure handle is somehow invalid, fall back
                % to the modal questdlg so we still respect the
                % user choice.
                answer = questdlg('Reset will undo all your curation. Are you sure?', ...
                    'Confirm Reset', 'Yes', 'No', 'No');
            end
            if ~strcmpi(answer, 'Yes')
                return;
            end
        end

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
        visPageGrp = currentPageVisibleGroups();
        for k = visPageGrp(:)'
            if isgraphics(groupSelectCheckboxes(k))
                groupSelectCheckboxes(k).Value = false;
            end
            if isgraphics(groupVisCheckboxes(k))
                groupVisCheckboxes(k).Value = false;
            end
        end




        selected_order = [];
        unitIsolation = repmat({'NA'},[numGroups,1]);
        % Recompute firing rate / ISI on the restored labels so the
        % Group rows show real numbers again (any earlier NaN values
        % from a removed-unit pass are wiped out here).
        preprocessSortedData();
        % Wipe every Auto-pair audit trail; every Auto-merge / drop /
        % strip has been undone by the label restore above.
        autoMergedPrimary    = [];
        autoMergedAbsorbed   = [];
        autoOverlapDroppedK  = [];
        autoOverlapTriggerJ  = [];
        autoCCGStripCont     = [];
        autoCCGStripClean    = [];
        lastAutoPair         = 0;
        if exist('uifAutoPair','var') && isvalid(uifAutoPair)
            uifAutoPair.Value  = '0';
        end
        if exist('uifAutoPair2','var') && isvalid(uifAutoPair2)
            uifAutoPair2.Value = '/0';
        end
        % Inline cache wipes (in addition to refresh* calls below).
        % Defensive against any edge case where refreshXxx's cell
        % reassignment fails to propagate to the parent workspace.
        xcorr_vals     = cell(numGroups, numGroups);
        isiCounts      = cell(numGroups, numGroups);
        isiCenters     = cell(numGroups, numGroups);
        isiViolations  = zeros(numGroups, numGroups);
        DenCenters     = cell(numGroups, numGroups);
        DenCounts      = cell(numGroups, numGroups);
        presence_ratio = zeros(numGroups, numGroups);
        refreshDensity();
        refreshISI();
        refreshCCG();
        refreshLabels();
        chkAllSel.Value = false;
        chkAllVis.Value = false;
        rdoNone.Value = true;
        updateChannelPlot();
        % Force plotOption to fully rebuild axChannels rather than just
        % deleting children -- otherwise a previously-merged unit's
        % stale tile / line objects can leave the new ACG looking flat.
        % lastMultiPlotOrder is reset too so plotMultiple recreates
        % its multiPlotPanels (otherwise those handles point to
        % children of the just-deleted axChannels and allchild()
        % errors).
        lastPlotted         = '';
        lastMultiPlotOrder  = [];
        selected_order_last = [];
        drawnow;
        plotOption();
    end

    function onSplit()
        % Split each ticked unit into N sub-clusters (set in the Split
        % options dropdown / num-clusters field next to the toolbar)
        % via the chosen method (K-means / GMM / Graph / UMAP+DBSCAN)
        % run on the loaded per-spike features + amplitude. Sub-clusters get
        % fresh labels = max(groupList,'omitnan') + 1, +2, ... ; the
        % highest-amplitude (highest-SNR) sub-cluster keeps the
        % original label so the user's "primary" identity is stable.
        % originalUnifiedLabels and originalGroupList are NOT
        % touched, so Reset / Undo unwind the split.
        %
        % While the routine is running we disable the Split button so
        % a second click can't enqueue an overlapping split. The
        % onCleanup hook re-enables it even if we error out.
        if exist('btnSplit','var') && isvalid(btnSplit)
            btnSplit.Enable = 'off';
            if isvalid(btnSplit) && isprop(btnSplit,'Text')
                btnSplit.Text = 'Splitting...';
            end
        end
        if exist('lblSave','var') && isvalid(lblSave)
            prevSaveLbl = lblSave.Text;
            lblSave.Text = 'Split in progress...';
        else
            prevSaveLbl = '';
        end
        drawnow;
        cleanupSplit = onCleanup(@() resetSplitButton(prevSaveLbl)); %#ok<NASGU>

        selectedIdx = find(chkBoxSelect);
        if isempty(selectedIdx)
            try
                uialert(parentFig, 'Tick at least one unit to split.', 'Split');
            catch
            end
            return;
        end
        if isempty(mappedData)
            if isfield(cfg,'inputFolder') && isfield(cfg,'outputFolder')
                mappedData = map_input_file(cfg.fullFilePath, cfg);
            else
                try
                    uialert(parentFig, 'Cannot split: input file is not accessible.', 'Split');
                catch
                end
                return;
            end
        end

        % Split clusters on per-spike features loaded from the result
        % (features.h5 + amplitude.h5). Both must align with spike_idx.
        haveFeat = isfield(sortedRes,'features') && ~isempty(sortedRes.features) && ...
                   size(sortedRes.features,1) == numel(sortedRes.spike_idx) && ...
                   isfield(sortedRes,'amplitude') && ~isempty(sortedRes.amplitude) && ...
                   numel(sortedRes.amplitude) == numel(sortedRes.spike_idx);
        if ~haveFeat
            try
                uialert(parentFig, ['Cannot split: per-spike features ' ...
                    '(features.h5 / amplitude.h5) are missing or do not match ' ...
                    'the spike count.'], 'Split');
            catch
            end
            return;
        end

        nSplitDone = 0;
        newSelected = [];
        for k = selectedIdx(:)'
            if isnan(groupList(k)),       continue; end
            if isnan(channelList(k)),     continue; end

            kLabel    = groupList(k);
            kRows     = find(sortedRes.unifiedLabels == kLabel);
            nSpikesK  = numel(kRows);
            if nSpikesK < 50, continue; end   % too few to split reliably

            % Features come straight from the sorter result: the first 3
            % loaded feature columns + per-spike amplitude, z-scored so
            % amplitude doesn't dominate. No waveform re-reading -- the
            % features already exist for every spike.
            nClusters = max(2, splitNumClusters);
            chCenter  = channelList(k);
            allRowsK  = kRows;

            nFeatLoad = min(3, size(sortedRes.features, 2));
            if nFeatLoad < 1, continue; end
            featLoad  = double(sortedRes.features(allRowsK, 1:nFeatLoad));
            ampLoad   = double(sortedRes.amplitude(allRowsK));
            feat_raw  = [featLoad, ampLoad(:)];
            feat_raw(~isfinite(feat_raw)) = 0;
            feat_mean = mean(feat_raw, 1);
            feat_std  = max(std(feat_raw, 0, 1), eps);
            feat      = (feat_raw - feat_mean) ./ feat_std;

            % Fit the clusterer on a subsample for the heavier methods;
            % the rest are assigned by nearest centroid in this same
            % feature space.
            sampleCap = min(nSpikesK, 2000);
            if nSpikesK > sampleCap
                sampleSel = sort(randperm(nSpikesK, sampleCap))';
            else
                sampleSel = (1:nSpikesK)';
            end
            featFit = feat(sampleSel, :);
            ampFit  = ampLoad(sampleSel);

            clustLabels = [];
            try
                switch lower(splitMethod)
                    case 'gmm'
                        gm = fitgmdist(featFit, nClusters, ...
                            'RegularizationValue', 0.01, ...
                            'Replicates', 10, ...
                            'Options', statset('MaxIter', 300), ...
                            'CovarianceType','full', ...
                            'SharedCovariance', false);
                        clustLabels = cluster(gm, featFit);
                    case {'graph','community','spectral'}
                        clustLabels = spectralcluster(featFit, nClusters);
                    case {'umap+dbscan','umap','dbscan'}
                        clustLabels = umapDbscanCluster(featFit, nClusters);
                    otherwise
                        clustLabels = kmeans(featFit, nClusters, ...
                            'Start','plus', ...
                            'Replicates', 10, ...
                            'Options', statset('MaxIter', 300));
                end
            catch
                try
                    clustLabels = kmeans(featFit, nClusters, ...
                        'Start','plus', 'Replicates', 10, ...
                        'Options', statset('MaxIter', 300));
                catch
                    continue;
                end
            end
            if isempty(clustLabels) || numel(unique(clustLabels)) < 2
                continue;
            end

            % Highest-amplitude cluster keeps the original label; the rest
            % become new units.
            uniqClusts   = unique(clustLabels(:))';
            nClustActual = numel(uniqClusts);
            ampClust     = zeros(nClustActual, 1);
            clustMeans   = zeros(nClustActual, size(feat, 2));
            for ci = 1:nClustActual
                mask = clustLabels == uniqClusts(ci);
                if ~any(mask), continue; end
                ampClust(ci)     = mean(ampFit(mask));
                clustMeans(ci,:) = mean(featFit(mask, :), 1, 'omitnan');
            end
            [~, ampOrder] = sort(abs(ampClust), 'descend');
            keepClust     = uniqClusts(ampOrder(1)); %#ok<NASGU>
            newClusts     = uniqClusts(ampOrder(2:end));

            % Assign every spike to its nearest cluster centroid.
            assignment            = nan(numel(allRowsK), 1);
            assignment(sampleSel) = clustLabels;
            unsampled = find(isnan(assignment));
            if ~isempty(unsampled)
                distU = pdist2(feat(unsampled, :), clustMeans);
                [~, nearest] = min(distU, [], 2);
                assignment(unsampled) = uniqClusts(nearest);
            end

            % Skip degenerate splits where every spike landed in one
            % cluster.
            survivingClusts = unique(assignment(:))';
            if numel(survivingClusts) < 2
                continue;
            end

            % --- Allocate one new unit slot per non-keeper cluster.
            for nc = newClusts(:)'
                if ~ismember(nc, survivingClusts), continue; end
                rowsForNewLabel = allRowsK(assignment == nc);
                if isempty(rowsForNewLabel), continue; end

                newLabel = max(groupList(~isnan(groupList))) + 1;
                newK     = numGroups + 1;

                groupList(newK,1)              = newLabel;
                channelList(newK,1)            = chCenter;
                mergedFlag(newK,1)             = false;
                unitIsolation{newK,1}          = unitIsolation{k};
                stable_length(newK,:)          = stable_length(k,:);
                mainPolarity(newK,1)           = mainPolarity(k);
                sidePolarity(newK,1)           = sidePolarity(k);
                detectblity(newK,1)            = detectblity(k);
                chkBoxSelect(newK,1)           = 1;
                chkBoxVis(newK,1)              = 1;
                preprocessed.firingRate(newK,1)    = 0;
                preprocessed.isiViolation(newK,1)  = 0;
                preprocessed.logFiringRate(newK,1) = log(eps);
                originalGroupList(newK,1)      = NaN;
                originalChannelList(newK,1)    = NaN;
                wfShape = size(sampleWaveform);
                sampleWaveform(newK,:,:)         = nan(1, wfShape(2), wfShape(3));
                originalSampleWaveform(newK,:,:) = nan(1, wfShape(2), wfShape(3));
                unitNotes{newK,1}             = sprintf('Split from G%g', kLabel);
                channelPlot(newK,:)           = chCenter + (-numChannelPlot:numChannelPlot);

                xcorr_vals(newK,:)     = cell(1, size(xcorr_vals,2));
                xcorr_vals(:,newK)     = cell(size(xcorr_vals,1), 1);
                isiCounts(newK,:)      = cell(1, size(isiCounts,2));
                isiCounts(:,newK)      = cell(size(isiCounts,1), 1);
                isiCenters(newK,:)     = cell(1, size(isiCenters,2));
                isiCenters(:,newK)     = cell(size(isiCenters,1), 1);
                isiViolations(newK,:)  = 0;
                isiViolations(:,newK)  = 0;
                DenCenters(newK,:)     = cell(1, size(DenCenters,2));
                DenCenters(:,newK)     = cell(size(DenCenters,1), 1);
                DenCounts(newK,:)      = cell(1, size(DenCounts,2));
                DenCounts(:,newK)      = cell(size(DenCounts,1), 1);
                presence_ratio(newK,:) = 0;
                presence_ratio(:,newK) = 0;

                numGroups   = newK;
                recomputeDisplayedGroups();

                sortedRes.unifiedLabels(rowsForNewLabel) = newLabel;
                sortedRes.channelNum(rowsForNewLabel)    = chCenter;

                nSplitDone = nSplitDone + 1;
                newSelected(end+1) = newK; %#ok<AGROW>
            end
        end

        if nSplitDone == 0
            try
                uialert(parentFig, ...
                    'No unit could be split (need >=50 spikes and a non-degenerate clustering).', ...
                    'Split');
            catch
            end
            return;
        end

        % Recompute per-unit metrics, refresh caches, rebuild page UI.
        try, preprocessSortedData(); catch, end
        try, refreshCCG();           catch, end
        try, refreshISI();           catch, end
        try, refreshDensity();       catch, end

        % Invalidate pair / unit walkers and the similarity matrix --
        % their stored indices reference the pre-split unit set, so
        % stepping through them now shows units with the wrong content.
        autoMergedPrimary    = [];
        autoMergedAbsorbed   = [];
        autoOverlapDroppedK  = [];
        autoOverlapTriggerJ  = [];
        autoCCGStripCont     = [];
        autoCCGStripClean    = [];
        lastAutoPair         = 0;
        rowGroup             = [];
        colGroup             = [];
        lastGroup            = 0;
        disimlarityScore     = sparse(numGroups, numGroups);
        ampSimilarity        = sparse(numGroups, numGroups);

        % Auto-tick the new sub-cluster slots so the user sees both
        % halves of every split unit side-by-side after the operation.
        for nk = newSelected
            if ~ismember(nk, selected_order)
                selected_order(end+1) = nk;
            end
        end

        % Rebuild page so the new rows / checkboxes / radios appear.
        currentPage = max(1, min(currentPage, totalPages));
        onChangePage(0);

        lastPlotted         = '';
        lastMultiPlotOrder  = [];
        selected_order_last = [];
        drawnow;
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

    function onEditNotes()
        % Edit the free-text note for the unit currently chosen by the
        % radio button (preferred, since it identifies a single unit
        % unambiguously). If no radio is selected, fall back to the
        % most recently ticked Select-checkbox unit.
        targetIdx = [];
        for i = 1:numel(groupRadioButtons)
            if isvalid(groupRadioButtons(i)) && isprop(groupRadioButtons(i),'Value') ...
                    && groupRadioButtons(i).Value
                targetIdx = i;
                break;
            end
        end
        if isempty(targetIdx)
            tickedIdx = find(chkBoxSelect, 1, 'last');
            if ~isempty(tickedIdx)
                targetIdx = tickedIdx;
            end
        end
        if isempty(targetIdx)
            try
                uialert(parentFig, ...
                    'Select a unit first (radio button or checkbox).', ...
                    'Notes');
            catch
            end
            return;
        end
        if numel(unitNotes) < targetIdx
            unitNotes(end+1:targetIdx) = {''};
        end
        existing = unitNotes{targetIdx};
        if isempty(existing), existing = ''; end
        prompt = sprintf('Note for unit G:%g (Ch %g):', ...
            groupList(targetIdx), channelList(targetIdx));
        try
            answer = inputdlg(prompt, 'Edit unit note', [3 70], {existing});
            if ~isempty(answer)
                unitNotes{targetIdx} = answer{1};
            end
        catch ME
            warning('Notes dialog failed: %s', ME.message);
        end
    end

    function onUndo()
        selectedIdx = find(chkBoxSelect);

        % If any of the ticked units are part of a merge (their
        % groupList value is shared with at least one other unit),
        % expand the undo set to include every constituent of that
        % merge. Otherwise undoing only the primary (or only the
        % absorbed) leaves the other side still wired into the
        % merged train, and metrics like ACG / ISI / CCG never come
        % back for the dissolved constituents.
        expandedIdx = selectedIdx(:);
        for k = selectedIdx(:)'
            gv = groupList(k);
            if ~isnan(gv)
                sameMerge = find(groupList == gv);
                if numel(sameMerge) > 1
                    expandedIdx = [expandedIdx; sameMerge(:)]; %#ok<AGROW>
                end
            end
        end
        expandedIdx = unique(expandedIdx);

        for k = expandedIdx'
            groupList(k) = originalGroupList(k);
            channelList(k) = originalChannelList(k);
            sampleWaveform(k,:,:) = originalSampleWaveform(k,:,:);
            sortedRes.unifiedLabels(originalUnifiedLabels==groupList(k)) = groupList(k);
            sortedRes.channelNum (originalUnifiedLabels==groupList(k)) = channelList(k);
            sortedRes.spike_idx(originalUnifiedLabels==groupList(k))   = originalSpikeIdx(originalUnifiedLabels==groupList(k));
            mergedFlag(k,1) = false;
            unitIsolation{k} = 'NA';
            stable_length(k,:) = [1, trialLength];

            % Auto-tick + auto-make-visible every dissolved
            % constituent so the user actually sees its ACG / ISI /
            % waveform after Undo. Without this, undoing just the
            % visible primary leaves the absorbed unit dissolved in
            % data but never re-rendered (selected_order still only
            % had the primary, so plotCCG never asks for the
            % absorbed's spike train).
            chkBoxSelect(k) = 1;
            chkBoxVis(k)    = 1;
            if ~ismember(k, selected_order)
                selected_order(end+1) = k; %#ok<AGROW>
            end
            if k <= numel(groupSelectCheckboxes) && ...
               isvalid(groupSelectCheckboxes(k)) && ...
               isprop(groupSelectCheckboxes(k),'Value')
                groupSelectCheckboxes(k).Value = true;
                groupVisCheckboxes(k).Value    = true;
                cu = getGroupColor(groupList(k));
                groupSelectCheckboxes(k).FontColor  = cu;
                groupSelectCheckboxes(k).FontWeight = 'Bold';
                groupVisCheckboxes(k).FontColor     = cu;
                groupVisCheckboxes(k).FontWeight    = 'Bold';
                if k <= numel(lblChannelHandles) && isvalid(lblChannelHandles(k))
                    lblChannelHandles(k).BackgroundColor = cu;
                    lblChannelHandles(k).FontColor       = bestTextColorFor(cu);
                end
            end
        end
        % Recompute the per-unit metrics on the restored labels so that
        % Group rows for previously-NaN units regain real ISI / rate.
        preprocessSortedData();

        % If any of the just-undone units are at the END of the array
        % AND were originally created by Split (originalGroupList is
        % NaN), trim them off entirely instead of leaving a ghost
        % NaN row hanging off the last page. We only trim trailing
        % slots so internal indices used by selected_order /
        % autoMergedPrimary etc. don't have to shift.
        trimmedAny = false;
        while numGroups > 0 ...
                && (numel(originalGroupList) >= numGroups) ...
                && isnan(groupList(numGroups)) ...
                && isnan(originalGroupList(numGroups))
            trimSplitSlot(numGroups);
            numGroups = numGroups - 1;
            trimmedAny = true;
        end
        if trimmedAny
            recomputeDisplayedGroups();
        end

        % channelList was just restored above for the dissolved units;
        % rebuild channelPlot inline so plotWaveforms / plotTrace can
        % resolve their footprints without going through
        % updateChannelPlot (which would force an early plotOption).
        channelPlot = channelList + (-numChannelPlot:numChannelPlot);
        % Inline cache wipes (in addition to refresh* calls below) so
        % the ACG / ISI / Density caches are guaranteed cleared even
        % if a refreshXxx variant is somehow not propagating its cell
        % reassignment up to the parent workspace.
        xcorr_vals     = cell(numGroups, numGroups);
        isiCounts      = cell(numGroups, numGroups);
        isiCenters     = cell(numGroups, numGroups);
        isiViolations  = zeros(numGroups, numGroups);
        DenCenters     = cell(numGroups, numGroups);
        DenCounts      = cell(numGroups, numGroups);
        presence_ratio = zeros(numGroups, numGroups);
        refreshDensity();
        refreshISI();
        refreshCCG();
        refreshLabels();
        if trimmedAny
            onChangePage(0);
        end
        % Force the next plotOption call to *fully* rebuild axChannels
        % (including the tiledlayout used by plotCCG / plotISI / etc.),
        % otherwise stale child objects from the pre-undo state can
        % leave the dissolved unit's ACG stuck on a flat row even when
        % the underlying xcorr_vals cache has been wiped. We also reset
        % selected_order_last so plotWaveforms re-snaps the y-window
        % and lastMultiPlotOrder so plotMultiple recreates its
        % multiPlotPanels (otherwise those handles point to children
        % of the just-deleted axChannels and allchild() errors).
        lastPlotted         = '';
        lastMultiPlotOrder  = [];
        selected_order_last = [];
        drawnow;
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
        % Limit only changes what's visible on the Density plot (the
        % shaded trim band) and what gets saved later. The Waveform /
        % Trace / CCG / ISI / Amplitude views don't reflect the trim,
        % so skip the heavy redraw on those plot types.
        if any(strcmpi(plotType, {'Density','Multiple'}))
            plotOption();
        end

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

    function onAutoCurate()
        % Surface a "running" status into the top Save label as soon as
        % Auto kicks off. The pipeline replaces the text with the merge
        % count when it finishes; this lets the user see Auto is still
        % working through Phases 1-5 without staring at a frozen GUI.
        if exist('lblSave','var') && isvalid(lblSave)
            lblSave.Text = 'Auto curation in progress...';
        end
        if exist('btnAutoCurate','var') && isvalid(btnAutoCurate)
            btnAutoCurate.Enable = 'off';
        end
        cleanupAuto = onCleanup(@() resetAutoButton());

        drawnow;

        % Reset to the original state first so each Auto run is
        % independent of any previous run (or manual edits). Skip
        % the confirmation dialog -- the user has just pressed
        % Auto, they don't need a second prompt.
        try, onReset(true); catch, end

        % Reset every Auto-pair audit trail so the Next/Prev walker on
        % the Auto panel starts at zero for whichever pair-type the
        % user has selected.
        autoMergedPrimary    = [];
        autoMergedAbsorbed   = [];
        autoOverlapDroppedK  = [];
        autoOverlapTriggerJ  = [];
        autoCCGStripCont     = [];
        autoCCGStripClean    = [];
        lastAutoPair         = 0;
        if exist('uifAutoPair','var') && isvalid(uifAutoPair)
            uifAutoPair.Value  = '0';
        end
        if exist('uifAutoPair2','var') && isvalid(uifAutoPair2)
            uifAutoPair2.Value = '/0';
        end

        % Automatic curation pipeline (removal -> cleaning -> merge,
        % with the necessary setup steps interleaved):
        %  1) Drop units that never fire (firingRate == 0).
        %  2) Classify SUA+ / SUA / MUA+ / MUA from a composite quality
        %     score that blends ISI %, ACG refractoriness ratio, and
        %     SNR (1 + detectability). Each component contributes
        %     independently so a unit can still be SUA when one factor
        %     is weak provided the others compensate.
        %  3) Run similarity analysis. yDistThr is forced to >=100 um
        %     for this run (a zero-distance setting would only compare
        %     same-y groups), then restored.
        %  4) Removal: drop units whose spike train is mostly a
        %     duplicate of a single neighbour's (>autoOverlapFrac
        %     coincidence within 0.5 ms) AND whose SNR < 1.5. Doing
        %     this BEFORE cleaning means the CCG-peak step doesn't
        %     waste effort on units we're about to throw away.
        %  5) Cleaning: per-pair CCG-peak cleanup. Pairs with a strong
        %     lag-0 peak strip coincident spikes from the noisier
        %     side; identity is preserved.
        %  6) Merge: walk candidate pairs and merge those that pass
        %     similarity, amplitude, PC, and merged-ISI gates.
        %  7) Apply the auto-limit to every surviving unit using the
        %     current Drop Rate field.

        % All gates come from the "Auto Curation" panel (closure
        % variables) so that Auto runs are independent of whatever the
        % user has typed in the Similarity Processing panel.
        thrIsiAbs    = autoIsiVal;
        thrIsiBudget = autoIsiBudget;
        thrPC        = autoPCVal;
        thrAmp       = autoAmpVal;
        thrSim       = autoSimVal;

        if strcmp(distanceEstType, 'XCorr')
            simBetter = @(s) s >= thrSim;
        else
            simBetter = @(s) s <= thrSim;
        end

        % Per-step counters for the end-of-Auto summary popup. The user
        % shouldn't have to read console output to know what changed.
        nAutoZeroDrops  = 0;     % zero-rate drops (Phase 1)
        nAutoCoincDrops = 0;     % coincidence-overlap drops (Phase 4)
        nAutoCCGStrips  = 0;     % CCG-peak coincident-spike strips (Phase 5)
        nAutoCCGDrops   = 0;     % CCG-peak full unit drops (Phase 5)

        % --- Phase 1: drop zero-firing-rate units ---------------------
        % Mark as removed but leave sortedRes.unifiedLabels alone so
        % the underlying spikes can still be visualised via effLabel.
        try
            for k = 1:numGroups
                if isnan(groupList(k)), continue; end
                if preprocessed.firingRate(k) == 0
                    groupList(k)      = NaN;
                    channelList(k)    = NaN;
                    unitIsolation{k}  = 'NA';
                    mergedFlag(k)     = true;
                    nAutoZeroDrops    = nAutoZeroDrops + 1;
                end
            end
        catch ME
            uialert(parentFig, sprintf('Phase 1 (drop zero-rate) failed: %s', ME.message), 'Auto Curation Error');
            return;
        end

        % --- Phase 2: composite-score isolation classification --------
        try
            classifyIsolationAll();
        catch ME
            uialert(parentFig, sprintf('Phase 2 (isolation classify) failed: %s', ME.message), 'Auto Curation Error');
            return;
        end

        % --- Phase 3: similarity analysis at autoDistVal um neighbourhood
        % Temporarily override yDistThr with the auto-curation Dist
        % value so onSimilarityEstimation uses the right radius, then
        % restore the user's previously typed value afterwards.
        prevYDist = yDistThr;
        if ~isempty(autoDistVal) && autoDistVal > 0
            yDistThr = autoDistVal;
        elseif yDistThr <= 0
            yDistThr = 100;
        end
        try
            onSimilarityEstimation();
        catch ME
            yDistThr = prevYDist;
            uialert(parentFig, sprintf('Phase 3 (similarity) failed: %s', ME.message), 'Auto Curation Error');
            return;
        end
        yDistThr = prevYDist;

        % --- Phase 4: removal (drop coincidence-noise units) ---------
        % For every non-SUA unit, count how many of its spikes fall
        % within 0.5 ms of ANY spike from a neighbour group whose
        % primary channel sits within autoDistVal um in y. If more
        % than autoOverlapFrac of the unit's spikes are coincident,
        % the unit is treated as a duplicate / shared-noise pickup
        % and set to NaN. Running removal BEFORE cleaning keeps the
        % CCG-peak phase from wasting cycles on units we're about
        % to throw away anyway.
        try
            coincSamples = round(0.5e-3 * cfg.samplingFrequency);
            coincFrac    = autoOverlapFrac;
            if ~isempty(autoDistVal) && autoDistVal > 0
                coincYUm = autoDistVal;
            elseif ~isempty(yDistThr) && yDistThr > 0
                coincYUm = yDistThr;
            else
                coincYUm = 200;
            end
            snrAll4 = 1 + detectblity;
            for k = 1:numGroups
                if isnan(groupList(k)),                              continue; end
                spk_k = sortedRes.spike_idx(sortedRes.unifiedLabels == groupList(k));
                if numel(spk_k) < 10, continue; end
                if isnan(channelList(k)), continue; end
                yk = ylocs(channelList(k));

                % Per-neighbour test: drop k as soon as ONE single
                % neighbour j alone accounts for > coincFrac of k's
                % spikes (within 0.5 ms) AND k itself is genuinely
                % low quality (SNR < 1.5). Without the SNR gate we
                % would NaN any unit that happens to share spikes
                % with a noisier neighbour -- including high-SNR
                % units the user explicitly wants to keep.
                shouldDrop = false;
                triggerJ   = NaN;
                threshold  = coincFrac * numel(spk_k);
                for j = 1:numGroups
                    if j == k || isnan(groupList(j)) || isnan(channelList(j)), continue; end
                    if abs(ylocs(channelList(j)) - yk) > coincYUm, continue; end
                    spk_j = sortedRes.spike_idx(sortedRes.unifiedLabels == groupList(j));
                    if isempty(spk_j), continue; end
                    [d_kj, ~] = nearest_distances(spk_k, spk_j);
                    if sum(d_kj <= coincSamples) > threshold
                        shouldDrop = true;
                        triggerJ   = j;
                        break;
                    end
                end

                snr_k = snrAll4(k);
                if isnan(snr_k), snr_k = 0; end
                if shouldDrop && snr_k < 1.5
                    groupList(k)     = NaN;
                    channelList(k)   = NaN;
                    unitIsolation{k} = 'NA';
                    mergedFlag(k)    = true;
                    nAutoCoincDrops  = nAutoCoincDrops + 1;
                    % Record (dropped unit, neighbour-that-caused-it)
                    % so the Auto-panel pair walker can revisit each
                    % overlap drop.
                    autoOverlapDroppedK(end+1,1) = k;          %#ok<AGROW>
                    autoOverlapTriggerJ(end+1,1) = triggerJ;   %#ok<AGROW>
                end
            end
            % counts/labels were rewritten; refresh metrics + UI labels
            preprocessSortedData();
            refreshLabels();
        catch ME
            uialert(parentFig, sprintf('Phase 4 (removal) failed: %s', ME.message), 'Auto Curation Error');
        end

        % --- Phase 5: cleaning (per-pair CCG-peak cleanup) -----------
        % Trigger criterion: the cross-correlogram has a peak AT lag 0
        % whose height is the GLOBAL MAX of the CCG and exceeds 2x the
        % 90th percentile of all the other bins. That's a stricter
        % "true peak at zero" test than a generic Poisson null --
        % structural temporal coupling at non-zero lags (real biology)
        % won't trip it because their global max is away from zero.
        %
        % CCG bins: 1 ms wide, ±100 ms half-window. We compute it via
        % a two-pointer pairwise-difference histogram so we don't have
        % to allocate a num_Samples-long binary vector per pair (the
        % old binary_xcorr approach).
        %
        % If triggered, pick the lower-SNR side as the contaminated
        % unit and strip its spikes that are coincident (within
        % 0.5 ms) with the cleaner unit by setting their
        % sortedRes.unifiedLabels rows to -1. The higher-SNR (cleaner)
        % unit is never touched, and identity is preserved (cleanup
        % never drops a class, only spikes).
        try
            coincSamples = round(0.5e-3 * cfg.samplingFrequency);
            if ~isempty(autoDistVal) && autoDistVal > 0
                coincYUm = autoDistVal;
            elseif ~isempty(yDistThr) && yDistThr > 0
                coincYUm = yDistThr;
            else
                coincYUm = 200;
            end
            snrAll = 1 + detectblity;

            % CCG parameters. ccgPeakRatio is wired to the Auto
            % panel's CCG field so the user can dial cleanup
            % strictness without touching the source.
            ccgBin       = round(1e-3 * cfg.samplingFrequency);  % 1 ms
            ccgMaxLag    = 100 * ccgBin;                          % ±100 ms
            ccgZeroTolBins = 1;                                    % "at zero" = ±1 ms
            ccgPeakRatio = autoCCGRatio;                          % peak > ratio * baseline median

            for k = 1:numGroups
                if isnan(groupList(k)) || isnan(channelList(k)), continue; end
                yk = ylocs(channelList(k));

                for j = (k+1):numGroups
                    if isnan(groupList(j)) || isnan(channelList(j)), continue; end
                    if groupList(k) == groupList(j), continue; end
                    if abs(ylocs(channelList(j)) - yk) > coincYUm, continue; end

                    spk_k_rows = find(sortedRes.unifiedLabels == groupList(k));
                    spk_j_rows = find(sortedRes.unifiedLabels == groupList(j));
                    spk_k = sortedRes.spike_idx(spk_k_rows);
                    spk_j = sortedRes.spike_idx(spk_j_rows);
                    if numel(spk_k) < 10 || numel(spk_j) < 10, continue; end

                    % Build the CCG histogram and run the peak test.
                    [ccg, zeroBinIdx] = pairCCG(spk_k, spk_j, ccgMaxLag, ccgBin);
                    if ~ccgPeakAtZero(ccg, zeroBinIdx, ccgZeroTolBins, ccgPeakRatio)
                        continue;
                    end

                    % Coincidence mask for the strip phase. Use
                    % nearest_distances to find each k spike's tightest
                    % j neighbour (and vice versa) within ±0.5 ms.
                    [d1_kj, d2_kj] = nearest_distances(spk_k, spk_j);
                    coincMask_k    = min(d1_kj, d2_kj) <= coincSamples;

                    % CCG peak triggered. Only clean similar-quality
                    % pairs: a large SNR gap belongs to overlap removal,
                    % not spike stripping. Within the pair, strip the
                    % lower firing-rate (fewer-spike) side -- the likely
                    % fragment. NaN SNR is treated as 0.
                    sk = snrAll(k); if isnan(sk), sk = 0; end
                    sj = snrAll(j); if isnan(sj), sj = 0; end
                    if abs(sk - sj) >= autoCCGSnrDiff, continue; end
                    if numel(spk_k) <= numel(spk_j)
                        contIdx    = k;
                        cleanIdx   = j;
                        cont_rows  = spk_k_rows;
                        cont_coinc = coincMask_k;
                        signedLag  = autoSignedNearestLag(spk_k, spk_j);
                    else
                        contIdx    = j;
                        cleanIdx   = k;
                        cont_rows  = spk_j_rows;
                        [d1_jk, d2_jk] = nearest_distances(spk_j, spk_k);
                        cont_coinc     = min(d1_jk, d2_jk) <= coincSamples;
                        signedLag  = autoSignedNearestLag(spk_j, spk_k);
                    end

                    if ~any(cont_coinc), continue; end

                    % --- Timing-consistency gate -----------------
                    % The coincident spikes' signed lags to the
                    % nearest clean spike must cluster tightly around
                    % their median. A real duplicate (same physical
                    % spike picked up by both units with a small
                    % sorter-induced offset) gives a sharp peak;
                    % random coincidences spread the lag distribution
                    % uniformly inside the +-0.5 ms window. The strip
                    % is applied only to the consistent-lag subset;
                    % if fewer than autoLagConsistencyFrac of the
                    % coincident lags fall inside +-autoLagTightMs of
                    % the median, skip the pair entirely.
                    tightSamples = max(1, round(autoLagTightMs * 1e-3 * cfg.samplingFrequency));
                    coincIdx     = find(cont_coinc);
                    coincLags    = signedLag(coincIdx);
                    coincLags    = coincLags(~isnan(coincLags));
                    if numel(coincLags) < 5
                        continue;
                    end
                    medLag           = median(coincLags);
                    % Median lag must also sit near zero; a non-zero
                    % median is a synaptic delay (real coupling), not a
                    % duplicate.
                    if abs(medLag) > tightSamples, continue; end
                    withinTight      = abs(coincLags - medLag) <= tightSamples;
                    consistencyFrac  = sum(withinTight) / numel(coincLags);
                    if consistencyFrac < autoLagConsistencyFrac
                        continue;
                    end
                    keepCoinc            = false(numel(signedLag), 1);
                    keepCoinc(coincIdx)  = withinTight;
                    cont_coinc           = cont_coinc & keepCoinc;
                    if ~any(cont_coinc), continue; end

                    % Strip-fraction cap: never gut the fragment.
                    if sum(cont_coinc) / numel(cont_rows) > autoCCGMaxStripFrac
                        continue;
                    end

                    % Strip the contaminated unit's coincident rows
                    % from sortedRes.unifiedLabels (set to -1).
                    sortedRes.unifiedLabels(cont_rows(cont_coinc)) = -1;
                    nAutoCCGStrips = nAutoCCGStrips + 1;
                    autoCCGStripCont(end+1,1)  = contIdx;  %#ok<AGROW>
                    autoCCGStripClean(end+1,1) = cleanIdx; %#ok<AGROW>

                    % Refresh the contaminated unit's ISI / firing-rate
                    % immediately so the GUI label and any downstream
                    % gate sees the cleaned numbers, not the pre-strip
                    % ones. (preprocessSortedData runs again at the
                    % end of the phase to catch every unit.) Cleaning
                    % only strips spikes; identity is preserved so the
                    % later merge phase sees the cleaned train.
                    spk_c_after = sortedRes.spike_idx(sortedRes.unifiedLabels == groupList(contIdx));
                    if numel(spk_c_after) >= 2
                        [~, ~, preprocessed.isiViolation(contIdx)] = ...
                            getISIViolations(spk_c_after, cfg.samplingFrequency, thresholdISI);
                    else
                        preprocessed.isiViolation(contIdx) = 0;
                    end
                    preprocessed.firingRate(contIdx)   = numel(spk_c_after) / max(trialLength, eps);
                    preprocessed.logFiringRate(contIdx) = log(max(preprocessed.firingRate(contIdx), eps));
                end
            end
            preprocessSortedData();
            refreshLabels();
        catch ME
            uialert(parentFig, sprintf('Phase 5 (cleaning) failed: %s', ME.message), 'Auto Curation Error');
        end

        % --- Phase 6: gated automatic merges -------------------------
        % Walk candidate pairs from the simMat that survived cleanup
        % and removal, and merge those that pass similarity, amplitude,
        % PC, and merged-ISI gates. simMat is keyed by group index;
        % NaN groupList entries are skipped so dropped units don't
        % participate even if they appear in simMat from Phase 3.
        nMerged = 0;
        try
            if ~isempty(disimlarityScore)
                simMat = full(disimlarityScore);
                ampMat = full(ampSimilarity);
                % Defensive: similarity matrices may be from a previous
                % numGroups (e.g. before a Split). If their shapes don't
                % match, snap both to a common dim and to current numGroups
                % so the indexing into groupList stays valid.
                if ~isequal(size(simMat), size(ampMat)) || ...
                        size(simMat,1) ~= size(simMat,2) || ...
                        size(simMat,1) ~= numGroups
                    nMin = min([size(simMat,1), size(simMat,2), ...
                                size(ampMat,1), size(ampMat,2), numGroups]);
                    if nMin < 1, nMin = 1; end
                    simMat = simMat(1:nMin, 1:nMin);
                    ampMat = ampMat(1:nMin, 1:nMin);
                end
                if any(simMat(:) ~= 0)
                    [rIdx, cIdx] = find(simBetter(simMat) & ampMat <= thrAmp);
                    keep = rIdx > cIdx;
                    rIdx = rIdx(keep); cIdx = cIdx(keep);
                    for p = 1:numel(rIdx)
                        gA = cIdx(p); gB = rIdx(p);
                        if gA > numel(groupList) || gB > numel(groupList), continue; end
                        if isnan(groupList(gA)) || isnan(groupList(gB)), continue; end
                        if groupList(gA) == groupList(gB), continue; end
                        if ~checkPCDistance(gA, gB, thrPC), continue; end
                        if ~checkMergedISI(gA, gB, thrIsiAbs, thrIsiBudget), continue; end

                        if preprocessed.firingRate(gA) >= preprocessed.firingRate(gB)
                            primary = gA; absorbed = gB;
                        else
                            primary = gB; absorbed = gA;
                        end
                        mergeLabel   = groupList(primary);
                        mergeChannel = channelList(primary);
                        absLab       = groupList(absorbed);

                        % Align absorbed spike times to primary frame
                        % using shared-channel XCorr lag. Done for every
                        % merge regardless of the similarity metric the
                        % user selected. originalSpikeIdx is untouched
                        % so Reset / Undo restore the original timing.
                        if ~isnan(mergeChannel)
                            mwA = localMeanWaveform(primary, mergeChannel);
                            mwB = localMeanWaveform(absorbed, mergeChannel);
                            if ~isempty(mwA) && ~isempty(mwB) && ...
                               numel(mwA) == numel(mwB) && numel(mwA) > 4
                                M2 = numel(mwA);
                                lagCap = max(1, round(M2/4));
                                [~, bestLag] = max_half_corr(mwA(:)', mwB(:)', 1, M2, lagCap, 0);
                                if isfinite(bestLag) && bestLag ~= 0
                                    shiftRows = find(sortedRes.unifiedLabels == absLab);
                                    shifted   = sortedRes.spike_idx(shiftRows) - bestLag;
                                    shifted(shifted < 1)            = 1;
                                    shifted(shifted > num_Samples)  = num_Samples;
                                    sortedRes.spike_idx(shiftRows)  = shifted;
                                end
                            end
                        end

                        % Find every unit currently sharing the
                        % absorbed label -- this is what fixes the
                        % transitive-merge bug. With chain merges like
                        % (3->5) then (3->9), the second merge would
                        % otherwise leave unit 5 stranded on its old
                        % label even though 3, 7, 9 all share the new
                        % one. Updating every member of the absorbed
                        % cluster keeps the whole transitive group on
                        % the same label.
                        absMembers = find(groupList == absLab);
                        sortedRes.unifiedLabels(sortedRes.unifiedLabels == absLab) = mergeLabel;
                        sortedRes.channelNum(sortedRes.unifiedLabels == absLab)   = mergeChannel;
                        groupList(absMembers)   = mergeLabel;
                        channelList(absMembers) = mergeChannel;
                        mergedFlag(absMembers)  = true;
                        nMerged = nMerged + 1;

                        % Record this merge for the Auto-panel pair
                        % walker. We store the (primary, absorbed)
                        % index pair so the user can step through every
                        % merge and inspect both constituents.
                        autoMergedPrimary(end+1,1)  = primary;  %#ok<AGROW>
                        autoMergedAbsorbed(end+1,1) = absorbed; %#ok<AGROW>
                    end
                end
            end
        catch ME
            uialert(parentFig, sprintf('Phase 6 (merge) failed: %s', ME.message), 'Auto Curation Error');
            return;
        end

        % Recompute per-unit metrics now that some labels collapsed.
        try
            preprocessSortedData();
            if nMerged > 0
                % Phase 2 already classified before any merges. Only
                % re-classify if Phase 6 actually changed any spike
                % trains; otherwise the existing labels are still
                % correct and we save a numGroups-sized loop.
                classifyIsolationAll();
            end
            % Keep channelPlot in sync with the post-merge channelList
            % (avoids calling updateChannelPlot which would also force
            % an early plotOption redraw mid-pipeline).
            channelPlot = channelList + [-numChannelPlot:numChannelPlot];
            refreshLabels();
            refreshDensity();
            refreshISI();
            refreshCCG();
        catch ME
            uialert(parentFig, sprintf('Recompute failed: %s', ME.message), 'Auto Curation Error');
        end

        % --- Phase 7: auto-limit each surviving unit (inline) --------
        % Compute presence-ratio density per unit and look for a
        % dropRate-fold rate change with findchangepts; trim
        % stable_length to the stable epoch. Done inline so we do not
        % trigger plotOption over hundreds of groups via onAutoCut.
        densityBinSize = max(1, round(trialLength/(2*ccgLag)));
        for k = 1:numGroups
            if isnan(groupList(k)), continue; end
            stable_length(k,:) = [1, trialLength];
            spk_k = unique(sortedRes.spike_idx( ...
                sortedRes.unifiedLabels == groupList(k)));
            if numel(spk_k) < 10, continue; end
            try
                [Dc_centers, Dc_counts, ~] = presenceRatio( ...
                    spk_k, cfg.samplingFrequency, trialLength, ...
                    densityBinSize, smoothN);
                % Do NOT store into DenCenters / DenCounts here -- those
                % cells are indexed by LABEL elsewhere, and caching them
                % under the group INDEX would later trip plotDensity
                % into reading a stale (zeroed) presence_ratio scalar.
                if numel(Dc_counts) >= 3
                    [sI, eI] = limitStableInterval(Dc_counts, dropRate);
                    xs = Dc_centers(sI); xe = Dc_centers(eI);
                    span = max(Dc_centers(end), eps);
                    stable_length(k,1) = max(1, min(trialLength, ...
                        trialLength * (xs - Dc_centers(1)) / span));
                    stable_length(k,2) = max(1, min(trialLength, ...
                        trialLength * xe / span));
                end
            catch
                % skip units whose density cannot be computed
            end
        end

        plotOption();

        if exist('lblSave','var') && isvalid(lblSave)
            lblSave.Text = sprintf('Auto curation: %d merges', nMerged);
        end
        fprintf('Auto curation: %d pairs merged.\n', nMerged);

        % Surface the audit count for whichever pair-type the user
        % currently has selected so the walker is ready to browse it.
        if exist('uifAutoPair2','var') && isvalid(uifAutoPair2)
            [aArr, ~] = currentAutoPairList();
            uifAutoPair2.Value = sprintf('/%d', numel(aArr));
        end

        % Plain-language summary popup for the user. No "Phase X"
        % jargon -- just what changed, ordered to match the
        % removal -> cleaning -> merge pipeline.
        summaryMsg = sprintf([ ...
            'Auto curation finished.\n\n', ...
            '  Removed (zero rate):    %d\n', ...
            '  Removed (overlap):      %d\n', ...
            '  Cleaned (CCG peak):     %d spike-strips\n', ...
            '  Merged pairs:           %d\n'], ...
            nAutoZeroDrops, nAutoCoincDrops, nAutoCCGStrips, ...
            nMerged);
        try
            uialert(parentFig, summaryMsg, 'Auto curation complete', ...
                'Icon', 'success');
        catch
            disp(summaryMsg);
        end
    end

    function resetAutoButton()
        % Re-enable the Auto button when onAutoCurate finishes (or
        % errors out). Wired via onCleanup so we never leave the
        % button stuck in disabled state.
        if exist('btnAutoCurate','var') && isvalid(btnAutoCurate)
            btnAutoCurate.Enable = 'on';
        end
    end

    function resetSplitButton(prevSaveLbl)
        % Re-enable the Split button when onSplit finishes / errors.
        if exist('btnSplit','var') && isvalid(btnSplit)
            btnSplit.Enable = 'on';
            if isprop(btnSplit,'Text')
                btnSplit.Text = 'Split';
            end
        end
        if exist('lblSave','var') && isvalid(lblSave) && nargin > 0
            try
                lblSave.Text = prevSaveLbl;
            catch
            end
        end
    end

    function [ccg, zeroBinIdx] = pairCCG(spkA, spkB, maxLagSamples, binSamples)
        % Build a 1-D cross-correlogram histogram of (spkB - spkA)
        % differences inside ±maxLagSamples, binned at binSamples
        % resolution. Uses a two-pointer sweep so total work is
        % O(nA + nB + #pairs-in-range) rather than O(num_Samples).
        spkA = sort(spkA(:));
        spkB = sort(spkB(:));
        nA   = numel(spkA);
        nB   = numel(spkB);
        halfBins   = floor(maxLagSamples / binSamples);
        nBins      = 2*halfBins + 1;
        zeroBinIdx = halfBins + 1;
        ccg        = zeros(nBins, 1);
        if nA == 0 || nB == 0, return; end

        % Pre-compute every pairwise difference within range using
        % two pointers. For each k spike we keep a moving window
        % [j_lo, j_hi] over spkB.
        b_lo = 1; b_hi = 0;
        diffsCell = cell(nA, 1);
        for i = 1:nA
            t = spkA(i);
            while b_lo <= nB && spkB(b_lo) < t - maxLagSamples
                b_lo = b_lo + 1;
            end
            while b_hi < nB && spkB(b_hi + 1) <= t + maxLagSamples
                b_hi = b_hi + 1;
            end
            if b_lo <= b_hi
                diffsCell{i} = spkB(b_lo:b_hi) - t;
            end
        end
        if all(cellfun(@isempty, diffsCell)), return; end
        allDiffs = vertcat(diffsCell{:});
        if isempty(allDiffs), return; end
        binIdx = round(double(allDiffs) / binSamples) + zeroBinIdx;
        keep   = binIdx >= 1 & binIdx <= nBins;
        binIdx = binIdx(keep);
        if isempty(binIdx), return; end
        ccg = accumarray(binIdx(:), 1, [nBins, 1]);
    end

    function tf = ccgPeakAtZero(ccg, zeroBinIdx, zeroTolBins, ratio)
        % True iff ALL of these hold:
        %   * the CCG has a peak in the lag-zero window
        %     (zeroBinIdx +/- zeroTolBins -- so with the default
        %     zeroTolBins=1, that's bins -1, 0, +1 in lag units)
        %     whose value is the (joint) global max of the whole
        %     CCG. Computing it as "max-in-window must equal
        %     max-anywhere" handles ties cleanly: if two non-zero
        %     bins are tied for the top, MATLAB's max(...) would
        %     only point at the first one, so a tie split between
        %     a zero-window bin and an off-zero bin used to slip
        %     past undetected -- now both have to be present in
        %     the window for the test to pass.
        %   * that peak is greater than `ratio` times the baseline,
        %     where baseline = MEDIAN of the *non-peak* bins (the
        %     rest of the CCG outside the lag-zero window). Using
        %     the median instead of the mean makes the baseline
        %     robust to a few high off-zero bins (real biological
        %     coupling, edge artefacts) -- a couple of outlier bins
        %     no longer drag the baseline up far enough to mask
        %     a genuine duplicate peak.
        % If the baseline is all-zero we don't trigger -- there's
        % nothing to compare the peak against, so conservatively
        % treat that as "not a real peak".
        tf = false;
        if isempty(ccg) || all(ccg == 0), return; end
        nB         = numel(ccg);
        peakRange  = max(1, zeroBinIdx-zeroTolBins) : min(nB, zeroBinIdx+zeroTolBins);
        if isempty(peakRange), return; end
        peakValue  = max(ccg(peakRange));
        offRange   = ccg;
        offRange(peakRange) = [];
        if isempty(offRange)
            return;     % CCG is entirely inside the peak window
        end
        % Off-window must NOT have a strictly larger bin -- if it
        % does, the actual global max sits away from zero and this
        % isn't a "true peak at zero" pair.
        if any(offRange > peakValue), return; end
        if ~any(offRange > 0)
            return;     % no baseline activity -> can't compare
        end
        baseline = median(offRange);
        % Sparse CCGs (most off-zero bins are zero) end up with
        % median = 0, which would make the test trivially true for
        % any non-zero peak. Fall back to the mean of non-zero
        % off-window bins so the baseline still reflects the
        % typical non-zero rate.
        if baseline <= 0
            nz = offRange(offRange > 0);
            if isempty(nz), return; end
            baseline = mean(nz);
            if baseline <= 0, return; end
        end
        tf = peakValue > ratio * baseline;
    end

    function lag = autoSignedNearestLag(A, B)
        % For each element of A, signed lag (samples) to its nearest
        % neighbour in B: lag(i) = B_nearest - A(i). Positive means
        % the B spike comes after the A spike. NaN when B is empty.
        % Used by the Phase 5 timing-consistency gate so we can tell
        % whether the coincident spikes cluster around a specific
        % sub-ms offset (real duplicate) or are scattered (random).
        A = A(:); B = B(:);
        nA = numel(A); nB = numel(B);
        lag = nan(nA, 1);
        if nA == 0 || nB == 0, return; end
        [Bs, ~]    = sort(B);
        [As, srtA] = sort(A);
        j = 1;
        for ii = 1:nA
            while j < nB && Bs(j+1) < As(ii)
                j = j + 1;
            end
            bestDiff = Bs(j) - As(ii);
            bestAbs  = abs(bestDiff);
            if j+1 <= nB
                d2 = Bs(j+1) - As(ii);
                if abs(d2) < bestAbs
                    bestDiff = d2; bestAbs = abs(d2);
                end
            end
            if j-1 >= 1
                d2 = Bs(j-1) - As(ii);
                if abs(d2) < bestAbs
                    bestDiff = d2; bestAbs = abs(d2);
                end
            end
            lag(srtA(ii)) = bestDiff;
        end
    end

    function trimSplitSlot(idx)
        % Remove the unit at index `idx` from every per-unit array.
        % Caller guarantees idx is at the end of the array (so other
        % indices don't have to shift) and that the slot is a
        % dissolved Split-created unit (groupList = NaN AND
        % originalGroupList = NaN). Pure structural cleanup --
        % does NOT touch sortedRes.* (the spike rows are restored
        % independently via originalUnifiedLabels in onUndo).
        if idx < 1, return; end

        % Drop scalar / vector arrays (length numGroups).
        if numel(groupList)        >= idx, groupList(idx)        = []; end
        if numel(channelList)      >= idx, channelList(idx)      = []; end
        if numel(mergedFlag)       >= idx, mergedFlag(idx)       = []; end
        if numel(unitIsolation)    >= idx, unitIsolation(idx)    = []; end
        if numel(mainPolarity)     >= idx, mainPolarity(idx)     = []; end
        if numel(sidePolarity)     >= idx, sidePolarity(idx)     = []; end
        if numel(detectblity)      >= idx, detectblity(idx)      = []; end
        if numel(chkBoxSelect)     >= idx, chkBoxSelect(idx)     = []; end
        if numel(chkBoxVis)        >= idx, chkBoxVis(idx)        = []; end
        if numel(unitNotes)        >= idx, unitNotes(idx)        = []; end
        if numel(originalGroupList)   >= idx, originalGroupList(idx)   = []; end
        if numel(originalChannelList) >= idx, originalChannelList(idx) = []; end

        % Per-unit struct arrays inside `preprocessed`.
        if isfield(preprocessed,'firingRate') && numel(preprocessed.firingRate) >= idx
            preprocessed.firingRate(idx) = [];
        end
        if isfield(preprocessed,'isiViolation') && numel(preprocessed.isiViolation) >= idx
            preprocessed.isiViolation(idx) = [];
        end
        if isfield(preprocessed,'logFiringRate') && numel(preprocessed.logFiringRate) >= idx
            preprocessed.logFiringRate(idx) = [];
        end

        % 2-D rows.
        if size(stable_length,1)            >= idx, stable_length(idx,:)          = []; end
        if size(channelPlot,1)              >= idx, channelPlot(idx,:)            = []; end
        if size(sampleWaveform,1)           >= idx, sampleWaveform(idx,:,:)       = []; end
        if size(originalSampleWaveform,1)   >= idx, originalSampleWaveform(idx,:,:) = []; end

        % NxN caches: drop both the row and the column at idx.
        if size(xcorr_vals,1)    >= idx, xcorr_vals(idx,:)    = []; end
        if size(xcorr_vals,2)    >= idx, xcorr_vals(:,idx)    = []; end
        if size(isiCounts,1)     >= idx, isiCounts(idx,:)     = []; end
        if size(isiCounts,2)     >= idx, isiCounts(:,idx)     = []; end
        if size(isiCenters,1)    >= idx, isiCenters(idx,:)    = []; end
        if size(isiCenters,2)    >= idx, isiCenters(:,idx)    = []; end
        if size(isiViolations,1) >= idx, isiViolations(idx,:) = []; end
        if size(isiViolations,2) >= idx, isiViolations(:,idx) = []; end
        if size(DenCenters,1)    >= idx, DenCenters(idx,:)    = []; end
        if size(DenCenters,2)    >= idx, DenCenters(:,idx)    = []; end
        if size(DenCounts,1)     >= idx, DenCounts(idx,:)     = []; end
        if size(DenCounts,2)     >= idx, DenCounts(:,idx)     = []; end
        if size(presence_ratio,1)>= idx, presence_ratio(idx,:)= []; end
        if size(presence_ratio,2)>= idx, presence_ratio(:,idx)= []; end

        % Drop / clean up UI handles for this slot. Nested functions
        % share the parent's workspace so we can mutate the handle
        % arrays directly.
        if numel(groupSelectCheckboxes) >= idx
            if isvalid(groupSelectCheckboxes(idx))
                delete(groupSelectCheckboxes(idx));
            end
            groupSelectCheckboxes(idx) = [];
        end
        if numel(groupVisCheckboxes) >= idx
            if isvalid(groupVisCheckboxes(idx))
                delete(groupVisCheckboxes(idx));
            end
            groupVisCheckboxes(idx) = [];
        end
        if numel(groupRadioButtons) >= idx
            if isvalid(groupRadioButtons(idx))
                delete(groupRadioButtons(idx));
            end
            groupRadioButtons(idx) = [];
        end
        if exist('lblGroupHandles','var') && numel(lblGroupHandles) >= idx
            if isvalid(lblGroupHandles(idx)), delete(lblGroupHandles(idx)); end
            lblGroupHandles(idx) = [];
        end
        if exist('lblIdHandles','var') && numel(lblIdHandles) >= idx
            if isvalid(lblIdHandles(idx)), delete(lblIdHandles(idx)); end
            lblIdHandles(idx) = [];
        end
        if exist('lblChannelHandles','var') && numel(lblChannelHandles) >= idx
            if isvalid(lblChannelHandles(idx)), delete(lblChannelHandles(idx)); end
            lblChannelHandles(idx) = [];
        end
        if exist('lblIsolation','var') && numel(lblIsolation) >= idx
            if isvalid(lblIsolation(idx)), delete(lblIsolation(idx)); end
            lblIsolation(idx) = [];
        end
        if exist('lblRate','var') && numel(lblRate) >= idx
            if isvalid(lblRate(idx)), delete(lblRate(idx)); end
            lblRate(idx) = [];
        end
        if exist('lblDetectablity','var') && numel(lblDetectablity) >= idx
            if isvalid(lblDetectablity(idx)), delete(lblDetectablity(idx)); end
            lblDetectablity(idx) = [];
        end
        if exist('lblISIViolation','var') && numel(lblISIViolation) >= idx
            if isvalid(lblISIViolation(idx)), delete(lblISIViolation(idx)); end
            lblISIViolation(idx) = [];
        end

        % selected_order: drop idx (no shifting; we only trim trailing).
        selected_order(selected_order == idx) = [];

        % All Auto-pair audit trails: drop entries that touched idx.
        % Any of the four lists may be non-empty after an Auto run.
        if ~isempty(autoMergedPrimary)
            keep = (autoMergedPrimary ~= idx) & (autoMergedAbsorbed ~= idx);
            autoMergedPrimary  = autoMergedPrimary(keep);
            autoMergedAbsorbed = autoMergedAbsorbed(keep);
        end
        if ~isempty(autoOverlapDroppedK)
            keep = (autoOverlapDroppedK ~= idx) & (autoOverlapTriggerJ ~= idx);
            autoOverlapDroppedK = autoOverlapDroppedK(keep);
            autoOverlapTriggerJ = autoOverlapTriggerJ(keep);
        end
        if ~isempty(autoCCGStripCont)
            keep = (autoCCGStripCont ~= idx) & (autoCCGStripClean ~= idx);
            autoCCGStripCont  = autoCCGStripCont(keep);
            autoCCGStripClean = autoCCGStripClean(keep);
        end
    end

    function classifyIsolationAll()
        % Composite ISI / ACG / SNR score; weights 0.4 / 0.4 / 0.2 so
        % a clean ISI+ACG can earn SUA even with modest SNR.
        snrAll = 1 + detectblity;
        for k = 1:numGroups
            if isnan(groupList(k))
                unitIsolation{k} = 'NA';
                continue;
            end
            isi      = preprocessed.isiViolation(k);
            sr       = snrAll(k);
            acgRatio = acgRefractoryRatio(k);
            if isnan(isi),      isi = 2.0;   end
            if isnan(sr),       sr  = 1.0;   end
            if isnan(acgRatio), acgRatio = 1; end

            sIsi  = max(0, 1 - isi / 2.0);
            sSnr  = max(0, min(1, (sr - 1) / 4));
            sAcg  = max(0, 1 - acgRatio / 0.10);
            score = 0.4*sIsi + 0.4*sAcg + 0.2*sSnr;

            if score >= 0.85 && isi < 0.5 && acgRatio < 0.05
                unitIsolation{k} = 'SUA+';
            elseif score >= 0.65 && isi < 1.0
                unitIsolation{k} = 'SUA';
            elseif score >= 0.45
                unitIsolation{k} = 'MUA+';
            else
                unitIsolation{k} = 'MUA';
            end
        end
    end


    function [startIdx, endIdx] = limitStableInterval(d, dropFactor)
        % Identical heuristic to onAutoCut's findStableInterval:
        % findchangepts -> if 1 cp, trim the side whose mean rate is
        % more than dropFactor lower than the other; if 2+ cps, walk
        % them and tighten startIdx / endIdx whenever the same drop
        % criterion is met (also confirmed against the next-segment
        % mean to avoid trimming a transient dip). Without this
        % multi-cp branch, Phase 5 silently never trims units that
        % have a drop+recovery+drop shape.
        n = numel(d);
        startIdx = 1; endIdx = n;
        if n < 3, return; end
        cp = findchangepts(d, 'Statistic', 'rms', 'MinThreshold', 5);
        if isscalar(cp)
            if mean(d(1:cp)) < mean(d(cp+1:end)) / dropFactor
                startIdx = cp + 1; endIdx = n;
            elseif mean(d(1:cp)) / dropFactor > mean(d(cp+1:end))
                startIdx = 1; endIdx = cp;
            end
        elseif numel(cp) >= 2
            for ci = 1:numel(cp)
                if mean(d(1:cp(ci))) < mean(d(cp(ci)+1:end))/dropFactor && ...
                   mean(d(1:cp(ci))) < mean(d(cp(ci)+1:min(2*cp(ci),n)))/dropFactor
                    startIdx = cp(ci)+1;
                elseif mean(d(startIdx:cp(ci)))/dropFactor > mean(d(cp(ci)+1:end)) && ...
                       mean(d(startIdx:cp(ci)))/dropFactor > mean(d(cp(ci)+1:min(2*cp(ci),n)))
                    endIdx = cp(ci);
                    break
                end
            end
        end
    end

    function setOverlapVal(val)
        if ~isempty(val) && isfinite(val) && val >= 0 && val <= 1
            autoOverlapFrac = val;
        end
    end

    function setAutoDistVal(val)
        if ~isempty(val) && isfinite(val) && val >= 0
            autoDistVal = val;
        end
    end

    function setAutoSimVal(val)
        if ~isempty(val) && isfinite(val) && val >= 0 && val <= 1
            autoSimVal = val;
        end
    end

    function setAutoAmpVal(val)
        if ~isempty(val) && isfinite(val) && val >= 0 && val <= 1
            autoAmpVal = val;
        end
    end

    function setAutoIsiVal(val)
        if ~isempty(val) && isfinite(val) && val >= 0
            autoIsiVal = val;
        end
    end

    function setAutoPCVal(val)
        if ~isempty(val) && isfinite(val) && val >= 0
            autoPCVal = val;
        end
    end

    function setAutoBudgetVal(val)
        if ~isempty(val) && isfinite(val) && val >= 1
            autoIsiBudget = val;
        end
    end

    function setAutoCCGRatio(val)
        % Auto-cleanup peak/baseline-median ratio. Clamped to
        % [0.5, 10] -- below 0.5 the test would trigger on the
        % median itself (i.e. nothing); above 10 even the most
        % obvious shared-spike pairs would slip past.
        if ~isempty(val) && isfinite(val)
            autoCCGRatio = max(0.5, min(10, val));
        end
    end

    function setSplitMethod(val)
        if ~isempty(val) && (ischar(val) || isstring(val))
            splitMethod = char(val);
        end
    end

    function setSplitNumClusters(val)
        if ~isempty(val) && isfinite(val) && val >= 2
            splitNumClusters = max(2, round(val));
        end
    end

    function lbl = effLabel(k)
        % Effective label of group k for plotting / spike lookup. When
        % a unit is marked as removed (groupList(k) = NaN) we still
        % want to be able to visualise it, so fall back to the unit's
        % pre-curation label in originalGroupList.
        if k >= 1 && k <= numGroups && ~isnan(groupList(k))
            lbl = groupList(k);
        elseif k >= 1 && k <= numel(originalGroupList)
            lbl = originalGroupList(k);
        else
            lbl = NaN;
        end
    end

    function lbls = effLabelVec(ks)
        % Vectorised version of effLabel — avoids the arrayfun /
        % function-call overhead when resolving the effective labels of
        % many indices at once (called every plotWaveforms / plotCCG /
        % plotISI redraw for every selected unit).
        ks   = ks(:);
        lbls = nan(size(ks));
        valid = ks >= 1 & ks <= numGroups;
        if any(valid)
            kv         = ks(valid);
            gl         = groupList(kv);
            useOrig    = isnan(gl);
            % primary path: groupList value when not NaN
            lbls(valid) = gl;
            % NaN-marked units fall back to the pre-curation label
            if any(useOrig)
                kfall = kv(useOrig);
                inOrig = kfall >= 1 & kfall <= numel(originalGroupList);
                if any(inOrig)
                    out = originalGroupList(kfall(inOrig));
                    % map back into lbls
                    validIdx = find(valid);
                    fallIdx  = validIdx(useOrig);
                    fallIdx  = fallIdx(inOrig);
                    lbls(fallIdx) = out;
                end
            end
        end
    end

    function r = acgRefractoryRatio(k)
        % Refractory ratio: count of ISI < 1.5 ms divided by count in
        % the 5-25 ms window. Smaller is better (cleaner refractory
        % hole). Returns 1 (worst) when there are too few spikes or
        % the baseline window is empty, so under-sampled units do not
        % spuriously look refractory.
        r = 1;
        if isnan(groupList(k)), return; end
        spk = sortedRes.spike_idx(sortedRes.unifiedLabels == groupList(k));
        if numel(spk) < 50, return; end
        d = diff(sort(double(spk))) / cfg.samplingFrequency * 1000;
        if isempty(d), return; end
        nearZero = sum(d < 1.5);
        baseline = sum(d >= 5 & d <= 25);
        if baseline < 5, return; end
        r = nearZero / baseline;
    end

    function ok = checkPCDistance(gA, gB, thr)
        % Project each group's mean waveform onto the shared channel's
        % PCA basis and compare PC1/PC2. Returns false only when both
        % groups have a valid projection AND their normalised distance
        % exceeds thr; empty PCA / out-of-footprint cases pass through.
        %
        % When the user-selected similarity metric is XCorr, mwB is
        % first shifted by the lag that maximises the cross-correlation
        % with mwA before projecting. This way the PC distance reflects
        % shape mismatch only, not a sub-sample timing offset that the
        % XCorr similarity score has already accounted for. KL / Bhatta
        % already operate on per-spike PCA scores, so no realignment
        % is applied for those metrics.
        ok = true;
        ch = channelList(gA);
        if ch < 1 || ch > numChannels, return; end
        if isempty(PCA{ch}) || ~isfield(PCA{ch},'coeff') || isempty(PCA{ch}.coeff)
            return;
        end
        mwA = localMeanWaveform(gA, ch);
        mwB = localMeanWaveform(gB, ch);
        if isempty(mwA) || isempty(mwB), return; end

        if strcmp(distanceEstType, 'XCorr')
            mwA_row = mwA(:)';
            mwB_row = mwB(:)';
            M2 = numel(mwA_row);
            if M2 == numel(mwB_row) && M2 > 4
                maxLag = max(1, round(M2/4));
                [~, bestLag] = max_half_corr(mwA_row, mwB_row, 1, M2, maxLag, 0);
                if isfinite(bestLag) && bestLag ~= 0
                    % bestLag is "B leads A by bestLag samples"; shift
                    % B by -bestLag to bring it into A's frame.
                    mwB = circshift(mwB_row, -bestLag);
                    mwB = mwB(:);
                end
            end
        end

        pcA = (mwA(:)' - PCA{ch}.mu) * PCA{ch}.coeff;
        pcB = (mwB(:)' - PCA{ch}.mu) * PCA{ch}.coeff;
        scale = max(PCA{ch}.max(:));
        if scale <= 0, scale = 1; end
        ok = norm(pcA - pcB) / scale <= thr;
    end

    function mw = localMeanWaveform(g, ch)
        % Mean waveform of group g on global channel ch (empty if ch is
        % outside g's local footprint stored in sampleWaveform).
        mw = [];
        if isnan(channelList(g)), return; end
        localIdx = ch - channelList(g) + numChannelPlot + 1;
        if localIdx < 1 || localIdx > 2*numChannelPlot + 1, return; end
        nLocal = size(sampleWaveform, 2);
        if localIdx > nLocal, return; end
        mw = squeeze(sampleWaveform(g, localIdx, :));
    end

    function ok = checkMergedISI(gA, gB, isiAbs, isiBudget)
        % Decide whether merging gA into gB (or vice versa) is safe based
        % on ISI-violation behaviour. Three independent gates:
        %
        %   1) Hard cap on the merged train: isiM <= isiAbs.
        %   2) Per-parent cap: each parent's own ISI must already be
        %      reasonable (<= isiAbs * isiBudget). This catches the
        %      "small but contaminated unit dissolves into a large clean
        %      one" failure mode -- the merged ISI looks fine because
        %      the large clean parent dominates the spike count, so the
        %      merge appears valid even though the small parent is
        %      clearly bad.
        %   3) Size-weighted budget: the merged ISI must not exceed
        %      isiBudget * weighted-average parent ISI, where the
        %      weights are the parent spike counts. The merged ISI of
        %      two independent point processes is approximately this
        %      weighted average, so the budget is naturally tighter
        %      when one parent is much smaller. The budget is also
        %      tightened toward 1.0 (no margin) as the size disparity
        %      grows.
        labA = groupList(gA); labB = groupList(gB);
        spkA = sortedRes.spike_idx(sortedRes.unifiedLabels == labA);
        spkB = sortedRes.spike_idx(sortedRes.unifiedLabels == labB);
        NA = numel(spkA); NB = numel(spkB);
        spkM = sort([spkA(:); spkB(:)]);
        if numel(spkM) < 2
            ok = true;
            return;
        end
        [~,~,isiM] = getISIViolations(spkM, cfg.samplingFrequency, thresholdISI);
        isiA = preprocessed.isiViolation(gA);
        isiB = preprocessed.isiViolation(gB);

        % Per-parent cap.
        parentCap = isiAbs * max(isiBudget, 1);

        % Spike-count-weighted average of parent ISIs.
        denom        = max(NA + NB, 1);
        weightedIsi  = (NA*isiA + NB*isiB) / denom;

        % Each parent's ISI must sit close to the size-weighted
        % average (within 0.5 %). This blocks merges where the two
        % parents have wildly different ISI profiles -- a strong sign
        % that they aren't drawn from the same source even when the
        % combined train looks fine.
        weightedDevCap = 0.5;
        devA = abs(isiA - weightedIsi);
        devB = abs(isiB - weightedIsi);

        % Size-disparity factor in [0, 1]: 1 when NA == NB, ~0 when
        % one is wildly larger. Used to shrink the effective budget
        % toward 1.0 when sizes are very imbalanced.
        sizeFactor   = 2 * min(NA, NB) / denom;
        effBudget    = 1 + (isiBudget - 1) * sizeFactor;

        budgetCap    = max(effBudget * weightedIsi, 1e-6);

        ok = (isiM   <= isiAbs)         && ...
             (isiA   <= parentCap)       && ...
             (isiB   <= parentCap)       && ...
             (devA   <= weightedDevCap)  && ...
             (devB   <= weightedDevCap)  && ...
             (isiM   <= budgetCap);
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

        % Pick one representative group index per unique surviving
        % label, skipping NaN entries (units marked as removed). The
        % older code indexed channelList / sampleWaveform / ... by
        % the LABEL VALUE itself, which (a) only worked when labels
        % happened to equal indices and (b) crashes on NaN.
        [uniqLabels_curated, repIdx] = unique(groupList);
        keepRow                      = ~isnan(uniqLabels_curated);
        uniqLabels_curated           = uniqLabels_curated(keepRow);
        repIdx                       = repIdx(keepRow);

        curatedSamples.unifiedLabels = uniqLabels_curated;
        curatedSamples.channelNum    = channelList(repIdx);
        curatedSamples.waveform      = sampleWaveform(repIdx,:,:);
        curatedSamples.unitIsolation = unitIsolation(repIdx);
        curatedSamples.validInterval = stable_length(repIdx,:);
        curatedSamples.spikePolarity = mainPolarity(repIdx,:);
        if exist('unitNotes','var')
            % Per-unit free-text notes (added in this version of the
            % GUI). Saved alongside the curated payload so they
            % survive across sessions.
            try
                curatedSamples.notes = unitNotes(repIdx);
            catch
            end
        end
        curatedSamples.original.unifiedLabels = originalUnifiedLabels;
        curatedSamples.original.channelNum    = originalChannelNum;
        curatedSamples.original.waveform      = originalSampleWaveform;
        curatedSamples.original.unitIsolation = unitIsolation;

        % --- Full session state (resume on next launch) -----------
        % Everything the GUI needs to stand back up exactly where
        % the user left off: per-unit arrays (length = numGroups,
        % including dropped/merged units in their original slots)
        % and the post-curation per-spike arrays (with -1 markers
        % for cleanup-stripped spikes intact). On the next launch,
        % the constructor checks for curated_sample.mat in this
        % folder and applies session.* on top of the freshly-loaded
        % raw sort data.
        curatedSamples.session.groupList     = groupList;
        curatedSamples.session.channelList   = channelList;
        curatedSamples.session.unitIsolation = unitIsolation;
        curatedSamples.session.stable_length = stable_length;
        curatedSamples.session.mainPolarity  = mainPolarity;
        curatedSamples.session.sidePolarity  = sidePolarity;
        curatedSamples.session.mergedFlag    = mergedFlag;
        curatedSamples.session.unitNotes     = unitNotes;
        curatedSamples.session.detectblity   = detectblity;
        curatedSamples.session.numGroups     = numGroups;
        curatedSamples.session.spikeLabels   = sortedRes.unifiedLabels;
        curatedSamples.session.spikeIdx      = sortedRes.spike_idx;
        curatedSamples.session.spikeChannels = sortedRes.channelNum;

        filename = fullfile(outputFolder,'curated_sample.mat');
        save(filename,'curatedSamples', '-v7.3');

        % Also write a plain-text CSV summarising every surviving unit
        % so the metrics can be opened directly in Excel / pandas /
        % awk without having to load the .mat. Columns follow the
        % same per-unit order as curated_sample.mat.
        try
            csvFile = fullfile(outputFolder, 'curated_metrics.csv');
            nKept   = numel(repIdx);
            ratesC  = preprocessed.firingRate(repIdx);
            isiC    = preprocessed.isiViolation(repIdx);
            snrC    = 1 + detectblity(repIdx);
            polC    = mainPolarity(repIdx);
            stableA = stable_length(repIdx,1);
            stableB = stable_length(repIdx,2);
            isoStr  = unitIsolation(repIdx);
            if exist('unitNotes','var')
                notesC = unitNotes(repIdx);
            else
                notesC = repmat({''}, nKept, 1);
            end
            % Build the table column by column to avoid issues if any
            % of the inputs has unexpected shape.
            T = table( ...
                uniqLabels_curated(:),     ...
                channelList(repIdx),       ...
                ratesC(:),                 ...
                isiC(:),                   ...
                snrC(:),                   ...
                isoStr(:),                 ...
                polC(:),                   ...
                stableA(:),                ...
                stableB(:),                ...
                notesC(:),                 ...
                'VariableNames', { ...
                    'Label','Channel','FiringRate_Hz', ...
                    'ISI_violation_pct','SNR','Isolation', ...
                    'Polarity','StableStart','StableEnd','Notes'});
            writetable(T, csvFile);
        catch ME
            warning('CSV export failed: %s', ME.message);
        end

        % Checkpoint user prefs alongside the curated outputs so the
        % next launch on this dataset remembers the panel tuning.
        % saveKiaPrefs() writes to RES_Sorted/kiaSort_curate_prefs.mat
        % and is best-effort (never throws), so a pref-write failure
        % doesn't mask the data-save success message above.
        saveKiaPrefs();

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

    function tf = plotUsesCosmetic()
        % Plot types whose appearance depends on ampScale / lineWidth /
        % alphaLevel. Cosmetic-only callbacks early-out when we are
        % showing a plot type that ignores them (CCG / ISI / Density),
        % saving a wasted full redraw.
        tf = any(strcmpi(plotType, {'Waveform','Trace','Amplitude', ...
                                    'PCs','Features','Multiple'}));
    end

    function tf = wfFastPathReady()
        % True iff plotWaveforms has populated wfPlotCache and every
        % cached line handle is still valid (not deleted by a later
        % redraw of a different plot type).
        tf = strcmpi(plotType,'Waveform') && ~isempty(wfPlotCache) ...
             && isfield(wfPlotCache,'lines') && ~isempty(wfPlotCache.lines);
        if tf
            for ii = 1:numel(wfPlotCache.lines)
                hs = wfPlotCache.lines{ii};
                if isempty(hs) || ~all(isvalid(hs))
                    tf = false;
                    return;
                end
            end
        end
    end

    function updateScale(sVal)
        ampScale = sVal;
        if wfFastPathReady()
            % Fast path: rescale the existing line YData in place.
            % set(handles, {'YData'}, cellOfRows) is a single MATLAB
            % call per channel-group, instead of one call per line.
            for ii = 1:numel(wfPlotCache.lines)
                hs = wfPlotCache.lines{ii};
                us = wfPlotCache.unscaledY{ii};
                yo = wfPlotCache.yOff(ii);
                yScaled = us * sVal + yo;
                nRows = size(yScaled, 1);
                if nRows == 1
                    set(hs(1), 'YData', yScaled);
                else
                    nUse = min(nRows, numel(hs));
                    if nUse > 0
                        ydCells = mat2cell(yScaled(1:nUse,:), ones(nUse,1), size(yScaled,2));
                        set(hs(1:nUse), {'YData'}, ydCells);
                    end
                end
            end
        elseif plotUsesCosmetic()
            plotOption();
        end
    end

    function updateLineWidth(sVal)
        if sVal == 0
            sVal = eps;
        end
        lineWidth = sVal;
        if wfFastPathReady()
            for ii = 1:numel(wfPlotCache.lines)
                hs   = wfPlotCache.lines{ii};
                lw   = sVal * (1 + double(wfPlotCache.isMain(ii)));   % main = 2×
                set(hs, 'LineWidth', lw);
            end
        elseif plotUsesCosmetic()
            plotOption();
        end
    end

    function updateAlphaLevel(sVal)
        if sVal == 0
            sVal = eps;
        end
        alphaLevel = sVal;
        if wfFastPathReady()
            for ii = 1:numel(wfPlotCache.lines)
                hs  = wfPlotCache.lines{ii};
                rgb = wfPlotCache.rgb(ii,:);
                set(hs, 'Color', [rgb, sVal]);
            end
        elseif plotUsesCosmetic()
            plotOption();
        end
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
        chanFilterActive = chkInclusion.Value || chkExclusion.Value;
        if chkInclusion.Value
            lastUnit = 0;
            incChan = isiInc & rateInc & snrInc & polarityInc;
            selectUnits = find(incChan);

            visibleGroups = currentPageVisibleGroups();
            chkBoxSelect(:) = 0;
            chkBoxSelect(selectUnits) = 1;
            if chkTogVis.Value
                chkBoxVis(:) = 0;
                chkBoxVis(selectUnits) = 1;
                selected_order = selectUnits';
            else
                selected_order = [];
            end
            
            for i = visibleGroups(:)'
                if isgraphics(groupSelectCheckboxes(i))
                    groupSelectCheckboxes(i).Value = incChan(i);
                end
                if chkTogVis.Value && isgraphics(groupVisCheckboxes(i))
                    groupVisCheckboxes(i).Value = incChan(i);
                end
            end
            unitNum = sprintf('/%d',length(selectUnits));
            changeuifUnitAll(unitNum);
            changeuifUnit(num2str(lastUnit),1);
            
        elseif ~chkExclusion.Value
            selected_order = [];
            onAllSelectChecked(false);
            onAllVisChecked(false);
            selectUnits = displayedGroups;
            lastUnit = 0;
            unitNum = sprintf('/%d',length(selectUnits));
            changeuifUnitAll(unitNum);
            changeuifUnit(num2str(lastUnit),1);

        end

        if chkExclusion.Value
            lastUnit = 0;
            exChan = ~isiInc | ~rateInc | ~snrInc & polarityInc;
            selectUnits = find(exChan);

            visibleGroups = currentPageVisibleGroups();
            chkBoxSelect(:) = 0;
            chkBoxSelect(selectUnits) = 1;
            if chkTogVis.Value
                chkBoxVis(:) = 0;
                chkBoxVis(selectUnits) = 1;
            end
            for i = visibleGroups(:)'
                if isgraphics(groupSelectCheckboxes(i))
                    groupSelectCheckboxes(i).Value = exChan(i);
                end
                if chkTogVis.Value && isgraphics(groupVisCheckboxes(i))
                    groupVisCheckboxes(i).Value = exChan(i);
                    if groupVisCheckboxes(i).Value
                        selected_order = [selected_order, i];
                    else
                        selected_order(selected_order==i) = [];
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
            selectUnits = displayedGroups;
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

        % Skip past hidden (zero-rate / dropped) units when navigating.
        step = sign(val);
        if step ~= 0
            while thisUnit >= 1 && thisUnit <= length(selectUnits) ...
                    && ~ismember(selectUnits(thisUnit), displayedGroups)
                thisUnit = thisUnit + step;
            end
        end

        if thisUnit >= 1 && thisUnit <= length(selectUnits)
            thisCol = selectUnits(thisUnit);

            % Only change page if the target unit isn't already on
            % the current page. Page index uses position inside
            % displayedGroups so it matches the filtered render.
            visNow = currentPageVisibleGroups();
            if ~ismember(thisCol, visNow)
                posInDisplay = find(displayedGroups == thisCol, 1);
                if isempty(posInDisplay)
                    lastUnit = 0;
                    if val ~= 0
                        changeuifUnit(num2str(lastUnit), 1);
                        changeuifUnitAll(sprintf('/%d', length(selectUnits)));
                    end
                    plotOption();
                    return;
                end
                targetPage = ceil(posInDisplay / PAGE_SIZE);
                if targetPage ~= currentPage
                    currentPage = targetPage;
                    onChangePage(0);
                end
            end

            % Make sure the unit is now visible
            if isgraphics(groupVisCheckboxes(thisCol))
                groupVisCheckboxes(thisCol).Value = true;
                if isgraphics(groupSelectCheckboxes(thisCol))
                    groupSelectCheckboxes(thisCol).Value = true;
                end
                chkBoxSelect(thisCol) = 1;
                chkBoxVis(thisCol) = 1;
                selected_order = [thisCol];
                if isgraphics(groupRadioButtons(thisCol))
                    groupRadioButtons(thisCol).Value = true;
                end
                c = getGroupColor(groupList(thisCol));
                groupVisCheckboxes(thisCol).FontColor = c;
                groupVisCheckboxes(thisCol).FontWeight ='Bold';
                if isgraphics(groupSelectCheckboxes(thisCol))
                    groupSelectCheckboxes(thisCol).FontColor = c;
                    groupSelectCheckboxes(thisCol).FontWeight ='Bold';
                end
                if isgraphics(lblChannelHandles(thisCol))
                    lblChannelHandles(thisCol).BackgroundColor = c;
                    lblChannelHandles(thisCol).FontColor = bestTextColorFor(c);
                end
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

            % Land on the page of the HIGHER-firing-rate member, but
            % fall back to whichever side is actually in displayedGroups
            % (the other one may have been hidden by the zero-rate
            % filter).
            if preprocessed.firingRate(thisCol) >= preprocessed.firingRate(thisRow)
                anchor = thisCol;  alt = thisRow;
            else
                anchor = thisRow;  alt = thisCol;
            end
            if ~ismember(anchor, displayedGroups) && ismember(alt, displayedGroups)
                anchor = alt;
            end
            visNow = currentPageVisibleGroups();
            if ~ismember(anchor, visNow)
                posAnchor = find(displayedGroups == anchor, 1);
                if ~isempty(posAnchor)
                    targetPage = ceil(posAnchor / PAGE_SIZE);
                    if targetPage ~= currentPage
                        currentPage = targetPage;
                        onChangePage(0);
                    end
                end
            end

            chkBoxVis(thisCol) = 1;
            chkBoxVis(thisRow) = 1;
            chkBoxSelect(thisCol) = 1;
            chkBoxSelect(thisRow) = 1;
            selected_order = [thisCol,thisRow];

            cColRow = getGroupColor(groupList(thisCol));
            if isgraphics(groupVisCheckboxes(thisCol))
                groupVisCheckboxes(thisCol).Value     = true;
                groupVisCheckboxes(thisCol).FontColor = cColRow;
                groupVisCheckboxes(thisCol).FontWeight = 'Bold';
            end
            if isgraphics(groupSelectCheckboxes(thisCol))
                groupSelectCheckboxes(thisCol).Value     = true;
                groupSelectCheckboxes(thisCol).FontColor = cColRow;
                groupSelectCheckboxes(thisCol).FontWeight = 'Bold';
            end
            if isgraphics(lblChannelHandles(thisCol))
                lblChannelHandles(thisCol).BackgroundColor = cColRow;
                lblChannelHandles(thisCol).FontColor       = bestTextColorFor(cColRow);
            end

            cRow = getGroupColor(groupList(thisRow));
            if isgraphics(groupVisCheckboxes(thisRow))
                groupVisCheckboxes(thisRow).Value     = true;
                groupVisCheckboxes(thisRow).FontColor = cRow;
                groupVisCheckboxes(thisRow).FontWeight = 'Bold';
            end
            if isgraphics(groupSelectCheckboxes(thisRow))
                groupSelectCheckboxes(thisRow).Value     = true;
                groupSelectCheckboxes(thisRow).FontColor = cRow;
                groupSelectCheckboxes(thisRow).FontWeight = 'Bold';
            end
            if isgraphics(lblChannelHandles(thisRow))
                lblChannelHandles(thisRow).BackgroundColor = cRow;
                lblChannelHandles(thisRow).FontColor       = bestTextColorFor(cRow);
            end

            selected_order = [thisCol,thisRow];
            if preprocessed.firingRate(thisCol) > preprocessed.firingRate(thisRow)
                if isgraphics(groupRadioButtons(thisCol))
                    groupRadioButtons(thisCol).Value = true;
                end
            else
                if isgraphics(groupRadioButtons(thisRow))
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

    function [aArr, bArr] = currentAutoPairList()
        % Resolve the audit arrays for whichever target the user
        % picked from the Auto-panel dropdown. Two return shapes:
        %   PAIR-mode: aArr and bArr are equal-length index vectors
        %              and the walker highlights both per click.
        %   UNIT-mode: bArr is empty; aArr is a list of unit
        %              indices and the walker highlights one per
        %              click.
        switch autoPairType
            case 'Overlap removed'
                aArr = autoOverlapDroppedK;
                bArr = autoOverlapTriggerJ;
            case 'CCG cleaned'
                % Legacy walker target; the dropdown no longer
                % exposes this, but a hot-reload could still set
                % autoPairType to it -- handle gracefully.
                aArr = autoCCGStripCont;
                bArr = autoCCGStripClean;
            case {'SUA+','SUA','MUA+','MUA'}
                % Unit-class walker. Iterate every surviving unit
                % whose unitIsolation matches the dropdown value.
                aArr = unitsByIsolation(autoPairType);
                bArr = [];
            case 'NaN'
                % Dropped / unclassified units: anything with NaN
                % groupList OR an explicit 'NA' isolation flag. We
                % unique() because zero-rate dropped units may carry
                % both markers and we don't want duplicate stops.
                isoNA = strcmp(unitIsolation, 'NA');
                if numel(isoNA) ~= numel(groupList)
                    isoNA = false(size(groupList));
                end
                aArr = unique(find(isnan(groupList(:)) | isoNA(:)));
                bArr = [];
            otherwise   % 'Merged'
                aArr = autoMergedPrimary;
                bArr = autoMergedAbsorbed;
        end
    end

    function idx = unitsByIsolation(isoStr)
        % Return the indices of every surviving unit whose
        % unitIsolation cell exactly matches isoStr. Skips NaN
        % groupList slots so the unit-class walker only stops at
        % units that still exist on the page.
        idx = [];
        if isempty(unitIsolation), return; end
        n = numel(unitIsolation);
        for ii = 1:n
            if ii > numel(groupList) || isnan(groupList(ii)), continue; end
            if strcmp(unitIsolation{ii}, isoStr)
                idx(end+1,1) = ii; %#ok<AGROW>
            end
        end
    end

    function tf = isUnitModePairType(typeStr)
        % True for dropdown items that walk a flat list of units
        % rather than (left,right) pairs.
        tf = any(strcmp(typeStr, {'SUA+','SUA','MUA+','MUA','NaN'}));
    end

    function onSelectMergedPairs(direction)
        % Walker over the targets the user picked from the dropdown.
        % Pair-mode targets (Merged / Overlap removed) keep their
        % (primary, absorbed) shape and the walker highlights BOTH
        % units. Unit-mode targets (SUA+ / SUA / MUA+ / MUA / NaN)
        % carry a flat list of unit indices and the walker
        % highlights ONE unit per click.
        [aArr, bArr] = currentAutoPairList();
        nPairs = numel(aArr);
        if nPairs == 0
            % Nothing to walk; just reset the counter UI and bail out.
            if exist('uifAutoPair','var') && isvalid(uifAutoPair)
                uifAutoPair.Value = '0';
            end
            if exist('uifAutoPair2','var') && isvalid(uifAutoPair2)
                uifAutoPair2.Value = '/0';
            end
            return;
        end

        thisIdx = lastAutoPair + direction;
        if thisIdx < 1 || thisIdx > nPairs
            % Wrap-style behaviour matches the Similarity walker: if we
            % run off either end, snap back to a valid range and stop
            % stepping further.
            thisIdx = max(1, min(nPairs, thisIdx));
        end

        unitMode  = isempty(bArr) || numel(bArr) ~= numel(aArr);
        primaryIdx = aArr(thisIdx);
        if unitMode
            absorbedIdx = [];
        else
            absorbedIdx = bArr(thisIdx);
        end

        % Anchor the page to whichever side of the pair is currently
        % displayable. If primary is dropped (NaN groupList -> hidden by
        % the zero-rate filter), use absorbed (the trigger neighbour);
        % otherwise use primary.
        anchorIdx = primaryIdx;
        if ~ismember(primaryIdx, displayedGroups)
            if ~unitMode && ismember(absorbedIdx, displayedGroups)
                anchorIdx = absorbedIdx;
            end
        end
        visNow = currentPageVisibleGroups();
        if ~ismember(anchorIdx, visNow)
            posAnchor = find(displayedGroups == anchorIdx, 1);
            if ~isempty(posAnchor)
                targetPage = ceil(posAnchor / PAGE_SIZE);
                if targetPage ~= currentPage
                    currentPage = targetPage;
                    onChangePage(0);
                end
            end
        end

        selected_order        = [];
        onAllSelectChecked(false);
        onAllVisChecked(false);

        chkBoxVis(primaryIdx)     = 1;
        chkBoxSelect(primaryIdx)  = 1;
        if ~unitMode
            chkBoxVis(absorbedIdx)    = 1;
            chkBoxSelect(absorbedIdx) = 1;
            selected_order            = [primaryIdx, absorbedIdx];
        else
            selected_order            = primaryIdx;
        end

        cP = getGroupColor(groupList(primaryIdx));
        if isgraphics(groupVisCheckboxes(primaryIdx))
            groupVisCheckboxes(primaryIdx).Value     = true;
            groupVisCheckboxes(primaryIdx).FontColor = cP;
            groupVisCheckboxes(primaryIdx).FontWeight = 'Bold';
        end
        if isgraphics(groupSelectCheckboxes(primaryIdx))
            groupSelectCheckboxes(primaryIdx).Value     = true;
            groupSelectCheckboxes(primaryIdx).FontColor = cP;
            groupSelectCheckboxes(primaryIdx).FontWeight = 'Bold';
        end
        if isgraphics(lblChannelHandles(primaryIdx))
            lblChannelHandles(primaryIdx).BackgroundColor = cP;
            lblChannelHandles(primaryIdx).FontColor       = bestTextColorFor(cP);
        end
        if ~unitMode
            cA = getGroupColor(groupList(absorbedIdx));
            if isgraphics(groupVisCheckboxes(absorbedIdx))
                groupVisCheckboxes(absorbedIdx).Value     = true;
                groupVisCheckboxes(absorbedIdx).FontColor = cA;
                groupVisCheckboxes(absorbedIdx).FontWeight = 'Bold';
            end
            if isgraphics(groupSelectCheckboxes(absorbedIdx))
                groupSelectCheckboxes(absorbedIdx).Value     = true;
                groupSelectCheckboxes(absorbedIdx).FontColor = cA;
                groupSelectCheckboxes(absorbedIdx).FontWeight = 'Bold';
            end
            if isgraphics(lblChannelHandles(absorbedIdx))
                lblChannelHandles(absorbedIdx).BackgroundColor = cA;
                lblChannelHandles(absorbedIdx).FontColor       = bestTextColorFor(cA);
            end
        end

        if isgraphics(groupRadioButtons(primaryIdx))
            groupRadioButtons(primaryIdx).Value = true;
        end

        lastAutoPair = thisIdx;
        if exist('uifAutoPair','var') && isvalid(uifAutoPair)
            uifAutoPair.Value  = num2str(thisIdx);
        end
        if exist('uifAutoPair2','var') && isvalid(uifAutoPair2)
            uifAutoPair2.Value = sprintf('/%d', nPairs);
        end
        plotOption();
    end

    function changeAutoPair(val)
        % Editable counter on the Auto panel: jump to a specific pair
        % by typing its index. Mirrors the parsing from changeuifPar.
        % Honors the dropdown's current Pair-type selection.
        [aArr, ~] = currentAutoPairList();
        nPairs = numel(aArr);
        if nPairs == 0, return; end
        if isnumeric(val)
            target = val;
        else
            target = str2double(val);
        end
        if isnan(target), return; end
        target = round(target);
        target = max(1, min(nPairs, target));
        lastAutoPair = target - 1;        % onSelectMergedPairs(+1) lands on `target`
        onSelectMergedPairs(1);
    end

    function setAutoPairType(val)
        % Switch which list the Auto Next/Prev walker traverses.
        % Resets the walker state, refreshes the /total counter, and
        % retitles the counter label so the user can see at a glance
        % whether the walker is in pair-mode (#Pair) or unit-mode
        % (#Unit).
        autoPairType = val;
        lastAutoPair = 0;
        if exist('uifAutoPair','var') && isvalid(uifAutoPair)
            uifAutoPair.Value = '0';
        end
        if exist('uifAutoPair2','var') && isvalid(uifAutoPair2)
            [aArr, ~] = currentAutoPairList();
            uifAutoPair2.Value = sprintf('/%d', numel(aArr));
        end
        if exist('lblAutoPair','var') && isvalid(lblAutoPair)
            if isUnitModePairType(val)
                lblAutoPair.Text = 'Unit #:';
            else
                lblAutoPair.Text = 'Pair #:';
            end
        end
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
        visPageGrp = currentPageVisibleGroups();
        for k = visPageGrp(:)'
            % isgraphics catches both placeholders AND deleted handles
            % (isvalid returns true for placeholders, so it's not enough
            % here -- Phase 4 of Auto Curation can change displayedGroups
            % without re-rendering the page).
            if ~isgraphics(lblGroupHandles(k))
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
            if isgraphics(lblChannelHandles(k))
                lblChannelHandles(k).Text = txt;
            end

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
            if isgraphics(lblISIViolation(k))
                lblISIViolation(k).Text = txt;
            end

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
            case 'Amp Dist'
                multipleCheckbox('off');
                plotAmpDistribution();
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

    function plotTable(panelOption)
        % Sortable per-unit summary table. Renders into the parent
        % panel passed as panelOption -- in normal use this is the
        % standalone-table window's content panel (openTableWindow).
        % The empty-arg path lower down is legacy from when Table
        % used to be a plot type and rendered into axChannels; it's
        % still here as a safe fallback if anyone calls plotTable
        % without an argument. Click any column header to sort
        % (built-in uitable behaviour via ColumnSortable); click
        % the same header again to flip the direction. Notes
        % column is editable in place.
        if nargin < 1
            panelOption = [];
        end
        if isempty(panelOption)
            if strcmp(lastPlotted,'Table')
                delete(allchild(axChannels));
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uipanel(bottomRightGrid,'BorderType','none',...
                    'BackgroundColor',figColor);
                axChannels.Layout.Column = [2 3];
                axChannels.Layout.Row    = [1 3];
                applyColorScheme(axChannels, figColor);
            end
            plotAxis = axChannels;
            lastPlotted = plotType;
        else
            plotAxis = panelOption;
            try, delete(allchild(plotAxis)); catch, end
        end

        % Build per-row data. Every unit appears -- including ones
        % the user has removed (groupList(i) == NaN). Removed rows
        % get ISO = "REM" so they stay searchable / sortable,
        % and paintSelectedTableRows draws them in italic grey so
        % they read as visually inactive without leaving the
        % table. Restoring via Undo / Reset re-fills their cells
        % on the next refresh.
        snrAllT = 1 + detectblity;
        nLive   = numGroups;
        removedMask = isnan(groupList(:));
        % Contiguous Group ID = rank among the shown (nonzero) groups,
        % matching the Groups-panel ID column. NaN for hidden/removed.
        gidAll = nan(numGroups, 1);
        if ~isempty(displayedGroups)
            gidAll(displayedGroups) = (1:numel(displayedGroups)).';
        end

        labels   = double(groupList(:));
        channels = double(channelList(:));
        rates    = double(preprocessed.firingRate(:));
        isis     = double(preprocessed.isiViolation(:));
        snrs     = double(snrAllT(:));
        stables  = (double(stable_length(:,2)) - double(stable_length(:,1))) ...
                   / max(trialLength,eps) * 100;
        isolats  = strings(nLive, 1);
        notesC   = strings(nLive, 1);
        for ii = 1:nLive
            if removedMask(ii)
                isolats(ii) = "REM";
            elseif ii <= numel(unitIsolation) && (ischar(unitIsolation{ii}) || isstring(unitIsolation{ii}))
                isolats(ii) = string(unitIsolation{ii});
            else
                isolats(ii) = "NA";
            end
            if ii <= numel(unitNotes) && (ischar(unitNotes{ii}) || isstring(unitNotes{ii}))
                notesC(ii) = string(unitNotes{ii});
            end
        end

        % Single-line headers so the unit label sits inline with
        % the metric name (compact, no "(Hz)" floating below).
        % Target column is absent: clicking a row sets the merge
        % target radio in the main GUI; MATLAB's own focused-cell
        % highlight shows where the cursor is.
        colNames = { 'ID', 'G', 'Ch', 'Rate (Hz)', 'ISI (%)', ...
                     'SNR', 'ISO', 'Stable', 'Notes' };

        % Notes (col 9) is the only editable column.
        nCols        = 9;
        editableMask = false(1, nCols);
        editableMask(9) = true;

        % We use a CELL-ARRAY Data so ColumnFormat = 'bank' is
        % actually honoured. uitable silently ignores ColumnFormat
        % when Data is a MATLAB `table`, which is why earlier
        % versions of this code showed full double precision in
        % Rate / ISI / SNR / Stable. Cells are still typed (numeric
        % vs char) so each column header still sorts by its native
        % type -- the 'bank' format only changes display, not sort.
        rows = cell(nLive, nCols);
        for ii = 1:nLive
            rows{ii, 1} = gidAll(ii);
            rows{ii, 2} = labels(ii);
            rows{ii, 3} = channels(ii);
            rows{ii, 4} = rates(ii);
            rows{ii, 5} = isis(ii);
            rows{ii, 6} = snrs(ii);
            rows{ii, 7} = char(isolats(ii));
            rows{ii, 8} = stables(ii);
            rows{ii, 9} = char(notesC(ii));
        end

        % 'fit' shrinks each column to its widest value so the
        % table reads compact. Notes is given an explicit cap so
        % a single long note can't blow the whole table wider
        % than the figure. 'bank' = fixed-point with two decimals
        % (e.g. 12.50, 0.07).
        colWidths  = repmat({'fit'}, 1, nCols);
        colWidths{9} = 180;     % Notes capped to keep the table compact
        colFormats = {'numeric','numeric','numeric', ...
                      'bank','bank','bank', ...
                      'char','bank','char'};

        try
            tbl = uitable(plotAxis, ...
                'Units','normalized','Position',[0 0 1 1], ...
                'Data', rows, ...
                'ColumnName', colNames, ...
                'ColumnWidth', colWidths, ...
                'ColumnFormat', colFormats, ...
                'ColumnEditable', editableMask, ...
                'ColumnSortable', true(1,nCols), ...
                'RowName', [], ...
                'CellSelectionCallback', @onTableCellSelect, ...
                'CellEditCallback',      @onTableCellEdit);
            applyColorScheme(tbl, figColor);
            tableHandle = tbl;
            paintSelectedTableRows();
        catch ME
            warning('Table render failed: %s', ME.message);
        end
    end

    function paintSelectedTableRows()
        % Repaint per-row background colours on the most recent
        % table view so selected units stand out using their group
        % colour (with auto-chosen text colour for contrast).
        % Called both at table render and on every selection toggle.
        if isempty(tableHandle) || ~isgraphics(tableHandle), return; end
        tbl = tableHandle;
        try
            % Wipe any previous styling so deselected rows revert
            % to the default theme.
            removeStyle(tbl);
        catch
            % Older MATLAB versions may not support removeStyle on
            % a freshly-created table -- ignore.
        end

        % Map each displayed row back to its unit array index (col 1 is
        % the contiguous Group ID, col 2 the label). Survives a column-
        % header sort, which reorders the displayed rows.
        nRowsT = size(tbl.Data, 1);
        idxCol = nan(nRowsT, 1);
        for r = 1:nRowsT
            idxCol(r) = tableUnitIdx(tbl.Data{r, 1}, tbl.Data{r, 2});
        end

        % Removed units carry ISO == "REM" (col 7); paint those rows
        % italic-grey. Selection styling below can override that on any
        % row that's both removed and still chkBoxSelect-ticked.
        try
            removedStyle = uistyle('FontAngle','italic', ...
                'FontColor',[0.55 0.55 0.55]);
            for r = 1:nRowsT
                isoV = tbl.Data{r, 7};
                if (ischar(isoV) || isstring(isoV)) && strcmp(char(isoV),'REM')
                    addStyle(tbl, removedStyle, 'row', r);
                end
            end
        catch
            % uistyle/addStyle not available -> just skip the
            % grey-out; the ISO="REM" tag still flags the row.
        end

        selUnits = find(chkBoxSelect(:));
        for ii = 1:numel(selUnits)
            k = selUnits(ii);
            if k > numel(groupList) || isnan(groupList(k)), continue; end
            r = find(idxCol == k, 1);
            if isempty(r), continue; end
            try
                c   = getGroupColor(groupList(k));
                tc  = bestTextColorFor(c);
                stl = uistyle('BackgroundColor', c, 'FontColor', tc);
                addStyle(tbl, stl, 'row', r);
            catch
                % uistyle / addStyle missing -> fall through, no
                % colour highlight on this MATLAB release.
            end
        end
    end

    function openTableWindow()
        % Refresh the main page first so the table opens on the latest
        % curation state (Group IDs, metrics, removed units).
        try, recomputeDisplayedGroups(); catch, end
        try, onChangePage(0);            catch, end
        % Pop the curation table out into its own uifigure. Because
        % it lives in a separate window, the per-unit selection
        % handlers and the multi-plot redraw never touch it -- so
        % the user's column sort, scroll position and cell focus
        % all survive across selection changes. Singleton: a
        % second click brings the existing window to the front
        % rather than spawning a duplicate.
        %
        % The window also carries:
        %   * A Spike Group Controller toolbar (Remove / Merge /
        %     Split / Limit / Isolation / Undo / Reset / Notes)
        %     so the user can curate from inside the table window
        %     without switching back to the main GUI.
        %   * The same Ctrl+letter keyboard shortcuts as the main
        %     figure (R/M/T/K/U/D/L/E/H/N/P/S/[ ]/Up/Down/...).
        %   * Plain Space toggles the focused unit's Vis/Select.
        %
        % WindowStyle is set to 'alwaysontop' so the table stays
        % above the main GUI until the user closes it. MATLAB
        % releases that don't accept that value just fall through
        % to a normal window via the try/catch.
        if ~isempty(tableFigure) && isgraphics(tableFigure)
            try, figure(tableFigure); catch, end
            return;
        end
        % Initial width is sized to the table's content. The eight
        % numeric/short columns auto-fit and the Notes column is
        % capped at 180 px (see plotTable), so the table itself
        % needs roughly 600 px. The figure remains user-resizable
        % if longer notes warrant a wider window.
        tableFigure = uifigure( ...
            'Name','Curation Table', ...
            'Position',[200 200 620 660], ...
            'Color',figColor);
        applyColorScheme(tableFigure, figColor);
        tableFigure.WindowKeyPressFcn = @keyPressHandler;
        try
            tableFigure.WindowStyle = 'alwaysontop';
        catch
            % MATLAB version doesn't support 'alwaysontop' for
            % uifigure (added in R2023b). Fall back to a normal
            % window -- the user can still bring it forward by
            % clicking it.
        end

        winGrid = uigridlayout(tableFigure,[3 1], ...
            'RowHeight',{'fit','fit','1x'}, ...
            'ColumnWidth',{'1x'}, ...
            'Padding',[6 6 6 6], ...
            'RowSpacing',4);
        applyColorScheme(winGrid, figColor);

        % --- Row 1: window controls (Refresh / Help / Close) -----
        toolRow = uigridlayout(winGrid,[1 4], ...
            'ColumnWidth',{'fit','fit','fit','1x'}, ...
            'RowHeight',{'fit'}, ...
            'Padding',[0 0 0 0], ...
            'ColumnSpacing',6);
        toolRow.Layout.Row = 1;
        applyColorScheme(toolRow, figColor);

        btnRefresh = uibutton(toolRow,'Text','Refresh', ...
            'Tooltip','Reload the table from the latest curation state.', ...
            'ButtonPushedFcn',@(btn,ev) refreshTableWindow());
        btnRefresh.Layout.Column = 1;
        applyColorScheme(btnRefresh, figColor);

        btnHelp = uibutton(toolRow,'Text','Help', ...
            'Tooltip','How to select units, set the target, sort, and read columns.', ...
            'ButtonPushedFcn',@(btn,ev) showTableHelp());
        btnHelp.Layout.Column = 2;
        applyColorScheme(btnHelp, figColor);

        btnClose = uibutton(toolRow,'Text','Close', ...
            'Tooltip','Close the table window.', ...
            'ButtonPushedFcn',@(btn,ev) delete(tableFigure));
        btnClose.Layout.Column = 3;
        applyColorScheme(btnClose, figColor);

        % --- Row 2: Spike Group Controller buttons ---------------
        % Same actions as the main GUI's top button row, wired to
        % the same handlers so toggles / undos / resets all flow
        % through the existing logic. After each one we refresh
        % the table data (data may have changed) and repaint
        % selection rows.
        ctrlRow = uigridlayout(winGrid,[1 9], ...
            'ColumnWidth',repmat({'1x'}, 1, 9), ...
            'RowHeight',{'fit'}, ...
            'Padding',[0 0 0 0], ...
            'ColumnSpacing',4);
        ctrlRow.Layout.Row = 2;
        applyColorScheme(ctrlRow, figColor);

        btnTblRemove = uibutton(ctrlRow,'Text','Remove', ...
            'Tooltip','Remove the selected units. Ctrl+D', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onRemove));
        btnTblRemove.Layout.Column = 1;
        applyColorScheme(btnTblRemove, figColor);

        btnTblMerge = uibutton(ctrlRow,'Text','Merge', ...
            'Tooltip','Merge selected units into the radio-button choice. Ctrl+M', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onMerge));
        btnTblMerge.Layout.Column = 2;
        applyColorScheme(btnTblMerge, figColor);

        btnTblSplit = uibutton(ctrlRow,'Text','Split', ...
            'Tooltip','Split the selected unit(s). Ctrl+T', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onSplit));
        btnTblSplit.Layout.Column = 3;
        applyColorScheme(btnTblSplit, figColor);

        btnTblLimit = uibutton(ctrlRow,'Text','Limit', ...
            'Tooltip','Trim selected units to their stable interval. Ctrl+L', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onAutoCut));
        btnTblLimit.Layout.Column = 4;
        applyColorScheme(btnTblLimit, figColor);

        btnTblIso = uibutton(ctrlRow,'Text','Isolation', ...
            'Tooltip','Cycle isolation label (SUA+/SUA/MUA+/MUA). Ctrl+K', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onMUA));
        btnTblIso.Layout.Column = 5;
        applyColorScheme(btnTblIso, figColor);

        btnTblUndo = uibutton(ctrlRow,'Text','Undo', ...
            'Tooltip','Revert selected units. Ctrl+U', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onUndo));
        btnTblUndo.Layout.Column = 6;
        applyColorScheme(btnTblUndo, figColor);

        btnTblReset = uibutton(ctrlRow,'Text','Reset', ...
            'Tooltip','Reset all edits. Ctrl+R', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onReset));
        btnTblReset.Layout.Column = 7;
        applyColorScheme(btnTblReset, figColor);

        btnTblNotes = uibutton(ctrlRow,'Text','Notes', ...
            'Tooltip','Edit the unit''s note. Ctrl+E', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onEditNotes));
        btnTblNotes.Layout.Column = 8;
        applyColorScheme(btnTblNotes, figColor);

        btnTblDeselect = uibutton(ctrlRow,'Text','Deselect', ...
            'Tooltip','Untick all selected units. Ctrl+H', ...
            'ButtonPushedFcn',@(b,e) tblWindowAct(@onDeselect));
        btnTblDeselect.Layout.Column = 9;
        applyColorScheme(btnTblDeselect, figColor);

        % --- Row 3: the table itself -----------------------------
        tableContentPanel = uipanel(winGrid, ...
            'BorderType','none','BackgroundColor',figColor);
        tableContentPanel.Layout.Row = 3;
        applyColorScheme(tableContentPanel, figColor);

        % Render the table into the panel. plotTable accepts a
        % parent panel as its argument, sets tableHandle for the
        % paintSelectedTableRows hook, and wires onTableCellSelect.
        plotTable(tableContentPanel);

        % When the user closes the figure (X button), clear our
        % closure handles so the next "Open table" click spawns a
        % fresh window instead of trying to reuse a dead handle.
        tableFigure.CloseRequestFcn = @(src,evt) onTableFigureClose(src);
    end

    function tblWindowAct(actionFcn)
        % Run a Spike Group Controller action from inside the
        % table window, then bring the table data + row colours
        % back into sync. Most actions mutate groupList /
        % unitIsolation / sortedRes, so we re-render the table
        % rather than try to patch in place. Selection state lives
        % in chkBoxSelect, which the actions update directly.
        try, actionFcn(); catch ME, fprintf('table-window action failed: %s\n', ME.message); end
        try, refreshTableWindow(); catch, end
    end

    function onTableFigureClose(src)
        try
            if ~isempty(tableHandle) && isgraphics(tableHandle) && ...
                    isequal(ancestor(tableHandle,'figure'), src)
                tableHandle = [];
            end
        catch
        end
        tableContentPanel = [];
        tableFigure       = [];
        try, delete(src); catch, end
    end

    function refreshTableWindow()
        % Re-render the standalone table from the current curation
        % state (groupList, channelList, isolation, notes, ...).
        % Note: this rebuild does reset the user's column sort -- it
        % only fires when the user explicitly clicks Refresh, so
        % that's the explicit trade-off (data is fresher, sort
        % gets cleared).
        if ~isempty(tableContentPanel) && isgraphics(tableContentPanel)
            plotTable(tableContentPanel);
        end
    end

    function showTableHelp()
        % Pop a small help dialog that explains the curation table.
        % Anchored on tableFigure when it's alive, otherwise
        % parentFig, so it always sits on top of the right window.
        if ~isempty(tableFigure) && isgraphics(tableFigure)
            anchor = tableFigure;
        else
            anchor = parentFig;
        end
        msg = sprintf([ ...
            'Curation Table\n' ...
            '\n' ...
            'Select a unit (toggle Vis + Select):\n' ...
            '  Click a row to focus it, then press SPACE.\n' ...
            '  Selected rows are tinted with the unit''s group colour.\n' ...
            '  Press SPACE again on the same row to deselect it.\n' ...
            '\n' ...
            'Set the merge target (radio button equivalent):\n' ...
            '  Just click any cell in the row -- the row your cursor\n' ...
            '  is on becomes the target. Merge into target via the\n' ...
            '  Merge button or Ctrl+M.\n' ...
            '\n' ...
            'Edit a note:\n' ...
            '  Double-click the Notes cell, type your text, then\n' ...
            '  press Enter or click outside the cell to commit.\n' ...
            '  Notes are saved with the rest of the curation when\n' ...
            '  you press Ctrl+S (Save).\n' ...
            '\n' ...
            'Sort:\n' ...
            '  Click any column header to sort by that column.\n' ...
            '  Click again to flip ascending/descending.\n' ...
            '  Click another header to sort by it instead.\n' ...
            '  Numeric columns sort numerically (not as strings).\n' ...
            '  The Refresh button and curation actions rebuild the\n' ...
            '  table and reset the sort.\n' ...
            '\n' ...
            'Columns:\n' ...
            '  G        Internal unit handle (group index).\n' ...
            '  ID       Unit label after merges.\n' ...
            '  Ch       Main channel.\n' ...
            '  Rate     Mean firing rate (Hz).\n' ...
            '  ISI      Refractory-period violation rate (%%).\n' ...
            '  SNR      Signal-to-noise ratio.\n' ...
            '  ISO      Isolation tier (SUA+/SUA/MUA+/MUA/NA).\n' ...
            '  Stable   Stable-window length as a percentage of\n' ...
            '           the recording.\n' ...
            '  Notes    Free-text notes; double-click to edit.\n' ...
            '\n' ...
            'Numeric columns are shown to two decimals; sort still\n' ...
            'uses native numeric ordering.\n' ...
            '\n' ...
            'Removed units: italic grey, ISO = "REM".\n' ...
            'Bring them back with Ctrl+U (Undo) or Ctrl+R (Reset).\n' ...
            '\n' ...
            'Keyboard (with the table window focused):\n' ...
            '  SPACE                  toggle Vis+Select on focused row\n' ...
            '  Ctrl+D                 Remove\n' ...
            '  Ctrl+M                 Merge into target\n' ...
            '  Ctrl+T                 Split\n' ...
            '  Ctrl+L                 Limit (auto-trim)\n' ...
            '  Ctrl+K                 cycle Isolation tier\n' ...
            '  Ctrl+U                 Undo selected\n' ...
            '  Ctrl+R                 Reset all edits\n' ...
            '  Ctrl+H                 Deselect all\n' ...
            '  Ctrl+E                 edit Notes (radio-selected unit)\n' ...
            '  Ctrl+S                 Save curation\n' ...
            '  Ctrl+P / Ctrl+N        previous / next unit\n' ...
            '  Ctrl+Up / Ctrl+Down    Auto walker\n' ...
            '  Ctrl+[ / Ctrl+]        Similarity walker\n' ...
            '  Ctrl+Left / Ctrl+Right page back / forward\n' ...
            '  Ctrl+/                 show full shortcuts list\n' ...
            ]);
        try
            uialert(anchor, msg, 'Curation Table Help', 'Icon', 'info');
        catch
            % Fallback for older MATLAB releases without uialert.
            msgbox(msg, 'Curation Table Help', 'help');
        end
    end

    function uidx = tableUnitIdx(gidVal, labVal)
        % Map a curation-table row back to its unit array index. Column 1
        % is the contiguous Group ID, column 2 the group label. The label
        % is the stable per-unit identifier, so resolve by it first (this
        % stays correct even if the page changed while the table was open);
        % fall back to the Group ID. NaN when neither resolves (e.g. a
        % removed unit, whose label is NaN).
        uidx = NaN;
        if isnumeric(labVal) && isscalar(labVal) && isfinite(labVal)
            f = find(groupList(:) == labVal, 1);
            if ~isempty(f), uidx = f; return; end
        end
        if isnumeric(gidVal) && isscalar(gidVal) && isfinite(gidVal) ...
                && gidVal >= 1 && gidVal <= numel(displayedGroups)
            uidx = displayedGroups(round(gidVal));
        end
    end

    function onTableCellSelect(src, evt)
        % Click handler for the curation table. The model:
        %   * Clicking ANY cell sets that row's unit as the merge
        %     target (mirrors the radio button in the main GUI).
        %     We don't touch src.Data here -- assigning into the
        %     Data property re-applies the underlying row order
        %     and would silently undo a column-header sort, so
        %     setMergeTarget alone is the entire effect of a
        %     click. MATLAB's native focused-cell highlight
        %     shows the user which row their cursor sits on.
        %   * Vis / Select toggling is reserved for the SPACE key
        %     (see toggleTableFocusedSelection), so a sort or a
        %     casual click never flips a unit's selection.
        try
            if isempty(evt.Indices), return; end
            row = evt.Indices(1);
            if row < 1, return; end
            if isa(src.Data, 'table')
                if row > height(src.Data), return; end
            else
                if row > size(src.Data, 1), return; end
            end
            unitIdx = tableUnitIdx(src.Data{row, 1}, src.Data{row, 2});
            if isempty(unitIdx) || ~isfinite(unitIdx) || ...
                    unitIdx < 1 || unitIdx > numGroups
                return;
            end
            setMergeTarget(unitIdx);
        catch
        end
    end

    function toggleUnitSelect(unitIdx)
        % Flip Vis / Select for unitIdx and route the radio + page.
        % Shared between Idx-column click and Space-key handler.
        if isempty(unitIdx) || ~isfinite(unitIdx) || ...
                unitIdx < 1 || unitIdx > numGroups
            return;
        end
        if chkBoxSelect(unitIdx)
            chkBoxSelect(unitIdx) = 0;
            chkBoxVis(unitIdx)    = 0;
            selected_order(selected_order == unitIdx) = [];
            if unitIdx <= numel(groupSelectCheckboxes) && ...
                    isvalid(groupSelectCheckboxes(unitIdx)) && ...
                    isprop(groupSelectCheckboxes(unitIdx),'Value')
                groupSelectCheckboxes(unitIdx).Value = false;
                groupVisCheckboxes(unitIdx).Value    = false;
            end
        else
            chkBoxSelect(unitIdx) = 1;
            chkBoxVis(unitIdx)    = 1;
            if ~ismember(unitIdx, selected_order)
                selected_order(end+1) = unitIdx;
            end
            visNow = currentPageVisibleGroups();
            if ~ismember(unitIdx, visNow)
                posUnit = find(displayedGroups == unitIdx, 1);
                if isempty(posUnit), posUnit = 1; end
                tgtPage = ceil(posUnit / PAGE_SIZE);
                if tgtPage ~= currentPage
                    currentPage = tgtPage;
                    onChangePage(0);
                end
            end
            if unitIdx <= numel(groupSelectCheckboxes) && ...
                    isvalid(groupSelectCheckboxes(unitIdx)) && ...
                    isprop(groupSelectCheckboxes(unitIdx),'Value')
                groupSelectCheckboxes(unitIdx).Value = true;
                groupVisCheckboxes(unitIdx).Value    = true;
            end
            if unitIdx <= numel(groupRadioButtons) && ...
                    isvalid(groupRadioButtons(unitIdx)) && ...
                    isprop(groupRadioButtons(unitIdx),'Value')
                groupRadioButtons(unitIdx).Value = true;
            end
        end
    end

    function setMergeTarget(unitIdx)
        % Set the merge-target radio for unitIdx. The radio buttons
        % live in a uibuttongroup so MATLAB clears the others
        % automatically when we set one.
        if isempty(unitIdx) || ~isfinite(unitIdx) || ...
                unitIdx < 1 || unitIdx > numGroups
            return;
        end
        if unitIdx <= numel(groupRadioButtons) && ...
                isgraphics(groupRadioButtons(unitIdx)) && ...
                isvalid(groupRadioButtons(unitIdx)) && ...
                isprop(groupRadioButtons(unitIdx),'Value')
            groupRadioButtons(unitIdx).Value = true;
        end
    end

    function onTableCellEdit(src, evt)
        % Notes column is editable in-place; write through to
        % unitNotes so the change persists into Save. The unit is
        % resolved from the row's Group ID / label (tableUnitIdx) so a
        % post-sort reorder routes the edit to the correct unit.
        try
            if isempty(evt.Indices), return; end
            row = evt.Indices(1);
            col = evt.Indices(2);
            % Notes is always the last column (col 9 here).
            if isa(src.Data, 'table')
                nCols   = width(src.Data);
            else
                nCols   = size(src.Data, 2);
            end
            unitIdx = tableUnitIdx(src.Data{row, 1}, src.Data{row, 2});
            if col ~= nCols, return; end
            if isempty(unitIdx) || ~isfinite(unitIdx) || ...
                    unitIdx < 1 || unitIdx > numGroups
                return;
            end
            if numel(unitNotes) < unitIdx
                unitNotes(end+1:unitIdx) = {''};
            end
            unitNotes{unitIdx} = char(evt.NewData);
        catch
        end
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
                        multiPlotPanels(k).Title = 'Amp Distribution';
                        plotAmpDistribution(multiPlotPanels(k))
                    case 7
                        multiPlotPanels(k).Title = 'Channel PC';
                        plotChannelPCA(multiPlotPanels(k))
                    case 8
                        multiPlotPanels(k).Title = 'Trace';
                        plot_on_Signal(multiPlotPanels(k))
                    case 9
                        multiPlotPanels(k).Title = 'Features';
                        plotFeatures(multiPlotPanels(k))
                end
                addMultiPlotGripHandles(k);
            end
        end

        lastMultiPlotOrder = multiplePlotOrder;
        reScaledFlag = false;
    end

    function addMultiPlotGripHandles(k)
        % Two small overlay buttons live on each multi-plot cell so the
        % user can rearrange and resize cells directly. The swap button
        % (top-right) marks the cell on first click and swaps with the
        % next click. The resize button (bottom-right) enters drag mode
        % so the next mouse motion drives gridMultiple's row height /
        % column width for that cell, committed on the next click.
        % SizeChangedFcn keeps both buttons pinned to the corners
        % whenever the cell or the figure resizes -- and we chain any
        % pre-existing SizeChangedFcn so we don't clobber other layout
        % handlers the plot helpers may have installed.
        if ~isvalid(multiPlotPanels(k)), return; end
        panel = multiPlotPanels(k);
        try
            sz = panel.InnerPosition(3:4);
        catch
            sz = [160 80];
        end
        if isempty(sz) || any(sz <= 0), sz = [160 80]; end

        swapColor = [0.55 0.55 0.55];
        if multiSwapPendingK == k
            swapColor = [1 0.55 0];   % highlight when this cell is the pending source
        end

        swapBtn = uibutton(panel,'Text',char(8644),...
            'FontSize',8,...
            'BackgroundColor',swapColor,...
            'Tooltip','Click to mark for swap; click another cell''s grip to swap.',...
            'ButtonPushedFcn',@(~,~) onMultiSwapClicked(k));

        resizeBtn = uibutton(panel,'Text',char(8690),...
            'FontSize',8,...
            'BackgroundColor',[0.55 0.55 0.55],...
            'Tooltip','Click to start drag-resize, move mouse, click again to commit.',...
            'ButtonPushedFcn',@(~,~) onMultiResizeClicked(k));

        repositionMultiGripHandles(panel, swapBtn, resizeBtn);

        prevFcn = panel.SizeChangedFcn;
        panel.SizeChangedFcn = @(src,evt) chainedMultiGripResize(src, evt, prevFcn, swapBtn, resizeBtn);
    end

    function chainedMultiGripResize(panel, evt, prevFcn, swapBtn, resizeBtn)
        repositionMultiGripHandles(panel, swapBtn, resizeBtn);
        if ~isempty(prevFcn)
            try
                if isa(prevFcn, 'function_handle')
                    prevFcn(panel, evt);
                end
            catch
            end
        end
    end

    function repositionMultiGripHandles(panel, swapBtn, resizeBtn)
        if ~isvalid(panel), return; end
        try
            sz = panel.InnerPosition(3:4);
        catch
            sz = [160 80];
        end
        if isempty(sz) || any(sz <= 0), sz = [160 80]; end
        btnW = 14; btnH = 14;     % 20% smaller than the original 18 px
        if isvalid(swapBtn)
            swapBtn.Position = [max(1,sz(1)-btnW-2), max(1,sz(2)-btnH-2), btnW, btnH];
        end
        if isvalid(resizeBtn)
            resizeBtn.Position = [max(1,sz(1)-btnW-2), 2, btnW, btnH];
        end
    end

    function onMultiSwapClicked(k)
        if multiSwapPendingK == 0
            multiSwapPendingK = k;
            if isvalid(multiPlotPanels(k))
                multiPlotPanels(k).BorderColor = [1 0.55 0];
                multiPlotPanels(k).BorderType  = 'line';
            end
        elseif multiSwapPendingK == k
            if isvalid(multiPlotPanels(k))
                multiPlotPanels(k).BorderColor = 1-figColor;
            end
            multiSwapPendingK = 0;
        else
            a = multiSwapPendingK;
            b = k;
            multiSwapPendingK = 0;
            swapMultiPanels(a, b);
        end
    end

    function swapMultiPanels(a, b)
        if a < 1 || b < 1 || a == b, return; end
        if a > numel(multiplePlotOrder) || b > numel(multiplePlotOrder), return; end
        tmp = multiplePlotOrder(a);
        multiplePlotOrder(a) = multiplePlotOrder(b);
        multiplePlotOrder(b) = tmp;
        % Force a clean rebuild so the panels regenerate in the new
        % order; the lastMultiPlotOrder reset is the existing rebuild
        % trigger inside plotMultiple.
        lastMultiPlotOrder = [];
        plotMultiple();
    end

    function onMultiResizeClicked(k)
        if multiResizeDrag.active
            endMultiResizeDrag();
            return;
        end
        if isempty(gridMultiple) || ~isvalid(gridMultiple), return; end
        if ~isvalid(multiPlotPanels(k)), return; end

        pnl = multiPlotPanels(k);
        cellPx = getpixelposition(pnl);

        multiResizeDrag.active     = true;
        multiResizeDrag.k          = k;
        multiResizeDrag.startMouse = parentFig.CurrentPoint;
        multiResizeDrag.cellW      = cellPx(3);
        multiResizeDrag.cellH      = cellPx(4);
        multiResizeDrag.rowIdx     = pnl.Layout.Row(1);
        multiResizeDrag.colIdx     = pnl.Layout.Column(1);
        multiResizeDrag.startRH    = gridMultiple.RowHeight;
        multiResizeDrag.startCW    = gridMultiple.ColumnWidth;
        multiResizeDrag.origMotionFcn = parentFig.WindowButtonMotionFcn;
        multiResizeDrag.origDownFcn   = parentFig.WindowButtonDownFcn;

        parentFig.WindowButtonMotionFcn = @(~,~) onMultiResizeMotion();
        parentFig.WindowButtonDownFcn   = @(~,~) endMultiResizeDrag();
        try, parentFig.Pointer = 'fleur'; catch, end
    end

    function onMultiResizeMotion()
        if ~multiResizeDrag.active || isempty(gridMultiple) || ~isvalid(gridMultiple)
            return;
        end
        cur = parentFig.CurrentPoint;
        deltaX = cur(1) - multiResizeDrag.startMouse(1);
        deltaY = cur(2) - multiResizeDrag.startMouse(2);
        % In figure coords y grows upward, so dragging the bottom-right
        % grip "down" (negative deltaY) shrinks height, "up" grows it.
        newW = max(60, multiResizeDrag.cellW + deltaX);
        newH = max(60, multiResizeDrag.cellH - deltaY);
        rowH = multiResizeDrag.startRH;
        colW = multiResizeDrag.startCW;
        if multiResizeDrag.rowIdx >= 1 && multiResizeDrag.rowIdx <= numel(rowH)
            rowH{multiResizeDrag.rowIdx} = newH;
        end
        if multiResizeDrag.colIdx >= 1 && multiResizeDrag.colIdx <= numel(colW)
            colW{multiResizeDrag.colIdx} = newW;
        end
        try
            gridMultiple.RowHeight   = rowH;
            gridMultiple.ColumnWidth = colW;
        catch
        end
    end

    function endMultiResizeDrag()
        if ~multiResizeDrag.active, return; end
        multiResizeDrag.active = false;
        try, parentFig.WindowButtonMotionFcn = multiResizeDrag.origMotionFcn; catch, end
        try, parentFig.WindowButtonDownFcn   = multiResizeDrag.origDownFcn;   catch, end
        try, parentFig.Pointer = 'arrow'; catch, end

        % Re-run plotMultiple so each cell's plot regenerates against
        % its new pixel dimensions. The rebuild guard inside
        % plotMultiple is keyed on lastMultiPlotOrder + reScaledFlag,
        % both unchanged here, so gridMultiple is preserved (our pixel
        % RowHeight/ColumnWidth survive) while the per-cell dispatch
        % loop still fires and redraws axes/legends to fit.
        try, plotMultiple(); catch, end
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
            % colorMapAll is built once at load with the original
            % numGroups; Split can grow numGroups beyond that, so
            % wrap on the colour map's actual length, not numGroups.
            nMap = size(colorMapAll,1);
            if nMap < 1
                c = [0.6 0.6 0.6];
                return;
            end
            idx = mod(round(lbl)-1, nMap) + 1;
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

        % Restrict yDataRange to channels actually covered by the
        % selected neurons' waveform footprints. Falls back to the
        % slider window if nothing is selected.
        if ~isempty(selected_order)
            selChans = unique(channelPlot(selected_order, :));
            selChans = selChans(selChans >= 1 & selChans <= sy);
            yDataRange = sort(selChans(:))';
            yDataRange = yDataRange(chan_wave_inclusion(yDataRange));
        else
            yDataRange = max(1, yWindowStart) : min(yWindowEnd, sy);
            yDataRange = yDataRange(chan_wave_inclusion(yDataRange));
        end

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

        % Feature plot spans the FULL recording on its time (X) axis -
        % a small visible window (e.g. xStepDropdown = 0.1 s) is great
        % for the Trace view but useless for the Feature view, where
        % drift / quality-over-time is the whole point. We still keep
        % the y-window slider for the channel (Z) axis.
        xDataRange = round(1:num_Samples);
        yDataRange = max(1, yWindowStart) : min(yWindowEnd, numChannels);
        yDataRange = yDataRange(chan_wave_inclusion(yDataRange));

        % Optimized spike finding
        inRange_mask = sortedRes.spike_idx > xDataRange(1) & sortedRes.spike_idx < xDataRange(end);
        inRange_idx = find(inRange_mask);

        unifiedLabels = sortedRes.unifiedLabels(inRange_idx);
        spike_idx = sortedRes.spike_idx(inRange_idx);
        spike_features = sortedRes.features(inRange_idx,:);
        % Pad to at least 2 feature columns. Some configs (e.g. when the
        % sorter only kept PC1) leave a single-column features matrix,
        % and the scatter3 below indexes column 2 unconditionally.
        if size(spike_features, 2) < 2
            spike_features = [spike_features, zeros(size(spike_features,1), 1)];
        end
        norm_spike_features = mapminmax(spike_features', -.5, .5)';
        if size(norm_spike_features, 2) < 2
            norm_spike_features = [norm_spike_features, zeros(size(norm_spike_features,1), 1)];
        end

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
                    if isnan(channelList(k)) && k >= 1 && k <= numel(originalChannelList) ...
                            && ~isnan(originalChannelList(k))
                        chBase = originalChannelList(k);
                    else
                        chBase = channelList(k);
                    end
                    if isnan(chBase), chBase = 1; end
                    scatter3(plotAxis, plot_spike_idx/cfg.samplingFrequency, ...
                        norm_spike_features(scatterIdx,1), ...
                        chBase + norm_spike_features(scatterIdx,2),...
                        25,'markerfacecolor',c,'markeredgecolor','none')
                end
            end
        end

        hold(plotAxis,'off')


        axis(plotAxis,'tight',[1/cfg.samplingFrequency trialLength -.5 .5 yWindowStart-.5 yWindowEnd+.5 ]);
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
            % NaN-marked units have channelPlot = NaN. Fall back to the
            % unit's pre-curation channel so the y-window auto-scroll
            % still lands on the right band when the user re-ticks an
            % undone / removed group.
            if all(isnan(chVals)) && lastSel >= 1 && lastSel <= numel(originalChannelList) ...
                    && ~isnan(originalChannelList(lastSel))
                chVals = originalChannelList(lastSel) + (-numChannelPlot:numChannelPlot);
            end
            minVal  = min(chVals) - 1;
            maxVal  = max(chVals) + 1;
            if ~isempty(minVal) && ~isnan(minVal) && ~isnan(maxVal)
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
        % originalUnifiedLabels is preserved at load time, so even when
        % merging has rewritten sortedRes.unifiedLabels we can still
        % look up each constituent unit's original spike train. This is
        % what lets us plot every merged-into constituent as a separate
        % waveform (same colour, same slot).
        if exist('originalUnifiedLabels','var') && numel(originalUnifiedLabels) == numel(sortedRes.unifiedLabels)
            origLabels = originalUnifiedLabels(inRange);
        else
            origLabels = sortedRes.unifiedLabels(inRange);
        end %#ok<NODEF>
        % Spikes that the Auto cleanup (Phase 4c) flagged as
        % contaminated are now sortedRes.unifiedLabels == -1. Exclude
        % them from waveform display so the user sees the cleaned
        % unit, while originalUnifiedLabels still keeps them for
        % Reset / Undo.
        currLabels = sortedRes.unifiedLabels(inRange);
        notMarked  = currLabels ~= -1;
        [~, spike_idx] = ismember(spike_idx, xDataRange);

        % Cache both the raw viewport read AND its bandpass-filtered
        % counterpart so that cosmetic redraws (scale / alpha / line
        % width sliders, group selection toggles, ...) reuse the heavy
        % filtfilt result instead of recomputing it on every plotOption.
        % Names are prefixed wf* to avoid clashing with the parent
        % function's `cachedFilteredData` persistent (which is owned by
        % the keypress / Trace path).
        persistent wfCachedRaw wfCachedFiltered wfLastXRange wfLastYRange wfLastBandpass
        rangeChanged = isempty(wfCachedRaw) || ...
                       ~isequal(wfLastXRange, xDataRange) || ...
                       ~isequal(wfLastYRange, yDataRange);
        if rangeChanged
            wfCachedRaw      = double(mappedData.Data.data(channel_mapping(yDataRange), xDataRange));
            wfCachedFiltered = [];
            wfLastXRange     = xDataRange;
            wfLastYRange     = yDataRange;
        end

        if strcmp(plotFilterType,'filtered')
            if isempty(wfCachedFiltered) || ~isequal(wfLastBandpass, cfg.bandpass)
                wfCachedFiltered = bandpass_filter_GUI(wfCachedRaw', ...
                    cfg.bandpass, cfg.samplingFrequency)';
                wfLastBandpass = cfg.bandpass;
            end
            inputData = wfCachedFiltered;
        else
            inputData = wfCachedRaw;
        end

        if isempty(maxSignal)
            sampleSize = min(size(mappedData.Data.data,2), 10000);
            sampleIdx  = randsample(size(mappedData.Data.data,2), sampleSize);
            maxSignal  = max(abs(double(mappedData.Data.data(channel_mapping(chan_wave_inclusion), sampleIdx))), [], 'all');
        end
        normIn = inputData ./ maxSignal;

        chanPlotList = [];
        hold(plotAxis,'on');

        % Reset the fast-path cache for this redraw. We only enable the
        % fast path in single-plot mode (plotType == 'Waveform'); the
        % Multiple layout would need per-panel handle bookkeeping that
        % isn't worth the complexity for that view.
        canFastPath = ~strcmp(plotType,'Multiple');
        if canFastPath
            wfPlotCache = struct('lines', {{}}, 'unscaledY', {{}}, ...
                                 'yOff', [], 'isMain', logical([]), ...
                                 'rgb', zeros(0,3), 'axis', plotAxis);
        end

        % Auto-expand selected_order to include every group that merged
        % into the same effective label. This way, ticking the primary
        % of a merged unit also brings up its absorbed constituents'
        % waveforms as separate traces (same colour, same slot).
        expandedOrder = [];
        seenK = false(numGroups, 1);
        for kk = selected_order(:)'
            if kk < 1 || kk > numGroups || seenK(kk), continue; end
            seenK(kk) = true;
            expandedOrder(end+1,1) = kk;        %#ok<AGROW>
            ekk = effLabel(kk);
            if isnan(ekk), continue; end
            for kj = 1:numGroups
                if seenK(kj), continue; end
                if isequal(effLabel(kj), ekk)
                    seenK(kj) = true;
                    expandedOrder(end+1,1) = kj;%#ok<AGROW>
                end
            end
        end

        % Slot bookkeeping.
        %   Overlay OFF: each unit (including each constituent of a
        %   merged group) gets its own horizontal slot, so the user can
        %   eyeball whether a merge was reasonable side-by-side.
        %   Overlay ON: every unit collapses into a single slot and
        %   only colour distinguishes them.
        if isempty(expandedOrder)
            ekVec = [];
        else
            ekVec = effLabelVec(expandedOrder);
        end
        nUnits = numel(expandedOrder);
        if plotOverlay
            slotPerK  = zeros(nUnits, 1);
            uniqPlots = 1;
        else
            slotPerK  = (0:nUnits-1)';
            uniqPlots = max(nUnits, 1);
        end

        plotCount = 0;            % preserved for empty-plot bookkeeping
        minAx = inf; maxAx = -inf;
        max_YAx = -inf; min_YAx = inf;

        if any(spike_idx) && ~isempty(expandedOrder)
            for kIdx = 1:numel(expandedOrder)
                k    = expandedOrder(kIdx);
                ek   = ekVec(kIdx);
                slot = slotPerK(kIdx);

                % Spike-row selection for this unit:
                %   * Merged constituents: filter via originalUnifiedLabels
                %     == originalGroupList(k) so each unit's pre-merge
                %     spikes stay distinct, AND restrict to rows whose
                %     CURRENT label still equals groupList(k) so a
                %     subsequent Split that moved some spikes to a new
                %     label doesn't leak them back into the parent.
                %   * Split-created units: originalGroupList(k) is NaN
                %     by design (the unit didn't exist pre-curation),
                %     so the originalUnifiedLabels-based query returns
                %     nothing. Fall through to the current label.
                %   * Phase-4c stripped rows (unifiedLabels == -1) are
                %     dropped via notMarked.
                origLab        = originalGroupList(k);
                currLab        = groupList(k);
                if isnan(currLab)
                    plot_spike_idx = [];
                else
                    currLabels = sortedRes.unifiedLabels(inRange);
                    if isnan(origLab)
                        plot_spike_idx = spike_idx((currLabels == currLab) & notMarked);
                    else
                        plot_spike_idx = spike_idx( ...
                            (origLabels == origLab) & ...
                            (currLabels == currLab) & notMarked);
                    end
                end

                % Channel-footprint lookup. For NaN-marked (removed)
                % units channelList(k) is NaN, so channelPlot(k,:) is
                % all-NaN and ismember would discard everything --
                % the unit's gray waveform would never render. Fall
                % back to the unit's pre-curation channel so it can
                % still be inspected.
                if k >= 1 && k <= numGroups && ~isnan(channelList(k))
                    chBase = channelList(k);
                elseif k >= 1 && k <= numel(originalChannelList)
                    chBase = originalChannelList(k);
                else
                    chBase = NaN;
                end
                if isnan(chBase) || chBase < 1
                    continue;
                end
                chPlot = chBase + (-numChannelPlot:numChannelPlot);
                chPlot = chPlot(ismember(chPlot, yDataRange));
                if ~any(chPlot)
                    continue;
                end
                % Colour rule:
                %   * NaN-marked (removed) units always draw in a flat
                %     gray so it is obvious they have been deselected
                %     by Auto / manual Remove, even though the user
                %     can still inspect their waveform.
                %   * Otherwise colour by the effective label so every
                %     constituent of a merged group draws in the
                %     merge colour.
                if k >= 1 && k <= numGroups && isnan(groupList(k))
                    c = [0.55 0.55 0.55];
                else
                    c = getGroupColor(ek);
                end

                spike_indices = plot_spike_idx(:);
                numSpikes     = numel(spike_indices);
                if ~plotMean && numSpikes > numWaveforms
                    spike_indices = spike_indices(randperm(numSpikes, numWaveforms));
                    numSpikes     = numel(spike_indices);
                end

                if isempty(spike_indices)
                    all_indices = zeros(0, numel(spike_Xaxis));
                else
                    all_indices = spike_indices + spike_Xaxis;
                end
                numChPlots  = numel(chPlot);

                if ~plotMean
                    Waveforms = zeros(numSpikes, numChPlots, numel(spike_Xaxis));
                else
                    Waveforms = nan(numSpikes, numChPlots, numel(spike_Xaxis));
                end

                % Vectorised read of every in-window snippet for this
                % unit: slice the cached normIn array per channel in a
                % single call instead of touching each (spike, sample)
                % cell individually. Replaces a numChPlots × numSpikes
                % scalar inner loop.
                wfLen = numel(spike_Xaxis);
                validRow = all(all_indices > 0 & all_indices <= size(inputData,2), 2);
                if any(validRow)
                    idxValid = all_indices(validRow, :);
                    nValid   = size(idxValid,1);
                    for cIdx = 1:numChPlots
                        chan  = chPlot(cIdx);
                        rowIdx = find(yDataRange == chan, 1);
                        if isempty(rowIdx), continue; end
                        slab  = ampScale * normIn(rowIdx, idxValid);    % nValid × wfLen
                        Waveforms(validRow, cIdx, :) = reshape(slab, nValid, 1, wfLen);
                    end
                end

                % Fallback for low-rate units: when the visible window
                % contains fewer in-window spikes than the user asked
                % for, pull additional spikes from anywhere in the
                % recording. We read each extra snippet directly from
                % the memory-mapped data (with a small filter pad to
                % avoid edge ringing) and append them to Waveforms.
                % Triggered only when needed, so the fast path on
                % high-rate units is unaffected.
                needed = max(0, numWaveforms - numSpikes);
                if ~plotMean && needed > 0
                    % Split-created units have origLab = NaN; fall back
                    % to the current sortedRes label so we can still
                    % find their spikes.
                    if isnan(origLab)
                        usedLocal = spike_idx(currLabels == currLab);
                        lookupLab = currLab;
                        useCurrent = true;
                    else
                        usedLocal = spike_idx(origLabels == origLab);
                        lookupLab = origLab;
                        useCurrent = false;
                    end
                    [extraWf, nAdded] = readExtraWaveforms(lookupLab, ...
                        usedLocal, needed, chPlot, xDataRange, useCurrent);
                    if nAdded > 0
                        Waveforms = cat(1, Waveforms, extraWf);
                        numSpikes = size(Waveforms, 1);
                    end
                end

                % Nothing to draw for this unit even after the fallback:
                % bookkeep the slot extents so other units stay aligned,
                % then move on.
                if size(Waveforms, 1) == 0
                    min_t = waveformXaxis + (uniqPlots/3) * min(xNorm_adj(chPlot));
                    min_t = min_t + 1.1 * slot * range(min_t);
                    max_t = waveformXaxis + (uniqPlots/3) * max(xNorm_adj(chPlot));
                    max_t = max_t + 1.1 * slot * range(max_t);
                    minAx = min(minAx, min(min_t));
                    maxAx = max(maxAx, max(max_t));
                    plotCount = plotCount + 1;
                    continue;
                end

                chanPlotList = [chanPlotList, chPlot];

                if plotMean
                    Waveforms = mean(Waveforms, 1, 'omitnan');
                end

                for cIdx = 1:numChPlots
                    chan = chPlot(cIdx);
                    t = waveformXaxis + (uniqPlots) * xNorm_adj(chan);
                    t = t + 1.1 * slot * range(t);
                    minAx = min(minAx, min(t));
                    maxAx = max(maxAx, max(t));
                    wfRaw   = squeeze(Waveforms(:, cIdx, :));
                    if size(Waveforms,1) == 1
                        wfRaw = wfRaw(:)';   % keep the row orientation for single-trace mean
                    end
                    yChan   = yNorm_adj(chan);
                    wf      = wfRaw + yChan;

                    min_YAx = min(min_YAx, yChan-1);
                    max_YAx = max(max_YAx, yChan+1);
                    isMainCh = (chan == channelList(k));
                    if isMainCh
                        hL = plot(plotAxis, t, wf', 'Color', [c, alphaLevel], 'LineWidth', lineWidth*2);
                    else
                        hL = plot(plotAxis, t, wf', 'Color', [c, alphaLevel], 'LineWidth', lineWidth);
                    end
                    if canFastPath
                        % wfRaw already has ampScale baked in; divide it
                        % out so the cache stores the unit-amplitude
                        % version, ready to be re-multiplied by a new
                        % ampScale on the fast path.
                        if ampScale ~= 0
                            unscaled = wfRaw / ampScale;
                        else
                            unscaled = wfRaw;
                        end
                        wfPlotCache.lines{end+1,1}     = hL;
                        wfPlotCache.unscaledY{end+1,1} = unscaled;
                        wfPlotCache.yOff(end+1,1)      = yChan;
                        wfPlotCache.isMain(end+1,1)    = isMainCh;
                        wfPlotCache.rgb(end+1,:)       = c(:)';
                    end
                    plot(plotAxis, [t(round(end/2)) t(round(end/2))], ...
                        [yChan-1 yChan+1], ...
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

    function [extraWf, nAdded] = readExtraWaveforms(origLab, alreadyUsedLocal, needed, chPlot, xDataRangeIn, useCurrent)
        % Read additional spike waveforms for a low-rate unit when the
        % current viewport does not contain enough spikes. Each chunk
        % is read with a filter pad so bandpass_filter_GUI does not
        % introduce edge ringing inside the spike window. Returns an
        % (nAdded × numel(chPlot) × numel(spike_Xaxis)) array.
        %
        % useCurrent (optional, default false): when true, search
        % sortedRes.unifiedLabels for origLab instead of
        % originalUnifiedLabels. Used for Split-created units, whose
        % originalUnifiedLabels never carried the new label.
        if nargin < 6, useCurrent = false; end
        nAdded  = 0;
        extraWf = zeros(0, numel(chPlot), numel(spike_Xaxis));
        if needed <= 0, return; end

        % All spikes for this unit (regardless of viewport). When
        % we resolve via originalUnifiedLabels we ALSO have to drop
        % rows whose current label is -1 -- those are spikes the
        % CCG cleanup phase stripped, and pulling them in here was
        % the reason a previously-cleaned unit looked "contaminated"
        % again after save+reload (the in-viewport pass was
        % filtering -1 correctly, but if the requested waveform
        % count exceeded the viewport, this helper backfilled
        % from anywhere in the recording -- including the rows
        % the cleanup had stripped).
        if useCurrent
            allSpkK = sortedRes.spike_idx(sortedRes.unifiedLabels == origLab);
        else
            allSpkK = sortedRes.spike_idx(originalUnifiedLabels == origLab & ...
                                          sortedRes.unifiedLabels ~= -1);
        end
        allSpkK = allSpkK(allSpkK > numSpikeSamples & ...
                          allSpkK <= num_Samples - numSpikeSamples);
        if isempty(allSpkK), return; end

        % Drop spikes already drawn from the cached viewport.
        usedGlobal = [];
        if ~isempty(alreadyUsedLocal) && ~isempty(xDataRangeIn)
            valid = alreadyUsedLocal > 0 & alreadyUsedLocal <= numel(xDataRangeIn);
            usedGlobal = xDataRangeIn(alreadyUsedLocal(valid));
        end
        candidates = setdiff(allSpkK(:), usedGlobal(:));
        if isempty(candidates), return; end

        nExtra = min(needed, numel(candidates));
        if numel(candidates) > nExtra
            extraSpikes = candidates(randperm(numel(candidates), nExtra));
        else
            extraSpikes = candidates;
        end

        filterPad = max(numSpikeSamples * 4, 256);
        extraWf   = nan(nExtra, numel(chPlot), numel(spike_Xaxis));
        for sIdx = 1:nExtra
            spk = extraSpikes(sIdx);
            readStart = max(1, spk - numSpikeSamples - filterPad);
            readEnd   = min(num_Samples, spk + numSpikeSamples + filterPad);
            chunk = double(mappedData.Data.data(channel_mapping(chPlot), readStart:readEnd));
            if strcmp(plotFilterType,'filtered')
                chunk = bandpass_filter_GUI(chunk', cfg.bandpass, cfg.samplingFrequency)';
            end
            centerIdx = spk - readStart + 1;
            snipIdx   = centerIdx + spike_Xaxis;
            if all(snipIdx >= 1 & snipIdx <= size(chunk,2))
                extraWf(sIdx, :, :) = ampScale * chunk(:, snipIdx) ./ maxSignal;
            end
        end
        nAdded = nExtra;
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


        plotList = unique(effLabelVec(selected_order));
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

        plotList = unique(effLabelVec(selected_order));
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


    function plotAmpDistribution(panelOption)
        % Per-unit RMS amplitude distribution. For each visible
        % unit we sample up to maxSpikes spikes, read each spike's
        % snippet on the unit's MAIN channel only, crop to a tight
        % ±0.5 ms window around the alignment sample, and take the
        % per-spike RMS (sqrt(mean(snip.^2))). RMS is the energy
        % of the spike inside the window, so it stays meaningful
        % even when a single sample is noisy -- the previous
        % max(abs(...)) variant got fooled by single-sample
        % transients and ranked clearly-bigger units below smaller
        % ones. RMS also remains polarity-agnostic so negative-
        % and positive-going spikes are summarised the same way.
        %
        % All reads use the curated label (sortedRes.unifiedLabels
        % == groupList(k)), so cleaned spikes (label set to -1) are
        % automatically excluded -- the visualisation always
        % reflects the post-curation spike train of each unit.
        %
        % We picked main-channel
        % over multi-channel RMS because:
        %   * It's the standard amplitude metric in spike sorters
        %     (Phy / KiloSort all gate on it).
        %   * Off-footprint channels mostly contribute noise, which
        %     would dilute a tight per-unit distribution.
        %   * One channel per spike => 1/N the memmap I/O of the
        %     full footprint read, so the plot stays interactive
        %     even with many visible units.
        %
        % Each unit's distribution is rendered as a kernel-density
        % patch with FaceAlpha 0.5 so overlapping curves stay
        % legible. Tight, single-mode = clean unit. Long tail or
        % bimodal = drift / contamination / merge candidate.
        if nargin < 1, panelOption = []; end

        % --- Set up axes (mirrors plotISI / plotDensity). --------
        if isempty(panelOption)
            if strcmp(lastPlotted, 'Amp Dist')
                delete(allchild(axChannels));
            else
                delete(allchild(axChannels));
                delete(axChannels);
                axChannels = uipanel(bottomRightGrid,'BorderType','none', ...
                    'BackgroundColor',figColor);
                axChannels.Layout.Column = [1 5];
                axChannels.Layout.Row    = [1 5];
                applyColorScheme(axChannels, figColor);
            end
            plotAxis = uiaxes(axChannels,'Units','normalized','Position',[0 0 1 1]);
            lastPlotted = plotType;
        else
            plotAxis = uiaxes(panelOption,'Units','normalized','Position',[0 0 1 1]);
        end
        applyColorScheme(plotAxis, figColor);
        hold(plotAxis,'on');
        title(plotAxis,'Per-spike RMS distribution','Color',1-figColor);
        % Match the density panel: no axis frame, ticks, or labels.
        % The colour-coded curves and title are enough to read the
        % plot, and a clean axis-less layout sits better in the
        % multi-plot grid.
        set(plotAxis,'Color',figColor, ...
            'XColor','none','YColor','none', ...
            'XTick',[],'YTick',[],'Box','off');
        xlabel(plotAxis,'');
        ylabel(plotAxis,'');

        visUnits = find(chkBoxVis);
        if isempty(visUnits)
            text(plotAxis, 0.5, 0.5, 'No visible units', ...
                'Units','normalized','HorizontalAlignment','center', ...
                'Color', 1-figColor);
            hold(plotAxis,'off');
            return;
        end

        if isempty(mappedData)
            if isfield(cfg,'inputFolder') && isfield(cfg,'outputFolder')
                try
                    mappedData = map_input_file(cfg.fullFilePath, cfg);
                catch
                    text(plotAxis,0.5,0.5,'Data file not accessible', ...
                        'Units','normalized','HorizontalAlignment','center', ...
                        'Color', 1-figColor);
                    hold(plotAxis,'off');
                    return;
                end
            else
                hold(plotAxis,'off');
                return;
            end
        end

        nDataChans   = size(mappedData.Data.data, 1);
        nSamplesData = size(mappedData.Data.data, 2);
        windowSize   = 2*numSpikeSamples + 1;
        relIdx       = -numSpikeSamples:numSpikeSamples;
        % Per-unit sample cap. Bumped to 1000 because the KDE
        % visibly wobbles below ~500 samples; the single-channel
        % fancy-indexed read scales linearly with sample count and
        % stays under ~10 ms per unit at 1000.
        maxSpikes    = 1000;
        % Tight per-spike window for the amplitude estimate. The
        % full 2*numSpikeSamples+1 snippet (~2 ms wide) is too
        % wide -- a noise transient or the tail of a coupled spike
        % at one of the edges drags max-min up and made well-
        % isolated big-amplitude units look smaller than they
        % actually are. We keep only the central +/-0.5 ms slice
        % around the spike's alignment sample, which is where the
        % action-potential peak/trough actually sits.
        centerIdx    = numSpikeSamples + 1;     % spike alignment sample
        halfTight    = round(0.0005 * cfg.samplingFrequency);
        tightStart   = max(1, centerIdx - halfTight);
        tightEnd     = min(windowSize, centerIdx + halfTight);

        unitAmps   = cell(numel(visUnits),1);
        unitColors = nan(numel(visUnits),3);
        for vi = 1:numel(visUnits)
            k = visUnits(vi);
            if k > numel(groupList)   || isnan(groupList(k)),   continue; end
            if k > numel(channelList) || isnan(channelList(k)), continue; end

            kLabel = groupList(k);
            kRows  = find(sortedRes.unifiedLabels == kLabel);
            if numel(kRows) < 5, continue; end

            spkPos = round(double(sortedRes.spike_idx(kRows)));
            edgeOK = spkPos > numSpikeSamples & ...
                     spkPos <= nSamplesData - numSpikeSamples;
            spkPos = spkPos(edgeOK);
            if numel(spkPos) < 5, continue; end

            % Random sample (sorted) so the sample reflects firing
            % across the whole recording, not just the start.
            nSamp = min(numel(spkPos), maxSpikes);
            if numel(spkPos) > nSamp
                rIdx    = randperm(numel(spkPos), nSamp);
                sampPos = sort(spkPos(rIdx));
            else
                sampPos = sort(spkPos);
            end

            % Main-channel only. Check the channel index and
            % memmap-row index for sanity before doing the read.
            chCenter = round(channelList(k));
            if chCenter < 1 || chCenter > numChannels, continue; end
            if chCenter > numel(channel_mapping),      continue; end
            chMap = channel_mapping(chCenter);
            if ~isfinite(chMap) || chMap < 1 || chMap > nDataChans
                continue;
            end
            chMap = round(chMap);

            % Single fancy-indexed read of every (mainCh, sample)
            % we need for this unit, reshaped so row i = spike i's
            % snippet. We build the index matrix as
            % windowSize x nSamp on purpose so MATLAB's column-
            % major flatten emits "spike1's window, then spike2's
            % window, ..." -- if we built it nSamp x windowSize the
            % flatten interleaves values across spikes (one offset
            % from every spike, then the next offset, ...) and the
            % reshape downstream packs values from many different
            % spikes into each "row", which silently scrambled the
            % RMS so amplitude no longer reflected unit size.
            try
                spikeIdxMat = relIdx(:) + sampPos(:)';      % windowSize x nSamp, each column = one spike's window
                allIdx      = spikeIdxMat(:)';              % spike-major flatten
                slab        = double(mappedData.Data.data(chMap, allIdx));
                slab2D      = reshape(slab, windowSize, numel(sampPos))';   % nSamp x windowSize
                tightSlab   = slab2D(:, tightStart:tightEnd);
                ampVals     = sqrt(mean(tightSlab.^2, 2));
            catch
                continue;
            end
            ampVals = ampVals(isfinite(ampVals));
            if numel(ampVals) < 5, continue; end

            unitAmps{vi}     = ampVals;
            unitColors(vi,:) = getGroupColor(kLabel);
        end

        nonEmpty = find(~cellfun(@isempty, unitAmps));
        if isempty(nonEmpty)
            text(plotAxis,0.5,0.5,'Not enough spikes for any visible unit', ...
                'Units','normalized','HorizontalAlignment','center', ...
                'Color', 1-figColor);
            hold(plotAxis,'off');
            return;
        end

        % Common x-grid covering the bulk of every distribution; we
        % use 1st/99th percentiles (with a 5% pad) so an isolated
        % giant outlier on one unit doesn't squish every other
        % curve against the y-axis.
        allValues = vertcat(unitAmps{nonEmpty});
        xMin = prctile(allValues,1);
        xMax = prctile(allValues,99);
        if ~isfinite(xMin) || ~isfinite(xMax) || xMax <= xMin
            xMin = min(allValues);
            xMax = max(allValues);
        end
        if xMax <= xMin, xMax = xMin + 1; end
        pad  = 0.05 * (xMax - xMin);
        xPts = linspace(xMin - pad, xMax + pad, 256);

        for vi = nonEmpty(:)'
            ampVals = unitAmps{vi};
            c       = unitColors(vi,:);
            try
                [f, xi] = ksdensity(ampVals, xPts);
            catch
                continue;
            end
            patch(plotAxis, [xi, fliplr(xi)], [f, zeros(1,numel(f))], ...
                c, 'FaceAlpha', 0.5, 'EdgeColor', c, ...
                'LineWidth', 1.0, 'EdgeAlpha', 0.85);
        end

        axis(plotAxis,'tight');
        box(plotAxis,'off');
        hold(plotAxis,'off');
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

        plotList = unique(effLabelVec(selected_order));
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
                                presenceRatio(spike_idx, cfg.samplingFrequency, trialLength, max(1, round(trialLength/(2*ccgLag))), smoothN);
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
                            % Hide axis lines, ticks and labels.
                            % The bars, title and the interactive
                            % imlines stay visible because they're
                            % drawn as separate graphics objects --
                            % only the axis frame is suppressed.
                            set(tiledAX{iDen,jDen}, 'XColor','none', 'YColor','none', ...
                                'XTick',[], 'YTick',[], 'Box','off');
                            xlabel(tiledAX{iDen,jDen},'');
                            ylabel(tiledAX{iDen,jDen},'');
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
                                scatter(plotAxis,xScatter,yScatter,lineWidth*100,'filled','MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerFaceAlpha',min(alphaLevel*2,1), 'MarkerEdgeAlpha', min(alphaLevel*2,1))
                            else
                                scatter(plotAxis,xScatter,yScatter,lineWidth*50,'filled','MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerFaceAlpha',min(alphaLevel*2,1), 'MarkerEdgeAlpha', min(alphaLevel*2,1))
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
                scatterIdx = find(unifiedLabels == effLabel(k));
                if length(scatterIdx) > 5E4
                    scatterIdx = scatterIdx(randperm(length(scatterIdx), 5E4));
                end
                plot_spike_idx = spike_idx(scatterIdx);
                if ~isempty(plot_spike_idx)
                    c = getGroupColor(groupList(k));
                    scatter(plotAxis, plot_spike_idx/cfg.samplingFrequency,classPolarity * abs(norm_spike_amplitude(scatterIdx)),...
                        lineWidth * 25,'markerfacecolor',c,'markeredgecolor','none','MarkerFaceAlpha',min(alphaLevel*2,1))
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
        if strcmpi(plotType, 'Trace')
            xStepDropdown.Value = '0.1';
        else
            xStepDropdown.Value = 'full';
        end
        updateXWindow(xSlider.Value, xStepDropdown.Value);
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

    function updatePlotOverlay(val)
        % Toggle between overlay (all groups in one column, distinguished
        % by colour) and side-by-side (each group in its own column).
        plotOverlay = logical(val);
        disp(['Overlay = ', mat2str(plotOverlay)]);
        plotOption();
    end


    function preprocessSortedData()
        % Vectorized preprocessing.
        % If a unit has been marked removed (groupList(i) = NaN), fall
        % back to the unit's pre-curation spike train so the displayed
        % ISI never becomes NaN -- the user still sees the original
        % unit quality even after a Remove / Auto sweep, and Reset /
        % Undo / a second Auto pass restore meaningful values.
        for i = 1:numGroups
            groupID_i = groupList(i);
            if isnan(groupID_i)
                origLab     = originalGroupList(i);
                spike_idx_i = sortedRes.spike_idx(originalUnifiedLabels == origLab);
            else
                spike_idx_i = sortedRes.spike_idx(sortedRes.unifiedLabels == groupID_i);
            end
            if isempty(spike_idx_i)
                preprocessed.isiViolation(i,1) = 0;
                preprocessed.firingRate(i,1)   = 0;
            else
                [~, ~, preprocessed.isiViolation(i,1)] = getISIViolations(spike_idx_i, cfg.samplingFrequency, thresholdISI);
                preprocessed.firingRate(i,1) = numel(spike_idx_i)./trialLength;
            end
        end
        preprocessed.logFiringRate = log(preprocessed.firingRate);
        recomputeDisplayedGroups();
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
                    % 200 spikes per group (doubled from 100) gives a
                    % more stable PCA / mean-waveform estimate at a
                    % small extra read cost. The stability is what
                    % feeds the similarity / amp-variance gates that
                    % the Auto pipeline relies on.
                    maxSpikesPerGroup = 200;
                    for jChan = 1:length(channels)
                        spike_idx_i = sortedRes.spike_idx(sortedRes.unifiedLabels == groupList(channels(jChan)));
                        numSamples = min(length(spike_idx_i),maxSpikesPerGroup);
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

    function saveKiaPrefs()
        % Build the prefs struct from current closure state and write
        % it next to the curated output (cfg.outputFolder/RES_Sorted).
        % Called from onSave (so tuning is checkpointed every time
        % the user hits Save) and from onParentClose (so a dirty
        % session that hasn't yet been Saved still gets its prefs
        % captured). Never throws -- prefs are best-effort.
        try
            prefs = struct();
            prefs.plotType       = plotType;
            prefs.plotFilterType = plotFilterType;
            prefs.numWaveforms   = numWaveforms;
            prefs.plotMean       = plotMean;
            prefs.plotOverlay    = plotOverlay;
            prefs.ccgLag         = ccgLag;
            prefs.smoothN        = smoothN;
            prefs.thresholdISI   = thresholdISI;
            prefs.ampScale       = ampScale;
            prefs.lineWidth      = lineWidth;
            prefs.alphaLevel     = alphaLevel;
            prefs.exportType     = exportType;
            % Auto-curation panel state.
            prefs.autoDistVal     = autoDistVal;
            prefs.autoSimVal      = autoSimVal;
            prefs.autoAmpVal      = autoAmpVal;
            prefs.autoIsiVal      = autoIsiVal;
            prefs.autoOverlapFrac = autoOverlapFrac;
            prefs.autoIsiBudget   = autoIsiBudget;
            prefs.autoPCVal       = autoPCVal;
            prefs.autoCCGRatio    = autoCCGRatio;
            prefs.autoLagTightMs           = autoLagTightMs;
            prefs.autoLagConsistencyFrac   = autoLagConsistencyFrac;
            prefs.autoPairType    = autoPairType;

            saveDir = projectPrefsDir;
            if isempty(saveDir) || ~ischar(saveDir)
                % Fallback if cfg.outputFolder went missing somehow.
                saveDir = getKiaPrefsDir();
            end
            if ~exist(saveDir, 'dir'), mkdir(saveDir); end
            save(fullfile(saveDir, 'kiaSort_curate_prefs.mat'), 'prefs');
        catch
            % Pref save failure is non-fatal.
        end
    end

    function onParentClose(src, evt, prevFcn) %#ok<INUSD>
        % Persist user prefs when the figure closes, then chain to any
        % previously-installed CloseRequestFcn (or the default delete).
        saveKiaPrefs();
        % Chain to any previous CloseRequestFcn so we don't strip
        % MATLAB's default close behaviour.
        try
            if ~isempty(prevFcn) && isa(prevFcn,'function_handle')
                feval(prevFcn, src, evt);
            else
                delete(src);
            end
        catch
            try, delete(src); catch, end
        end
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

function d = getKiaPrefsDir()
% Cross-platform location for the curate-GUI's persistent prefs file.
% Falls back to the temp dir if the home directory is somehow not
% writable (e.g. read-only mount).
try
    if ispc
        d = fullfile(getenv('USERPROFILE'), '.kiaSort');
    else
        d = fullfile(getenv('HOME'), '.kiaSort');
    end
    if isempty(d) || strcmp(d, '.kiaSort')
        d = fullfile(tempdir, 'kiaSort');
    end
catch
    d = fullfile(tempdir, 'kiaSort');
end
end


function lbls = umapDbscanCluster(X, targetK)
% Embed X with UMAP (Python package) when available, else PCA, then
% DBSCAN forced to exactly targetK clusters (the most-similar clusters
% are merged when DBSCAN finds more). Falls back to k-means only if
% DBSCAN cannot reach targetK at all.
nComp = min(2, size(X,2));
Y = [];
try
    Y = pythonUMAP(X, nComp);
catch
    Y = [];
end
if isempty(Y) || size(Y,1) ~= size(X,1) || any(~isfinite(Y(:)))
    try
        [~, Y] = pca(X, 'NumComponents', nComp);
    catch
        Y = X;
    end
end
lbls = dbscanToTargetK(Y, targetK);
if numel(unique(lbls)) ~= targetK
    try
        lbls = kmeans(X, targetK, 'Start','plus', 'Replicates', 5, ...
            'Options', statset('MaxIter', 300));
    catch
    end
end
end


function lbls = dbscanToTargetK(Y, targetK)
% DBSCAN epsilon sweep that aims for at least targetK clusters, then
% merges the most-similar clusters down to exactly targetK. Noise (-1)
% is folded into the nearest cluster centroid first. Returns labels
% 1..targetK (or fewer when DBSCAN never reaches targetK).
n = size(Y,1);
minpts = max(5, 2*size(Y,2));
minpts = min(minpts, max(2, n-1));
try
    eps0 = estimate_dbscan_par(Y, minpts);
catch
    eps0 = [];
end
if isempty(eps0) || ~isfinite(eps0) || eps0 <= 0
    dd = pdist(Y(1:min(n,500), :));
    eps0 = median(dd(dd > 0));
    if isempty(eps0) || ~isfinite(eps0) || eps0 <= 0, eps0 = 1; end
end
% Prefer the fewest clusters that is still >= targetK (least merging);
% otherwise the most clusters DBSCAN could produce.
factors  = [1.6 1.3 1.1 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.15];
bestGE = []; bestGEn = inf;
bestAny = []; bestAnyN = 0;
for f = factors
    lb = dbscan(Y, eps0*f, minpts);
    nC = numel(unique(lb(lb > 0)));
    if nC < 1, continue; end
    if nC >= targetK && nC < bestGEn
        bestGEn = nC; bestGE = lb;
    end
    if nC > bestAnyN
        bestAnyN = nC; bestAny = lb;
    end
end
if ~isempty(bestGE)
    lbls = bestGE;
elseif ~isempty(bestAny)
    lbls = bestAny;
else
    lbls = ones(n,1);
end
lbls = foldNoiseToNearest(Y, lbls);
lbls = mergeClustersToK(Y, lbls, targetK);
uc = unique(lbls);
remap = zeros(max(uc), 1);
remap(uc) = 1:numel(uc);
lbls = remap(lbls);
end


function lbls = foldNoiseToNearest(Y, lbls)
% Assign DBSCAN noise points (label <= 0) to the nearest cluster centroid.
noise = lbls <= 0;
if all(noise)
    lbls(:) = 1;
    return;
end
if any(noise)
    uc   = unique(lbls(~noise));
    cent = zeros(numel(uc), size(Y,2));
    for ii = 1:numel(uc)
        cent(ii,:) = mean(Y(lbls == uc(ii), :), 1);
    end
    dn = pdist2(Y(noise,:), cent);
    [~, nn] = min(dn, [], 2);
    lbls(noise) = uc(nn);
end
end


function lbls = mergeClustersToK(Y, lbls, K)
% Merge the two nearest cluster centroids repeatedly until exactly K
% clusters remain (no-op when the count is already <= K).
uc = unique(lbls);
while numel(uc) > K
    cent = zeros(numel(uc), size(Y,2));
    for ii = 1:numel(uc)
        cent(ii,:) = mean(Y(lbls == uc(ii), :), 1);
    end
    DD = squareform(pdist(cent));
    DD(1:size(DD,1)+1:end) = inf;
    [~, idx] = min(DD(:));
    [a, b]   = ind2sub(size(DD), idx);
    lbls(lbls == uc(b)) = uc(a);
    uc = unique(lbls);
end
end