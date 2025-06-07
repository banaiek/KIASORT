function channel_inclusion = channel_selection_gui(cfg, channel_mapping, channel_inclusion, parentPanel, figColor, parentFig)



if isfield(cfg,'outputFolder')
    fileDir = fullfile(cfg.outputFolder,'RES_Samples');
    filename = fullfile(fileDir,'channel_info.mat');
    if exist(filename,"file")
        load(filename, 'channel_inclusion', 'thresh_MAD');
    end
end

numChannels = min(cfg.numChannels, length(channel_mapping));

if isempty(channel_mapping)
    channel_mapping = [1:numChannels]';
end

if isempty(channel_inclusion)
    channel_inclusion = true(numChannels, 1);    
end
thresh_MAD = nan(numChannels, 1);
channel_inclusion_plot = channel_inclusion;
trialLength = 100;
mappedData = [];

% left panels
leftPanel = uipanel(parentPanel, ...
    'Title','Channel Controller', ...
    'BackgroundColor',figColor, ...
    'Scrollable','on');
leftPanel.Layout.Column = [1 2];
leftPanel.Layout.Row    = [2 3];
applyColorScheme(leftPanel, figColor);


leftGrid = uigridlayout(leftPanel, [2 1], ...
    'RowHeight',{'fit','1x'}, ...
    'ColumnWidth',{'fit'}, ...
    'Padding',10, ...
    'RowSpacing',10);
applyColorScheme(leftGrid, figColor);

% top left panel and its buttons
topPanel = uipanel(leftGrid, 'BorderType','none');
topPanel.Layout.Row = 1;
topPanel.Layout.Column = 1;
applyColorScheme(topPanel, figColor);

topButtonLayout = uigridlayout(topPanel,[1 2], ...
    'ColumnWidth',{'fit'}, ...
    'ColumnSpacing',7.5, ...
    'Padding',[0 0 0 0]);
applyColorScheme(topButtonLayout, figColor);


btnSave = uibutton(topButtonLayout,'Text','Save',...
    'ButtonPushedFcn',@(btn,ev)onSave());
btnSave.Layout.Column = 1;
applyColorScheme(btnSave, figColor);



%% Channel check boxes and radio panels
checkPanel = uipanel(leftGrid, ...
    'Title','Channel List', ...
    'BackgroundColor',figColor,...
    'Scrollable','on');
checkPanel.Layout.Row = 2;
checkPanel.Layout.Column = 1;
applyColorScheme(checkPanel, figColor);

totalRows = numChannels + 1;  % +1 for the "All" row


mainCheckGrid = uigridlayout(checkPanel, [totalRows, 3], ...
    'RowHeight', repmat({'fit'},1,totalRows), ...
    'ColumnWidth', repmat({'fit'},1,3), ...
    'RowSpacing',2, ...
    'Scrollable','on',...
    'ColumnSpacing',5, ...
    'Padding',[5 5 5 5]);
applyColorScheme(mainCheckGrid, figColor);

% qrrays for UI handles
lblChannelIDs       = gobjects(numChannels,1);
lblChannelHandles     = gobjects(numChannels,1);
channelCheckboxes    = gobjects(numChannels,1);

% all selection row
lblAll = uilabel(mainCheckGrid,'Text','ID #:',...
    'FontWeight','bold');
lblAll.Layout.Row = 1;
lblAll.Layout.Column = 1;
applyColorScheme(lblAll, figColor);

lblChannel = uilabel(mainCheckGrid,'Text','Channel #:');
lblChannel.Layout.Row = 1;
lblChannel.Layout.Column = 2;
applyColorScheme(lblChannel, figColor);

chkAllChannel = uicheckbox(mainCheckGrid,'Text','Include All',...
    'FontWeight','bold',...
    'Value',true,...
    'ValueChangedFcn',@(cb,ev)onAllChanChecked(cb.Value));
chkAllChannel.Layout.Row = 1;
chkAllChannel.Layout.Column = 3;
applyColorScheme(chkAllChannel, figColor);

%% Channel rows
for i = 1:numChannels
    rowIdx = i + 1;

    % channel list
    lblChannelIDs(i) = uilabel(mainCheckGrid, ...
        'Text', sprintf('Ch: %d', channel_mapping(i)),...
        'FontColor',1-figColor,...
        'FontWeight','bold');
    lblChannelIDs(i).Layout.Row = rowIdx;
    lblChannelIDs(i).Layout.Column = 1;

    % channel list
    lblChannelHandles(i) = uilabel(mainCheckGrid, ...
        'Text', sprintf('Channel %d', i),...
        'FontColor',1-figColor,...
        'FontWeight','bold');
    lblChannelHandles(i).Layout.Row = rowIdx;
    lblChannelHandles(i).Layout.Column = 2;



    % visualization checkboxes
    channelCheckboxes(i) = uicheckbox(mainCheckGrid,...
        'Value',true, 'Text','Include',...
        'FontColor',1-figColor,...
        'ValueChangedFcn',@(cb,ev)onChanCheckboxChanged(i));
    channelCheckboxes(i).Layout.Row = rowIdx;
    channelCheckboxes(i).Layout.Column = 3;
end
refreshLabels();

%% Right plot and control panel
topRightPanel = uipanel(parentPanel,...
    'Title','Visualization Control',...
    'BackgroundColor',figColor);
topRightPanel.Layout.Column = 3;
topRightPanel.Layout.Row    = [1 2];
applyColorScheme(topRightPanel, figColor);

bottomRightPanel = uipanel(parentPanel,...
    'Title','Plots',...
    'BackgroundColor',figColor);
bottomRightPanel.Layout.Column = 3;
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

plotButtonGrid = uigridlayout(plotButtonPanel,[2 7],...
    'ColumnWidth',{'1x','1x','1x','1x','1x','1x','2x'},...
    'RowHeight',{'fit','fit'},...
    'Padding',[0 0 0 0],...
    'RowSpacing', 1,...
    'ColumnSpacing',20);
applyColorScheme(plotButtonGrid, figColor);


plotFilterTypeDD = uidropdown(plotButtonGrid,...
    'Items',{'filtered','raw'},...
    'Value','filtered',...
    'ValueChangedFcn',@(dd,ev)updatePlotFilterType(dd.Value));
plotFilterTypeDD.Layout.Row = 2;
plotFilterTypeDD.Layout.Column = 1;
applyColorScheme(plotFilterTypeDD, figColor);
plotFilterType = 'filtered';

ampSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 1,'MajorTicksMode','auto', 'Limits',[0 10]);
ampSlider.MajorTicks = [0 10];
ampSlider.Layout.Row = 2;
ampSlider.Layout.Column = 2;
applyColorScheme(ampSlider, figColor);
ampSlider.ValueChangedFcn      = @(sld,ev)updateScale(sld.Value);

lineWidthSlider = uislider(plotButtonGrid,'Orientation','horizontal','Value', 1,'MajorTicksMode','auto', 'Limits',[0 5]);
lineWidthSlider.MajorTicks = [0 5];
lineWidthSlider.Layout.Row = 2;
lineWidthSlider.Layout.Column = 3;
applyColorScheme(lineWidthSlider, figColor);
lineWidthSlider.ValueChangedFcn      = @(sld,ev)updateLineWidth(sld.Value);
lineWidth = 1;

hueSlider = uislider(plotButtonGrid,'Value', .1, 'Limits',[0 1]);
hueSlider.MajorTicks = [0 1];
hueSlider.Layout.Row = 2;
hueSlider.Layout.Column = 4;
applyColorScheme(hueSlider, figColor);
hueSlider.ValueChangedFcn      = @(sld,ev)updatePlotLineHue(sld.Value);
plotLineHue = hueSlider.Value;


satSlider = uislider(plotButtonGrid,'Value', .1, 'Limits',[0 1]);
satSlider.MajorTicks = [0 1];
satSlider.Layout.Row = 2;
satSlider.Layout.Column = 5;
applyColorScheme(satSlider, figColor);
satSlider.ValueChangedFcn      = @(sld,ev)updatePlotLineSat(sld.Value);
plotLineSat = satSlider.Value;



valueSlider = uislider(plotButtonGrid,'Value', .7, 'Limits',[0 1]);
valueSlider.MajorTicks = [0 1];
valueSlider.Layout.Row = 2;
valueSlider.Layout.Column = 6;
applyColorScheme(valueSlider, figColor);
valueSlider.ValueChangedFcn      = @(sld,ev)updatePlotLineValue(sld.Value);
plotLineValue = valueSlider.Value;


% labels for buttons

lblFilterType = uilabel(plotButtonGrid,'Text','Data Type:',...
    'HorizontalAlignment','right');
lblFilterType.Layout.Row = 1;
lblFilterType.Layout.Column = 1;
applyColorScheme(lblFilterType, figColor);

lblXAmpScale = uilabel(plotButtonGrid,'Text','Amp. Scale:',...
    'HorizontalAlignment','right');
lblXAmpScale.Layout.Row = 1;
lblXAmpScale.Layout.Column = 2;
applyColorScheme(lblXAmpScale, figColor);

lblXLineWidth = uilabel(plotButtonGrid,'Text','Line Width:',...
    'HorizontalAlignment','right');
lblXLineWidth.Layout.Row = 1;
lblXLineWidth.Layout.Column = 3;
applyColorScheme(lblXLineWidth, figColor);

lblHue = uilabel(plotButtonGrid,'Text','Hue:',...
    'HorizontalAlignment','right');
lblHue.Layout.Row = 1;
lblHue.Layout.Column = 4;
applyColorScheme(lblHue, figColor);


lblsat = uilabel(plotButtonGrid,'Text','Saturation:',...
    'HorizontalAlignment','right');
lblsat.Layout.Row = 1;
lblsat.Layout.Column = 5;
applyColorScheme(lblsat, figColor);


lblValue = uilabel(plotButtonGrid,'Text','Value:',...
    'HorizontalAlignment','right');
lblValue.Layout.Row = 1;
lblValue.Layout.Column = 6;
applyColorScheme(lblValue, figColor);


lblMAD = uilabel(plotButtonGrid,'Text','MAD Thresh:',...
    'FontWeight','bold','HorizontalAlignment','Right');
lblMAD.Layout.Row=1;
lblMAD.Layout.Column = 7;
applyColorScheme(lblMAD, figColor);

uifMAD = uieditfield(plotButtonGrid,'numeric','Value',1,...
    'ValueChangedFcn',@(ui,ev)setMADVal(ui.Value,0));
uifMAD.Layout.Row=2;
uifMAD.Layout.Column=7;
applyColorScheme(uifMAD, figColor);

sliderMAD = uislider(plotButtonGrid,'Value', 1, 'Limits',[0 20]);
sliderMAD.MajorTicks = [0 20];
sliderMAD.Layout.Row = [1 2];
sliderMAD.Layout.Column = 8;
applyColorScheme(sliderMAD, figColor);
sliderMAD.ValueChangedFcn = @(sld,ev)setMADVal(sld.Value,1);
thresh_factor = 1;

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
    'Items',{'0.001','0.01','0.1','1','2','5','10','100'},...
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


%% Defualt run
plotChannels();


%% functions
    function refreshLabels
       for k = 1:numChannels
            channelCheckboxes(k).Value = channel_inclusion_plot(k);
            if channelCheckboxes(k).Value
                channelCheckboxes(k).Text = 'Include';
                channelCheckboxes(k).FontColor = 1-figColor;
                lblChannelIDs(k).FontColor = 1-figColor;
                lblChannelHandles(k).FontColor = 1-figColor;
                drawnow;
            else
                channelCheckboxes(k).Text = 'Exclude';
                channelCheckboxes(k).FontColor = [.8 .1 .1];
                lblChannelIDs(k).FontColor = [.8 .1 .1];
                lblChannelHandles(k).FontColor = [.8 .1 .1];
                drawnow;
            end
       end
       
       if any(~channel_inclusion_plot)
           chkAllChannel.Value = false;
       end
    end

    function onAllChanChecked(val)
        for k = 1:cfg.numChannels
            channelCheckboxes(k).Value = val;
            if val
                channelCheckboxes(k).Text = 'Include';
                channelCheckboxes(k).FontColor = 1-figColor;
                lblChannelIDs(k).FontColor = 1-figColor;
                lblChannelHandles(k).FontColor = 1-figColor;
                drawnow;
            else
                channelCheckboxes(k).Text = 'Exclude';
                channelCheckboxes(k).FontColor = [.8 .1 .1];
                lblChannelIDs(k).FontColor = [.8 .1 .1];
                lblChannelHandles(k).FontColor = [.8 .1 .1];
                drawnow;
            end
            channel_inclusion_plot(k) = val;
        end
        plotChannels();
    end


    function onChanCheckboxChanged(idx)
        channel_inclusion_plot(idx) = channelCheckboxes(idx).Value;
        if channelCheckboxes(idx).Value
            channelCheckboxes(idx).Text = 'Include';
            channelCheckboxes(idx).FontColor = 1-figColor;
             lblChannelIDs(idx).FontColor = 1-figColor;
             lblChannelHandles(idx).FontColor = 1-figColor;
             drawnow;
        else
            channelCheckboxes(idx).Text = 'Exclude';
            channelCheckboxes(idx).FontColor = [.8 .1 .1];
            lblChannelIDs(idx).FontColor = [.8 .1 .1];
            lblChannelHandles(idx).FontColor = [.8 .1 .1];
            drawnow;
        end
        plotChannels();
    end

    function onSave()
        channel_inclusion = channel_inclusion_plot;
        setappdata(parentFig, 'channel_inclusion', channel_inclusion);
        setappdata(parentFig, 'thresh_MAD', thresh_factor);
        uialert(parentFig,'Channel inclusion updated','Success','Icon','success');
    end

    function updatePlotLineSat(sldVal)
        plotLineSat = sldVal;
        disp(['Plot line saturation set to ', num2str(plotLineSat)]);
        plotChannels();
    end

    function updatePlotLineHue(sldVal)
        plotLineHue = sldVal;
        disp(['Plot line hue set to ', num2str(plotLineHue)]);
        plotChannels();
    end

    function updatePlotLineValue(sldVal)
        plotLineValue = sldVal;
        disp(['Plot line hue set to ', num2str(plotLineValue)]);
        plotChannels();
    end

    function updateScale(sVal)
        ampScale = sVal;
        plotChannels();
    end

    function updateLineWidth(sVal)
        if sVal == 0
            sVal = eps;
        end
        lineWidth = sVal;
        plotChannels();
    end

    
    function setMADVal(val,type)
        if type == 0
            sliderMAD.Value = val;
            drawnow;
        else
            uifMAD.Value = val;
            drawnow;
        end
        thresh_factor = val;
        plotChannels()
    end

    function updateXWindow(sVal, stepStr)
        if strcmp(stepStr,'full')
            xStepVal = 100;
        else
            xStepVal = str2double(stepStr);
        end
        xWindowStart = round(sVal);
        xWindowEnd   = xWindowStart + xStepVal;
        if xWindowEnd > xSlider.Limits(2)
            xWindowEnd   = xSlider.Limits(2);
            xWindowStart = xWindowEnd - xStepVal;
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
        plotChannels();
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
        plotChannels();
    end




    function plotChannels()
        disp('--- Plot filtered data function called ---');
        if isempty(mappedData)
            if isfield(cfg, 'inputFolder') && isfield(cfg, 'outputFolder')
                mappedData = map_input_file(cfg.fullFilePath, cfg);
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


        cla(axChannels,'reset');

        if any(isnan(thresh_MAD))
            fs = cfg.samplingFrequency;
            maxSampleSpan       = cfg.maxSampleSpan * 60;
            total_duration      = size(mappedData.Data.data,2) / fs;
            sampleSpan          = min(maxSampleSpan, total_duration);
            qualityCheckLength  = 5 * fs;
            qualityCheck_idx =  round(sampleSpan * fs/2);            
            qualityCheckSignal = double(mappedData.Data.data(channel_mapping, qualityCheck_idx:qualityCheck_idx + qualityCheckLength-1));
            if cfg.denoising
                qualityCheckSignal = removeSharedNoise_percentile(qualityCheckSignal, cfg.noisePrctile, cfg.noiseCorr);
            end
            if isfield(cfg, 'commonRef')
                switch cfg.commonRef
                    case 'median'
                        qualityCheckSignal = qualityCheckSignal - median(qualityCheckSignal,'omitmissing');
                    case 'mean'
                        qualityCheckSignal = qualityCheckSignal - mean(qualityCheckSignal,'omitmissing');
                    case 'none'
                end
            end
            % [qualityOut] = kiaSort_signal_quality_check(qualityCheckSignal, cfg);
            [qualityOut] = kiaSort_signal_quality_check(qualityCheckSignal, cfg);
            scale_factor = qualityOut.scale_factor;
            thresh_MAD   = qualityOut.thresh_MAD;
            channel_inclusion_plot = scale_factor' < cfg.badChannel_factor & channel_inclusion;
            refreshLabels;
        end


        [sy, sx] = size(mappedData.Data.data);

        xDataRange = round(xWindowStart * cfg.samplingFrequency)+1 : min(round(xWindowEnd * cfg.samplingFrequency), sx);
        yDataRange = max(1, yWindowStart) : min(yWindowEnd, sy);
        yDataRange = yDataRange(channel_inclusion_plot(yDataRange));
        
        if ~isempty(yDataRange)
            inputData = double(mappedData.Data.data(:, xDataRange));

            if strcmp(plotFilterType,'filtered')
                inputData = bandpass_filter_GUI(inputData', cfg.bandpass, cfg.samplingFrequency)';
            end

            if cfg.denoising
                inputData = removeSharedNoise_percentile(inputData, cfg.noisePrctile, cfg.noiseCorr);
            end
            inputData = inputData(channel_mapping(yDataRange),:);

            normIn = inputData./max(abs(inputData),[],'all');
            thresh_MAD_in = thresh_MAD(yDataRange)./max(abs(inputData),[],'all');
            size(normIn)
            size((yDataRange' + thresh_MAD_in * ampScale)')
            plot(axChannels, xDataRange/cfg.samplingFrequency, (yDataRange' + normIn * ampScale)', 'LineWidth', lineWidth/2);                        
            axChannels = plot_grouped_lines(numel(yDataRange), 0, axChannels, plotLineHue, plotLineSat, plotLineValue);
            hold(axChannels,"on")
            plot(axChannels, [xDataRange(1)/cfg.samplingFrequency xDataRange(end)/cfg.samplingFrequency], [(yDataRange' - thresh_factor * thresh_MAD_in * ampScale)' ;(yDataRange' - thresh_factor * thresh_MAD_in * ampScale)'], 'LineWidth', lineWidth,'color',[.8 .5 .1]);
            plot(axChannels, [xDataRange(1)/cfg.samplingFrequency xDataRange(end)/cfg.samplingFrequency], [(yDataRange' + thresh_factor * thresh_MAD_in * ampScale)' ;(yDataRange' + thresh_factor * thresh_MAD_in * ampScale)'], 'LineWidth', lineWidth,'color',[.1 .5 .8]);
            hold(axChannels,"off")
            disp(['Time range: ', num2str(xWindowStart), ' to ', num2str(xWindowEnd)]);
            disp(['Channel range: ', num2str(yDataRange(1)), ' to ', num2str(yDataRange(end))]);
        end

            title(axChannels,'Channel Explorer','Color','white');
            xlabel(axChannels, 'Time (s)');
            ylabel(axChannels, 'Channel #');
            set(axChannels,'tickdir','out');
            axChannels.Color   = figColor;
            axChannels.XColor  = 1 - figColor;
            axChannels.YColor  = 1 - figColor;
            box(axChannels,'off')
            axis(axChannels,[xWindowStart xWindowEnd yWindowStart-.5 yWindowEnd+.5]);
            yTicksIntervals = round((yWindowEnd-yWindowStart)/10);
            yTicksIntervals = max(yTicksIntervals,1);
            set(axChannels,'ytick',[yWindowStart:yTicksIntervals:yWindowEnd])       
    end

    function updatePlotFilterType(newVal)
        plotFilterType = newVal;
        disp(['Plot filter type set to ', plotFilterType]);
        plotChannels();
    end


end

