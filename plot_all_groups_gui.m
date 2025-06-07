function plot_all_groups_gui(cfg, parentPanel, figColor, parentFig)
    

    sampleRes = load_sorted_samples_gui(cfg.outputFolder);
    groupList   = sampleRes.crossChannelStats.unified_labels.label;  
    channelList = sampleRes.crossChannelStats.unified_labels.channelID;
    channelPlot = channelList + [-cfg.num_channel_extract : cfg.num_channel_extract];
    waveforms   = sampleRes.crossChannelStats.unified_labels.meanWaveforms;
    numGroups   = length(groupList);
    halfSpikeWaveDur = cfg.spikeDuration/2;
    waveformXaxis = -cfg.spikeDuration/2 : 1000/(cfg.samplingFrequency) : cfg.spikeDuration/2;
    xlocs = sampleRes.channel_info.channel_locations(:,1);
    ylocs = sampleRes.channel_info.channel_locations(:,2);    
    xNorm_adj = []; 
    yNorm_adj = [];
    normalizeTimeAmp();

    numChannels = min(cfg.numChannels, length(sampleRes.channel_info.channel_mapping));

    chanLabel = cell(1,numChannels);
    for i = 1 : numChannels
        chanLabel{i,1} = sprintf('Ch %d',i);
    end

    originalGroupList = groupList;        % For reset

    colorMapAll = generateDistinctColors(numGroups);

    mainLayout = uigridlayout(parentPanel,[1 2], ...
        'ColumnWidth',{400,'1x'}, ...
        'RowHeight',{'1x'}, ...
        'Padding',10, ...
        'ColumnSpacing',10);
    applyColorScheme(mainLayout, figColor);

    % left panels
    leftPanel = uipanel(mainLayout, ...
        'Title','Group Control', ...
        'BackgroundColor',figColor, ...
        'Scrollable','on');
    leftPanel.Layout.Column = 1;
    leftPanel.Layout.Row    = 1;
    applyColorScheme(leftPanel, figColor);


    leftGrid = uigridlayout(leftPanel, [2 1], ...
        'RowHeight',{'fit','1x'}, ...
        'ColumnWidth',{'1x'}, ...
        'Padding',10, ...
        'RowSpacing',10);
    applyColorScheme(leftGrid, figColor);

    % top left panel and its buttons
    topPanel = uipanel(leftGrid, 'BorderType','none');
    topPanel.Layout.Row = 1;
    topPanel.Layout.Column = 1;
    applyColorScheme(topPanel, figColor);

    topButtonLayout = uigridlayout(topPanel,[1 4], ...
        'ColumnWidth',{'fit','fit','fit','fit'}, ...
        'ColumnSpacing',10, ...
        'Padding',[0 0 0 0]);
    applyColorScheme(topButtonLayout, figColor);

    btnRemove = uibutton(topButtonLayout,'Text','Remove',...
        'ButtonPushedFcn',@(btn,ev)onRemove());
    btnRemove.Layout.Column = 1;
    applyColorScheme(btnRemove, figColor);

    btnMerge = uibutton(topButtonLayout,'Text','Merge',...
        'ButtonPushedFcn',@(btn,ev)onMerge());
    btnMerge.Layout.Column = 2;
    applyColorScheme(btnMerge, figColor);

    btnReset = uibutton(topButtonLayout,'Text','Reset',...
        'ButtonPushedFcn',@(btn,ev)onReset());
    btnReset.Layout.Column = 3;
    applyColorScheme(btnReset, figColor);


    %% Groups check boxes and radio panels
    checkPanel = uipanel(leftGrid, ...
        'Title','Unified Groups', ...
        'BackgroundColor',figColor,...
        'Scrollable','on');
    checkPanel.Layout.Row = 2;
    checkPanel.Layout.Column = 1;
    applyColorScheme(checkPanel, figColor);

    totalRows = numGroups + 1;  % +1 for the "All" row

    
    mainCheckGrid = uigridlayout(checkPanel, [totalRows, 5], ...
        'RowHeight', repmat({'fit'},1,totalRows), ...
        'ColumnWidth', {'fit','fit','fit','fit','fit'}, ...
        'RowSpacing',2, ...
        'Scrollable','on',...
        'ColumnSpacing',5, ...
        'Padding',[5 5 5 5]);
    applyColorScheme(mainCheckGrid, figColor);

    % qrrays for UI handles
    lblGroupHandles       = gobjects(numGroups,1);
    lblChannelHandles     = gobjects(numGroups,1);
    groupVisCheckboxes    = gobjects(numGroups,1);
    groupSelectCheckboxes = gobjects(numGroups,1);
    groupRadioButtons     = gobjects(numGroups,1);

    radioGroupBG = uibuttongroup(mainCheckGrid, ...
        'Title','', ...
        'Scrollable','on');
    radioGroupBG.Layout.Row = [2 totalRows];
    radioGroupBG.Layout.Column = 5;
    applyColorScheme(radioGroupBG, figColor);
    rowPixelHeight = 20;
    rowOffsetTop   = 5;   

    lblMerge = uilabel(mainCheckGrid,'Text','Merge To:',...
        'FontWeight','bold');
    lblMerge.Layout.Row = 1;
    lblMerge.Layout.Column = 5;
    applyColorScheme(lblMerge, figColor);

    % all selection row
    lblAll = uilabel(mainCheckGrid,'Text','T #:',...
        'FontWeight','bold');
    lblAll.Layout.Row = 1;
    lblAll.Layout.Column = 1;
    applyColorScheme(lblAll, figColor);

    lblDummy = uilabel(mainCheckGrid,'Text','Channel #:');
    lblDummy.Layout.Row = 1;
    lblDummy.Layout.Column = 2;
    applyColorScheme(lblDummy, figColor);

    chkAllVis = uicheckbox(mainCheckGrid,'Text','Vis. All',...
        'FontWeight','bold',...
        'Value',true,...
        'ValueChangedFcn',@(cb,ev)onAllVisChecked(cb.Value));
    chkAllVis.Layout.Row = 1;
    chkAllVis.Layout.Column = 3;
    applyColorScheme(chkAllVis, figColor);

    chkAllSel = uicheckbox(mainCheckGrid,'Text','Select All',...
        'FontWeight','bold',...
        'ValueChangedFcn',@(cb,ev)onAllSelectChecked(cb.Value));
    chkAllSel.Layout.Row = 1;
    chkAllSel.Layout.Column = 4;
    applyColorScheme(chkAllSel, figColor);
    
    yRow1 = (totalRows-1)*rowPixelHeight + rowOffsetTop;  % totalRows-1 from the bottom
    rdoNone = uiradiobutton(radioGroupBG, 'Text','none',...
        'FontWeight','bold','fontColor',1-figColor,...
        'Position',[10 yRow1 80 22]);
    rdoNone.Value = true; % default selection

    %% group rows
    for i = 1:numGroups
        rowIdx = i + 1; 
        
        % group list
        lblGroupHandles(i) = uilabel(mainCheckGrid, ...
            'Text', sprintf('G: %d', groupList(i)),...
            'FontWeight','bold');
        lblGroupHandles(i).Layout.Row = rowIdx;
        lblGroupHandles(i).Layout.Column = 1;

        % channel list
        lblChannelHandles(i) = uilabel(mainCheckGrid, ...
            'Text', sprintf('Channel %d', channelList(i)),...
            'FontWeight','bold');
        lblChannelHandles(i).Layout.Row = rowIdx;
        lblChannelHandles(i).Layout.Column = 2;

        % visualization checkboxes
        groupVisCheckboxes(i) = uicheckbox(mainCheckGrid,...
            'Value',true, 'Text','Vis.',...
            'ValueChangedFcn',@(cb,ev)onVisCheckboxChanged(i));
        groupVisCheckboxes(i).Layout.Row = rowIdx;
        groupVisCheckboxes(i).Layout.Column = 3;

        % selection checkboxes
        groupSelectCheckboxes(i) = uicheckbox(mainCheckGrid,...
            'Text','Select','Value',false);
        groupSelectCheckboxes(i).Layout.Row = rowIdx;
        groupSelectCheckboxes(i).Layout.Column = 4;

        % radio buttons
        yPos = (totalRows - rowIdx)*rowPixelHeight + rowOffsetTop;
        groupRadioButtons(i) = uiradiobutton(radioGroupBG,...
            'Text', sprintf('G%d',i),...
            'FontWeight','bold',...
            'Position',[10, yPos, 80, 22]);
    end

    %% Plot panel

    rightPanel = uipanel(mainLayout, ...
        'Title','Group Plot', ...
        'BackgroundColor',figColor);
    rightPanel.Layout.Column = 2;
    rightPanel.Layout.Row    = 1;
    applyColorScheme(rightPanel, figColor);

    rightGrid = uigridlayout(rightPanel, [10 3], ...
        'RowHeight',{'fit','fit','fit','fit','fit','fit','fit','fit','1x','10x'},...
        'ColumnWidth',{'1x', '10x', '1x'},...
        'Padding',[2 2 2 2],...
        'RowSpacing',5);
    applyColorScheme(rightGrid, figColor);
    
    
    channelValues = 2.^(0:10);
    tickIntIdx = nearest(channelValues, numChannels/10);
    ySlider = uislider(rightGrid,'Orientation','vertical','Value', 1,'MajorTicksMode','auto', 'Limits',[1 numChannels]);        
    ySlider.MajorTicks = [0 : channelValues(tickIntIdx): numChannels];
    ySlider.Layout.Row = [9 10];
    ySlider.Layout.Column = 1;
    applyColorScheme(ySlider, figColor);


    channelValues(channelValues > numChannels) = [];
    channelValues = arrayfun(@num2str, channelValues, 'UniformOutput', false);
    yStepDropdown = uidropdown(rightGrid,...
        'Items',[{'full'}, channelValues],...
        'Value','full');
    yStepDropdown.Layout.Row = 8;
    yStepDropdown.Layout.Column = 1;
    applyColorScheme(yStepDropdown, figColor);

    lblXStepSize = uilabel(rightGrid,'Text','# Channels:',...
        'HorizontalAlignment','right');
    lblXStepSize.Layout.Row = 7;
    lblXStepSize.Layout.Column = 1;
    applyColorScheme(lblXStepSize, figColor);

    ySlider.ValueChangedFcn      = @(sld,ev)updateYWindow(sld.Value, yStepDropdown.Value);
    yStepDropdown.ValueChangedFcn= @(dd,ev)updateYWindow(ySlider.Value, dd.Value);


    lblAmpSclae = uilabel(rightGrid,'Text','Amp. Scale:',...
        'HorizontalAlignment','right');
    lblAmpSclae.Layout.Row = 1;
    lblAmpSclae.Layout.Column = 1;
    applyColorScheme(lblAmpSclae, figColor);
    
    ampSlider = uislider(rightGrid,'Orientation','horizontal','Value', 1,'MajorTicksMode','auto', 'Limits',[1 10]);        
    ampSlider.MajorTicks = [1 10];
    ampSlider.Layout.Row = 2;
    ampSlider.Layout.Column = 1;
    applyColorScheme(ampSlider, figColor);
    ampSlider.ValueChangedFcn      = @(sld,ev)updateScale(sld.Value);

    lblLineWidth = uilabel(rightGrid,'Text','Line width:',...
        'HorizontalAlignment','right');
    lblLineWidth.Layout.Row = 3;
    lblLineWidth.Layout.Column = 1;
    applyColorScheme(lblLineWidth, figColor);

    lineSlider = uislider(rightGrid,'Orientation','horizontal','Value', 1,'MajorTicksMode','auto', 'Limits',[0 5]);        
    lineSlider.MajorTicks = [0 5];
    lineSlider.Layout.Row = 4;
    lineSlider.Layout.Column = 1;
    applyColorScheme(lineSlider, figColor);
    lineSlider.ValueChangedFcn      = @(sld,ev)updateLineWidth(sld.Value);



  
    lblchType = uilabel(rightGrid,'Text','Channel Show:',...
        'HorizontalAlignment','right');
    lblchType.Layout.Row = 5;
    lblchType.Layout.Column = 1;
    applyColorScheme(lblchType, figColor);

    chTypeDropdown = uidropdown(rightGrid,...
        'Items',{'All' , 'Main'},...
        'Value','All');
    chTypeDropdown.Layout.Row = 6;
    chTypeDropdown.Layout.Column = 1;
    applyColorScheme(chTypeDropdown, figColor);
    chTypeDropdown.ValueChangedFcn      = @(sld,ev)updateChType(sld.Value);

    AllChannelPlot = true;






    % Initialize X/Y window
    yWindowStart = 1;  
    yWindowEnd   = numChannels; 
    yStepVal     = numChannels;
    ampScale     = 1;
    lineWidth    = 1;

    axChannels = uiaxes(rightGrid);
    axChannels.Layout.Row = [1 10];
    axChannels.Layout.Column = [2 3];
    axChannels.Color   = [0.2 0.2 0.2];
    axChannels.XColor  = [1 1 1];
    axChannels.YColor  = [1 1 1];
    title(axChannels,'Channels','Color','white');
    xlabel(axChannels,'Time','Color','white');
    ylabel(axChannels,'Amplitude','Color','white');

    %% Add scroller and channel number to visualize

    %% =============== INITIAL REFRESH ===============
    refreshLabels();
    updatePlot();

    %% =============== CALLBACKS ===============

    function onAllVisChecked(val)
        for k = 1:numGroups
            groupVisCheckboxes(k).Value = val;
        end
        updatePlot();
    end

    function onAllSelectChecked(val)
        for k = 1:numGroups
            groupSelectCheckboxes(k).Value = val;
        end
    end

    function onVisCheckboxChanged(idx)
        updatePlot();
    end

    function onRemove()
        % For each group that is selected => set label = NaN => color => gray
        for k = 1:numGroups
            if groupSelectCheckboxes(k).Value
                groupList(k) = NaN;
            end
        end
        refreshLabels();
        updatePlot();
    end

    function onMerge()
        % Find which radio is selected
        selIdx = find([groupRadioButtons.Value] == 1, 1);
        if isempty(selIdx)
            % means "none" or no radio
            uialert(parentFig,'No valid merge group selected.','Merge Warning');
            return;
        end
        % The user must also have that group selected for merging => else "Inconsistent"
        if ~groupSelectCheckboxes(selIdx).Value
            uialert(parentFig,'Inconsistent merging groups. Merging should be to one selected sample.','Merge Error');
            return;
        end

        mergeLabel = groupList(selIdx);
        if isnan(mergeLabel)
            uialert(parentFig,'Cannot merge into a removed (NaN) sample.','Merge Error');
            return;
        end

        % Perform the merge => all selected become the same label
        for k = 1:numGroups
            if groupSelectCheckboxes(k).Value
                groupList(k) = mergeLabel;
            end
        end
        refreshLabels();
        updatePlot();
    end

    function onReset()
        groupList = originalGroupList;

        % Uncheck all selection
        for k = 1:numGroups
            groupSelectCheckboxes(k).Value = false;
        end
        chkAllSel.Value = false;

        % Set all vis = true
        for k = 1:numGroups
            groupVisCheckboxes(k).Value = true;
        end
        chkAllVis.Value = false;

        % Radio => none
        rdoNone.Value = true;

        refreshLabels();
        updatePlot();
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
        updatePlot();
    end

    function updateScale(sVal)
        ampScale = sVal;
        updatePlot();
    end

    function updateLineWidth(sVal)        
        lineWidth = sVal + eps;
        updatePlot();
    end

    function updateChType(sVal)

        if strcmp(sVal,'All') 
            AllChannelPlot = true;
        else
            AllChannelPlot = false;
        end
        updatePlot();
    end


    %% =============== HELPER FUNCTIONS ===============
    function refreshLabels()
        % Update group labels & colors
        for k = 1:numGroups
            lblVal = groupList(k);
            if isnan(lblVal)
                txt = 'G: NaN';
            else
                txt = sprintf('G: %g', lblVal);
            end
            lblGroupHandles(k).Text = txt;

            c = getGroupColor(lblVal);
            % We'll color the label's background if supported
            if isprop(lblGroupHandles(k), 'BackgroundColor')
                lblGroupHandles(k).BackgroundColor = c;
            end
            if isprop(lblGroupHandles(k), 'FontColor')
                lblGroupHandles(k).FontColor = bestTextColorFor(c);
            end

            % Channel label
            if isprop(lblChannelHandles(k), 'BackgroundColor')
                lblChannelHandles(k).BackgroundColor = c;
            end
            if isprop(lblChannelHandles(k), 'FontColor')
                lblChannelHandles(k).FontColor = bestTextColorFor(c);
            end

            % Visualization checkbox (often no BackgroundColor property, so skip it)
            if isprop(groupVisCheckboxes(k), 'FontColor')
                groupVisCheckboxes(k).FontColor = 1-figColor;
            end

            % Selection checkbox
            if isprop(groupSelectCheckboxes(k), 'FontColor')
                groupSelectCheckboxes(k).FontColor = 1-figColor;
            end

            % Radio button
            if isprop(groupRadioButtons(k), 'FontColor')
                groupRadioButtons(k).FontColor = c;
            end
        end
    end


    function updatePlot()
        cla(axChannels,'reset');
        hold(axChannels,'on');
        axChannels.Color  = [0.2 0.2 0.2];
        axChannels.XColor = [1 1 1];
        axChannels.YColor = [1 1 1];
        chanPlotList = [];
        for k = 1:numGroups
            if groupVisCheckboxes(k).Value                
                middleChannel = 1 + cfg.num_channel_extract;
                middleChannelIdx = channelList(k);
                chPlot = channelPlot(k,:);                
                validCh = (chPlot >= 1 & chPlot <= numChannels & ...
                    chPlot >= yWindowStart & chPlot < yWindowEnd);
                chPlot = chPlot(validCh);               
                c = getGroupColor(groupList(k));
                waveformsPlot = squeeze(waveforms(k, validCh, :));
                middleWaveform = squeeze(waveforms(k, middleChannel, :));
                waveformsPlot(mean(waveformsPlot, 2) == 0, :) = nan;
                waveformsPlot = ampScale * waveformsPlot./max(abs(waveforms),[],'all');
                middleWaveform= ampScale * middleWaveform./max(abs(waveforms),[],'all');
                if AllChannelPlot
                plot(axChannels, (waveformXaxis + xNorm_adj(chPlot))',(waveformsPlot+yNorm_adj(chPlot))', 'Color', c, 'LineWidth', lineWidth)
                end
                plot(axChannels, (waveformXaxis + xNorm_adj(middleChannelIdx))',(middleWaveform+yNorm_adj(middleChannelIdx))',...
                    'Color', c, 'LineWidth', 1.5 * lineWidth)           
                chanPlotList = [chanPlotList, chPlot];
            end
        end
        chanPlotList = unique(chanPlotList);        
        text(axChannels, xNorm_adj(chanPlotList)+waveformXaxis(1), yNorm_adj(chanPlotList), chanLabel(chanPlotList),...
            'color',1-figColor,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        yTicksIntervals = round((yWindowEnd-yWindowStart)/10);
        yTicksIntervals = max(yTicksIntervals,1);
        hold(axChannels,'off');
        title(axChannels,'Channels','Color','white');
        axis(axChannels,[min(xNorm_adj)-halfSpikeWaveDur max(xNorm_adj)+halfSpikeWaveDur  min(yNorm_adj(yWindowStart))-1 max(yNorm_adj(yWindowEnd))+1 ]);
        box(axChannels,'off')
        axis(axChannels,'off')
       
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


end



function colorMap = generateDistinctColors(N)
    % GENERATEDISTINCTCOLORS(N)
    % Returns Nx3 "vibrant" colors, spaced out using a golden-ratio approach.
    % For large N, consider advanced methods (e.g., "distinguishable_colors" on FEX).

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
            % h = mod(h + goldenRatio, 1);
            rgb_color = hsv2rgb([h, s ,v]);
            brightness = sum(rgb_color .* [0.299, 0.587, 0.114]);
        end
        colorMap(i,:) = rgb_color; 
    end
end



function tc = bestTextColorFor(bg)
    
    brightness = sum(bg .* [0.299, 0.587, 0.114]); % weighted luminosity
    if brightness < 0.5
        tc = [1 1 1]; % white text
    else
        tc = [0 0 0]; % black text
    end
end