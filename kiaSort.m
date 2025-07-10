function kiaSort

    %%  layout structure
    cfg = kiaSort_main_configs();              % user config
    hp  = sorting_hyperparameters_in();     % hyperparameter struct

    channel_mapping = [];
    channel_inclusion = [];
    channel_locations = [];
    mappedData        = [];

    
    choiceMap = struct();
    choiceMap.dataType             = {'int16','int32','float32','unit16'};
    choiceMap.method               = {'direct','indirect'};
    choiceMap.commonRef            = {'median','mean','none'};
    choiceMap.modelType            = {'mlp','cnn','svm'};

    hpChoiceMap = struct();
    hpChoiceMap.mlpLossFcn           = {'mse','mae'};
    hpChoiceMap.cnnLossFcn           = {'dlconv','crossentropy'};
    hpChoiceMap.executionEnvironment = {'auto','gpu','cpu'};
    hpChoiceMap.KernelFunction       = {'polynomial','rbf','linear','auto'};
    hpChoiceMap.KernelScale          = {'auto','manual'};  

    if ~isfield(cfg,'numChannels')
        cfg.numChannels = 16; 
    end

    executionControl = 'running';

    %% Main figure
    figColor = [0.15 0.15 0.15]; 
    tabColors = [.8 .5 .5];
    gridColors = [0.05 0.05 0.05];
    fig = uifigure('Name','KIASort: Knowledge-Integrated Automated Spike Sorting',...
                   'Color',figColor,...
                   'Position',[.9 .9 1 .9].*get(0,'ScreenSize'),'Scrollable','on');

    mainTabGroup = uitabgroup(fig,'Units','normalized','Position',[0 0 1 1]);
    applyColorScheme(mainTabGroup, figColor);

    % 6 primary tabs
    tab1 = uitab(mainTabGroup, 'Title','Main');
    applyColorScheme(tab1, tabColors);
    tab2 = uitab(mainTabGroup, 'Title','Channels');
    applyColorScheme(tab2, tabColors);
    tab3 = uitab(mainTabGroup, 'Title','Curation');
    applyColorScheme(tab3, tabColors);
    tab4 = uitab(mainTabGroup, 'Title','Clusters');
    applyColorScheme(tab4, tabColors);
    tab5 = uitab(mainTabGroup, 'Title','Help');
    applyColorScheme(tab5, tabColors);
    tab6 = uitab(mainTabGroup, 'Title','About');
    applyColorScheme(tab6, tabColors);

    %% Tab 1
    subTabGroup1 = uitabgroup(tab1,'Units','normalized','Position',[0 0 1 1]);
    applyColorScheme(subTabGroup1, figColor);

    configSubTab = uitab(subTabGroup1,'Title','Run');
    applyColorScheme(configSubTab, tabColors);

    configExtTab  = uitab(subTabGroup1,'Title','Extended Config');
    applyColorScheme(configExtTab, tabColors);

    hyperSubTab  = uitab(subTabGroup1,'Title','Hyperparameters');
    applyColorScheme(hyperSubTab, tabColors);

    

    %% Tab 1: Sub-Tab (a) Run

    gridTab1_L = uigridlayout(configSubTab,[4 2],...
        'ColumnWidth',{'1.5x','3.5x'},...
        'RowHeight',{'1.45x','1x','1.5x','9.5x'},...
        'Padding',10,...
        'ColumnSpacing',10,'Scrollable','on');
    applyColorScheme(gridTab1_L, gridColors);

        % Left panel for config fields & folder selection

    topLeftPanel = uipanel(gridTab1_L,...
        'Title','Run KIA Sort',...
        'Scrollable','on',...
        'BackgroundColor',figColor);
    topLeftPanel.Layout.Column = 1;
    topLeftPanel.Layout.Row    = [1 3];
    applyColorScheme(topLeftPanel, figColor);

    bottomLeftPanel = uipanel(gridTab1_L,...
        'Title','Sorting Config',...
        'Scrollable','on',...
        'BackgroundColor',figColor);
    bottomLeftPanel.Layout.Column = 1;
    bottomLeftPanel.Layout.Row    = 4;
    applyColorScheme(bottomLeftPanel, figColor);

    topRightPanel = uipanel(gridTab1_L,...
        'Title','Visualization Control',...
        'BackgroundColor',figColor);
    topRightPanel.Layout.Column = 2;
    topRightPanel.Layout.Row    = 1;
    applyColorScheme(topRightPanel, figColor);

    bottomRightPanel = uipanel(gridTab1_L,...
        'Title','Plots',...
        'BackgroundColor',figColor);
    bottomRightPanel.Layout.Column = 2;
    bottomRightPanel.Layout.Row    = [2 4];
    applyColorScheme(bottomRightPanel, figColor);


    topLeftGrid = uigridlayout(topLeftPanel,...
        'RowHeight', {'1.5x','fit','fit','.75x','.75x','.75x'},...
        'ColumnWidth',{'fit','fit','fit','1x','1x','1x'},...
        'RowSpacing',10,...
        'ColumnSpacing',5,...
        'Padding',[10 5 10 5]);
    applyColorScheme(topLeftGrid, figColor);


    % Row1: push buttons

    if exist('Run_KIASORT.png','file')
        btnRun = uibutton(topLeftGrid, 'Icon','Run_KIASORT.png', 'Text','Run & Sort',...
            'ButtonPushedFcn',@(btn,ev)onRunAndSort());
    else
        btnRun = uibutton(topLeftGrid, 'Text','Run & Sort',...
            'ButtonPushedFcn',@(btn,ev)onRunAndSort());
    end
    btnRun.Layout.Row = 1;  
    btnRun.Layout.Column = [1 6];    
    applyColorScheme(btnRun, figColor);

    

    if exist('SortD.png','file')
    btnStop = uibutton(topLeftGrid, 'Text', '', 'Icon', 'stop.png',...
        'ButtonPushedFcn', @(btn,ev)onStop());
    else
    btnStop = uibutton(topLeftGrid, 'Text', 'Stop',...
        'ButtonPushedFcn', @(btn,ev)onStop());
    end
    btnStop.Layout.Row = 2;
    btnStop.Layout.Column = 1;
    applyColorScheme(btnStop, figColor);


    if exist('SortD.png','file')
    btnPause = uibutton(topLeftGrid, 'Text', '', 'Icon', 'pause.png',...
        'ButtonPushedFcn', @(btn,ev)onPause());
    else
    btnPause = uibutton(topLeftGrid, 'Text', 'Pause',...
        'ButtonPushedFcn', @(btn,ev)onPause());
    end
    btnPause.Layout.Row = 2;
    btnPause.Layout.Column = 2;
    applyColorScheme(btnPause, figColor);   

    if exist('SortD.png','file')
        btnResume = uibutton(topLeftGrid, 'Text', '', 'Icon', 'resume.png',...
            'ButtonPushedFcn', @(btn,ev)onResume());
    else
        btnResume = uibutton(topLeftGrid, 'Text', 'Resume',...
            'ButtonPushedFcn', @(btn,ev)onResume());
    end
    btnResume.Layout.Row = 2;
    btnResume.Layout.Column = 3;
    applyColorScheme(btnResume, figColor);



    btnExtract = uibutton(topLeftGrid, 'Text','Extract Samples',...
            'ButtonPushedFcn',@(btn,ev)onExtractSamples()); 
    btnExtract.Layout.Row = 2;  
    btnExtract.Layout.Column = 4;
    applyColorScheme(btnExtract, figColor);

    btnSortT = uibutton(topLeftGrid, 'Text','Sort Samples',...
            'ButtonPushedFcn',@(btn,ev)onSortGroups());
    btnSortT.Layout.Row = 2;  
    btnSortT.Layout.Column = 5;
    applyColorScheme(btnSortT, figColor);

    btnSortD = uibutton(topLeftGrid, 'Text','Sort Data',...
            'ButtonPushedFcn',@(btn,ev)onSortData());
    btnSortD.Layout.Row = 2;  
    btnSortD.Layout.Column = 6;
    applyColorScheme(btnSortD, figColor);
    
    

    % Row3: progress bar and label
    progressBar = uigauge(topLeftGrid,"linear",...
        'Limits', [0 100], ...
        'ScaleColors', figColor,...
        'ScaleColorLimits',[0 100],...
        'Value', 0, ...       
        'MajorTicks', 0:25:100 ...
        );
    progressBar.Layout.Row = 3;  
    progressBar.Layout.Column = [5 6];
    applyColorScheme(progressBar, figColor);

    progressLabel = uilabel(topLeftGrid,...
        'HorizontalAlignment','Right',...
        'Text', 'No Sorting in Progress');
    progressLabel.Layout.Row = 3;  
    progressLabel.Layout.Column = [1 4];
    applyColorScheme(progressLabel, figColor);

    % Row 4-6: folder/file selection
    rowFolder1 = 4;
    rowFolder2 = 5;
    rowMatFile = 6;

    btnFolder1 = uibutton(topLeftGrid,'Text','Input data file',...
        'ButtonPushedFcn',@(btn,ev)loadDataFile());
    btnFolder1.Layout.Row = rowFolder1;
    btnFolder1.Layout.Column = [1 3];
    applyColorScheme(btnFolder1, figColor);

    lblFolder1 = uilabel(topLeftGrid,'Text','Select the data input ".dat" or ".bin" file.',...
        'HorizontalAlignment','left');
    lblFolder1.Layout.Row = rowFolder1;
    lblFolder1.Layout.Column = [4 6];
    applyColorScheme(lblFolder1, figColor);

    btnFolder2 = uibutton(topLeftGrid,'Text','Output Path',...
        'ButtonPushedFcn',@(btn,ev)selectFolder(2));
    btnFolder2.Layout.Row = rowFolder2;
    btnFolder2.Layout.Column = [1 3];
    applyColorScheme(btnFolder2, figColor);

    lblFolder2 = uilabel(topLeftGrid,'Text','Select a folder to save the results.',...
        'HorizontalAlignment','left');
    lblFolder2.Layout.Row = rowFolder2;
    lblFolder2.Layout.Column = [4 6];
    applyColorScheme(lblFolder2, figColor);

    btnLoadMat = uibutton(topLeftGrid,'Text','Channel Map File',...
        'ButtonPushedFcn',@(btn,ev)loadMatFile());
    btnLoadMat.Layout.Row = rowMatFile;
    btnLoadMat.Layout.Column = [1 3];
    applyColorScheme(btnLoadMat, figColor);

    lblMatFile = uilabel(topLeftGrid,'Text','Load ".mat" channel map file.',...
        'HorizontalAlignment','left');
    lblMatFile.Layout.Row = rowMatFile;
    lblMatFile.Layout.Column = [4 6];
    applyColorScheme(lblMatFile, figColor);

    % bottom left panel

    cfgFields = fieldnames(cfg);
    totalRows = 10+ceil(numel(cfgFields)/2);

    bottomleftGrid = uigridlayout(bottomLeftPanel,...
        'RowHeight', repmat({'fit'},1,totalRows),...
        'ColumnWidth',{'1.25x','0.75x','1.25x','0.75x'},...
        'RowSpacing',10,...
        'ColumnSpacing',10,...
        'Scrollable','on',...
        'Padding',[10 10 10 10]);
    applyColorScheme(bottomleftGrid, figColor);
    
    rowIndex = 4;  
    colIndex = 1;
    for iF = 1:numel(cfgFields)-1

        if rowIndex == 4
            textLbl = sprintf('                _______________________________ Data Config __________________________ \n');
            lbl = uilabel(bottomleftGrid, 'Text',textLbl,...
                'HorizontalAlignment','Left');
            lbl.Layout.Row    = rowIndex;
            lbl.Layout.Column = [1 4];
            applyColorScheme(lbl, figColor);
            rowIndex=rowIndex+1;
        end

        fName = cfgFields{iF};
        fVal  = cfg.(fName);

        lbl = uilabel(bottomleftGrid, 'Text', cfg.varName.(cfgFields{iF}),...
            'HorizontalAlignment','right');
        lbl.Layout.Row    = rowIndex;
        lbl.Layout.Column = colIndex;
        applyColorScheme(lbl, figColor);

        if isfield(choiceMap, fName)
            ddItems = choiceMap.(fName);
            ddVal   = char(fVal);
            dd = uidropdown(bottomleftGrid,...
                'Items', ddItems,...
                'Value', ddVal,...
                'ValueChangedFcn',@(ui,ev)updateDropdown(fName, ui.Value));
            dd.Layout.Row    = rowIndex;
            dd.Layout.Column = colIndex+1;
            applyColorScheme(dd, figColor);

        elseif islogical(fVal)
            cb = uicheckbox(bottomleftGrid, 'Value',fVal, 'Text','',...
                'ValueChangedFcn',@(ui,ev)updateBoolean(fName, ui.Value));
            cb.Layout.Row    = rowIndex;
            cb.Layout.Column = colIndex+1;
            applyColorScheme(cb, figColor);

        elseif isnumeric(fVal) && numel(fVal)==2
            subP = uipanel(bottomleftGrid,'BorderType','none');
            subP.Layout.Row    = rowIndex;
            subP.Layout.Column = colIndex+1;
            applyColorScheme(subP, figColor);

            miniG = uigridlayout(subP,[1 2],'ColumnWidth',{'1x','1x'},...
                'Padding',[0 0 0 0]);
            applyColorScheme(miniG, figColor);

            e1 = uieditfield(miniG,'numeric','Value',fVal(1),...
                'ValueChangedFcn',@(ui,ev)setArrayField(fName,1,ui.Value));
            e1.Layout.Row=1; 
            e1.Layout.Column=1;
            applyColorScheme(e1, figColor);

            e2 = uieditfield(miniG,'numeric','Value',fVal(2),...
                'ValueChangedFcn',@(ui,ev)setArrayField(fName,2,ui.Value));
            e2.Layout.Row=1; 
            e2.Layout.Column=2;
            applyColorScheme(e2, figColor);

        elseif isnumeric(fVal)
            numF = uieditfield(bottomleftGrid,'numeric','Value',fVal,...
                'ValueChangedFcn',@(ui,ev)updateNumericField(fName, ui.Value));
            numF.Layout.Row = rowIndex;
            numF.Layout.Column = colIndex+1;
            applyColorScheme(numF, figColor);

        elseif isstring(fVal) || ischar(fVal)
            txt = uieditfield(bottomleftGrid,'text','Value',char(fVal),...
                'ValueChangedFcn',@(ui,ev)updateTextOrNumeric(fName, ui.Value));
            txt.Layout.Row = rowIndex;
            txt.Layout.Column = colIndex+1;
            applyColorScheme(txt, figColor);
        else
            txt = uieditfield(bottomleftGrid,'text','Value',mat2str(fVal),...
                'ValueChangedFcn',@(ui,ev)setField(fName, ui.Value));
            txt.Layout.Row = rowIndex;
            txt.Layout.Column = colIndex+1;
            applyColorScheme(txt, figColor);
        end

        if colIndex==1
            colIndex=3;
        else
            colIndex=1;
            rowIndex=rowIndex+1;
        end

        if colIndex==1 && rowIndex == 7
            
            textLbl = sprintf(' \n                _______________________________ Denoising __________________________\n');
            lbl = uilabel(bottomleftGrid, 'Text',textLbl,...
                'HorizontalAlignment','Left');
            lbl.Layout.Row    = rowIndex;
            lbl.Layout.Column = [1 4];
            applyColorScheme(lbl, figColor);

            colIndex=1;
            rowIndex=rowIndex+1;

        elseif colIndex==1 && rowIndex == 10
            % rowIndex=rowIndex+1;
            textLbl = sprintf(' \n                 _____________________________ Sorting Config _______________________\n');
            lbl = uilabel(bottomleftGrid, 'Text', textLbl,...
                'HorizontalAlignment','Left');
            lbl.Layout.Row    = rowIndex;
            lbl.Layout.Column = [1 4];
            applyColorScheme(lbl, figColor);

            colIndex=1;
            rowIndex=rowIndex+1;

        elseif colIndex==3 && rowIndex == 13
            rowIndex=rowIndex+1;
            textLbl = sprintf(' \n                _______________________________ Options __________________________\n');
            lbl = uilabel(bottomleftGrid, 'Text',textLbl,...
                'HorizontalAlignment','Left');
            lbl.Layout.Row    = rowIndex;
            lbl.Layout.Column = [1 4];
            applyColorScheme(lbl, figColor);

            colIndex=1;
            rowIndex=rowIndex+1;

        end

    end

    % Top right panel: plotting controller
    topRightGrid = uigridlayout(topRightPanel,...
        'RowHeight',{'fit'},...
        'ColumnWidth',{'1x'},...        
        'Padding',10);
    applyColorScheme(topRightGrid, figColor);

    bottomRightGrid = uigridlayout(bottomRightPanel,[3 3],...
        'RowHeight',{'1x','8x','1x'},...
        'ColumnWidth',{'1x','10x','1x'},...
        'Padding',10);
    applyColorScheme(bottomRightGrid, figColor);


    % Row1: Buttons for Plot raw data, Plot filtered data, and dropdown
    plotButtonPanel = uipanel(topRightGrid,'BorderType','none','BackgroundColor',figColor);
    plotButtonPanel.Layout.Row = 1;
    plotButtonPanel.Layout.Column = 1;
    applyColorScheme(plotButtonPanel, figColor);

    plotButtonGrid = uigridlayout(plotButtonPanel,[2 8],...
        'ColumnWidth',{'3x','1x', '1x', '1x', '1x', '1x', '1x', '1x','1x'},...
        'RowHeight',{'fit','fit'},...
        'Padding',[0 0 0 0],...
        'RowSpacing', 1,...
        'ColumnSpacing',20);
    applyColorScheme(plotButtonGrid, figColor);

    btnPlotData = uibutton(plotButtonGrid,'Text','Plot data',...
        'ButtonPushedFcn',@(btn,ev)onPlotData());
    btnPlotData.Layout.Row = [1 2];
    btnPlotData.Layout.Column = 1;
    applyColorScheme(btnPlotData, figColor);
    current_cfg = cfg;
    replotVal = false;

    plotTypeDD = uidropdown(plotButtonGrid,...
        'Items',{'images','lines'},...
        'Value','lines',...
        'ValueChangedFcn',@(dd,ev)updatePlotType(dd.Value));
    plotTypeDD.Layout.Row = 2;
    plotTypeDD.Layout.Column = 2;
    applyColorScheme(plotTypeDD, figColor);
    lastPlotted = false;
    currentPlotType = 'lines';

    plotFilterTypeDD = uidropdown(plotButtonGrid,...
        'Items',{'filtered','raw'},...
        'Value','filtered',...
        'ValueChangedFcn',@(dd,ev)updatePlotFilterType(dd.Value));
    plotFilterTypeDD.Layout.Row = 2;
    plotFilterTypeDD.Layout.Column = 3;
    applyColorScheme(plotFilterTypeDD, figColor);
    plotFilterType = 'filtered';
   
    numChannels = cfg.numChannels;
    groupValues = [0, 2.^(0:10)];
    groupValues(groupValues > cfg.numChannels) = [];
    groupValues = arrayfun(@num2str, groupValues, 'UniformOutput', false);
    plotGroupLineColor = uidropdown(plotButtonGrid,...
        'Items', groupValues,...
        'Value', '0',...
        'ValueChangedFcn',@(dd,ev)updatePlotLineColor(dd.Value));
    plotGroupLineColor.Layout.Row = 2;
    plotGroupLineColor.Layout.Column = 4;
    applyColorScheme(plotGroupLineColor, figColor);
    plotLineColor = '0';

    ampSlider = uislider(plotButtonGrid,'Value', 1, 'Limits',[0 10]);        
    ampSlider.MajorTicks = [0 10];
    ampSlider.Layout.Row = 2;
    ampSlider.Layout.Column = 5;
    applyColorScheme(ampSlider, figColor);
    ampSlider.ValueChangedFcn      = @(sld,ev)updatePlotLineScale(sld.Value);
    plotLineScale = 1;
   
    hueSlider = uislider(plotButtonGrid,'Value', .1, 'Limits',[0 1]);        
    hueSlider.MajorTicks = [0 1];
    hueSlider.Layout.Row = 2;
    hueSlider.Layout.Column = 6;
    applyColorScheme(hueSlider, figColor);
    hueSlider.ValueChangedFcn      = @(sld,ev)updatePlotLineHue(sld.Value);
    plotLineHue = hueSlider.Value;

    satSlider = uislider(plotButtonGrid,'Value', .1, 'Limits',[0 1]);        
    satSlider.MajorTicks = [0 1];
    satSlider.Layout.Row = 2;
    satSlider.Layout.Column = 7;
    applyColorScheme(satSlider, figColor);
    satSlider.ValueChangedFcn      = @(sld,ev)updatePlotLineSat(sld.Value);
    plotLineSat = satSlider.Value;

    valueSlider = uislider(plotButtonGrid,'Value', .8, 'Limits',[0 1]);        
    valueSlider.MajorTicks = [0 1];
    valueSlider.Layout.Row = 2;
    valueSlider.Layout.Column = 8;
    applyColorScheme(valueSlider, figColor);
    valueSlider.ValueChangedFcn      = @(sld,ev)updatePlotLineValue(sld.Value);
    plotLineValue = valueSlider.Value;

    btnExportFig = uibutton(plotButtonGrid,'Text','Export Fig.',...
        'ButtonPushedFcn',@(btn,ev)onExportFig());
    btnExportFig.Layout.Row = 2;
    btnExportFig.Layout.Column = 9;
    applyColorScheme(btnExportFig, figColor);



    % labels for buttons
    lblplotType = uilabel(plotButtonGrid,'Text','Plot Type:',...
        'HorizontalAlignment','right');
    lblplotType.Layout.Row = 1;
    lblplotType.Layout.Column = 2;
    applyColorScheme(lblplotType, figColor);

    lblFilterType = uilabel(plotButtonGrid,'Text','Data Type:',...
        'HorizontalAlignment','right');
    lblFilterType.Layout.Row = 1;
    lblFilterType.Layout.Column = 3;
    applyColorScheme(lblFilterType, figColor);

    lblLineGroup = uilabel(plotButtonGrid,'Text','# Line Grouping:',...
        'HorizontalAlignment','right');
    lblLineGroup.Layout.Row = 1;
    lblLineGroup.Layout.Column = 4;
    applyColorScheme(lblLineGroup, figColor);

    lblLineMag = uilabel(plotButtonGrid,'Text','Amp Scale:',...
        'HorizontalAlignment','right');
    lblLineMag.Layout.Row = 1;
    lblLineMag.Layout.Column = 5;
    applyColorScheme(lblLineMag, figColor);

    lblHue = uilabel(plotButtonGrid,'Text','Hue:',...
        'HorizontalAlignment','right');
    lblHue.Layout.Row = 1;
    lblHue.Layout.Column = 6;
    applyColorScheme(lblHue, figColor);

    lblsat = uilabel(plotButtonGrid,'Text','Saturation:',...
        'HorizontalAlignment','right');
    lblsat.Layout.Row = 1;
    lblsat.Layout.Column = 7;
    applyColorScheme(lblsat, figColor);

    lblValue = uilabel(plotButtonGrid,'Text','Value:',...
        'HorizontalAlignment','right');
    lblValue.Layout.Row = 1;
    lblValue.Layout.Column = 8;
    applyColorScheme(lblValue, figColor);

    exportFormatDD = uidropdown(plotButtonGrid,...
        'Items',{'Image','Vector'},...
        'Value','Image',...
        'ValueChangedFcn',@(dd,ev)updateExportFormat(dd.Value));
    exportFormatDD.Layout.Row = 1;
    exportFormatDD.Layout.Column = 9;
    applyColorScheme(exportFormatDD, figColor);
    exportType = 'Image';
    
    % bottom right panel: Time-axis slider
    xPanel = uipanel(bottomRightGrid,'BackgroundColor',figColor,'BorderColor',figColor);
    xPanel.Layout.Row = 3;
    xPanel.Layout.Column = [2 3];
    applyColorScheme(xPanel, figColor);

    xGrid = uigridlayout(xPanel,[3 4],...
        'ColumnWidth',{'1x','8x','1x','1x'},...
        'RowHeight',{'1x', '1x', '1x'},...
        'Padding',[5 5 5 5],...
        'ColumnSpacing',10);
    applyColorScheme(xGrid, figColor);

    trialLength = 100;
    xSlider = uislider(xGrid,'Value',0,...
        'MajorTicksMode','auto');
    xSlider.Layout.Row = [1 3];
    xSlider.Layout.Column = [1 3];    
    applyColorScheme(xSlider, figColor);

    lblXStepSize = uilabel(xGrid,'Text','Time window (s):',...
        'HorizontalAlignment','right');
    lblXStepSize.Layout.Row = 1;
    lblXStepSize.Layout.Column = 4;
    applyColorScheme(lblXStepSize, figColor);

    xStepDropdown = uidropdown(xGrid,...
        'Items',{'0.001','0.01','0.1','1','2','5','10','100'},...
        'Value','1');
    xStepDropdown.Layout.Row = [2 3];
    xStepDropdown.Layout.Column = 4;
    applyColorScheme(xStepDropdown, figColor);

    xSlider.ValueChangedFcn     = @(sld,ev)updateXWindow(sld.Value, xStepDropdown.Value);
    xStepDropdown.ValueChangedFcn = @(dd,ev)updateXWindow(xSlider.Value, dd.Value);

    % bottom right panel: channel slider and figure axis
    bottomRightGrid_L = uigridlayout(bottomRightGrid,[3 2],...
        'RowHeight',{'1x','8x','1x'},...
        'ColumnWidth',{70,'1x'},... 
        'Padding',[5 5 5 5],...
        'ColumnSpacing',5);
    bottomRightGrid_L.Layout.Row = [1 3];
    bottomRightGrid_L.Layout.Column = 1;
    applyColorScheme(bottomRightGrid_L, figColor);

    yPanel = uipanel(bottomRightGrid_L,'BackgroundColor',figColor,'BorderColor',figColor);
    yPanel.Layout.Row = [1 3];
    yPanel.Layout.Column = 1;
    applyColorScheme(yPanel, figColor);

    yGrid2 = uigridlayout(yPanel,[4 1],...
        'RowHeight',{'1x', '1x' ,'20x', '1x'},...
        'ColumnWidth',{'fit'},...
        'Padding',[2 2 2 2],...
        'RowSpacing',5);
    applyColorScheme(yGrid2, figColor);

    ySlider = uislider(yGrid2,'Orientation','vertical','Value', 1,'MajorTicksMode','auto', 'Limits',[1 cfg.numChannels]);        
    ySlider.MajorTicks = 0 : 10: cfg.numChannels;
    ySlider.Layout.Row = [3 4];
    ySlider.Layout.Column = 1;
    applyColorScheme(ySlider, figColor);

    channelValues = 2.^(0:10);
    channelValues(channelValues > cfg.numChannels) = [];
    channelValues = arrayfun(@num2str, channelValues, 'UniformOutput', false);
    yStepDropdown = uidropdown(yGrid2,...
        'Items',[{'full'}, channelValues],...
        'Value','full');
    yStepDropdown.Layout.Row = 2;
    yStepDropdown.Layout.Column = 1;
    applyColorScheme(yStepDropdown, figColor);

    lblYStepSize = uilabel(yGrid2,'Text','# Channels:',...
        'HorizontalAlignment','right');
    lblYStepSize.Layout.Row = 1;
    lblYStepSize.Layout.Column = 1;
    applyColorScheme(lblYStepSize, figColor);

    ySlider.ValueChangedFcn      = @(sld,ev)updateYWindow(sld.Value, yStepDropdown.Value);
    yStepDropdown.ValueChangedFcn= @(dd,ev)updateYWindow(ySlider.Value, dd.Value);

    axConfig = uiaxes(bottomRightGrid);
    axConfig.Layout.Row = [1 2];
    axConfig.Layout.Column = [2 3];
    axConfig.Color   = [0.2 0.2 0.2];
    axConfig.XColor  = 1-figColor;
    axConfig.YColor  = 1-figColor;
    axConfig.GridColor = [0.7 0.7 0.7];
    axConfig.TickDir   =  "out";
    title(axConfig,'Raw/Filtered Data','Color','white');
    xlabel(axConfig,'Time (s)','Color','white');
    ylabel(axConfig,'Channel #','Color','white');

    % Initialize X/Y window
    xWindowStart = 0;  
    xWindowEnd   = 1;  
    xStepVal     = 1;
    yWindowStart = 0;  
    yWindowEnd   = cfg.numChannels; 
    yStepVal     = cfg.numChannels;

%% Tab1 - Sub-Tab (b) Extended Config 
cfg = kiaSort_extended_configs(cfg);
exCfgPanel = uipanel(configExtTab, ...
    'Title','Extended Config', ...
    'Scrollable','on', ...
    'BackgroundColor',figColor, ...
    'Units','normalized','Position',[0 0 1 1]);
applyColorScheme(exCfgPanel, figColor);

exCfgFields = fieldnames(cfg);
totalExRows = ceil(numel(exCfgFields)/2);
exCfgGrid = uigridlayout(exCfgPanel, ...
    'RowHeight', repmat({'fit'},1,totalExRows), ...
    'ColumnWidth',{'1.5x','0.5x','1.5x','0.5x'}, ...
    'RowSpacing',20, ...
    'ColumnSpacing',8, ...
    'Scrollable','on', ...
    'Padding',[10 10 10 10]);
applyColorScheme(exCfgGrid, figColor);

  rowIndex = 5;  
    colIndex = 1;
    for iF = numel(cfgFields) : numel(exCfgFields) -1
        fName = exCfgFields{iF};
        fVal  = cfg.(fName);

        lbl = uilabel(exCfgGrid, 'Text', cfg.varName.(exCfgFields{iF}),...
            'HorizontalAlignment','right');
        lbl.Layout.Row    = rowIndex;
        lbl.Layout.Column = colIndex;
        applyColorScheme(lbl, figColor);

        if isfield(choiceMap, fName)
            ddItems = choiceMap.(fName);
            ddVal   = char(fVal);
            dd = uidropdown(exCfgGrid,...
                'Items', ddItems,...
                'Value', ddVal,...
                'ValueChangedFcn',@(ui,ev)updateDropdown(fName, ui.Value));
            dd.Layout.Row    = rowIndex;
            dd.Layout.Column = colIndex+1;
            applyColorScheme(dd, figColor);

        elseif islogical(fVal)
            cb = uicheckbox(exCfgGrid, 'Value',fVal, 'Text','',...
                'ValueChangedFcn',@(ui,ev)updateBoolean(fName, ui.Value));
            cb.Layout.Row    = rowIndex;
            cb.Layout.Column = colIndex+1;
            applyColorScheme(cb, figColor);

        elseif isnumeric(fVal) && numel(fVal)==2
            subP = uipanel(exCfgGrid,'BorderType','none');
            subP.Layout.Row    = rowIndex;
            subP.Layout.Column = colIndex+1;
            applyColorScheme(subP, figColor);

            miniG = uigridlayout(subP,[1 2],'ColumnWidth',{'1x','1x'},...
                'Padding',[0 0 0 0]);
            applyColorScheme(miniG, figColor);

            e1 = uieditfield(miniG,'numeric','Value',fVal(1),...
                'ValueChangedFcn',@(ui,ev)setArrayField(fName,1,ui.Value));
            e1.Layout.Row=1; 
            e1.Layout.Column=1;
            applyColorScheme(e1, figColor);

            e2 = uieditfield(miniG,'numeric','Value',fVal(2),...
                'ValueChangedFcn',@(ui,ev)setArrayField(fName,2,ui.Value));
            e2.Layout.Row=1; 
            e2.Layout.Column=2;
            applyColorScheme(e2, figColor);

        elseif isnumeric(fVal)
            numF = uieditfield(exCfgGrid,'numeric','Value',fVal,...
                'ValueChangedFcn',@(ui,ev)updateNumericField(fName, ui.Value));
            numF.Layout.Row = rowIndex;
            numF.Layout.Column = colIndex+1;
            applyColorScheme(numF, figColor);

        elseif isstring(fVal) || ischar(fVal)
            txt = uieditfield(exCfgGrid,'text','Value',char(fVal),...
                'ValueChangedFcn',@(ui,ev)updateTextOrNumeric(fName, ui.Value));
            txt.Layout.Row = rowIndex;
            txt.Layout.Column = colIndex+1;
            applyColorScheme(txt, figColor);
        else
            txt = uieditfield(exCfgGrid,'text','Value',mat2str(fVal),...
                'ValueChangedFcn',@(ui,ev)setField(fName, ui.Value));
            txt.Layout.Row = rowIndex;
            txt.Layout.Column = colIndex+1;
            applyColorScheme(txt, figColor);
        end

        if colIndex==1
            colIndex=3;
        else
            colIndex=1;
            rowIndex=rowIndex+1;
        end
    end
cfg = kiaSort_hidden_configs(cfg);

    %% Tab1 - Sub-Tab (c) Hyperparameters 
hpPanel = uipanel(hyperSubTab, ...
    'Title','Hyperparameters', ...
    'Scrollable','on', ...
    'BackgroundColor',figColor, ...
    'Units','normalized','Position',[0 0 1 1]);
applyColorScheme(hpPanel, figColor);

hpFields = fieldnames(hp);
totalHpRows = ceil(numel(hpFields)/2);
hpGrid = uigridlayout(hpPanel, ...
    'RowHeight', repmat({'fit'},1,totalHpRows), ...
    'ColumnWidth',{'1.5x','0.5x','1.5x','0.5x'}, ...
    'RowSpacing',20, ...
    'ColumnSpacing',8, ...
    'Scrollable','on', ...
    'Padding',[10 10 10 10]);
applyColorScheme(hpGrid, figColor);

dynamicPanels = struct();
dynamicEdits  = struct();

rIdx = 1;
cIdx = 1;
for iF = 1:numel(hpFields)-1
    fName = hpFields{iF};
    fVal  = hp.(fName);

    lblHP = uilabel(hpGrid, 'Text', hp.varName.(fName), ...
        'HorizontalAlignment','right');
    lblHP.Layout.Row    = rIdx;
    lblHP.Layout.Column = cIdx;
    applyColorScheme(lblHP, figColor);

    if isfield(hpChoiceMap, fName)
        ddItems = hpChoiceMap.(fName);
        ddVal   = convertToStringForDropdown(fVal);
        if ~ismember(ddVal, ddItems)
            ddVal = ddItems{1};
        end
        ddHP = uidropdown(hpGrid, ...
            'Items', ddItems, ...
            'Value', ddVal, ...
            'ValueChangedFcn', @(ui,ev)updateHPDropdown(fName, ui.Value));
        ddHP.Layout.Row    = rIdx;
        ddHP.Layout.Column = cIdx+1;
        applyColorScheme(ddHP, figColor);

    elseif islogical(fVal)
        cbHP = uicheckbox(hpGrid, 'Value', fVal, 'Text', '', ...
            'ValueChangedFcn', @(ui,ev)updateHPBoolean(fName, ui.Value));
        cbHP.Layout.Row    = rIdx;
        cbHP.Layout.Column = cIdx+1;
        applyColorScheme(cbHP, figColor);

    elseif isnumeric(fVal) && numel(fVal) > 1
        subPH = uipanel(hpGrid, 'BorderType', 'none');
        subPH.Layout.Row    = rIdx;
        subPH.Layout.Column = cIdx+1;
        applyColorScheme(subPH, figColor);

        arrVal = fVal(:)';
        miniGH = uigridlayout(subPH, [1 numel(arrVal)], ...
            'ColumnWidth', repmat({'1x'},1,numel(arrVal)), ...
            'Padding', [0 0 0 0], ...
            'ColumnSpacing', 5);
        applyColorScheme(miniGH, figColor);

        if any(strcmp(fName, {'hiddenLayerSize','numFilters','dropoutRate'}))
            dynamicPanels.(fName) = miniGH;
        end

        for idx = 1:numel(arrVal)
            eX = uieditfield(miniGH, 'numeric', 'Value', arrVal(idx), ...
                'ValueChangedFcn', @(ui,ev)setHpArrayField(fName, idx, ui.Value));
            eX.Layout.Row = 1; 
            eX.Layout.Column = idx;
            applyColorScheme(eX, figColor);
            if any(strcmp(fName, {'hiddenLayerSize','numFilters','dropoutRate'}))
                dynamicEdits.(fName)(idx) = eX;
            end
        end

    elseif isnumeric(fVal)
        % Scalar numeric field.
        numHP = uieditfield(hpGrid, 'numeric', 'Value', fVal, ...
            'ValueChangedFcn', @(ui,ev)updateHPNumeric(fName, ui.Value));
        numHP.Layout.Row = rIdx;
        numHP.Layout.Column = cIdx+1;
        applyColorScheme(numHP, figColor);

    elseif isstring(fVal) || ischar(fVal)
        txtHP = uieditfield(hpGrid, 'text', 'Value', char(fVal), ...
            'ValueChangedFcn', @(ui,ev)updateHPTextOrNumeric(fName, ui.Value));
        txtHP.Layout.Row = rIdx;
        txtHP.Layout.Column = cIdx+1;
        applyColorScheme(txtHP, figColor);

    else
        txtHP = uieditfield(hpGrid, 'text', 'Value', mat2str(fVal), ...
            'ValueChangedFcn', @(ui,ev)setHPField(fName, ui.Value));
        txtHP.Layout.Row = rIdx;
        txtHP.Layout.Column = cIdx+1;
        applyColorScheme(txtHP, figColor);
    end

    if cIdx == 1
        cIdx = 3;
    else
        cIdx = 1;
        rIdx = rIdx + 1;
    end
end

function updateHPNumeric(fName, newVal)
    hp.(fName) = newVal;   
    if strcmp(fName, 'numHidLayer')
        updateDynamicArrayField('hiddenLayerSize', newVal);
        updateDynamicArrayField('dropoutRate', newVal);
    elseif strcmp(fName, 'numConvBlocks')
        updateDynamicArrayField('numFilters', newVal);
        updateDynamicArrayField('dropoutRate', newVal);
    end
end

function updateDynamicArrayField(fieldName, newCount)
    if isfield(dynamicPanels, fieldName)
        container = dynamicPanels.(fieldName);
        delete(container.Children);
        
        oldVals = hp.(fieldName);
        if numel(oldVals) >= newCount
            newVals = oldVals(1:newCount);
        else
            newVals = [oldVals, zeros(1, newCount - numel(oldVals))];
        end
        hp.(fieldName) = newVals;
        
        newGrid = uigridlayout(container.Parent, [1 newCount], ...
            'ColumnWidth', repmat({'1x'},1,newCount), ...
            'Padding', [0 0 0 0], 'ColumnSpacing', 5);
        dynamicPanels.(fieldName) = newGrid; 
        dynamicEdits.(fieldName) = [];  % reset
        applyColorScheme(newGrid, figColor);

        for iCount = 1:newCount
            eX = uieditfield(newGrid, 'numeric', 'Value', newVals(iCount), ...
                'ValueChangedFcn', @(ui,ev)setHpArrayField(fieldName, iCount, ui.Value));
            eX.Layout.Row = 1;
            eX.Layout.Column = iCount;            
            dynamicEdits.(fieldName)(iCount) = eX;
            applyColorScheme(eX, figColor);
        end
    end
end
 
    %% Tab2: Channels

    gridTab2 = uigridlayout(tab2,[3 3],...
        'RowHeight',{'fit','fit','1x'},...
        'ColumnWidth',{'fit','fit','1x'},...
        'Padding',10,...
        'RowSpacing',10);
    applyColorScheme(gridTab2, figColor);


    btnLoadChannels = uibutton(gridTab2, ...
        'Text','Plot Channels', ...
        'ButtonPushedFcn',@(btn,ev)onLoadChannels());
    btnLoadChannels.Layout.Row = 1;
    btnLoadChannels.Layout.Column = 1;
    applyColorScheme(btnLoadChannels, figColor);

    

    %% Tab3: Curation

    gridTab3 = uigridlayout(tab3,[3 6],...
        'RowHeight',{'fit','fit','1x'},...
        'ColumnWidth',{'fit','fit','fit','fit','1x','4x'},...
        'Padding',10,...
        'RowSpacing',10);
    applyColorScheme(gridTab3, figColor);

    btnLoadSorting = uibutton(gridTab3, ...
        'Text','Plot Sorting', ...
        'ButtonPushedFcn',@(btn,ev)onPlotSorting());
    btnLoadSorting.Layout.Row = 1;
    btnLoadSorting.Layout.Column = 1;
    applyColorScheme(btnLoadSorting, figColor);

    btnAltRes = uibutton(gridTab3,'Text','Alt Res Folder',...
        'ButtonPushedFcn',@(btn,ev)onSelectAltResFolder(),'Tooltip',...
        sprintf('Select an alternative result folder if \n you have used Sort Only option'));
    btnAltRes.Layout.Column = 2;
    btnAltRes.Layout.Row = 1;
    applyColorScheme(btnAltRes, figColor);

    btnClear = uibutton(gridTab3,'Text','Clear',...
        'ButtonPushedFcn',@(btn,ev)onClearSorting());
    btnClear.Layout.Column = 3;
    btnClear.Layout.Row = 1;
    applyColorScheme(btnClear, figColor);

        %% Tab4:  Clusters

    gridTab4 = uigridlayout(tab4,[2 1],...
        'RowHeight',{'fit','1x'},...
        'Padding',10,...
        'RowSpacing',10);
    applyColorScheme(gridTab4, figColor);

    topGroupPanel = uipanel(gridTab4,'BorderType','none','BackgroundColor',figColor);
    topGroupPanel.Layout.Row = 1;
    topGroupPanel.Layout.Column = 1;
    applyColorScheme(topGroupPanel, figColor);

    topGroupGrid = uigridlayout(topGroupPanel,[1 5],...
        'ColumnWidth',{'fit','fit','fit','fit','fit'},...
        'RowHeight',{'fit'},...
        'Padding',[0 0 0 0],...
        'ColumnSpacing',5);
    applyColorScheme(topGroupGrid, figColor);

    % Button "Plot Single Channel"
    btnPlotSingle = uibutton(topGroupGrid,'Text','Plot Single Channel',...
        'ButtonPushedFcn',@(btn,ev)onPlotSingleChannel());
    btnPlotSingle.Layout.Row = 1;
    btnPlotSingle.Layout.Column = 1;
    applyColorScheme(btnPlotSingle, figColor);

    % Button "Plot All Groups"
    btnPlotAll = uibutton(topGroupGrid,'Text','Plot All Groups',...
        'ButtonPushedFcn',@(btn,ev)onPlotAllGroups());
    btnPlotAll.Layout.Row = 1;
    btnPlotAll.Layout.Column = 2;
    applyColorScheme(btnPlotAll, figColor);

        btnPlotAll = uibutton(topGroupGrid,'Text','Clear',...
        'ButtonPushedFcn',@(btn,ev)onClearGroupPlots());
    btnPlotAll.Layout.Row = 1;
    btnPlotAll.Layout.Column = 3;
    applyColorScheme(btnPlotAll, figColor);

    % panel for visualizing the group plots
    groupsPlotPanel = uipanel(gridTab4,'Title','Groups Visualization',...
        'BackgroundColor',figColor);
    groupsPlotPanel.Layout.Row = 2;
    groupsPlotPanel.Layout.Column = 1;
    applyColorScheme(groupsPlotPanel, figColor);


%% Tab 5: Help
helpSplitLayout = uigridlayout(tab5, [1 2], ...
    'ColumnWidth', {'1x','2x'}, 'RowHeight', {'1x'}, ...
    'Padding', 10, 'ColumnSpacing', 10);
applyColorScheme(helpSplitLayout, figColor);

helpTree = uitree(helpSplitLayout);
helpTree.Layout.Row = 1;
helpTree.Layout.Column = 1;
helpTree.BackgroundColor = figColor;
helpTree.FontColor = 1 - figColor;

% Create the root node
rootNode = uitreenode(helpTree, 'Text', 'KIASort GUI Help');
applyColorScheme(rootNode, figColor);

%%  MAIN TAB 
nodeMain = uitreenode(rootNode, 'Text', 'Main Tab');

%% Run & Sort Panel (Main Config) 
nodeRunSort = uitreenode(nodeMain, 'Text', 'Run & Sort Panel');
applyColorScheme(nodeRunSort, figColor);

nodeField = uitreenode(nodeRunSort, 'Text', 'Data Format'); % was 'dataType'
nodeField.UserData = 'Data format (e.g., int16, int32, float32, uint16).';

nodeField = uitreenode(nodeRunSort, 'Text', 'Num. Channels'); % was 'numChannels'
nodeField.UserData = 'Total number of channels.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Sampling Freq (Hz)'); % was 'samplingFrequency'
nodeField.UserData = 'Sampling frequency in Hz.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Common Ref. Type'); % was 'commonRef'
nodeField.UserData = 'Type of common referencing (e.g., median).';

nodeField = uitreenode(nodeRunSort, 'Text', 'Denoising'); % was 'commonRef'
nodeField.UserData = 'Denoising the neural recording by subtracting shared correlation accross channnels.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Percentile Chan. Denoising'); % was 'commonRef'
nodeField.UserData = 'Percent channels to extract shared correlation (default: median)';

nodeField = uitreenode(nodeRunSort, 'Text', 'Noise Corr Thr.'); % was 'commonRef'
nodeField.UserData = 'Threshold for correlation of shared noise between channels to be subtracted.';

nodeField = uitreenode(nodeRunSort, 'Text', 'High Freq. Noise'); % was 'commonRef'
nodeField.UserData = 'Denoising the neural recording by subtracting shared high frequency noise only.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Use GPU'); % was 'useGPU'
nodeField.UserData = 'Enable GPU acceleration if available.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Parallel Proc'); % was 'parallelProcessing'
nodeField.UserData = 'Enable parallel processing if supported.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Half Channels Extract'); % was 'num_channel_extract'
nodeField.UserData = 'Number of channels to extract on either side of a spike.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Bandpass Range (Hz)'); % was 'bandpass'
nodeField.UserData = 'Frequency range (Hz) used for filtering.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Min Thr.'); % combining 'min_threshold' and 'max_threshold'
nodeField.UserData = 'Median deviation detection thresholds. It is best between 7-9 for most datasets, but can be explored in the channel window.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Num. Sample Chunks'); % was 'numSampleChunks'
nodeField.UserData = 'Number of chunks to extract.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Sort Only'); % was 'minRate'
nodeField.UserData = 'True if you have trained the classifier on anothe data set, and want to run the sorting on a new dataset.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Sort Chunk Duration'); % was 'sortingChunkDuration'
nodeField.UserData = 'Duration (s) for sorting each chunk of signal.';

nodeField = uitreenode(nodeRunSort, 'Text', 'Save Sorted Waveforms'); % was 'extractWaveform'
nodeField.UserData = '(true/false) Extract and save spike waveforms after sorting.';

%% Extended Config 
nodeExtConfig = uitreenode(nodeMain, 'Text', 'Extended Config');
applyColorScheme(nodeExtConfig, figColor);

nodeField = uitreenode(nodeExtConfig, 'Text', 'Num. Parallel Workers'); % was 'numParallelWorker'
nodeField.UserData = 'Number of parallel workers (or ''auto'').';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Bad Ch. Factor'); % was 'badChannel_factor'
nodeField.UserData = 'Factor to flag channels with low spike power.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Noisy Ch. Factor'); % was 'noisyChannel_factor'
nodeField.UserData = 'Factor to flag noisy channels.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Min Clust Rate'); % was 'maxClusteringPoints'
nodeField.UserData = 'Clustering rate limits (dbScan).';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Max Clust Rate'); % was 'maxClusteringPoints'
nodeField.UserData = 'Clustering rate limits (dbScan).';

nodeField = uitreenode(nodeExtConfig, 'Text', '# UMAP Components'); % was 'umapNComp'
nodeField.UserData = 'Number of UMAP components.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Test Fraction'); % was 'testFraction'
nodeField.UserData = 'Fraction of data for testing the classifier.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Spike Dur. (ms)'); % was 'spikeDuration'
nodeField.UserData = 'Duration (ms) of a spike used for classification.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Spike Distance (ms)'); % was 'spikeDistance'
nodeField.UserData = 'Minimum time (ms) between spikes.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Sample Chunk Dur. (s)'); % was 'sampleChunkDuration'
nodeField.UserData = 'Duration (s) for each sample chunk.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Max Sample Span (min)'); % was 'maxSampleSpan'
nodeField.UserData = 'Maximum time span (min) for sample extraction.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Sorting: Model Type'); % was 'modelType'
nodeField.UserData = 'Sorting model type (e.g., svm, mlp, cnn, etc.). SVM is fastest, most accurate, and recommended.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Sorting: Use PCA'); % was 'usePCA'
nodeField.UserData = 'Whether to use PCA during training (true/false).';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Sorting: # PCA Comp.'); % was 'nPCAcomp'
nodeField.UserData = 'Number of PCA components.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Clustering Spike Dur'); % was 'clusteringSpikeDuration'
nodeField.UserData = 'Spike duration (ms) for clustering.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'High Corr Thr'); % was 'corrThresholdHigh'
nodeField.UserData = 'High correlation threshold for merging clusters.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Realign Corr Thr'); % was 'corrThresholdMisAlign'
nodeField.UserData = 'Correlation threshold for spike realignment.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Spike Density Corr Thr'); % was 'corrThresholdSpkDensity'
nodeField.UserData = 'Negative correlation threshold for spike density merging.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Amp. Var. Thr'); % was 'ampVarThreshold'
nodeField.UserData = 'Maximum allowed amplitude variance for merging groups.';

nodeField = uitreenode(nodeExtConfig, 'Text', 'Amp. Var. Thr Xcorr'); % was 'ampVarThreshold'
nodeField.UserData = 'Maximum allowed amplitude variance for merging misaligned groups.';

%% Hyperparameters
nodeHyper = uitreenode(nodeMain, 'Text', 'Hyperparameters');
applyColorScheme(nodeHyper, figColor);

% MLP Parameters
nodeField = uitreenode(nodeHyper, 'Text', 'MLP: Number of Hidden Layers'); % was 'numHidLayer'
nodeField.UserData = 'Number of hidden layers in the MLP.';

nodeField = uitreenode(nodeHyper, 'Text', 'MLP: Hidden Layer Sizes'); % was 'hiddenLayerSize'
nodeField.UserData = 'Vector specifying neurons in each hidden layer.';

nodeField = uitreenode(nodeHyper, 'Text', 'MLP: Loss Function'); % was 'mlpLossFcn'
nodeField.UserData = 'Loss function for MLP (mse or mae).';

% CNN Parameters
nodeField = uitreenode(nodeHyper, 'Text', 'CNN: Number of Conv Blocks'); % was 'numConvBlocks'
nodeField.UserData = 'Number of convolutional blocks in the CNN.';

nodeField = uitreenode(nodeHyper, 'Text', 'CNN: Filter Size'); % was 'filterSize'
nodeField.UserData = 'Size of the CNN filters.';

nodeField = uitreenode(nodeHyper, 'Text', 'CNN: Filters per Block'); % was 'numFilters'
nodeField.UserData = 'Vector of filters per block.';

nodeField = uitreenode(nodeHyper, 'Text', 'CNN: Pooling Size'); % was 'poolSize'
nodeField.UserData = 'Pooling size in the CNN.';

nodeField = uitreenode(nodeHyper, 'Text', 'CNN: Loss Function'); % was 'cnnLossFcn'
nodeField.UserData = 'Loss function for CNN (crossentropy or dlconv).';

nodeField = uitreenode(nodeHyper, 'Text', 'CNN/MLP: Dropout Rates'); % was 'dropoutRate'
nodeField.UserData = 'Dropout rates for each CNN layer (also used in MLP).';

% Training Hyperparameters
nodeField = uitreenode(nodeHyper, 'Text', 'Training: Max Epochs'); % was 'maxEpochs'
nodeField.UserData = 'Maximum training epochs.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: Mini Batch Size'); % was 'miniBatchSize'
nodeField.UserData = 'Mini-batch size for training.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: Initial Learning Rate'); % was 'initialLearningRate'
nodeField.UserData = 'Starting learning rate.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: Learning Rate Drop Factor'); % was 'learningRateDropFactor'
nodeField.UserData = 'Factor by which the learning rate is reduced.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: Learning Rate Drop Period'); % was 'learningRateDropPeriod'
nodeField.UserData = 'Epoch interval for reducing the learning rate.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: L2 Regularization'); % was 'l2Reg'
nodeField.UserData = 'L2 regularization strength.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: Gradient Threshold'); % was 'gradientThreshold'
nodeField.UserData = 'Threshold for gradient clipping.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: Validation Ratio'); % was 'validationRatio'
nodeField.UserData = 'Fraction of data used for validation.';

nodeField = uitreenode(nodeHyper, 'Text', 'Training: Execution Environment'); % was 'executionEnvironment'
nodeField.UserData = 'Execution environment (auto, cpu, or gpu).';


% SVM Parameters
nodeField = uitreenode(nodeHyper, 'Text', 'SVM: Kernel Function'); % was 'KernelFunction'
nodeField.UserData = 'SVM kernel type (e.g., polynomial).';

nodeField = uitreenode(nodeHyper, 'Text', 'SVM: Polynomial Order'); % was 'PolynomialOrder'
nodeField.UserData = 'Order of the polynomial kernel.';

nodeField = uitreenode(nodeHyper, 'Text', 'SVM: Box Constraint'); % was 'BoxConstraint'
nodeField.UserData = 'SVM box constraint parameter.';

nodeField = uitreenode(nodeHyper, 'Text', 'SVM: Kernel Scale'); % was 'KernelScale'
nodeField.UserData = 'Scale parameter for the SVM kernel.';

% Miscellaneous
nodeField = uitreenode(nodeHyper, 'Text', 'Use Parallel Processing'); % was 'UseParallel'
nodeField.UserData = 'Enable parallel processing during training.';

%%  CHANNELS TAB 
nodeChannels = uitreenode(rootNode, 'Text', 'Channels Tab');
applyColorScheme(nodeChannels, figColor);

nodeField = uitreenode(nodeChannels, 'Text', 'Plot Channels');
nodeField.UserData = 'Load a channel mapping file and plots the data. If samples are extracted, the check mark for excluded channels is off.  Users can mannually also select channels to be included/excluded and save the data for the sorting.';

nodeField = uitreenode(nodeChannels, 'Text', 'Save botton');
nodeField.UserData = 'Save the modified channel inclusion info.';

nodeField = uitreenode(nodeChannels, 'Text', 'MAD Threshold');
nodeField.UserData = 'Threshold of Median Average Deviation, which can be changed through shifting the slider.';


%%  Clusters TAB 
nodeGroups = uitreenode(rootNode, 'Text', 'Clusters Tab');
applyColorScheme(nodeGroups, figColor);

nodeField = uitreenode(nodeGroups, 'Text', 'Plot Single Channel');
nodeField.UserData = 'View sorted samples for each individual channel. In the waveform window, the main channel (the one that a group is going to be sorted on) is shown in thicker line';

nodeField = uitreenode(nodeGroups, 'Text', 'Plot All Groups');
nodeField.UserData = 'View an overview of grouped sample sorting results across all channels. This tab shows the mean waveform for each group, and can be used for sanity check and exploration of unified sorted sample groups.';

%%  Curation TAB
nodeSorting = uitreenode(rootNode, 'Text', 'Curation Tab');
applyColorScheme(nodeSorting, figColor);

%% Spike Group Controller Panel 
nodeSpikeGroup = uitreenode(nodeSorting, 'Text', 'Spike Group Controller');
applyColorScheme(nodeSpikeGroup, figColor);

% Group Management fields
nodeField = uitreenode(nodeSpikeGroup, 'Text', 'Remove');
nodeField.UserData = 'Marks selected groups as removed (sets label and channel to NaN).';

nodeField = uitreenode(nodeSpikeGroup, 'Text', 'Merge');
nodeField.UserData = 'Merges selected groups into a target group (selected via its radio button in "Merge To").';

nodeField = uitreenode(nodeSpikeGroup, 'Text', 'MUA/SUA');
nodeField.UserData = 'Toggles the classification (SUA vs. MUA) for selected groups.';

nodeField = uitreenode(nodeSpikeGroup, 'Text', 'Undo');
nodeField.UserData = 'Restores selected groups to their original settings.';

nodeField = uitreenode(nodeSpikeGroup, 'Text', 'Reset');
nodeField.UserData = 'Resets all groups to the original sorted state and clears selections.';

nodeField = uitreenode(nodeSpikeGroup, 'Text', 'Limit');
nodeField.UserData = 'Automatically selects the stable areas of spiking activity.';

nodeField = uitreenode(nodeSpikeGroup, 'Text', 'Drop rate');
nodeField.UserData = 'Drop in the spike rate threshold for exclusion.';

nodeField = uitreenode(nodeSpikeGroup, 'Text', 'Realign (ms) Slider');
nodeField.UserData = 'Adjusts spike timing by shifting spike indices and circularly shifting waveforms.';

%% Similarity Processing Panel 
nodeSimilarity = uitreenode(nodeSpikeGroup, 'Text', 'Similarity Processing Panel');
applyColorScheme(nodeSimilarity, figColor);

nodeField = uitreenode(nodeSimilarity, 'Text', 'Process Button');
nodeField.UserData = 'Samples waveforms, computes similarity scores (using XCorr, KL-div, or Bhattacharyya).';

nodeField = uitreenode(nodeSimilarity, 'Text', 'Metric Dropdown');
nodeField.UserData = 'Choose the similarity metric to compare different groups.';

nodeField = uitreenode(nodeSimilarity, 'Text', 'Similarity Field');
nodeField.UserData = 'Set the threshold for similarity between groups. A value between 0 to 1. In XCorr, it is the 1-max(XCorr) value.';

nodeField = uitreenode(nodeSimilarity, 'Text', 'Amp. Var. Field');
nodeField.UserData = 'Set the maximum acceptable amplitude variance between groups. A value between 0 to 1.';

nodeField = uitreenode(nodeSimilarity, 'Text', 'Distance Field & Slider');
nodeField.UserData = 'Set the spatial threshold (in µm) between channels. This means it only looks for similarity between channels closer than this distance in depth. If no coordinates are not provided, then it considers channel numbers.';

nodeField = uitreenode(nodeSimilarity, 'Text', 'Next/Previous Buttons');
nodeField.UserData = 'Navigate through group pairs that meet the similarity criteria.';

%% Group Selection Panel 
nodeGroupSelection = uitreenode(nodeSpikeGroup, 'Text', 'Group Selection Panel');
applyColorScheme(nodeGroupSelection, figColor);

nodeField = uitreenode(nodeGroupSelection, 'Text', 'Firing Rate');
nodeField.UserData = 'Use the numeric field and slider to set the minimum firing rate.';

nodeField = uitreenode(nodeGroupSelection, 'Text', 'ISI %');
nodeField.UserData = 'Use the numeric field and slider to set the maximum allowed ISI violations (%).';

nodeField = uitreenode(nodeGroupSelection, 'Text', 'SNR');
nodeField.UserData = 'Use the numeric field and slider to set the minimum signal-to-noise ratio.';

nodeField = uitreenode(nodeGroupSelection, 'Text', 'Polarity Dropdown');
nodeField.UserData = 'Filter groups by spike polarity (All (all spikes), Main Neg. (main channel is negatively polarized), Main Pos. (main channel is positively polarized), All Pos. (all channels are positively polarized) , 1 Chan Neg. (at least one channel is negatively polarized)).';

nodeField = uitreenode(nodeGroupSelection, 'Text', 'Vis Toggle');
nodeField.UserData = 'Checkbox to auto-update visualization based on criteria.';

nodeField = uitreenode(nodeGroupSelection, 'Text', 'Inclusion/Exclusion Checkboxes');
nodeField.UserData = 'Select groups that meet (or do not meet) the criteria.';

nodeField = uitreenode(nodeGroupSelection, 'Text', 'Next/Previous Buttons');
nodeField.UserData = 'Navigate through each group that meet the Inclusion/Exclusion criteria.';

%% Unified Groups Panel 
nodeUnifiedGroups = uitreenode(nodeSpikeGroup, 'Text', 'Unified Groups Panel');
applyColorScheme(nodeUnifiedGroups, figColor);

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'Group #');
nodeField.UserData = 'Unique label for each group.';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'Channel #');
nodeField.UserData = 'Channel number associated with the group.';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'Isolation');
nodeField.UserData = 'Classification as SUA+, SUA, MUA+, and MUA. (default not assigned NA.) ';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'Firing Rate');
nodeField.UserData = 'Calculated firing rate.';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'SNR');
nodeField.UserData = 'Signal-to-noise ratio (e.g., 1 + detectability).';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'ISI %');
nodeField.UserData = 'Percentage of inter-spike interval violations.';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'Vis. Checkbox');
nodeField.UserData = 'Include/exclude the group from plots.';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'Select Checkbox');
nodeField.UserData = 'Select the group for operations (e.g., merging).';

nodeField = uitreenode(nodeUnifiedGroups, 'Text', 'Merge Target Radio');
nodeField.UserData = 'Radio button to choose a target group for merging.';

%% Visualization Control Panel
nodeVizControl = uitreenode(nodeSorting, 'Text', 'Visualization Control Panel');
applyColorScheme(nodeVizControl, figColor);

nodeField = uitreenode(nodeVizControl, 'Text', 'Plot Type Dropdown');
nodeField.UserData = 'Choose view (Amplitude, Waveform, CCG, ISI, Density, PC_channel, Trace, Features, Multiple).';

nodeField = uitreenode(nodeVizControl, 'Text', 'Data Type Dropdown');
nodeField.UserData = 'Select filtered (bandpass) or raw data.';

nodeField = uitreenode(nodeVizControl, 'Text', '# Channel Show Dropdown');
nodeField.UserData = 'Set the number of neighboring channels to display (or use default).';

nodeField = uitreenode(nodeVizControl, 'Text', '# Spikes Dropdown');
nodeField.UserData = 'Set the maximum number of spike waveforms to display.';

nodeField = uitreenode(nodeVizControl, 'Text', 'Waveform Duration Dropdown');
nodeField.UserData = 'Set spike waveform duration (ms).';

nodeField = uitreenode(nodeVizControl, 'Text', 'Waveform Type Dropdown');
nodeField.UserData = 'Choose between individual or mean waveform display.';

nodeField = uitreenode(nodeVizControl, 'Text', 'Bins Dropdown');
nodeField.UserData = 'For CCG plots, set number of bins (lags).';

nodeField = uitreenode(nodeVizControl, 'Text', 'Smoothness Dropdown');
nodeField.UserData = 'Adjust smoothing for CCG curves.';

nodeField = uitreenode(nodeVizControl, 'Text', 'ISI Threshold Dropdown');
nodeField.UserData = 'Set threshold (ms) for ISI violation detection.';

nodeField = uitreenode(nodeVizControl, 'Text', 'Amplitude Scale Slider');
nodeField.UserData = 'Adjust vertical scaling of waveforms.';

nodeField = uitreenode(nodeVizControl, 'Text', 'Line Width Slider');
nodeField.UserData = 'Set thickness of plotted lines.';

nodeField = uitreenode(nodeVizControl, 'Text', 'Export Fig. Button');
nodeField.UserData = 'Export the current figure as an image or vector graphic.';

nodeField = uitreenode(nodeVizControl, 'Text', 'Export Format Dropdown');
nodeField.UserData = 'Select the export format (Image or Vector).';

%% Plots Panel
nodePlots = uitreenode(nodeSorting, 'Text', 'Plots Panel');
applyColorScheme(nodePlots, figColor);

% X-Axis Controls
nodeField = uitreenode(nodePlots, 'Text', 'Time Window Slider');
nodeField.UserData = 'Adjust the starting time (s) of the plot.';

nodeField = uitreenode(nodePlots, 'Text', 'Time Window Dropdown');
nodeField.UserData = 'Set the width (s) of the time window.';

% Y-Axis Controls
nodeField = uitreenode(nodePlots, 'Text', '# Channels Dropdown');
nodeField.UserData = 'Set the vertical span (number of channels) to display (or select ''full'').';

nodeField = uitreenode(nodePlots, 'Text', 'Channel Slider');
nodeField.UserData = 'Adjust the starting channel for the display.';

% Main Axes Description
nodeField = uitreenode(nodePlots, 'Text', 'Main Axes');
nodeField.UserData = ['Displays the visualization according to the selected Plot Type:' newline ...
    '• Amplitude: Peak amplitude of the spike on the main channel.' newline ...
    '• Trace: Line plots of raw/filtered signal with overlaid spike waveforms.' newline ...
    '• Features: 3D scatter plot of PCA-normalized spike features.' newline ...
    '• Waveform: Detailed display of individual or mean spike waveforms.' newline ...
    '• CCG: Tiled cross-/auto-correlograms comparing groups.' newline ...
    '• ISI: Tiled histograms of inter-spike intervals with violation percentages.' newline ...
    '• Density: Tiled density plots showing spike occurrence and presence ratios.' newline ...
    '• ChannelPCs: PC1 and PC2 of each group on each channel.' newline ...
    '• Multiple: Plotting multiple plot types in one window. The size and selection order of each plot can be changed.'];

expand(rootNode);

helpDetails = uitextarea(helpSplitLayout, 'Editable', 'off', ...
    'Value', {'Select a topic from the tree to view details.'});
helpDetails.Layout.Row = 1;
helpDetails.Layout.Column = 2;
helpDetails.BackgroundColor = figColor;
helpDetails.FontColor = 1 - figColor;
helpDetails.WordWrap = 'on';

% Tree Selection Callback
helpTree.SelectionChangedFcn = @(src, event) updateHelpDetails(event.SelectedNodes, helpDetails);

function updateHelpDetails(selectedNodes, helpDetails)
    if isempty(selectedNodes)
        helpDetails.Value = {'Select a topic from the tree.'};
        return;
    end
    node = selectedNodes(1);
    if ~isempty(node.UserData)
        helpText = node.UserData;
    else
        helpText = {node.Text};
    end
    helpDetails.Value = helpText;
    drawnow;
end


%% Tab 6: About

aboutGrid = uigridlayout(tab6, [2 1], ...
    'RowHeight', {'1x','1.75x'}, ...
    'Padding', 10, 'RowSpacing', 10);
applyColorScheme(aboutGrid, figColor);

websitePanel = uipanel(aboutGrid, 'BorderColor','none','BorderType','none');
websitePanel.Layout.Row = 1;
websitePanel.Layout.Column = 1;
applyColorScheme(websitePanel, figColor);

websiteGrid = uigridlayout(websitePanel,[4,16],...
    'RowHeight', {'1x','4x','.5x','1x'}, ...dd
    'Padding', 10, 'RowSpacing', 10, 'ColumnSpacing',1);
applyColorScheme(websiteGrid, figColor);

lblTitle = uilabel(websiteGrid, ...
    'Text', 'About KIASort:', ...    
    'FontSize', 14, ...
    'HorizontalAlignment', 'left', ...
    'FontWeight','bold',...
    'FontColor', 1 - figColor);  
lblTitle.Layout.Row = 1;
lblTitle.Layout.Column = 3;
applyColorScheme(lblTitle, figColor);

imglogo = uiimage(websiteGrid, 'ImageSource', 'Run_KIASORT.png');
imglogo.Layout.Row = 1;
imglogo.Layout.Column = [1 2];
applyColorScheme(imglogo, figColor);

txtDeveloperInfo = uitextarea(websiteGrid, ...
    'Value', {['KIASort is a spike sorting platform developed by Kianoush Banaie Boroujeni (Kia) ' ...
    ['designed to overcome challenges associated with heterogenous drift patterns, and nonlinear waveform variabilities in high-density neural recordings. KIASort integrates knowledge from channel-specific classifiers trained on clustered sample data.' ...
    ' It evaluates channel quality, automatically excludes noisy recordings, and identifies spike classes using a hybrid dimensionality reduction approach combining PCA with nonlinear UMAP embeddings. '],...
'Conventional one-dimensional drift correction methods cannot address heterogeneous neuron-specific drift and nonlinear waveform changes. KIASORT geometry-free, per-neuron tracking approach overcomes these limitations without assumptions about cluster shape or probe motion. ',...
['The multi-tab graphical user interface (GUI) provides high degree of control over each stage, from initial configuration and channel inspection through detailed clustering analyses and classification performance assessments.' ...
' Post-sorting, users can easily curate and validate spike sorting quality through advanced '],...
'visualization tools implemented, including cross-correlograms, inter-spike interval (ISI) histograms, density analyses, similarity detection, and waveform and feature explorations. ',...
'KIASort is a flexible, user-friendly, and geometry-free spike sorting approach, which can be uses for high density laminar recordings as well as different probe geometry.'
]}, ...
    'Editable', 'off', ...
    'FontSize', 14, ...
    'FontColor', 1 - figColor, ...
    'WordWrap', 'on');  
txtDeveloperInfo.Layout.Row = 2;
txtDeveloperInfo.Layout.Column = [3 14];
applyColorScheme(txtDeveloperInfo, figColor);

lblWebsite = uilabel(websiteGrid, ...
    'Text', 'Website:', ...
    'FontSize', 14, ...
    'HorizontalAlignment', 'left', ...
    'FontColor', 1 - figColor);  
lblWebsite.Layout.Row = 3;
lblWebsite.Layout.Column = 6;
applyColorScheme(lblWebsite, figColor);

btnLink = uibutton(websiteGrid, ...
    'Text', 'Get KIASort', ...
    'FontSize', 14, ...
    'FontColor', 'blue', ...     
    'ButtonPushedFcn', @(btn, event) web('https://github.com/banaiek/KIASORT', '-browser'));
btnLink.Layout.Row = 4;
btnLink.Layout.Column = [7 10];
applyColorScheme(btnLink, figColor);

imagesPanel = uipanel(aboutGrid, 'BorderColor','none','BorderType','none');
imagesPanel.Layout.Row = 2;
imagesPanel.Layout.Column = 1;
applyColorScheme(imagesPanel, figColor);

imagesLayout = uigridlayout(imagesPanel, [1 2], ...
    'ColumnWidth', {'1x','1x'}, ...
    'Padding', 10, 'ColumnSpacing', 10);
applyColorScheme(imagesLayout, figColor);

img1 = uiimage(imagesLayout, 'ImageSource', 'Schematic.png');
img1.Layout.Row = 1;
img1.Layout.Column = 1;
applyColorScheme(img1, figColor);

img2 = uiimage(imagesLayout, 'ImageSource', 'Sorting_Sample.png');
img2.Layout.Row = 1;
img2.Layout.Column = 2;
applyColorScheme(img2, figColor);


%%  Callbacks
    % cfg update callbacks
    function updateDropdown(fieldName, newVal)
        cfg.(fieldName) = string(newVal);
        disp([fieldName, ' updated to ', char(cfg.(fieldName))]);
    end

    function updateBoolean(fieldName, val)
        cfg.(fieldName) = val;
        disp([fieldName, ' updated to ', mat2str(val)]);
        if cfg.denoising ~= current_cfg.denoising
            current_cfg.denoising = cfg.denoising;
            if lastPlotted
                onPlotData();
            end
        end

        if cfg.extremeNoise ~= current_cfg.extremeNoise
            current_cfg.extremeNoise = cfg.extremeNoise;
            if lastPlotted
                onPlotData();
            end
        end
    end

    function setArrayField(fieldName, idx, val)
        oldVal = cfg.(fieldName);
        oldVal(idx) = val;
        cfg.(fieldName) = oldVal;
        disp([fieldName, ' updated to ', mat2str(cfg.(fieldName))]);

        if any(cfg.bandpass ~= current_cfg.bandpass)
            current_cfg.bandpass = cfg.bandpass;
            if lastPlotted
                onPlotData();
            end
        end

        if any(cfg.numChannels ~= current_cfg.numChannels)
            current_cfg.numChannels = cfg.numChannels;
            numChannels =  min(cfg.numChannels, length(channel_mapping));
            if lastPlotted
                onPlotData();
            end
        end

        if any(cfg.noisePrctile ~= current_cfg.noisePrctile)
            current_cfg.noisePrctile = cfg.noisePrctile;
            if lastPlotted
                onPlotData();
            end
        end

        if any(cfg.noiseCorr ~= current_cfg.noiseCorr)
            current_cfg.noiseCorr = cfg.noiseCorr;
            if lastPlotted
                onPlotData();
            end
        end        

    end

    function updateNumericField(fieldName, val)
        cfg.(fieldName) = val;
        disp([fieldName, ' updated to ', num2str(val)]);
        
       if cfg.numChannels ~= numChannels
           channelValues = 2.^(0:10);
           channelValues(channelValues > cfg.numChannels) = [];
           channelValues = arrayfun(@num2str, channelValues, 'UniformOutput', false);
           yStepDropdown.Items = [{'full'} channelValues];
           plotGroupLineColor.Items = [{'0'} , channelValues];
           if ~isempty(channel_mapping)
               ySlider.Limits = [1 min(cfg.numChannels, length(channel_mapping))];
           else
               ySlider.Limits = [1 cfg.numChannels];
           end
           ticks = 5:5:100;
           interval_idx = nearest(ticks, floor(cfg.numChannels/10));
           ySlider.MajorTicks = 0:ticks(interval_idx):cfg.numChannels;
           numChannels = cfg.numChannels;
           drawnow;
           if lastPlotted
               % replotVal = false;
               onPlotData();
           end
           
        end

    end

    function updateTextOrNumeric(fieldName, rawVal)
        tmp = str2double(rawVal);
        if ~isnan(tmp)
            cfg.(fieldName) = tmp;
            disp([fieldName, ' updated to numeric: ', num2str(tmp)]);
        else
            cfg.(fieldName) = string(rawVal);
            disp([fieldName, ' updated to string: ', rawVal]);
        end
    end

    function setField(fieldName, rawVal)
        cfg.(fieldName) = rawVal;
        disp([fieldName, ' updated to ', mat2str(cfg.(fieldName))]);
    end

    % hp update callbacks
    function ddVal = convertToStringForDropdown(valIn)
        if isnumeric(valIn)
            ddVal = num2str(valIn);
        else
            ddVal = char(valIn);
        end
    end

    % time channel slider callbacks

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
        replotIfActive();
    end

    function updateYWindow(sVal, stepStr)
        if strcmp(stepStr,'full')
            yStepVal = min(cfg.numChannels, length(channel_mapping));             
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
        replotIfActive();
    end

    function replotIfActive()
        if lastPlotted
                onPlotData();                            
        end
    end

    %  Buttons in left panel (folder selection, load .mat)
    function selectFolder(folderNumber)
        
        figure(fig);  
        drawnow;
        oldStyle = fig.WindowStyle;
        fig.WindowStyle = 'modal';
        drawnow;

        folder = uigetdir();
        fig.WindowStyle = oldStyle;
        figure(fig);
        drawnow;

        if folder==0, return; end
        switch folderNumber
            case 1
                cfg.inputFolder = folder;
                lblFolder1.Text = folder;
                disp(['Folder 1 selected: ', folder]);
            case 2
                cfg.outputFolder = folder;
                lblFolder2.Text = folder;
                disp(['Folder 2 selected: ', folder]);
        end
    end

    function loadDataFile()
        figure(fig);
        oldStyle = fig.WindowStyle;
        fig.WindowStyle = 'modal';
        drawnow;
        [file, path] = uigetfile({...
            '*.dat;*.bin', 'All .dat & .bin files (*.dat, *.bin)'; ...
            '*.dat','Data File (*.dat)'; ...
            '*.bin','Binary File (*.bin)'; ...
            },...
            'Select Input Data File');
        fig.WindowStyle = oldStyle;
        figure(fig);
        drawnow;

        if isequal(file,0), return; end
        cfg.inputFolder = path;
        fullPath = fullfile(path, file);
        cfg.fullFilePath = fullPath;
        lblFolder1.Text  = fullPath;
        disp(['Input file: ', fullPath]);
        replotVal = 0;
    end

    function loadMatFile()
        figure(fig);
        oldStyle = fig.WindowStyle;
        fig.WindowStyle = 'modal';
        drawnow;
        [file, path] = uigetfile('*.mat','Select Channel-Mapping MAT-file');
        fig.WindowStyle = oldStyle;
        figure(fig);
        drawnow;

        if isequal(file,0), return; end
        fullPath = fullfile(path, file);
        cfg.channel_info = fullPath;
        lblMatFile.Text  = fullPath;
        disp(['Loaded MAT file: ', fullPath]);
        [channel_mapping, channel_locations, channel_inclusion] = load_channel_map(fullPath, cfg);
        replotVal = 0;
        setappdata(fig, 'channel_inclusion', channel_inclusion);
        ySlider.Limits = [1 min(cfg.numChannels, length(channel_mapping))];
        updateYWindow(yWindowStart, 'full');        
    end

    %  calling  sorting modules   

    function onRunAndSort()
        if ~isfield(cfg,'outputFolder')
            if isfield(cfg,'inputFolder')
            cfg.outputFolder = cfg.inputFolder;
            end
        end
        executionControl = 'Running';
        channel_inclusion = getappdata(fig, 'channel_inclusion');
        progressBar.Value = 0;
        progressLabel.Text = 'Sorting in progress. Extraction: 0%';
        progressBar.ScaleColors = [.9 .9 .9];
        progressBar.ScaleColorLimits = [0 100];
        drawnow;
        progressFcn = @(pct, msg) updateProgress(pct, msg);
        try
           kiaSort_main_extract_sample_data(cfg.fullFilePath, cfg.outputFolder, cfg, 'channel_mapping',channel_mapping,...
                'channel_inclusion', channel_inclusion, 'channel_locations', channel_locations,...
                'progressfcn', progressFcn);
            progressBar.Value = 100;
            progressLabel.Text = 'All samples are extracted.';
        catch ME
            progressLabel.Text = sprintf('Error: %s', ME.message);
            uialert(fig, ['Error: ' ME.message], 'Error');
            error('Extracting samples stopped.');
        end
        drawnow;

        progressBar.Value = 0;
        if cfg.parallelProcessing
            progressBar.ScaleColors = [.9 .9 .9];
            progressLabel.Text = 'Sorting samples in parallel. No progress report';
        else
            progressLabel.Text = 'Sorting in progress. Sorting samples: 0%';
            progressBar.ScaleColors = [.9 .9 .9];
            progressBar.ScaleColorLimits = [0 100];
        end
        drawnow;

        try
            kiaSort_main_sort_samples(cfg.outputFolder, cfg, hp,...
                'progressfcn', progressFcn);
            progressBar.Value = 100;
            progressLabel.Text = 'All samples are sorted';
        catch ME
            progressLabel.Text = sprintf('Error: %s', ME.message);
            uialert(fig, ['Error: ' ME.message], 'Error');
            error('Sorting samples stopped.');
        end
        drawnow;
        progressBar.Value = 0;
        progressLabel.Text = 'Sorting in progress. Sorting data: 0%';
        progressBar.ScaleColors = [.9 .9 .9];
        progressBar.ScaleColorLimits = [0 100];
        drawnow;

        try
            kiaSort_main_sortData(cfg.fullFilePath, cfg.outputFolder, cfg,...
                'progressfcn', progressFcn);
            progressBar.Value = 100;
            progressLabel.Text = 'All chunks are sorted';
            uialert(fig, 'Sorting Complete!', 'Success','Icon','success');
        catch ME
            progressLabel.Text = sprintf('Error: %s', ME.message);
            uialert(fig, ['Error: ' ME.message], 'Error');
            error('Sorting stopped.');

        end
        drawnow;
    end


    function onExtractSamples() 
        if ~isfield(cfg,'outputFolder')
            if isfield(cfg,'inputFolder')
                cfg.outputFolder = cfg.inputFolder;
            end
        end
        executionControl = 'Running';
        channel_inclusion = getappdata(fig, 'channel_inclusion');
        progressBar.Value = 0;
        progressLabel.Text = 'Extraction: 0%';
        progressBar.ScaleColors = [.9 .9 .9];    
        progressBar.ScaleColorLimits = [0 100];
        drawnow;
        progressFcn = @(pct, msg) updateProgress(pct, msg);
        try
            kiaSort_main_extract_sample_data(cfg.fullFilePath, cfg.outputFolder, cfg, 'channel_mapping',channel_mapping,...
                'channel_inclusion', channel_inclusion, 'channel_locations', channel_locations,...
                'progressfcn', progressFcn);
            progressBar.Value = 100;
            progressLabel.Text = 'All samples are extracted.';
            uialert(fig, 'Extraction Complete!', 'Success','Icon','success');
        catch ME
            progressLabel.Text = sprintf('Error: %s', ME.message);
            uialert(fig, ['Error: ' ME.message], 'Error');
            error('Extracting samples stopped.');
        end
    end

    function onSortGroups()
        if ~isfield(cfg,'outputFolder')
            if isfield(cfg,'inputFolder')
                cfg.outputFolder = cfg.inputFolder;
            end
        end
        executionControl = 'Running';
        progressBar.Value = 0;
        if cfg.parallelProcessing
            progressLabel.Text = 'Sorting samples in parallel. No progress report';
        else
        progressLabel.Text = 'Sorting samples: 0%';
        end      
        progressBar.ScaleColors = [.9 .9 .9];
        progressBar.ScaleColorLimits = [0 100];
        drawnow;
        progressFcn = @(pct, msg) updateProgress(pct, msg);
        try
            kiaSort_main_sort_samples(cfg.outputFolder, cfg, hp,...
                'progressfcn', progressFcn);
            progressBar.Value = 100;
            progressLabel.Text = 'All samples are sorted';
            uialert(fig, 'Sorting Samples Complete!', 'Success', 'Icon','success');
        catch ME
            progressLabel.Text = sprintf('Error: %s', ME.message);
            uialert(fig, ['Error: ' ME.message], 'Error');
            error('Sorting samples stopped.');
        end
        drawnow;
    end

    function onSortData()
        if ~isfield(cfg,'outputFolder')
            if isfield(cfg,'inputFolder')
            cfg.outputFolder = cfg.inputFolder;
            end
        end
        executionControl = 'Running';
        progressBar.Value = 0;
        progressLabel.Text = 'Sorting in progress: 0%';
        progressBar.ScaleColors = [.9 .9 .9];
        progressBar.ScaleColorLimits = [0 100];
        drawnow;
        progressFcn = @(pct, msg) updateProgress(pct, msg);
        try
            kiaSort_main_sortData(cfg.fullFilePath, cfg.outputFolder, cfg,...
                'progressfcn', progressFcn);
            progressBar.Value = 100;
            progressLabel.Text = 'All chunks are sorted';
            uialert(fig, 'Sorting Complete!', 'Success', 'Icon','success');
        catch ME
            progressLabel.Text = sprintf('Error: %s', ME.message);
            uialert(fig, ['Error: ' ME.message], 'Error');
            error('Sorting stopped.');
        end
        drawnow;
    end

    function updateProgress(pct, msg)
        
        if strcmp(executionControl, 'stopped')
            error('Execution stopped by user.');
        end
        while strcmp(executionControl, 'paused')
            pause(0.1);
            drawnow;
        end
        progressBar.ScaleColors =[.1 .8 .2; .9 .9 .9];
        pct = max(0, min(100, pct * 100));
        progressBar.Value = pct;
        progressLabel.Text = msg;
        if pct > 0 && pct < 100
            progressBar.ScaleColorLimits = [0 pct; pct 100];
        elseif pct == 100
            progressBar.ScaleColors = [.1 .8 .2];
            progressBar.ScaleColorLimits = [0 100];
        end

        drawnow;
    end

    % plotting callbacks

    function updatePlotType(newVal)
        currentPlotType = newVal;
        disp(['Plot type set to ', currentPlotType]);
        replotIfActive()
    end

    function updatePlotFilterType(newVal)
        plotFilterType = newVal;
        disp(['Plot filter type set to ', plotFilterType]);
        replotIfActive()
    end

    function updatePlotLineScale(newVal)
        plotLineScale = newVal;
        disp(['Plot line scale set to ', plotLineScale]);
        replotIfActive()
    end

    function updatePlotLineColor(newVal)
        plotLineColor = newVal;
        disp(['Plot group lines set to ', plotLineColor]);
        replotIfActive()
    end
    
    function updatePlotLineSat(sldVal)
        plotLineSat = sldVal;
        disp(['Plot line saturation set to ', num2str(plotLineSat)]);
        replotIfActive()
    end

    function updatePlotLineHue(sldVal)
        plotLineHue = sldVal;
        disp(['Plot line hue set to ', num2str(plotLineHue)]);
        replotIfActive()
    end

    function updatePlotLineValue(sldVal)
        plotLineValue = sldVal;
        disp(['Plot line hue set to ', num2str(plotLineValue)]);
        replotIfActive()
    end


    function onPlotData()
        disp('--- Plot filtered data function called ---');        
        if  ~replotVal || isempty(mappedData)
            if isfield(cfg, 'fullFilePath') 
                replotVal = true;                
                mappedData = map_input_file(cfg.fullFilePath, cfg);
                trialLength = size(mappedData.Data.data,2)/cfg.samplingFrequency;
                xSlider.Limits = [0 trialLength];
                roundTicks = 25:25:trialLength;
                nearestRoundTick = nearest(roundTicks,floor(trialLength/10));
                tickInterval = roundTicks(nearestRoundTick);
                xSlider.MajorTicks = 0:tickInterval:trialLength;
                updateXWindow(0,xStepDropdown.Value);
                disp(' Data loaded ');
            else
                disp('Data loading failed');
            end
        end

        disp(['Plot Type = ', currentPlotType]);


        cla(axConfig,'reset');
        lastPlotted = true;
        d_plotLineColor = str2double(plotLineColor);

        [sy, sx] = size(mappedData.Data.data);

        xDataRange = round(xWindowStart * cfg.samplingFrequency)+1 : min(round(xWindowEnd * cfg.samplingFrequency), sx);
        yDataRange = max(1, yWindowStart) : min(yWindowEnd, sy);
        
        
        inputData = double(mappedData.Data.data(:, xDataRange));
        

        if strcmp(plotFilterType,'filtered')        
        inputData = bandpass_filter_GUI(inputData', cfg.bandpass, cfg.samplingFrequency)';  
        end
        
        inputData = inputData(channel_mapping(yDataRange),:);
        
        if cfg.extremeNoise
            tic
            inputData = remove_res_shared_noise(inputData, cfg);
            toc

        end

        if cfg.denoising
            tic            
            inputData = remove_shared_noise(inputData, cfg);
            toc
        end


        if strcmp(currentPlotType,'images')
            normInImg = -inputData;
            imagesc(axConfig, xDataRange/cfg.samplingFrequency, yDataRange , normInImg,[-median(range(inputData')/2) median(range(inputData')/2)]);
            axis(axConfig,'xy')
            cmap = polar_colormap([plotLineHue, plotLineSat, plotLineValue]);
            colormap(axConfig,cmap);
            colorbar(axConfig,"eastoutside","Color",1-figColor,"Box","off",'TickDirection','out')
            title(axConfig,'Image plots','Color','white');
        else
            normIn = inputData./ max(abs(inputData),[],'all');                        
            plot(axConfig, xDataRange/cfg.samplingFrequency, (yDataRange' + normIn * plotLineScale)');                
            axConfig = plot_grouped_lines(numel(yDataRange), d_plotLineColor, axConfig, plotLineHue, plotLineSat, plotLineValue);

            title(axConfig,'Line plots','Color','white');
            xlabel(axConfig, 'Time (s)');
            ylabel(axConfig, 'Channel #');
        end
        xlabel(axConfig, 'Time (s)');
        ylabel(axConfig, 'Channel #');
        set(axConfig,'tickdir','out');
        axConfig.Color   = [0.1 0.1 0.1];
        axConfig.XColor  = [1 1 1];
        axConfig.YColor  = [1 1 1];
        axis(axConfig,[xWindowStart xWindowEnd yWindowStart-.5 yWindowEnd+.5]);

        disp(['Time range: ', num2str(xWindowStart), ' to ', num2str(xWindowEnd)]);
        disp(['Channel range: ', num2str(yDataRange(1)), ' to ', num2str(yDataRange(end))]);
    end

    function updateExportFormat(ddVal)
        exportType = ddVal;
    end


    function onExportFig()
        saveFigPath = fullfile(cfg.outputFolder, 'exported_figures');
        if ~exist("saveFigPath",'dir')
            mkdir(saveFigPath);
        end
        files = dir(fullfile(saveFigPath, 'exported_Fig_Signal_*.eps'));
        if isempty(files)
            nextNum = 1;
        else
            nums = cellfun(@(s) sscanf(s, 'exported_Fig_Signal_%d.eps'), {files.name});
            nextNum = max(nums) + 1;
        end
        filename = fullfile(saveFigPath, sprintf('exported_Fig_Signal_%d.eps', nextNum));
        exportgraphics(axConfig, filename,'BackgroundColor','white','ContentType', exportType)
    end

    function onStop()
        executionControl = 'stopped';
        disp('Execution stopped.');
    end

    function onPause()
        executionControl = 'paused';
        disp('Execution paused.');
    end

    function onResume()
        executionControl = 'running';
        disp('Execution resumed.');
    end


%% subtab Hyperparameters
    function updateHPDropdown(fieldName, newVal)
        tmp = str2double(newVal);
        if ~isnan(tmp)
            hp.(fieldName) = tmp;
        else
            hp.(fieldName) = string(newVal);
        end
        disp(['[HP] ',fieldName, ' updated to ', mat2str(hp.(fieldName))]);
    end

    function updateHPBoolean(fieldName, val)
        hp.(fieldName) = val;
        disp(['[HP] ',fieldName, ' updated to ', mat2str(val)]);
    end

    function setHpArrayField(fieldName, idx, val)
        arrOld = hp.(fieldName);
        arrOld(idx) = val;
        hp.(fieldName) = arrOld;
        disp(['[HP] ',fieldName, ' updated to ', mat2str(hp.(fieldName))]);
    end

    function updateHPTextOrNumeric(fieldName, rawVal)
        tmp = str2double(rawVal);
        if ~isnan(tmp)
            hp.(fieldName) = tmp;
            disp(['[HP] ',fieldName,' updated to numeric: ',num2str(tmp)]);
        else
            hp.(fieldName) = string(rawVal);
            disp(['[HP] ',fieldName,' updated to string: ',rawVal]);
        end
    end

    function setHPField(fieldName, rawVal)
        hp.(fieldName) = rawVal;
        disp(['[HP] ',fieldName,' updated to ',mat2str(hp.(fieldName))]);
    end


    %%   Channel Tab: callbacks

    function onLoadChannels()
        channel_inclusion = getappdata(fig, 'channel_inclusion');
        allChildren = gridTab2.Children;
        delete(allChildren(allChildren ~= btnLoadChannels));        
        if isfield(cfg,"outputFolder") || ~isempty(channel_inclusion)
            channel_inclusion = channel_selection_gui(cfg, channel_mapping, channel_inclusion, gridTab2, figColor, fig);
        else
            uialert(fig, ['Error: ' 'Path unspecified or no channel inclusion data loaded.'], 'Error');
        end
        channel_inclusion = getappdata(fig, 'channel_inclusion');
    end


    %%   Clusters Tab: callbacks

    function onPlotSingleChannel()
        if isfield(cfg,"outputFolder")
            delete(allchild(groupsPlotPanel));
            plot_single_channel_samples(cfg, groupsPlotPanel, figColor);
        else
            uialert(fig, ['Error: ' 'Path unspecified.'], 'Error');
        end
    end

    function onPlotAllGroups()
        if isfield(cfg,"outputFolder")
            delete(allchild(groupsPlotPanel));
            plot_all_groups_gui(cfg, groupsPlotPanel, figColor, fig);
        else
            uialert(fig, ['Error: ' 'Path unspecified.'], 'Error');
        end
    end

    function onClearGroupPlots()
        if isfield(cfg,"outputFolder")
            delete(allchild(groupsPlotPanel));
        else
            uialert(fig, ['Error: ' 'Path unspecified.'], 'Error');
        end
    end


    %% Curation Tab callbacks
    function onPlotSorting()
        allChildren = gridTab3.Children;
        delete(setdiff(allChildren, [btnLoadSorting, btnClear, btnAltRes]));
        if isfield(cfg,"outputFolder")
            kiaSort_curate_results(cfg, gridTab3, figColor, fig);
        else
            uialert(fig, ['Error: ' 'Path unspecified.'], 'Error');
        end
    end

    function onSelectAltResFolder()
        figure(fig);
        drawnow;
        oldStyle = fig.WindowStyle;
        fig.WindowStyle = 'modal';
        drawnow;

        folder = uigetdir();
        fig.WindowStyle = oldStyle;
        figure(fig);
        drawnow;

        if folder==0, return; end
        cfg.altResFolder = folder;
        disp(['Result Folder selected: ', folder]);
    end

    function onClearSorting()
        channel_inclusion = getappdata(fig, 'channel_inclusion');
        cfg.altResFolder = [];
        allChildren = gridTab3.Children;
        delete(setdiff(allChildren, [btnLoadSorting, btnClear, btnAltRes]));
    end

end 

    function applyColorScheme(h, parentBg)
        brightness = mean(parentBg);
        if brightness < 0.5
            textColor = [1 1 1];
        else
            textColor = [0 0 0];
        end
        if isprop(h,'BackgroundColor')
            h.BackgroundColor = parentBg;
        end
        if isprop(h,'FontColor')
            h.FontColor = textColor;
        end
        if isprop(h,'ForegroundColor')
            h.ForegroundColor = textColor;
        end
        if isprop(h,'Children')
            kids = h.Children;
            for k = 1:numel(kids)
                applyColorScheme(kids(k), parentBg);
            end
        end
    end

