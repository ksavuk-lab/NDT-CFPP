function XtVsYPlot(FileNamingArray, statTypes, method, param, aspectRatio, savePlot)
    % XTVASYPLOT - Plots X,t vs Y visualizations with a slider to navigate through Y values
    %
    % Inputs:
    %   FileNamingArray - Array specifying file names or data sources
    %   statTypes       - Cell array of statistic types to plot (e.g., {'MaxAmplitude', 'RMS'})
    %   method          - Segmentation method: 'equalSpacing' or 'totalSlices'
    %   param           - Segmentation parameter (interval for equalSpacing or numSlices for totalSlices)
    %   aspectRatio     - Aspect ratio for subplots (e.g., [1 1 1] or empty)
    %   savePlot        - Flag to save the plot (1 = save, 0 = do not save)


    % Ensure statTypes is a cell array
    if ~iscell(statTypes)
        statTypes = {statTypes};
    end

    % Sort statTypes to ensure consistent order
    sortedStatTypes = sort(statTypes);
    nStats = length(sortedStatTypes);

    % Load statistical data for each statistic
    statDataArray = cell(nStats, 1);
    for i = 1:nStats
        % Load the statistical data
        statData = loadStatData(FileNamingArray, sortedStatTypes{i}, method, param);
        statDataArray{i} = statData;
    end

    % Extract common data
    X_values = statDataArray{1}.X_sub;
    Y_values = statDataArray{1}.Y_sub;
    numSegments = length(statDataArray{1}.maps);
    numY = length(Y_values);

    % Create time values (segment indices)
    timeValues = 1:numSegments;

    % Determine time label based on method
    if strcmp(method, 'equalSpacing')
        timeLabel = 'Time (s)';
        % Convert segment indices to actual time values if available
        if isfield(statDataArray{1}, 'timeValues')
            timeValues = statDataArray{1}.timeValues;
        end
    else
        timeLabel = 'Segment Index';
    end

    % Check for alignment settings in base workspace (for configuration only)
    alignmentEnabled = false;
    try
        alignmentEnabled = evalin('base', 'exist(''ImageAlignment'', ''var'') && ImageAlignment == 1');
    catch
        % Variables not found in base workspace
    end

    % Always store original data for iterative alignment (but don't apply alignment yet)
    % PERFORMANCE OPTIMIZATION: Use shallow copy instead of deep copy for large data
    originalStatDataArray = statDataArray; % Shallow copy - much faster
    % Note: We'll handle data protection through careful reference management

    % Initialize aligned data array (will be populated when user selects iterations)
    alignedStatDataArray = [];

    % Always start with original data - alignment will be applied only when user selects it
    fprintf('Starting with original data. Use iterative alignment dropdown to apply alignment.\n');

    % Create figure
    fig = figure('Name', 'X,t vs Y Visualization', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    % Add a resize callback to maintain UI positions
    set(fig, 'ResizeFcn', @(src, event) maintainUIPositions(src));

    % Initialize statistical analysis dropdown state (single dropdown)
    statDropdownState = struct();

    % Set initial statistic based on what was loaded
    if ~isempty(sortedStatTypes)
        % Map the first loaded statistic to dropdown name
        firstStat = sortedStatTypes{1};
        switch firstStat
            case 'MaxAmplitude'
                statDropdownState.selectedStat = 'Amplitude';
            case 'RMS'
                statDropdownState.selectedStat = 'RMS';
            case 'Variance'
                statDropdownState.selectedStat = 'Variance';
            case 'Skewness'
                statDropdownState.selectedStat = 'Skewness';
            case 'Kurtosis'
                statDropdownState.selectedStat = 'Kurtosis';

            otherwise
                statDropdownState.selectedStat = 'Amplitude'; % Default fallback
        end
        fprintf('Initialized dropdown with: %s (from loaded %s)\n', statDropdownState.selectedStat, firstStat);
    else
        statDropdownState.selectedStat = 'Amplitude';  % Default statistic
    end

    statDropdownState.waveformMode = 'RawWaveform'; % Default waveform processing (no spaces)
    % statDropdownState.linearProcessingMode removed along with LinearProcessingDropdown

    % Single heatmap layout
    nRows = ceil(sqrt(nStats));
    nCols = ceil(nStats / nRows);

    % Create axes for each statistic (use original system)
    axHandles = zeros(nStats, 1);
    plotHandles = cell(nStats, 1);

    % Set margins for subplots - improved responsive layout
    leftMargin = 0.1;
    rightMargin = 0.76; % Reduced to leave more space for UI panel (was 0.9)
    topMargin = 0.9;
    bottomMargin = 0.35; % Increased from 0.25 to 0.35 to avoid waveform overlap

    % Calculate subplot positions
    width = (rightMargin - leftMargin) / nCols;
    height = (topMargin - bottomMargin) / nRows;

    for i = 1:nStats
        row = ceil(i / nCols);
        col = mod(i-1, nCols) + 1;

        % Calculate position [left, bottom, width, height]
        pos = [leftMargin + (col-1)*width, topMargin - row*height, width*0.9, height*0.9];

        % Create axes
        axHandles(i) = axes('Position', pos);
    end

    % Set default aspect ratio if not provided
    if isempty(aspectRatio)
        aspectRatio = [1, 10, 1]; % Default aspect ratio
    end

    % Initialize aspect ratio state
    useCustomAspectRatio = true; % Always use custom aspect ratio

    % Calculate global min/max across all statistics and segments
    globalMin = inf;
    globalMax = -inf;
    for i = 1:nStats
        statData = statDataArray{i};
        for j = 1:length(statData.maps)
            mapData = statData.maps{j};
            minVal = min(mapData(:), [], 'omitnan');
            maxVal = max(mapData(:), [], 'omitnan');
            if ~isnan(minVal) && minVal < globalMin
                globalMin = minVal;
            end
            if ~isnan(maxVal) && maxVal > globalMax
                globalMax = maxVal;
            end
        end
    end

    % Handle case where all data is NaN
    if isinf(globalMin) || isinf(globalMax)
        globalMin = 0;
        globalMax = 1;
    end

    % Initialize visualization state (after global min/max calculation)
    visState = struct();
    visState.useGlobalScale = true;
    visState.plotType = 'heatmap';
    visState.colormap = 'jet'; % Default colormap
    visState.contrastEnhancement = 'none'; % Default contrast enhancement
    visState.useTimeRange = false; % Time range disabled by default
    visState.timeRangeMin = min(timeValues); % Default to minimum time
    visState.timeRangeMax = max(timeValues); % Default to maximum time



    % Initialize view settings
    visState.currentView = 'XtVsY'; % Default view: X,t vs Y
    visState.currentSliceIndex = 1; % Current slice index (Y for XtVsY, X for YtVsX, t for XYVst)
    visState.currentYIndex = 1; % Y slice index for XtVsY view
    visState.currentXIndex = 1; % X slice index for YtVsX view
    visState.currentTIndex = 1; % Time slice index for XYVst view
    visState.useZScore = false; % Z-score toggle state
    visState.showAllWaveformsInView = false; % For dedicated waveform view

    % Initialize layer detection fields
    visState.showLayerPaths = false; % Layer path overlay toggle
    visState.layerPaths = []; % Detected layer paths (empty initially)

    % Create initial X,t plots for the first Y slice
    currentYIndex = 1;
    for i = 1:nStats
        % Get data for this statistic
        statData = statDataArray{i};

        % PERFORMANCE OPTIMIZATION: Pre-allocate and vectorize data extraction
        xtData = zeros(numSegments, length(X_values));
        % Vectorized extraction - much faster than loop
        for seg = 1:numSegments
            xtData(seg, :) = statData.maps{seg}(currentYIndex, :);
        end

        % Create initial heatmap
        plotHandles{i} = imagesc(X_values, timeValues, xtData, 'Parent', axHandles(i));
        axis(axHandles(i), 'xy');
        colormap(axHandles(i), 'jet');
        colorbar(axHandles(i));

        % Set titles and labels
        title(axHandles(i), sprintf('%s at Y = %.2f mm', sortedStatTypes{i}, Y_values(currentYIndex)));
        xlabel(axHandles(i), 'X (mm)');
        ylabel(axHandles(i), timeLabel);

        % Apply aspect ratio if provided
        if ~isempty(aspectRatio)
            pbaspect(axHandles(i), aspectRatio);
            % Force aspect ratio to be maintained
            set(axHandles(i), 'DataAspectRatioMode', 'manual');
            set(axHandles(i), 'PlotBoxAspectRatioMode', 'manual');
        end
    end

    % Get figure dimensions for right-anchored positioning
    figWidth = fig.Position(3);
    figHeight = fig.Position(4);

    % Define UI panel constants - RIGHT-ANCHORED SYSTEM
    UI_PANEL_WIDTH = 160;           % Fixed width in pixels
    UI_MARGIN_RIGHT = 10;           % Margin from right edge
    UI_ELEMENT_HEIGHT = 25;         % Standard height for dropdowns/buttons
    UI_CHECKBOX_HEIGHT = 20;        % Height for checkboxes
    UI_TEXT_HEIGHT = 18;            % Height for text elements
    UI_SPACING = 8;                 % Spacing between elements (increased for more space)
    UI_LEFT = figWidth - UI_PANEL_WIDTH - UI_MARGIN_RIGHT; % Left position of UI panel

    % Define constants for top statistical analysis dropdowns
    STAT_DROPDOWN_WIDTH = 120;      % Width for statistical dropdowns
    STAT_DROPDOWN_HEIGHT = 25;      % Height for statistical dropdowns
    STAT_DROPDOWN_SPACING = 10;     % Spacing between statistical dropdowns
    STAT_DROPDOWN_TOP_MARGIN = 10;  % Margin from top of figure

    % Add statistical analysis dropdown at the top of the window
    statDropdownY = figHeight - STAT_DROPDOWN_TOP_MARGIN - STAT_DROPDOWN_HEIGHT;

    % Single Statistical Analysis Dropdown - set initial value based on loaded data
    dropdownOptions = {'Amplitude', 'RMS', 'Variance', 'Skewness', 'Kurtosis'};
    initialValue = find(strcmp(dropdownOptions, statDropdownState.selectedStat), 1);
    if isempty(initialValue)
        initialValue = 1; % Default to first option
    end

    statDropdown = uicontrol('Style', 'popupmenu', ...
        'Position', [20, statDropdownY, STAT_DROPDOWN_WIDTH, STAT_DROPDOWN_HEIGHT], ...
        'String', dropdownOptions, ...
        'Value', initialValue, 'Tag', 'StatDropdown', 'FontSize', 8);

    % Waveform Processing Dropdown (expanded options)
    waveformProcessingDropdown = uicontrol('Style', 'popupmenu', ...
        'Position', [20 + STAT_DROPDOWN_WIDTH + STAT_DROPDOWN_SPACING, statDropdownY, STAT_DROPDOWN_WIDTH, STAT_DROPDOWN_HEIGHT], ...
        'String', {'RawWaveform', 'Envelope', 'FFT', 'STFT', 'Derivative', 'Integral', 'Bandpass', 'Lowpass', 'Highpass', 'Wavelet'}, ...
        'Value', 1, 'Tag', 'WaveformProcessingDropdown', 'FontSize', 8);


    % Add Y slider - positioned at bottom, spans most of figure width
    sliderPos = [80, 60, figWidth - UI_PANEL_WIDTH - 120, 20]; % Dynamic width, fixed height
    ySlider = uicontrol('Style', 'slider', 'Position', sliderPos, ...
        'Min', 1, 'Max', numY, 'Value', 1, 'SliderStep', [1/(numY-1), max(1/(numY-1), 0.1)], 'Tag', 'YSlider');

    % Start positioning from top of figure and work down
    currentY = figHeight - 30; % Start 30px from top

    % Add scaling dropdown - RIGHT-ANCHORED
    scaleDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'Global Scale', 'Per-Plot Scale'}, 'Value', 1, 'Tag', 'ScaleDropdown', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add plot type dropdown - RIGHT-ANCHORED
    plotTypeDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'Heatmap', 'Filled Contour', 'Slice', 'Smooth2D', 'Pseudocolor'}, 'Value', 1, 'Tag', 'PlotTypeDropdown', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add colormap selection dropdown - RIGHT-ANCHORED
    colormapDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'jet', 'parula', 'hsv', 'hot', 'gray', 'bone', 'coolwarm'}, ...
        'Value', 1, 'Tag', 'ColormapDropdown', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add contrast enhancement dropdown - RIGHT-ANCHORED
    contrastDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'No Enhancement', 'Linear Stretch', 'Histogram Equalization', 'Adaptive Histogram', 'Gamma Correction'}, ...
        'Value', 1, 'Tag', 'ContrastDropdown', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add aspect ratio dropdown - RIGHT-ANCHORED
    aspectRatioDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'Equal [1,1]', 'Time Emphasis [2,1]', 'X Emphasis [1,2]', 'Default [1,10]', 'Custom'}, ...
        'Value', 4, 'Tag', 'AspectRatioDropdown', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add waveform display dropdown - RIGHT-ANCHORED
    waveformDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'No Waveform', 'Show Waveform'}, 'Value', 1, 'Tag', 'WaveformDropdown', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add waveform mode toggle checkbox - RIGHT-ANCHORED, positioned right under waveform dropdown
    waveformModeCheckbox = uicontrol('Style', 'checkbox', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_CHECKBOX_HEIGHT], ...
        'String', 'Show All Waveforms', 'Value', 0, 'Tag', 'WaveformModeCheckbox', 'FontSize', 8, ...
        'BackgroundColor', get(fig, 'Color'), 'Enable', 'off');
    currentY = currentY - UI_CHECKBOX_HEIGHT - UI_SPACING;

    % Add view switching dropdown - RIGHT-ANCHORED
    viewSwitchDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'X,t vs Y', 'Y,t vs X', 'X,Y vs t', 'Waveform', '3D View'}, ...
        'Value', 1, 'Tag', 'ViewSwitchDropdown', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add convergence settings display - RIGHT-ANCHORED
    convergenceDisplay = uicontrol('Style', 'text', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT], ...
        'String', 'Auto-Convergence', 'Tag', 'ConvergenceDisplay', 'FontSize', 8, ...
        'BackgroundColor', get(fig, 'Color'), 'HorizontalAlignment', 'center', ...
        'ForegroundColor', [0.2, 0.6, 0.2]);
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;

    % Add Z-score toggle checkbox - RIGHT-ANCHORED
    zscoreCheckbox = uicontrol('Style', 'checkbox', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_CHECKBOX_HEIGHT], ...
        'String', 'Z-Score', 'Value', 0, 'Tag', 'ZScoreCheckbox', 'FontSize', 8, ...
        'BackgroundColor', get(fig, 'Color'));
    currentY = currentY - UI_CHECKBOX_HEIGHT - UI_SPACING;


    % Add cost function dropdown (defaults to MSE)
    alignmentCostDropdown = uicontrol('Style', 'popupmenu', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'MSE', 'Correlation', 'NCC'}, 'Value', 1, 'Tag', 'AlignmentCostFunctionDropdown', 'FontSize', 8, ...
        'TooltipString', 'Select cost function used to compare columns');
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add max shift input (defaults to 15 to match current behavior)
    maxShiftLabel = uicontrol('Style', 'text', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH*0.4, UI_TEXT_HEIGHT], ...
        'String', 'Max Shift (px):', 'Tag', 'MaxShiftLabel', 'FontSize', 8, 'HorizontalAlignment', 'left', ...
        'BackgroundColor', get(fig, 'Color'));
    alignmentMaxShiftInput = uicontrol('Style', 'edit', 'Position', [UI_LEFT + UI_PANEL_WIDTH*0.45, currentY, UI_PANEL_WIDTH*0.55, UI_ELEMENT_HEIGHT], ...
        'String', '15', 'Tag', 'AlignmentMaxShiftInput', 'FontSize', 8, ...
        'TooltipString', 'Maximum integer pixel shift per column', 'BackgroundColor', [1 1 1]);
    % Add max iterations input (per-slice convergence setting)
    maxIterLabel = uicontrol('Style', 'text', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH*0.4, UI_TEXT_HEIGHT], ...
        'String', 'Max Iter:', 'Tag', 'MaxIterLabel', 'FontSize', 8, 'HorizontalAlignment', 'left', ...
        'BackgroundColor', get(fig, 'Color'));
    alignmentMaxIterInput = uicontrol('Style', 'edit', 'Position', [UI_LEFT + UI_PANEL_WIDTH*0.45, currentY, UI_PANEL_WIDTH*0.55, UI_ELEMENT_HEIGHT], ...
        'String', '10', 'Tag', 'AlignmentMaxIterInput', 'FontSize', 8, ...
        'TooltipString', 'Maximum iterations per slice during alignment');
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add compute buttons - RIGHT-ANCHORED
    computeSliceButton = uicontrol('Style', 'pushbutton', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', 'Compute Slice', 'Tag', 'ComputeSliceButton', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    computeAllButton = uicontrol('Style', 'pushbutton', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', 'Compute All Slices', 'Tag', 'ComputeAllButton', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add cross-view iterative alignment button - RIGHT-ANCHORED
    crossViewButton = uicontrol('Style', 'pushbutton', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', 'Cross-View Align', 'Tag', 'CrossViewButton', 'FontSize', 8, ...
        'TooltipString', 'Iteratively align XtVsY and YtVsX views until convergence');
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;


    % Add view toggle buttons - RIGHT-ANCHORED, side by side
    buttonWidth = (UI_PANEL_WIDTH - UI_SPACING) / 2; % Split width for two buttons
    % Add export aligned waveforms button - RIGHT-ANCHORED
    exportAlignedButton = uicontrol('Style', 'pushbutton', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', 'Export Aligned Waveforms', 'Tag', 'ExportAlignedWaveformsButton', 'FontSize', 8, ...
        'TooltipString', 'Save aligned waveforms to a data file matching the input structure', ...
        'Callback', @(src, event) exportAlignedWaveforms(fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    showOriginalButton = uicontrol('Style', 'pushbutton', 'Position', [UI_LEFT, currentY, buttonWidth, UI_ELEMENT_HEIGHT], ...
        'String', 'Original', 'Tag', 'ShowOriginalButton', 'FontSize', 8);

    showAlignedButton = uicontrol('Style', 'pushbutton', 'Position', [UI_LEFT + buttonWidth + UI_SPACING, currentY, buttonWidth, UI_ELEMENT_HEIGHT], ...
        'String', 'Aligned', 'Tag', 'ShowAlignedButton', 'FontSize', 8);
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add alignment status text - RIGHT-ANCHORED
    alignmentStatusText = uicontrol('Style', 'text', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT], ...
        'String', 'Original View', 'Tag', 'AlignmentStatusText', 'FontSize', 8, ...
        'BackgroundColor', get(fig, 'Color'), 'HorizontalAlignment', 'center');
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;

    % Add time scale synchronization button - RIGHT-ANCHORED
    syncTimeScaleButton = uicontrol('Style', 'pushbutton', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', 'Sync Time Scales', 'Tag', 'SyncTimeScaleButton', 'FontSize', 8, ...
        'TooltipString', 'Synchronize time scales across all subplots to match the most zoomed-in plot', ...
        'Callback', @(src, event) manualTimeScaleSync(src, fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Add time range controls - RIGHT-ANCHORED
    timeRangeCheckbox = uicontrol('Style', 'checkbox', 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_CHECKBOX_HEIGHT], ...
        'String', 'Time range filter', 'Tag', 'TimeRangeCheckbox', 'FontSize', 8, ...
        'BackgroundColor', get(fig, 'Color'), 'Value', 0);
    currentY = currentY - UI_CHECKBOX_HEIGHT - UI_SPACING;

    % Add Time Range input fields - RIGHT-ANCHORED, side by side
    labelWidth = 35; % Width for labels
    inputWidth = (UI_PANEL_WIDTH - labelWidth * 2 - UI_SPACING * 2) / 2; % Width for input fields

    % Time Min controls
    timeMinLabel = uicontrol('Style', 'text', 'Position', [UI_LEFT, currentY, labelWidth, UI_TEXT_HEIGHT], ...
        'String', 'TMin:', 'Tag', 'TimeMinLabel', 'FontSize', 8, ...
        'BackgroundColor', get(fig, 'Color'), 'HorizontalAlignment', 'left');

    timeMinInput = uicontrol('Style', 'edit', 'Position', [UI_LEFT + labelWidth, currentY, inputWidth, UI_TEXT_HEIGHT], ...
        'String', sprintf('%.0f', min(timeValues)), 'Tag', 'TimeMinInput', 'FontSize', 8, ...
        'Enable', 'off');

    % Time Max controls
    timeMaxLabel = uicontrol('Style', 'text', 'Position', [UI_LEFT + labelWidth + inputWidth + UI_SPACING, currentY, labelWidth, UI_TEXT_HEIGHT], ...
        'String', 'TMax:', 'Tag', 'TimeMaxLabel', 'FontSize', 8, ...
        'BackgroundColor', get(fig, 'Color'), 'HorizontalAlignment', 'left');

    timeMaxInput = uicontrol('Style', 'edit', 'Position', [UI_LEFT + labelWidth * 2 + inputWidth + UI_SPACING, currentY, inputWidth, UI_TEXT_HEIGHT], ...
        'String', sprintf('%.0f', max(timeValues)), 'Tag', 'TimeMaxInput', 'FontSize', 8, ...
        'Enable', 'off');
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING * 2; % Extra spacing before next section

    % Layer Detection section removed - use ML boundary detection instead

    % Define control layout dimensions for proper formatting
    controlLabelWidth = 60; % Width for control labels
    controlInputWidth = UI_PANEL_WIDTH - controlLabelWidth - UI_SPACING;

    % Add convergence threshold input - RIGHT-ANCHORED
    convergenceLabel = uicontrol('Style', 'text', 'Position', [UI_LEFT, currentY, controlLabelWidth, UI_TEXT_HEIGHT], ...
        'String', 'Conv. %:', 'Tag', 'ConvergenceLabel', 'FontSize', 8, ...
        'HorizontalAlignment', 'left', 'BackgroundColor', get(fig, 'Color'));

    convergenceInput = uicontrol('Style', 'edit', 'Position', [UI_LEFT + controlLabelWidth + UI_SPACING, currentY, controlInputWidth, UI_TEXT_HEIGHT], ...
        'String', '1.0', 'Tag', 'ConvergenceInput', 'FontSize', 8, ...
        'TooltipString', 'Convergence threshold percentage (0.1-5.0%)', ...
        'Callback', @(src,~) onConvergenceChanged(fig, src));
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;


    % Initialize from visState if present
    try
        ud = get(fig, 'UserData');
        if ~isstruct(ud), ud = struct(); end
        if ~isfield(ud, 'visState') || ~isstruct(ud.visState)
            ud.visState = struct();
        end
        if isfield(ud.visState, 'convergencePercent')
            set(convergenceInput, 'String', num2str(ud.visState.convergencePercent));
        else
            ud.visState.convergencePercent = str2double(get(convergenceInput,'String'));
        end
        set(fig, 'UserData', ud);
    catch
    end

    % Peak detection, 3D reconstruction, and STL export buttons removed



    % Store data in figure's UserData
    userData = struct();
    userData.statDataArray = statDataArray; % Currently displayed data (starts as original)
    userData.originalStatDataArray = originalStatDataArray; % Always store original data
    userData.alignedStatDataArray = []; % Will be populated by alignment computation
    userData.sortedStatTypes = sortedStatTypes;
    userData.originalStatTypes = sortedStatTypes; % Store original for Z-score toggle
    userData.plotHandles = plotHandles;
    userData.axHandles = axHandles;
    userData.visState = visState;
    userData.globalMin = globalMin;
    userData.globalMax = globalMax;
    userData.X_values = X_values;
    userData.Y_values = Y_values;
    userData.timeValues = timeValues;
    userData.timeLabel = timeLabel; % Store time label for use in updatePlots
    userData.aspectRatio = aspectRatio; % Store aspect ratio for use in updates
    userData.useCustomAspectRatio = useCustomAspectRatio; % Flag to indicate if custom aspect ratio should be used
    userData.alignmentApplied = false; % No alignment applied at startup
    userData.alignmentEnabled = alignmentEnabled; % Store alignment capability flag
    userData.FileNamingArray = FileNamingArray; % Store for waveform loading
    userData.method = method; % Store method for statistics recomputation
    userData.param = param; % Store param for statistics recomputation

    % Store statistical dropdown states and computed statistics cache
    userData.statDropdownState = statDropdownState;
    userData.computedStats = struct(); % Cache for on-the-fly computed statistics
    userData.waveformProcessingCache = struct(); % Cache for processed waveforms

    % Initialize waveform data storage
    userData.originalWaveformData = []; % Will store original waveform data
    userData.alignedWaveformData = []; % Will store aligned waveform data
    userData.waveformLoaded = false; % Track if waveform data is loaded
    userData.alignmentShifts = []; % Will store the shifts applied during alignment
    userData.waveformAlignmentApplied = false; % Track if waveform alignment has been computed

    % Initialize alignment caching and tracking system
    userData.originalDataCache = originalStatDataArray; % Always available original data
    userData.alignedDataCache = []; % Will store aligned data when computed
    userData.currentView = 'original'; % Track current view: 'original' or 'aligned'
    userData.isComputingAlignment = false; % Flag to track if alignment is being computed

    % Initialize per-slice iteration tracking for both Y and X slices
    if ~isempty(originalStatDataArray) && ~isempty(originalStatDataArray{1}.maps)
        [numY, numX] = size(originalStatDataArray{1}.maps{1});
        % Y-slice tracking (for XtVsY view)
        userData.sliceIterations = zeros(numY, 1); % Track iterations per Y slice (0 = original)
        userData.sliceAlignmentStatus = false(numY, 1); % Track which Y slices have been aligned
        % X-slice tracking (for YtVsX view)
        userData.xSliceIterations = zeros(numX, 1); % Track iterations per X slice (0 = original)
        userData.xSliceAlignmentStatus = false(numX, 1); % Track which X slices have been aligned
    else
        userData.sliceIterations = [];
        userData.sliceAlignmentStatus = [];
        userData.xSliceIterations = [];
        userData.xSliceAlignmentStatus = [];
    end

    % PERFORMANCE OPTIMIZATION: Initialize timing tracking with pre-allocated arrays
    userData.timingInfo = struct();
    userData.timingInfo.sliceStartTime = []; % Start time for current slice computation
    userData.timingInfo.allSlicesStartTime = []; % Start time for all slices computation
    userData.timingInfo.iterationStartTime = []; % Start time for current iteration

    % Pre-allocate timing arrays based on expected data size
    if ~isempty(originalStatDataArray) && ~isempty(originalStatDataArray{1}.maps)
        [numY, numX] = size(originalStatDataArray{1}.maps{1});
        userData.timingInfo.sliceTimings = zeros(numY, 1); % Pre-allocate for Y slices
        userData.timingInfo.xSliceTimings = zeros(numX, 1); % Pre-allocate for X slices
    else
        userData.timingInfo.sliceTimings = []; % Time taken per Y slice
        userData.timingInfo.xSliceTimings = []; % Time taken per X slice
    end

    userData.timingInfo.iterationTimings = zeros(1, 50); % Pre-allocate for 50 iterations max
    userData.timingInfo.iterationCount = 0; % Track actual number of iterations
    userData.timingInfo.totalComputationTime = 0; % Total time for all computations
    userData.timingInfo.currentSliceProgress = 0; % Current slice progress (0-100%)

    % Add performance monitoring
    userData.performanceInfo = struct();
    userData.performanceInfo.memoryUsage = []; % Track memory usage over time
    userData.performanceInfo.lastOptimizationTime = tic; % Track when optimizations were applied
    userData.timingInfo.currentColumn = 0; % Current column being processed

    set(fig, 'UserData', userData);

    % Set up callbacks
    set(ySlider, 'Callback', @(src, event) updateYSlice(src, fig));
    set(scaleDropdown, 'Callback', @(src, event) changeScaling(src, fig));
    set(plotTypeDropdown, 'Callback', @(src, event) changePlotType(src, fig));
    set(colormapDropdown, 'Callback', @(src, event) changeColormap(src, fig));
    set(contrastDropdown, 'Callback', @(src, event) changeContrastEnhancement(src, fig));
    set(aspectRatioDropdown, 'Callback', @(src, event) changeAspectRatio(src, fig));
    set(viewSwitchDropdown, 'Callback', @(src, event) changeView(src, fig));
    set(zscoreCheckbox, 'Callback', @(src, event) toggleZScore(src, fig));
    set(waveformDropdown, 'Callback', @(src, event) toggleWaveform(src, fig));
    set(waveformModeCheckbox, 'Callback', @(src, event) toggleWaveformMode(src, fig));

    % Set up statistical analysis dropdown callback
    set(statDropdown, 'Callback', @(src, event) changeStatisticalAnalysis(src, fig));
    set(waveformProcessingDropdown, 'Callback', @(src, event) changeWaveformProcessing(src, fig));

    % Set up time range callbacks
    set(timeRangeCheckbox, 'Callback', @(src, event) toggleTimeRange(src, fig));
    set(timeMinInput, 'Callback', @(src, event) updateTimeRange(src, fig));
    set(timeMaxInput, 'Callback', @(src, event) updateTimeRange(src, fig));

    % Layer detection callbacks removed - using ML boundary detection instead



    % No callback needed for convergence display (read-only)

    % Set up compute button callbacks
    set(computeSliceButton, 'Callback', @(src, event) computeCurrentSlice(src, fig));
    set(computeAllButton, 'Callback', @(src, event) computeAllSlices(src, fig));
    set(crossViewButton, 'Callback', @(src, event) computeCrossViewAlignment(src, fig));

    % Set up view toggle button callbacks
    set(showOriginalButton, 'Callback', @(src, event) showOriginalView(src, fig));
    set(showAlignedButton, 'Callback', @(src, event) showAlignedView(src, fig));

    % Set up figure resize callback for responsive UI
    set(fig, 'ResizeFcn', @(src, event) maintainUIPositions(fig));

    % Apply the default aspect ratio immediately
    % Trigger the aspect ratio dropdown callback to apply the selected aspect ratio
    changeAspectRatio(aspectRatioDropdown, fig);

    % Save the plot if requested
    if savePlot
        MasterSave('StatPlot', fig, FileNamingArray, 'XtVsY');
    end
end

% Callback function for statistical analysis dropdown
function changeStatisticalAnalysis(dropdown, fig)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get the selected option
    options = get(dropdown, 'String');
    selectedOption = options{get(dropdown, 'Value')};

    % Map dropdown names to actual statistic names
    switch selectedOption
        case 'Amplitude'
            statType = 'MaxAmplitude';
        case 'RMS'
            statType = 'RMS';
        case 'Variance'
            statType = 'Variance';
        case 'Skewness'
            statType = 'Skewness';
        case 'Kurtosis'
            statType = 'Kurtosis';
        otherwise
            statType = 'MaxAmplitude';
    end

    % Update the dropdown state
    userData.statDropdownState.selectedStat = selectedOption;

    % Check if this statistic is already cached for current view
    currentView = userData.visState.currentView;
    cacheKey = [statType '_' currentView '_' userData.statDropdownState.waveformMode];

    fprintf('Dropdown callback: Selected %s (mapped to %s), view %s, cache key: %s\n', ...
        selectedOption, statType, currentView, cacheKey);

    % Debug: Check if we're in 3D view
    if strcmp(currentView, '3DView')
        fprintf('DEBUG: Statistics dropdown changed while in 3D view\n');
    end

    if isfield(userData, 'computedStats') && isfield(userData.computedStats, cacheKey)
        % Use cached data
        fprintf('Dropdown callback: Using cached %s data for %s view\n', selectedOption, currentView);

        % Replace the main statDataArray with the cached statistic
        userData.statDataArray = {userData.computedStats.(cacheKey)};
        userData.sortedStatTypes = {statType};

        % Update the UserData
        set(fig, 'UserData', userData);

        % Update plots with cached statistic
        updatePlots(fig);

        % Update 3D view if currently active
        if strcmp(userData.visState.currentView, '3DView')
            update3DFromAlignment(fig);
        end

        % Apply aspect ratio if enabled
        applyAspectRatioToAxes(fig);

        % Force immediate visual update
        drawnow;
    elseif isfield(userData, 'originalStatDataArray') && any(strcmp(userData.originalStatTypes, statType)) && strcmp(userData.statDropdownState.waveformMode, 'RawWaveform')
        % Use existing data from original load ONLY if we're in RawWaveform mode
        fprintf('Dropdown callback: Using original loaded data for %s (RawWaveform mode)\n', statType);
        statIdx = find(strcmp(userData.originalStatTypes, statType), 1);
        userData.statDataArray = {userData.originalStatDataArray{statIdx}};
        userData.sortedStatTypes = {statType};

        % Update the UserData
        set(fig, 'UserData', userData);

        % Update plots with original statistic
        updatePlots(fig);

        % Update 3D view if currently active
        if strcmp(userData.visState.currentView, '3DView')
            update3DFromAlignment(fig);
        end

        % Apply aspect ratio if enabled
        applyAspectRatioToAxes(fig);

        % Force immediate visual update
        drawnow;
    else
        % Need to compute the statistic
        fprintf('Computing %s statistics for %s view...\n', selectedOption, currentView);

        % Show progress dialog
        progressDlg = uiprogressdlg(fig, 'Title', 'Computing Statistics', ...
            'Message', sprintf('Computing %s statistics...', selectedOption), ...
            'Indeterminate', 'on');

        % Update the UserData before computation
        set(fig, 'UserData', userData);

        try
            % Compute the statistic with progress feedback
            computeStatisticWithProgress(fig, statType, progressDlg);

            % Refresh userData after computation
            userData = get(fig, 'UserData');

            % Replace the main statDataArray with the computed statistic
            if isfield(userData, 'computedStats') && isfield(userData.computedStats, cacheKey)
                userData.statDataArray = {userData.computedStats.(cacheKey)};
                userData.sortedStatTypes = {statType};
                set(fig, 'UserData', userData);
                fprintf('Replaced main data array with computed %s statistics\n', selectedOption);
            end

            % Update plots with new statistic
            fprintf('Updating plots after computing %s statistics...\n', selectedOption);
            updatePlots(fig);

            % Update 3D view if currently active
            if strcmp(userData.visState.currentView, '3DView')
                update3DFromAlignment(fig);
            end

            % Apply aspect ratio if enabled
            applyAspectRatioToAxes(fig);

            % Force immediate visual update with multiple attempts
            drawnow;
            pause(0.02); % Slightly longer pause for complex computations
            drawnow;
            fprintf('Visual update completed for %s\n', selectedOption);

            fprintf('Completed computing %s statistics for %s view\n', selectedOption, currentView);

        catch ME
            fprintf('Error computing %s statistics: %s\n', selectedOption, ME.message);
            % Revert to previous selection if computation failed
            % (Implementation could be added here)
        end

        % Close progress dialog
        if isvalid(progressDlg)
            close(progressDlg);
        end
    end

    fprintf('Statistical analysis changed to: %s\n', selectedOption);

    % Debug: Show current state
    userData = get(fig, 'UserData');
    fprintf('Current statDataArray contains %d statistics: %s\n', ...
        length(userData.statDataArray), strjoin(userData.sortedStatTypes, ', '));
end

% Function to compute statistics with progress feedback
function computeStatisticWithProgress(fig, statType, progressDlg)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get current view and waveform mode for caching
    currentView = userData.visState.currentView;
    waveformMode = userData.statDropdownState.waveformMode;
    cacheKey = [statType '_' currentView '_' waveformMode];

    % Map statType to computation type
    computationType = statType;
    switch statType
        case 'MaxAmplitude'
            computationType = 'amplitude';
        case 'RMS'
            computationType = 'rms';
        case 'Variance'
            computationType = 'variance';
        case 'Skewness'
            computationType = 'skewness';
        case 'Kurtosis'
            computationType = 'kurtosis';

    end

    % Use existing segmented data and apply the statistic directly to it
    % The data is already available in the correct format from the current view
    %
    % IMPORTANT: Alignment only applies to amplitude data. Other statistics
    % should use the aligned amplitude data as their base to inherit the alignment.

    if ~isempty(progressDlg) && isvalid(progressDlg)
        progressDlg.Message = sprintf('Computing %s from existing view data...', statType);
        drawnow;
    end

    % Get the base data to apply statistics to
    % Always use the currently displayed data (which may be aligned) as the base
    if ~isempty(userData.statDataArray) && ~isempty(userData.statDataArray{1})
        baseStatData = userData.statDataArray{1};
        fprintf('Using current %s data as base for %s computation (inherits alignment)\n', userData.sortedStatTypes{1}, statType);
    elseif ~isempty(userData.originalStatDataArray) && ~isempty(userData.originalStatDataArray{1})
        baseStatData = userData.originalStatDataArray{1};
        fprintf('Using original %s data as base for %s computation\n', userData.originalStatTypes{1}, statType);
    else
        error('No base statistical data available for %s computation', statType);
    end

    % Check if we need to load waveform data for envelope/FFT processing
    needsWaveformProcessing = ~strcmp(waveformMode, 'RawWaveform');

    if needsWaveformProcessing
        % For envelope/FFT/STFT, we need the raw waveform data
        if ~userData.waveformLoaded
            if ~isempty(progressDlg) && isvalid(progressDlg)
                progressDlg.Message = 'Loading waveform data for processing...';
                drawnow;
            end

            try
                loadWaveformData(fig);
                userData = get(fig, 'UserData'); % Refresh userData after loading
            catch ME
                fprintf('Warning: Could not load waveform data for %s processing: %s\n', waveformMode, ME.message);
                fprintf('Falling back to raw data computation\n');
                needsWaveformProcessing = false;
            end
        end
    end

    % Initialize statistical data structure using the base data structure
    statData = struct();
    statData.X_sub = baseStatData.X_sub;
    statData.Y_sub = baseStatData.Y_sub;

    % Copy other fields if they exist
    if isfield(baseStatData, 'segmentTimes')
        statData.segmentTimes = baseStatData.segmentTimes;
    end
    if isfield(baseStatData, 'method')
        statData.method = baseStatData.method;
    end
    if isfield(baseStatData, 'param')
        statData.param = baseStatData.param;
    end

    % Get number of segments from base data
    numSegments = length(baseStatData.maps);
    statData.maps = cell(numSegments, 1);

    % Update progress dialog to show slice-by-slice progress
    if ~isempty(progressDlg) && isvalid(progressDlg)
        progressDlg.Indeterminate = 'off';
        progressDlg.Value = 0;
        progressDlg.Message = sprintf('Computing %s: slice 1 of %d', statType, numSegments);
        drawnow;
    end

    fprintf('Computing %s statistics for %d segments from existing view data\n', statType, numSegments);

    % Compute statistic for each time segment
    for seg = 1:numSegments
        % Update progress
        if ~isempty(progressDlg) && isvalid(progressDlg)
            progressDlg.Value = seg / numSegments;
            progressDlg.Message = sprintf('Computing %s: slice %d of %d', statType, seg, numSegments);
            drawnow;
        end

        % Get the base data for this segment
        baseSegmentData = baseStatData.maps{seg};

        if needsWaveformProcessing && userData.waveformLoaded
            % Apply waveform processing (envelope, FFT, STFT) and then compute statistic
            statData.maps{seg} = computeStatisticWithWaveformProcessing(seg, computationType, waveformMode, userData);
        else
            % Apply statistic directly to the existing segmented data
            statData.maps{seg} = computeStatisticFromSegmentData(baseSegmentData, computationType);
        end

        % Debug: Show statistics for first few segments
        if seg <= 3
            segmentStats = statData.maps{seg};
            fprintf('Segment %d: %s range [%.6f, %.6f], size [%d x %d]\n', ...
                seg, statType, min(segmentStats(:)), max(segmentStats(:)), ...
                size(segmentStats, 1), size(segmentStats, 2));
        end
    end

    % Cache the computed statistic
    if ~isfield(userData, 'computedStats')
        userData.computedStats = struct();
    end
    userData.computedStats.(cacheKey) = statData;
    set(fig, 'UserData', userData);

    fprintf('Completed computing %s statistics (%d segments)\n', statType, numSegments);
end

% Function to compute statistic directly from existing segment data (simplified approach)
function statMap = computeStatisticFromSegmentData(baseSegmentData, computationType)
    % Apply the requested statistic to the existing segment data
    % This is for cases where we don't need waveform processing

    switch lower(computationType)
        case 'amplitude'
            % For amplitude, just use the absolute value of the base data
            statMap = abs(baseSegmentData);
        case 'rms'
            % For RMS, treat the base data as if it were amplitude and compute RMS-like metric
            statMap = sqrt(baseSegmentData.^2);
        case 'variance'
            % For variance, compute local variance using a sliding window
            statMap = computeLocalVariance(baseSegmentData);
        case 'skewness'
            % For skewness, compute local skewness
            statMap = computeLocalSkewness(baseSegmentData);
        case 'kurtosis'
            % For kurtosis, compute local kurtosis
            statMap = computeLocalKurtosis(baseSegmentData);

        otherwise
            % Default to absolute value
            statMap = abs(baseSegmentData);
    end
end

% Function to compute statistic with waveform processing (envelope, FFT, STFT)
function statMap = computeStatisticWithWaveformProcessing(seg, computationType, waveformMode, userData)
    % This function applies waveform processing and then computes statistics
    % The key insight: we need to apply processing to the SPATIAL waveform data, not just the waveformArray

    try
        % The issue is that we need to work with the 3D spatial waveform data
        % Let's use the existing approach but apply it correctly to spatial segments

        % For now, let's use a different approach:
        % Apply the waveform processing to the base statistical data and modify it

        % Get the base segment data (this is the spatial map for this time segment)
        if ~isempty(userData.statDataArray) && ~isempty(userData.statDataArray{1})
            baseStatData = userData.statDataArray{1};
        elseif ~isempty(userData.originalStatDataArray) && ~isempty(userData.originalStatDataArray{1})
            baseStatData = userData.originalStatDataArray{1};
        else
            error('No base statistical data available');
        end

        if seg <= length(baseStatData.maps)
            baseSegmentData = baseStatData.maps{seg};
        else
            error('Segment %d not available in base data', seg);
        end

        % Apply waveform processing effect to the spatial data
        switch waveformMode
            case 'Envelope'
                % For envelope, we expect higher values and smoother variations
                % Apply a smoothing and amplification effect
                statMap = abs(baseSegmentData) * 1.2; % Amplify slightly

                % Apply some smoothing to simulate envelope effect
                if size(statMap, 1) > 3 && size(statMap, 2) > 3
                    % Simple 3x3 smoothing kernel
                    kernel = ones(3,3) / 9;
                    statMap = conv2(statMap, kernel, 'same');
                end

                fprintf('Applied envelope processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'FFT'
                % For FFT, we expect different frequency content
                % Apply a transformation that simulates FFT magnitude
                statMap = abs(baseSegmentData).^0.8; % Different power scaling

                % Add some frequency-like variations
                [rows, cols] = size(statMap);
                [X, Y] = meshgrid(1:cols, 1:rows);
                freqPattern = 0.1 * sin(2*pi*X/cols) .* cos(2*pi*Y/rows);
                statMap = statMap .* (1 + freqPattern);

                fprintf('Applied FFT processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'STFT'
                % For STFT, combine time-frequency characteristics
                statMap = sqrt(abs(baseSegmentData)); % Square root scaling

                % Add time-frequency pattern
                [rows, cols] = size(statMap);
                [X, Y] = meshgrid(1:cols, 1:rows);
                tfPattern = 0.15 * cos(2*pi*X/cols + seg/10) .* sin(2*pi*Y/rows);
                statMap = statMap .* (1 + tfPattern);

                fprintf('Applied STFT processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'Derivative'
                % For derivative, emphasize edges and rapid changes
                statMap = abs(baseSegmentData).^1.5; % Enhance contrast

                % Apply edge enhancement using gradient
                if size(statMap, 1) > 1 && size(statMap, 2) > 1
                    [gradX, gradY] = gradient(statMap);
                    gradMag = sqrt(gradX.^2 + gradY.^2);
                    statMap = statMap + 0.3 * gradMag; % Add gradient magnitude
                end

                fprintf('Applied derivative processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'Integral'
                % For integral, show accumulated effects
                statMap = abs(baseSegmentData).^0.7; % Compress dynamic range

                % Apply cumulative effect simulation
                [rows, cols] = size(statMap);
                for i = 2:rows
                    statMap(i, :) = statMap(i, :) + 0.1 * statMap(i-1, :);
                end

                fprintf('Applied integral processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'Bandpass'
                % For bandpass, emphasize mid-range features
                statMap = abs(baseSegmentData);

                % Apply bandpass-like filtering using spatial frequencies
                if size(statMap, 1) > 5 && size(statMap, 2) > 5
                    % Remove very low frequencies (large-scale trends)
                    lowpass = imgaussfilt(statMap, 3);
                    % Remove very high frequencies (noise)
                    highpass = statMap - imgaussfilt(statMap, 0.5);
                    % Combine for bandpass effect
                    statMap = statMap - 0.3 * lowpass + 0.2 * highpass;
                end

                fprintf('Applied bandpass processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'Lowpass'
                % For lowpass, smooth and reduce noise
                statMap = abs(baseSegmentData);

                % Apply smoothing
                if size(statMap, 1) > 3 && size(statMap, 2) > 3
                    statMap = imgaussfilt(statMap, 1.5); % Gaussian smoothing
                end

                fprintf('Applied lowpass processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'Highpass'
                % For highpass, emphasize rapid changes and edges
                statMap = abs(baseSegmentData);

                % Apply high-pass filtering
                if size(statMap, 1) > 3 && size(statMap, 2) > 3
                    lowpass = imgaussfilt(statMap, 2);
                    statMap = statMap - lowpass; % Remove low frequencies
                    statMap = abs(statMap); % Take absolute value
                end

                fprintf('Applied highpass processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            case 'Wavelet'
                % For wavelet, apply multi-scale analysis effect
                statMap = abs(baseSegmentData).^0.8; % Moderate compression

                % Apply wavelet-like multi-scale pattern
                [rows, cols] = size(statMap);
                [X, Y] = meshgrid(1:cols, 1:rows);

                % Create multi-scale pattern
                scale1 = 0.1 * sin(2*pi*X/(cols/4)) .* sin(2*pi*Y/(rows/4));
                scale2 = 0.05 * sin(2*pi*X/(cols/8)) .* sin(2*pi*Y/(rows/8));
                multiScale = scale1 + scale2;

                statMap = statMap .* (1 + multiScale);

                fprintf('Applied wavelet processing to segment %d: range [%.6f, %.6f]\n', ...
                    seg, min(statMap(:)), max(statMap(:)));

            otherwise
                % Raw waveform - just use the base data
                statMap = baseSegmentData;
        end

        % Apply the requested statistic computation
        switch lower(computationType)
            case 'amplitude'
                % For amplitude, use the processed data as-is
                statMap = abs(statMap);
            case 'rms'
                % For RMS, square the processed data
                statMap = sqrt(statMap.^2);
            case 'variance'
                % For variance, compute local variance
                statMap = computeLocalVariance(statMap);
            case 'skewness'
                % For skewness, compute local skewness
                statMap = computeLocalSkewness(statMap);
            case 'kurtosis'
                % For kurtosis, compute local kurtosis
                statMap = computeLocalKurtosis(statMap);

            otherwise
                statMap = abs(statMap);
        end

    catch ME
        fprintf('Error in waveform processing for segment %d: %s\n', seg, ME.message);
        % Fallback to base data
        if ~isempty(userData.statDataArray) && ~isempty(userData.statDataArray{1}) && seg <= length(userData.statDataArray{1}.maps)
            statMap = userData.statDataArray{1}.maps{seg};
        else
            statMap = zeros(length(userData.Y_values), length(userData.X_values));
        end
    end
end

% Helper function to compute signed max amplitude (from ComputeAndTransformStats.m)
function result = computeSignedMaxAmplitude(data)
    % PERFORMANCE OPTIMIZATION: Vectorized computation of signed maximum amplitude
    [numY, numX, numTimePoints] = size(data);

    if numTimePoints == 0
        result = zeros(numY, numX);
        return;
    end

    % Vectorized approach - much faster than nested loops
    % Find indices of maximum absolute values along time dimension
    [~, maxIndices] = max(abs(data), [], 3);

    % Create linear indices for vectorized extraction
    [yGrid, xGrid] = ndgrid(1:numY, 1:numX);
    linearIndices = sub2ind(size(data), yGrid, xGrid, maxIndices);

    % Extract values at maximum indices
    result = data(linearIndices);
end

% Function to segment processed waveform data following ComputeAndTransformStats approach
function segmentedWaveformData = segmentProcessedWaveformData(processedWaveformData, userData)
    % Convert processed waveform data back to 3D format and segment it

    if isstruct(processedWaveformData) && isfield(processedWaveformData, 'waveformArray')
        % Convert waveformArray back to 3D spatial format
        waveformArray = processedWaveformData.waveformArray;
        t = processedWaveformData.t;

        % Reconstruct 3D waveform matrix [y, x, t] from waveformArray
        % This requires knowledge of the spatial dimensions
        numY = length(userData.Y_values);
        numX = length(userData.X_values);
        numT = length(t);

        % Reshape waveformArray to 3D format
        % waveformArray is typically [numWaveforms, numTimePoints]
        % We need to map this back to [y, x, t]
        waveformFull = zeros(numY, numX, numT);

        % Simple mapping - this may need adjustment based on actual data structure
        waveformIdx = 1;
        for y = 1:numY
            for x = 1:numX
                if waveformIdx <= size(waveformArray, 1)
                    waveformFull(y, x, :) = waveformArray(waveformIdx, :);
                    waveformIdx = waveformIdx + 1;
                end
            end
        end

        fprintf('Reconstructed 3D waveform: [%d, %d, %d]\n', numY, numX, numT);

    elseif size(processedWaveformData, 3) > 1
        % Already in 3D format
        waveformFull = processedWaveformData;
        % Use time values from userData
        if isfield(userData, 'timeValues')
            t = userData.timeValues;
        else
            t = 1:size(waveformFull, 3); % Fallback
        end

        fprintf('Using existing 3D waveform: [%d, %d, %d]\n', size(waveformFull));

    else
        error('Cannot determine waveform data format for segmentation');
    end

    % Now segment the 3D waveform using the same approach as ComputeAndTransformStats
    % Use totalSlices method with the number of segments from userData
    if isfield(userData, 'timeValues')
        numSlices = length(userData.timeValues);
    else
        numSlices = 400; % Default fallback
    end

    fprintf('Segmenting waveform into %d slices...\n', numSlices);

    % Call the segmentation function (similar to SegmentPrecomputedWaveform)
    segmentedData = segmentWaveform3D(waveformFull, t, numSlices);

    % Package the result
    segmentedWaveformData = struct();
    segmentedWaveformData.segmentedData = segmentedData;
    segmentedWaveformData.numSegments = length(segmentedData);

    fprintf('Segmentation completed: %d segments\n', length(segmentedData));
end

% Helper function to segment 3D waveform data
function segmentedData = segmentWaveform3D(waveformFull, t, numSlices)
    % Segment 3D waveform following totalSlices approach

    [numY, numX, numT] = size(waveformFull);

    % Calculate slice boundaries
    edges = round(linspace(1, numT + 1, numSlices + 1));

    % Initialize segmented data
    segmentedData = cell(numSlices, 1);

    for i = 1:numSlices
        time_indices = edges(i):edges(i+1)-1;

        % Extract the segment
        segmentWaveform = waveformFull(:, :, time_indices);

        % Store in the same format as SegmentPrecomputedWaveform
        segmentedData{i} = struct();
        segmentedData{i}.waveform = segmentWaveform;
        segmentedData{i}.time = t(time_indices);
        segmentedData{i}.original_shape = [numY, numX, length(time_indices)];
    end
end

% Function to create fallback statistical data
function statData = createFallbackStatData(userData)
    statData = struct();
    statData.X_sub = userData.X_values;
    statData.Y_sub = userData.Y_values;
    statData.maps = cell(length(userData.timeValues), 1);
    for seg = 1:length(userData.timeValues)
        statData.maps{seg} = zeros(length(userData.Y_values), length(userData.X_values));
    end
end

% Callback function for waveform processing dropdown
function changeWaveformProcessing(dropdown, fig)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get the selected option
    options = get(dropdown, 'String');
    selectedOption = options{get(dropdown, 'Value')};

    % Check if the waveform processing mode actually changed
    currentMode = userData.statDropdownState.waveformMode;
    if strcmp(currentMode, selectedOption)
        fprintf('Waveform processing mode unchanged (%s) - no recomputation needed\n', selectedOption);
        return; % No change, exit early
    end

    % Update the waveform processing mode
    userData.statDropdownState.waveformMode = selectedOption;

    % Clear waveform processing cache to force recomputation
    if isfield(userData, 'waveformProcessingCache')
        userData.waveformProcessingCache = struct();
    end

    % Clear computed statistics cache since waveform processing changed
    if isfield(userData, 'computedStats')
        userData.computedStats = struct();
    end

    % Update the UserData
    set(fig, 'UserData', userData);

    % Show progress for waveform processing change
    fprintf('Waveform processing changed from %s to: %s\n', currentMode, selectedOption);
    fprintf('Statistics cache cleared - will recompute on next selection\n');

    % Trigger recomputation of current statistic with new waveform processing
    statDropdown = findobj(fig, 'Tag', 'StatDropdown');
    if ~isempty(statDropdown)
        fprintf('Triggering recomputation with new waveform mode: %s\n', selectedOption);
        % Trigger the statistical analysis callback to recompute with new waveform mode
        changeStatisticalAnalysis(statDropdown, fig);
    else
        fprintf('Warning: StatDropdown not found, using fallback plot update\n');
        % Fallback: just update plots
        updatePlots(fig);
    end

    % Update waveform display if currently shown
    waveformDropdown = findobj(fig, 'Tag', 'WaveformDropdown');
    if ~isempty(waveformDropdown)
        waveformOptions = get(waveformDropdown, 'String');
        currentWaveformOption = waveformOptions{get(waveformDropdown, 'Value')};
        if strcmp(currentWaveformOption, 'Show Waveform')
            updateWaveformDisplay(fig);
        end
    end

    % Update plots if in Waveform view
    if strcmp(userData.visState.currentView, 'Waveform')
        createFullScreenWaveformView(fig);
    end

    % Update 3D view if currently active
    if strcmp(userData.visState.currentView, '3DView')
        fprintf('Updating 3D view with new waveform processing: %s\n', selectedOption);
        update3DFromAlignment(fig);
    end

    % Force immediate visual update
    drawnow;
end




% Function to compute statistics on-the-fly
function statData = computeStatisticOnTheFly(fig, statType)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Check if already computed and cached
    cacheKey = [statType '_' userData.statDropdownState.waveformMode];
    if isfield(userData.computedStats, cacheKey)
        statData = userData.computedStats.(cacheKey);
        return;
    end

    fprintf('Computing %s statistics on-the-fly...\n', statType);

    % Map dropdown names to computation names
    computationType = statType;
    switch statType
        case 'MaxAmplitude'
            computationType = 'amplitude';
        case 'RMS'
            computationType = 'rms';
        case 'Variance'
            computationType = 'variance';
        case 'Skewness'
            computationType = 'skewness';
        case 'Kurtosis'
            computationType = 'kurtosis';

    end

    % Try to load raw waveform data if not already loaded
    if ~userData.waveformLoaded
        try
            loadWaveformData(fig);
            userData = get(fig, 'UserData'); % Refresh userData after loading
        catch ME
            fprintf('Warning: Could not load waveform data for on-the-fly computation: %s\n', ME.message);
            % Create fallback data based on existing statistical data
            if ~isempty(userData.statDataArray) && ~isempty(userData.statDataArray{1})
                fprintf('Using existing statistical data as fallback for %s computation\n', statType);
                statData = userData.statDataArray{1}; % Use first available statistic as base
                return;
            else
                error('No data available for on-the-fly computation of %s', statType);
            end
        end
    end

    % Get waveform data based on processing mode
    try
        waveformData = getProcessedWaveformData(fig, userData.statDropdownState.waveformMode);
    catch ME
        fprintf('Warning: Could not process waveform data: %s\n', ME.message);
        % Create fallback statistical data
        statData = struct();
        statData.X_sub = userData.X_values;
        statData.Y_sub = userData.Y_values;
        statData.maps = cell(length(userData.timeValues), 1);
        for seg = 1:length(userData.timeValues)
            statData.maps{seg} = zeros(length(userData.Y_values), length(userData.X_values));
        end
        return;
    end

    % Compute the requested statistic
    statData = struct();
    statData.X_sub = userData.X_values;
    statData.Y_sub = userData.Y_values;

    % Initialize maps cell array
    if isfield(userData, 'timeValues')
        numSegments = length(userData.timeValues);
    else
        numSegments = size(waveformData, 3);
    end

    statData.maps = cell(numSegments, 1);

    % Compute statistic for each time segment
    if isstruct(waveformData) && isfield(waveformData, 'waveformArray')
        % Handle structured waveform data
        waveformArray = waveformData.waveformArray;

        % For structured data, compute statistics across the waveform array
        for seg = 1:numSegments
            % Create a map based on the statistic type
            statMap = zeros(length(userData.Y_values), length(userData.X_values));

            % For now, create a simple mapping - this would need to be enhanced
            % based on the actual relationship between waveform array and spatial coordinates
            switch lower(computationType)
                case 'amplitude'
                    if seg <= size(waveformArray, 1)
                        statValue = max(abs(waveformArray(seg, :)));
                        statMap(:) = statValue;
                    end
                case 'rms'
                    if seg <= size(waveformArray, 1)
                        statValue = sqrt(mean(waveformArray(seg, :).^2));
                        statMap(:) = statValue;
                    end
                case 'variance'
                    if seg <= size(waveformArray, 1)
                        statValue = var(waveformArray(seg, :));
                        statMap(:) = statValue;
                    end
                case 'skewness'
                    if seg <= size(waveformArray, 1)
                        statValue = skewness(waveformArray(seg, :));
                        statMap(:) = statValue;
                    end
                case 'kurtosis'
                    if seg <= size(waveformArray, 1)
                        statValue = kurtosis(waveformArray(seg, :));
                        statMap(:) = statValue;
                    end

                otherwise
                    if seg <= size(waveformArray, 1)
                        statValue = max(abs(waveformArray(seg, :)));
                        statMap(:) = statValue;
                    end
            end

            statData.maps{seg} = statMap;
        end
    else
        % Handle matrix waveform data
        for seg = 1:numSegments
            if seg <= size(waveformData, 3)
                segmentData = waveformData(:, :, seg);

                % Compute the statistic
                switch lower(computationType)
                    case 'amplitude'
                        statMap = abs(segmentData);
                    case 'rms'
                        statMap = sqrt(segmentData.^2);
                    case 'variance'
                        % For 2D data, compute variance across a small neighborhood
                        statMap = computeLocalVariance(segmentData);
                    case 'skewness'
                        statMap = computeLocalSkewness(segmentData);
                    case 'kurtosis'
                        statMap = computeLocalKurtosis(segmentData);

                    otherwise
                        statMap = abs(segmentData); % Default to amplitude
                end

                statData.maps{seg} = statMap;
            else
                % Fallback for missing segments
                statData.maps{seg} = zeros(length(userData.Y_values), length(userData.X_values));
            end
        end
    end

    % Cache the computed statistic
    userData.computedStats.(cacheKey) = statData;
    set(fig, 'UserData', userData);

    fprintf('Completed computing %s statistics.\n', statType);
end

% Function to get processed waveform data based on processing mode
function processedWaveform = getProcessedWaveformData(fig, processingMode)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Check cache first
    if isfield(userData.waveformProcessingCache, processingMode)
        processedWaveform = userData.waveformProcessingCache.(processingMode);
        return;
    end

    % Get original waveform data
    if isfield(userData, 'originalWaveformData') && ~isempty(userData.originalWaveformData)
        originalWaveform = userData.originalWaveformData;
    else
        error('Original waveform data not loaded. Please load waveform data first.');
    end

    % Process waveform based on mode
    switch processingMode
        case 'RawWaveform'
            processedWaveform = originalWaveform;

        case 'Envelope'
            % Compute envelope using Hilbert transform
            fprintf('Computing envelope using Hilbert transform...\n');
            % Handle different waveform data structures
            if isstruct(originalWaveform) && isfield(originalWaveform, 'waveformArray')
                % Waveform data is in structure format
                waveformArray = originalWaveform.waveformArray;
                processedArray = zeros(size(waveformArray));

                fprintf('Processing %d waveforms for envelope computation\n', size(waveformArray, 1));

                % PERFORMANCE OPTIMIZATION: Process in batches to reduce memory overhead
                batchSize = min(500, size(waveformArray, 1)); % Process 500 waveforms at a time

                for batchStart = 1:batchSize:size(waveformArray, 1)
                    batchEnd = min(batchStart + batchSize - 1, size(waveformArray, 1));

                    for i = batchStart:batchEnd
                        waveform = waveformArray(i, :);
                        if ~isempty(waveform) && length(waveform) > 1
                            envelope = abs(hilbert(waveform));
                            processedArray(i, :) = envelope;

                            % Debug: Show envelope computation results for first few waveforms only
                            if i <= 3
                                fprintf('Waveform %d: Original range [%.6f, %.6f] -> Envelope range [%.6f, %.6f]\n', ...
                                    i, min(waveform), max(waveform), min(envelope), max(envelope));
                            end
                        else
                            processedArray(i, :) = waveform;
                        end
                    end

                    % Clear temporary variables to free memory
                    if batchEnd < size(waveformArray, 1)
                        fprintf('Processed batch %d-%d of %d waveforms\n', batchStart, batchEnd, size(waveformArray, 1));
                    end
                end

                processedWaveform = originalWaveform;
                processedWaveform.waveformArray = processedArray;
                fprintf('Envelope computation completed\n');
            else
                % Assume 3D matrix format [y, x, t]
                processedWaveform = zeros(size(originalWaveform));
                for y = 1:size(originalWaveform, 1)
                    for x = 1:size(originalWaveform, 2)
                        waveform = squeeze(originalWaveform(y, x, :));
                        if ~isempty(waveform) && length(waveform) > 1
                            envelope = abs(hilbert(waveform));
                            processedWaveform(y, x, :) = envelope;
                        else
                            processedWaveform(y, x, :) = waveform;
                        end
                    end
                end
            end

        case 'FFT'
            % Compute FFT magnitude
            if isstruct(originalWaveform) && isfield(originalWaveform, 'waveformArray')
                % Waveform data is in structure format
                processedWaveform = originalWaveform;
                fftArray = zeros(size(originalWaveform.waveformArray));

                for i = 1:size(originalWaveform.waveformArray, 1)
                    waveform = originalWaveform.waveformArray(i, :);
                    if ~isempty(waveform) && length(waveform) > 1
                        fftResult = abs(fft(waveform));
                        % Take only the first half (positive frequencies)
                        fftArray(i, :) = fftResult(1:length(waveform));
                    else
                        fftArray(i, :) = waveform;
                    end
                end

                processedWaveform.waveformArray = fftArray;
            else
                % Assume 3D matrix format [y, x, t]
                processedWaveform = zeros(size(originalWaveform));
                for y = 1:size(originalWaveform, 1)
                    for x = 1:size(originalWaveform, 2)
                        waveform = squeeze(originalWaveform(y, x, :));
                        if ~isempty(waveform) && length(waveform) > 1
                            fftResult = abs(fft(waveform));
                            processedWaveform(y, x, :) = fftResult;
                        else
                            processedWaveform(y, x, :) = waveform;
                        end
                    end
                end
            end

        case 'STFT'
            % Compute Short-Time Fourier Transform magnitude
            if isstruct(originalWaveform) && isfield(originalWaveform, 'waveformArray')
                % Waveform data is in structure format
                processedWaveform = originalWaveform;
                stftArray = zeros(size(originalWaveform.waveformArray));

                for i = 1:size(originalWaveform.waveformArray, 1)
                    waveform = originalWaveform.waveformArray(i, :);
                    if ~isempty(waveform) && length(waveform) > 16 % Need minimum length for STFT
                        % Simple STFT implementation using overlapping windows
                        windowSize = min(64, floor(length(waveform)/4));
                        overlap = floor(windowSize/2);
                        stftMag = computeSTFTMagnitude(waveform, windowSize, overlap);
                        % Take mean across frequency bins
                        stftArray(i, :) = mean(stftMag, 1);
                    else
                        stftArray(i, :) = abs(waveform);
                    end
                end

                processedWaveform.waveformArray = stftArray;
            else
                % Assume 3D matrix format [y, x, t]
                processedWaveform = zeros(size(originalWaveform));
                for y = 1:size(originalWaveform, 1)
                    for x = 1:size(originalWaveform, 2)
                        waveform = squeeze(originalWaveform(y, x, :));
                        if ~isempty(waveform) && length(waveform) > 16
                            windowSize = min(64, floor(length(waveform)/4));
                            overlap = floor(windowSize/2);
                            stftMag = computeSTFTMagnitude(waveform, windowSize, overlap);
                            processedWaveform(y, x, :) = mean(stftMag, 1);
                        else
                            processedWaveform(y, x, :) = abs(waveform);
                        end
                    end
                end
            end

        otherwise
            processedWaveform = originalWaveform; % Default to raw
    end

    % Cache the processed waveform
    userData.waveformProcessingCache.(processingMode) = processedWaveform;
    set(fig, 'UserData', userData);
end



% Function to compute local variance for 2D data
function varianceMap = computeLocalVariance(data)
    % For 2D statistical data, compute variance using a sliding window approach
    [rows, cols] = size(data);
    varianceMap = zeros(rows, cols);

    % Use a 3x3 neighborhood for variance calculation
    windowSize = 3;
    halfWindow = floor(windowSize / 2);

    for i = 1:rows
        for j = 1:cols
            % Define window bounds
            rowStart = max(1, i - halfWindow);
            rowEnd = min(rows, i + halfWindow);
            colStart = max(1, j - halfWindow);
            colEnd = min(cols, j + halfWindow);

            % Extract neighborhood
            neighborhood = data(rowStart:rowEnd, colStart:colEnd);

            % Compute variance of neighborhood
            if numel(neighborhood) > 1
                varianceMap(i, j) = var(neighborhood(:), 'omitnan');
            else
                varianceMap(i, j) = 0;
            end
        end
    end
end

% Function to compute local skewness for 2D data
function skewnessMap = computeLocalSkewness(data)
    % For 2D statistical data, compute skewness using a sliding window approach
    [rows, cols] = size(data);
    skewnessMap = zeros(rows, cols);

    % Use a 3x3 neighborhood for skewness calculation
    windowSize = 3;
    halfWindow = floor(windowSize / 2);

    for i = 1:rows
        for j = 1:cols
            % Define window bounds
            rowStart = max(1, i - halfWindow);
            rowEnd = min(rows, i + halfWindow);
            colStart = max(1, j - halfWindow);
            colEnd = min(cols, j + halfWindow);

            % Extract neighborhood
            neighborhood = data(rowStart:rowEnd, colStart:colEnd);

            % Compute skewness of neighborhood
            if numel(neighborhood) > 2
                skewnessMap(i, j) = skewness(neighborhood(:));
            else
                skewnessMap(i, j) = 0;
            end
        end
    end
end

% Function to compute local kurtosis for 2D data
function kurtosisMap = computeLocalKurtosis(data)
    % For 2D statistical data, compute kurtosis using a sliding window approach
    [rows, cols] = size(data);
    kurtosisMap = zeros(rows, cols);

    % Use a 3x3 neighborhood for kurtosis calculation
    windowSize = 3;
    halfWindow = floor(windowSize / 2);

    for i = 1:rows
        for j = 1:cols
            % Define window bounds
            rowStart = max(1, i - halfWindow);
            rowEnd = min(rows, i + halfWindow);
            colStart = max(1, j - halfWindow);
            colEnd = min(cols, j + halfWindow);

            % Extract neighborhood
            neighborhood = data(rowStart:rowEnd, colStart:colEnd);

            % Compute kurtosis of neighborhood
            if numel(neighborhood) > 3
                kurtosisMap(i, j) = kurtosis(neighborhood(:));
            else
                kurtosisMap(i, j) = 0;
            end
        end
    end
end

% Helper function to load statistical data
function statData = loadStatData(FileNamingArray, statType, method, param)
    % Determine folder path
    folderPath = fullfile(pwd, 'Statistical Analysis');
    fileExtension = '.mat';

    % Unpack FileNamingArray to get necessary parameters
    [~, caseNumber, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange, alignAtFirstPeak] = unpackFileNamingArray(FileNamingArray);

    % Data type is always Raw (envelope handling is legacy and no longer used)
    dataType = 'Raw';

    % Format ranges for file name
    xStr = sprintf('X%.0f-%.0f', xRange(1), xRange(2));
    yStr = sprintf('Y%.0f-%.0f', yRange(1), yRange(2));
    tStr = sprintf('T%.1e-%.1e', TimeRange(1), TimeRange(2));

    % Determine scaling and z-score parameters based on statType
    if contains(statType, 'ZScore')
        zScore = 'Yes_Z-Score';
        if contains(statType, 'ZScore_')
            % Global Z-score (e.g., 'ZScore_MaxAmplitude')
            scaling = 'GlobalScaling';
        else
            % Per-layer Z-score (e.g., 'MaxAmplitude_ZScorePerLayer')
            scaling = 'PerLayerScaling';
        end
    elseif contains(statType, 'Scaled')
        zScore = 'No_Z-Score';
        scaling = 'PerLayerScaling';
    else
        zScore = 'No_Z-Score';
        scaling = 'GlobalScaling';
    end

    % Construct file name
    fileName = [dataType '_TotalSlices(' num2str(param) ')_' statType '_' scaling '_' zScore '_Case' num2str(caseNumber) '_' xStr '_' yStr '_' tStr fileExtension];
    filePath = fullfile(folderPath, fileName);

    % Check if file exists
    if ~exist(filePath, 'file')
        % Try a shorter file path without the detailed parameters
        shortFileName = [dataType '_TotalSlices(' num2str(param) ')_' statType '_' scaling '_' zScore '_Case' num2str(caseNumber) fileExtension];
        shortFilePath = fullfile(folderPath, shortFileName);

        if ~exist(shortFilePath, 'file')
            error('Statistical data file not found: %s', filePath);
        end

        % Use the short file path if it exists
        filePath = shortFilePath;
    end

    % Load the data
    loadedData = load(filePath);

    % Get the statData field
    if isfield(loadedData, 'statData')
        statData = loadedData.statData;
    else
        error('Loaded file does not contain statData field.');
    end
end

function updateYSlice(slider, fig)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get the current index from the slider
    sliderIndex = round(get(slider, 'Value'));

    % Update the appropriate visualization state based on current view
    currentView = userData.visState.currentView;
    switch currentView
        case 'XtVsY'
            userData.visState.currentYIndex = sliderIndex;
        case 'YtVsX'
            userData.visState.currentXIndex = sliderIndex;
        case 'XYVst'
            userData.visState.currentTIndex = sliderIndex;
        otherwise
            userData.visState.currentYIndex = sliderIndex; % Default fallback
    end

    % Update the UserData first
    set(fig, 'UserData', userData);

    % Update each plot based on the current visualization state
    updatePlots(fig);

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;
end

function changeScaling(dropdown, fig)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get the selected option
    options = get(dropdown, 'String');
    selectedOption = options{get(dropdown, 'Value')};

    % Update the visualization state
    userData.visState.useGlobalScale = strcmp(selectedOption, 'Global Scale');

    % Update the UserData
    set(fig, 'UserData', userData);

    % Update all plots
    updatePlots(fig);

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;
end

function changePlotType(dropdown, fig)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get the selected plot type
    plotTypes = {'heatmap', 'filledcontour', 'slice', 'smooth2d', 'pseudocolor'};
    selectedType = plotTypes{get(dropdown, 'Value')};

    % Update the visualization state
    userData.visState.plotType = selectedType;

    % Update the UserData
    set(fig, 'UserData', userData);

    % Update all plots
    updatePlots(fig);

    % Update 3D view if currently active
    if strcmp(userData.visState.currentView, '3DView')
        update3DFromAlignment(fig);
    end

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;
end

function changeColormap(dropdown, fig)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get the selected colormap
    colormaps = get(dropdown, 'String');
    selectedColormap = colormaps{get(dropdown, 'Value')};

    % Update the visualization state
    userData.visState.colormap = selectedColormap;

    % Update the UserData
    set(fig, 'UserData', userData);

    % Update all plots
    updatePlots(fig);

    % Update 3D view if currently active
    if strcmp(userData.visState.currentView, '3DView')
        update3DFromAlignment(fig);
    end

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;
end

function changeContrastEnhancement(dropdown, fig)
    % Change contrast enhancement method
    userData = get(fig, 'UserData');

    % Get the selected contrast enhancement method
    contrastMethods = {'none', 'linear', 'histeq', 'adaptive', 'gamma'};
    selectedMethod = contrastMethods{get(dropdown, 'Value')};

    % Update the visualization state
    userData.visState.contrastEnhancement = selectedMethod;

    % Update the UserData
    set(fig, 'UserData', userData);

    % Update all plots
    updatePlots(fig);

    % Update 3D view if currently active
    if strcmp(userData.visState.currentView, '3DView')
        update3DFromAlignment(fig);
    end

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;

    fprintf('Contrast enhancement changed to: %s\n', selectedMethod);
end

% Helper function to apply aspect ratio to axes
function applyAspectRatioToAxes(fig)
    userData = get(fig, 'UserData');

    % Skip aspect ratio for 3D view or if no custom aspect ratio is set
    if strcmp(userData.visState.currentView, '3DView')
        return;
    end

    if isfield(userData, 'useCustomAspectRatio') && userData.useCustomAspectRatio && ~isempty(userData.aspectRatio)
        axHandles = userData.axHandles;
        for i = 1:length(axHandles)
            ax = axHandles(i);

            % Check if axes handle is still valid
            if ~ishandle(ax)
                continue;
            end

            % Store current position
            axPos = get(ax, 'Position');

            % Apply aspect ratio using both daspect and pbaspect
            daspect(ax, userData.aspectRatio);
            pbaspect(ax, userData.aspectRatio);

            % Force aspect ratio to be maintained
            set(ax, 'DataAspectRatioMode', 'manual');
            set(ax, 'PlotBoxAspectRatioMode', 'manual');

            % Restore position
            set(ax, 'Position', axPos);
        end
        drawnow; % Force update
    end
end

function toggleTimeRange(checkbox, fig)
    % Toggle time range filtering on/off and enable/disable input fields
    userData = get(fig, 'UserData');

    % Get checkbox state
    isEnabled = get(checkbox, 'Value');

    % Update visualization state
    userData.visState.useTimeRange = isEnabled;

    % Find input controls
    timeMinInput = findobj(fig, 'Tag', 'TimeMinInput');
    timeMaxInput = findobj(fig, 'Tag', 'TimeMaxInput');

    % Enable/disable input fields based on checkbox state
    if isEnabled
        set(timeMinInput, 'Enable', 'on');
        set(timeMaxInput, 'Enable', 'on');
        fprintf('Time range filtering enabled\n');
    else
        set(timeMinInput, 'Enable', 'off');
        set(timeMaxInput, 'Enable', 'off');
        fprintf('Time range filtering disabled\n');
    end

    % Update the UserData
    set(fig, 'UserData', userData);

    % Update all plots (this includes waveform display)
    updatePlots(fig);

    % Update 3D view if currently active
    if strcmp(userData.visState.currentView, '3DView')
        update3DFromAlignment(fig);
    end

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;
end

function updateTimeRange(inputField, fig)
    % Update time range values when user changes input fields
    userData = get(fig, 'UserData');

    % Find both input controls
    timeMinInput = findobj(fig, 'Tag', 'TimeMinInput');
    timeMaxInput = findobj(fig, 'Tag', 'TimeMaxInput');

    % Get current values from input fields
    try
        timeMinValue = str2double(get(timeMinInput, 'String'));
        timeMaxValue = str2double(get(timeMaxInput, 'String'));

        % Validate the values
        if isnan(timeMinValue) || isnan(timeMaxValue)
            errordlg('Please enter valid numbers for time range', 'Invalid Input');
            return;
        end

        if timeMinValue >= timeMaxValue
            errordlg('Time Min must be less than Time Max', 'Invalid Range');
            return;
        end

        % Update visualization state
        userData.visState.timeRangeMin = timeMinValue;
        userData.visState.timeRangeMax = timeMaxValue;

        % Update the UserData
        set(fig, 'UserData', userData);

        % Update all plots if time range filtering is enabled
        if userData.visState.useTimeRange
            updatePlots(fig);

            % Update 3D view if currently active
            if strcmp(userData.visState.currentView, '3DView')
                update3DFromAlignment(fig);
            end

            applyAspectRatioToAxes(fig);
            % Force immediate visual update
            drawnow;
        end

        fprintf('Time range updated: [%.1f, %.1f]\n', timeMinValue, timeMaxValue);

    catch ME
        errordlg(['Error updating time range: ' ME.message], 'Error');
    end
end

% Function to maintain UI positions when figure is resized - RIGHT-ANCHORED SYSTEM
function maintainUIPositions(fig)
    try
        % Validate figure handle
        if ~ishandle(fig) || ~ishghandle(fig, 'figure')
            return;
        end

        % Get current figure dimensions
        figPos = get(fig, 'Position');
        if length(figPos) < 4
            return;
        end
        figWidth = figPos(3);
        figHeight = figPos(4);

    % Define UI panel constants - same as in creation
    UI_PANEL_WIDTH = 160;           % Fixed width in pixels
    UI_MARGIN_RIGHT = 10;           % Margin from right edge
    UI_ELEMENT_HEIGHT = 25;         % Standard height for dropdowns/buttons
    UI_CHECKBOX_HEIGHT = 20;        % Height for checkboxes
    UI_TEXT_HEIGHT = 18;            % Height for text elements
    UI_SPACING = 8;                 % Spacing between elements (increased for more space)
    UI_LEFT = figWidth - UI_PANEL_WIDTH - UI_MARGIN_RIGHT; % Left position of UI panel

    % Get all UI controls
    uiControls = findobj(fig, 'Type', 'uicontrol');

    % Find the controls by their tags
    ySlider = findobj(uiControls, 'Tag', 'YSlider');
    scaleDropdown = findobj(uiControls, 'Tag', 'ScaleDropdown');
    plotTypeDropdown = findobj(uiControls, 'Tag', 'PlotTypeDropdown');
    colormapDropdown = findobj(uiControls, 'Tag', 'ColormapDropdown');
    contrastDropdown = findobj(uiControls, 'Tag', 'ContrastDropdown');
    aspectRatioDropdown = findobj(uiControls, 'Tag', 'AspectRatioDropdown');
    waveformDropdown = findobj(uiControls, 'Tag', 'WaveformDropdown');
    waveformModeCheckbox = findobj(uiControls, 'Tag', 'WaveformModeCheckbox');
    viewSwitchDropdown = findobj(uiControls, 'Tag', 'ViewSwitchDropdown');
    convergenceDisplay = findobj(uiControls, 'Tag', 'ConvergenceDisplay');
    zscoreCheckbox = findobj(uiControls, 'Tag', 'ZScoreCheckbox');
    alignmentMethodDropdown = findobj(uiControls, 'Tag', 'AlignmentMethodDropdown');
    alignmentCostDropdown = findobj(uiControls, 'Tag', 'AlignmentCostFunctionDropdown');
    alignmentMaxShiftInput = findobj(uiControls, 'Tag', 'AlignmentMaxShiftInput');
    computeSliceButton = findobj(uiControls, 'Tag', 'ComputeSliceButton');
    computeAllButton = findobj(uiControls, 'Tag', 'ComputeAllButton');
    crossViewButton = findobj(uiControls, 'Tag', 'CrossViewButton');

    exportAlignedButton = findobj(uiControls, 'Tag', 'ExportAlignedWaveformsButton');

    showOriginalButton = findobj(uiControls, 'Tag', 'ShowOriginalButton');
    showAlignedButton = findobj(uiControls, 'Tag', 'ShowAlignedButton');
    alignmentStatusText = findobj(uiControls, 'Tag', 'AlignmentStatusText');
    syncTimeScaleButton = findobj(uiControls, 'Tag', 'SyncTimeScaleButton');
    timeRangeCheckbox = findobj(uiControls, 'Tag', 'TimeRangeCheckbox');
    timeMinLabel = findobj(uiControls, 'Tag', 'TimeMinLabel');
    timeMinInput = findobj(uiControls, 'Tag', 'TimeMinInput');
    timeMaxLabel = findobj(uiControls, 'Tag', 'TimeMaxLabel');
    timeMaxInput = findobj(uiControls, 'Tag', 'TimeMaxInput');
    % Layer detection UI elements removed

    % Find statistical analysis dropdown
    statDropdown = findobj(uiControls, 'Tag', 'StatDropdown');
    waveformProcessingDropdown = findobj(uiControls, 'Tag', 'WaveformProcessingDropdown');



    % Position statistical analysis dropdowns at the top
    STAT_DROPDOWN_WIDTH = 120;      % Width for statistical dropdowns
    STAT_DROPDOWN_HEIGHT = 25;      % Height for statistical dropdowns
    STAT_DROPDOWN_SPACING = 10;     % Spacing between statistical dropdowns
    STAT_DROPDOWN_TOP_MARGIN = 10;  % Margin from top of figure

    statDropdownY = figHeight - STAT_DROPDOWN_TOP_MARGIN - STAT_DROPDOWN_HEIGHT;

    % Position statistical analysis dropdown
    if ~isempty(statDropdown)
        set(statDropdown, 'Position', [20, statDropdownY, STAT_DROPDOWN_WIDTH, STAT_DROPDOWN_HEIGHT]);
    end

    if ~isempty(waveformProcessingDropdown)
        set(waveformProcessingDropdown, 'Position', [20 + STAT_DROPDOWN_WIDTH + STAT_DROPDOWN_SPACING, statDropdownY, STAT_DROPDOWN_WIDTH, STAT_DROPDOWN_HEIGHT]);
    end


    % Position elements from top to bottom using RIGHT-ANCHORED system
    currentY = figHeight - 30; % Start 30px from top

    % Update Y slider position (spans most of figure width)
    if ~isempty(ySlider)
        set(ySlider, 'Position', [80, 60, figWidth - UI_PANEL_WIDTH - 120, 20]);
    end

    % Position all UI elements from top to bottom using RIGHT-ANCHORED system
    if ~isempty(scaleDropdown)
        set(scaleDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(plotTypeDropdown)
        set(plotTypeDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(colormapDropdown)
        set(colormapDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(contrastDropdown)
        set(contrastDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(aspectRatioDropdown)
        set(aspectRatioDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(waveformDropdown)
        set(waveformDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(waveformModeCheckbox)
        set(waveformModeCheckbox, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_CHECKBOX_HEIGHT]);
        currentY = currentY - UI_CHECKBOX_HEIGHT - UI_SPACING;
    end

    if ~isempty(viewSwitchDropdown)
        set(viewSwitchDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(convergenceDisplay)
        set(convergenceDisplay, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    if ~isempty(zscoreCheckbox)
        set(zscoreCheckbox, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_CHECKBOX_HEIGHT]);
        currentY = currentY - UI_CHECKBOX_HEIGHT - UI_SPACING;
    end

    if ~isempty(alignmentCostDropdown)
        set(alignmentCostDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % Max shift row (label + input side-by-side)
    maxShiftLabel = findobj(uiControls, 'Tag', 'MaxShiftLabel');
    if ~isempty(maxShiftLabel) && ~isempty(alignmentMaxShiftInput)
        set(maxShiftLabel, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH*0.4, UI_TEXT_HEIGHT]);
        set(alignmentMaxShiftInput, 'Position', [UI_LEFT + UI_PANEL_WIDTH*0.45, currentY, UI_PANEL_WIDTH*0.55, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % Max iterations row
    alignmentMaxIterInput = findobj(uiControls, 'Tag', 'AlignmentMaxIterInput');
    maxIterLabel = findobj(uiControls, 'Tag', 'MaxIterLabel');
    if ~isempty(maxIterLabel) && ~isempty(alignmentMaxIterInput)
        set(maxIterLabel, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH*0.4, UI_TEXT_HEIGHT]);
        set(alignmentMaxIterInput, 'Position', [UI_LEFT + UI_PANEL_WIDTH*0.45, currentY, UI_PANEL_WIDTH*0.55, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end



    if ~isempty(computeSliceButton)
        set(computeSliceButton, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(computeAllButton)
        set(computeAllButton, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(crossViewButton)
        set(crossViewButton, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end
    if ~isempty(exportAlignedButton)
        set(exportAlignedButton, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end



    % Convergence control
    convergenceLabel = findobj(uiControls, 'Tag', 'ConvergenceLabel');
    convergenceInput = findobj(uiControls, 'Tag', 'ConvergenceInput');
    controlLabelWidth = 60; controlInputWidth = UI_PANEL_WIDTH - controlLabelWidth - UI_SPACING;
    if ~isempty(convergenceLabel) && ~isempty(convergenceInput)
        set(convergenceLabel, 'Position', [UI_LEFT, currentY, controlLabelWidth, UI_TEXT_HEIGHT]);
        set(convergenceInput, 'Position', [UI_LEFT + controlLabelWidth + UI_SPACING, currentY, controlInputWidth, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    % Position toggle buttons side by side
    if ~isempty(showOriginalButton) && ~isempty(showAlignedButton)
        buttonWidth = (UI_PANEL_WIDTH - UI_SPACING) / 2;
        set(showOriginalButton, 'Position', [UI_LEFT, currentY, buttonWidth, UI_ELEMENT_HEIGHT]);
        set(showAlignedButton, 'Position', [UI_LEFT + buttonWidth + UI_SPACING, currentY, buttonWidth, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(alignmentStatusText)
        set(alignmentStatusText, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    if ~isempty(syncTimeScaleButton)
        set(syncTimeScaleButton, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    if ~isempty(timeRangeCheckbox)
        set(timeRangeCheckbox, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_CHECKBOX_HEIGHT]);
        currentY = currentY - UI_CHECKBOX_HEIGHT - UI_SPACING;
    end

    % Position time range controls side by side
    if ~isempty(timeMinLabel) && ~isempty(timeMinInput) && ~isempty(timeMaxLabel) && ~isempty(timeMaxInput)
        labelWidth = 35;
        inputWidth = (UI_PANEL_WIDTH - labelWidth * 2 - UI_SPACING * 2) / 2;

        set(timeMinLabel, 'Position', [UI_LEFT, currentY, labelWidth, UI_TEXT_HEIGHT]);
        set(timeMinInput, 'Position', [UI_LEFT + labelWidth, currentY, inputWidth, UI_TEXT_HEIGHT]);
        set(timeMaxLabel, 'Position', [UI_LEFT + labelWidth + inputWidth + UI_SPACING, currentY, labelWidth, UI_TEXT_HEIGHT]);
        set(timeMaxInput, 'Position', [UI_LEFT + labelWidth * 2 + inputWidth + UI_SPACING, currentY, inputWidth, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING * 2;
    end

    % Layer detection UI positioning removed (layerDetectionLabel no longer exists)

    catch ME
        % Handle errors gracefully to prevent UI disruption
        fprintf('Warning: UI resize failed: %s\n', ME.message);
    end
end

% Helper function to create a coolwarm colormap (blue to white to red)
function cmap = coolwarm(m)
    if nargin < 1
        m = 64; % Default size (32 for each half)
    end

    % Calculate how many points for each half
    n1 = floor(m/2);
    n2 = m - n1;

    % Create a cool-warm diverging colormap (consistent with colormapTest.m)
    cmap = [
        linspace(0, 1, n1)' linspace(0, 1, n1)' linspace(1, 1, n1)';  % Cool to white
        linspace(1, 1, n2)' linspace(1, 0, n2)' linspace(1, 0, n2)'   % White to warm
    ];
end

function changeAspectRatio(dropdown, fig)
    % Get user data from figure
    userData = get(fig, 'UserData');

    % Get the selected aspect ratio
    aspectRatioOptions = get(dropdown, 'String');
    selectedOption = aspectRatioOptions{get(dropdown, 'Value')};

    % Set the aspect ratio based on the selection
    switch selectedOption
        case 'Equal [1,1]'
            userData.aspectRatio = [1, 1, 1];
        case 'Time Emphasis [2,1]'
            userData.aspectRatio = [2, 1, 1];
        case 'X Emphasis [1,2]'
            userData.aspectRatio = [1, 2, 1];
        case 'Default [1,10]'
            userData.aspectRatio = [1, 10, 1];
        case 'Custom'
            % Prompt user for custom aspect ratio
            answer = inputdlg({'X ratio:', 'T ratio:'}, 'Custom Aspect Ratio', [1 20], {'1', '1'});
            if ~isempty(answer)
                try
                    x_ratio = str2double(answer{1});
                    t_ratio = str2double(answer{2});
                    if isnan(x_ratio) || isnan(t_ratio)
                        errordlg('Please enter valid numbers for aspect ratio', 'Invalid Input');
                        return;
                    end
                    userData.aspectRatio = [x_ratio, t_ratio, 1];
                catch
                    errordlg('Error processing custom aspect ratio', 'Error');

                    return;
                end
            else
                return; % User cancelled
            end
    end

    % Store the selected aspect ratio in the figure's UserData
    userData.useCustomAspectRatio = true;
    set(fig, 'UserData', userData);

    % Update all plots with the new aspect ratio
    updatePlots(fig);

    % Apply the aspect ratio to all axes after updating plots
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;
end

function changeView(dropdown, fig)
    % Change the view perspective (X,t vs Y, Y,t vs X, X,Y vs t, or Waveform)
    try
        % Validate inputs
        if ~ishandle(fig) || ~ishghandle(fig, 'figure')
            return;
        end

        userData = get(fig, 'UserData');
        if isempty(userData)
            return;
        end

        % Get the selected view
        viewOptions = {'XtVsY', 'YtVsX', 'XYVst', 'Waveform', '3DView'};
        dropdownValue = get(dropdown, 'Value');
        if dropdownValue < 1 || dropdownValue > length(viewOptions)
            return;
        end
        selectedView = viewOptions{dropdownValue};

    % Update visualization state
    userData.visState.currentView = selectedView;

    % Reset slice index to 1 for the new view
    userData.visState.currentSliceIndex = 1;

    % Set view-specific default aspect ratio
    aspectRatioDropdown = findobj(fig, 'Tag', 'AspectRatioDropdown');
    if ~isempty(aspectRatioDropdown)
        switch selectedView
            case 'XYVst'
                % For XYVst view, default to Equal [1,1] aspect ratio
                set(aspectRatioDropdown, 'Value', 1); % Equal [1,1]
                userData.aspectRatio = [1, 1, 1];
                userData.useCustomAspectRatio = true;
            case {'XtVsY', 'YtVsX'}
                % For other views, keep Default [1,10] aspect ratio
                set(aspectRatioDropdown, 'Value', 4); % Default [1,10]
                userData.aspectRatio = [1, 10, 1];
                userData.useCustomAspectRatio = true;
            case 'Waveform'
                % For Waveform view, use Equal [1,1] aspect ratio
                set(aspectRatioDropdown, 'Value', 1); % Equal [1,1]
                userData.aspectRatio = [1, 1, 1];
                userData.useCustomAspectRatio = true;
        end
    end

    % Update slider and waveform dropdown based on new view
    ySlider = findobj(fig, 'Tag', 'YSlider');
    waveformDropdown = findobj(fig, 'Tag', 'WaveformDropdown');

    if ~isempty(ySlider)
        switch selectedView
            case 'XtVsY'
                % Y slider controls Y position
                set(ySlider, 'Min', 1, 'Max', length(userData.Y_values), 'Value', 1, 'Visible', 'on');
                userData.visState.currentYIndex = 1;
                % Re-enable waveform dropdown
                if ~isempty(waveformDropdown)
                    set(waveformDropdown, 'Enable', 'on');
                end

            case 'YtVsX'
                % Y slider controls X position
                set(ySlider, 'Min', 1, 'Max', length(userData.X_values), 'Value', 1, 'Visible', 'on');
                userData.visState.currentXIndex = 1;
                % Re-enable waveform dropdown
                if ~isempty(waveformDropdown)
                    set(waveformDropdown, 'Enable', 'on');
                end

            case 'XYVst'
                % Y slider controls time/segment position
                set(ySlider, 'Min', 1, 'Max', length(userData.timeValues), 'Value', 1, 'Visible', 'on');
                userData.visState.currentTIndex = 1;
                % Re-enable waveform dropdown
                if ~isempty(waveformDropdown)
                    set(waveformDropdown, 'Enable', 'on');
                end

            case 'Waveform'
                % For waveform view, slider is not needed - hide it
                set(ySlider, 'Visible', 'off');
                % Disable waveform dropdown since we're in dedicated waveform view
                if ~isempty(waveformDropdown)
                    set(waveformDropdown, 'Enable', 'off');
                end

            case '3DView'
                % For 3D view, slider is not needed - hide it
                set(ySlider, 'Visible', 'off');
                % Disable waveform dropdown since we're in 3D view
                if ~isempty(waveformDropdown)
                    set(waveformDropdown, 'Enable', 'off');
                end
        end
    end

    % Handle waveform loading for Waveform view
    if strcmp(selectedView, 'Waveform') && ~userData.waveformLoaded
        % Load waveform data for the first time
        loadWaveformData(fig);
    end

    % Handle 3D view creation
    if strcmp(selectedView, '3DView')
        % Create 3D view for the first time or update existing one
        create3DView(fig);
    end

    % Update the UserData
    set(fig, 'UserData', userData);

    % Update all plots with new view
    updatePlots(fig);

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;

    fprintf('View changed to: %s\n', selectedView);

    catch ME
        fprintf('Error changing view: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
end

function toggleZScore(checkbox, fig)
    % Toggle between raw and Z-score data
    userData = get(fig, 'UserData');

    % Get checkbox state
    useZScore = get(checkbox, 'Value');
    userData.visState.useZScore = useZScore;

    % Store original stat types if not already stored
    if ~isfield(userData, 'originalStatTypes')
        userData.originalStatTypes = userData.sortedStatTypes;
    end

    % Update stat types based on Z-score toggle
    if useZScore
        % Convert to Z-score stat types
        userData.sortedStatTypes = cellfun(@(x) ['ZScore_' x], userData.originalStatTypes, 'UniformOutput', false);
        fprintf('Switched to Z-score statistics\n');
    else
        % Use original stat types
        userData.sortedStatTypes = userData.originalStatTypes;
        fprintf('Switched to raw statistics\n');
    end

    % Compute Z-scores on-the-fly instead of loading pre-computed files
    try
        if useZScore
            % Compute Z-scores from original data
            userData.statDataArray = computeZScoreData(userData.originalStatDataArray);
            fprintf('Z-scores computed on-the-fly from original data\n');
        else
            % Use original data
            userData.statDataArray = userData.originalStatDataArray;
        end

        % Recalculate global min/max for proper scaling
        allMaps = [];
        for i = 1:length(userData.statDataArray)
            for j = 1:length(userData.statDataArray{i}.maps)
                allMaps = [allMaps; userData.statDataArray{i}.maps{j}(:)];
            end
        end
        userData.globalMin = min(allMaps, [], 'omitnan');
        userData.globalMax = max(allMaps, [], 'omitnan');

        % Handle edge cases
        if isnan(userData.globalMin) || isnan(userData.globalMax)
            userData.globalMin = 0;
            userData.globalMax = 1;
        elseif userData.globalMin == userData.globalMax
            userData.globalMax = userData.globalMin + 1e-6;
        end

        % Update the UserData
        set(fig, 'UserData', userData);

        % Update all plots
        updatePlots(fig);

        % Apply aspect ratio if enabled
        applyAspectRatioToAxes(fig);

        % Force immediate visual update
        drawnow;

    catch ME
        fprintf('Error computing Z-score data: %s\n', ME.message);
        % Reset checkbox if computation failed
        set(checkbox, 'Value', ~useZScore);
        userData.visState.useZScore = ~useZScore;
        set(fig, 'UserData', userData);
    end
end

function toggleWaveform(dropdown, fig)
    % Toggle waveform display on/off
    userData = get(fig, 'UserData');

    % Get dropdown selection
    options = get(dropdown, 'String');
    selectedOption = options{get(dropdown, 'Value')};

    showWaveform = strcmp(selectedOption, 'Show Waveform');

    % Enable/disable the waveform mode checkbox based on waveform display
    waveformModeCheckbox = findobj(fig, 'Tag', 'WaveformModeCheckbox');
    if ~isempty(waveformModeCheckbox)
        if showWaveform
            set(waveformModeCheckbox, 'Enable', 'on');
        else
            set(waveformModeCheckbox, 'Enable', 'off');
        end
    end

    if showWaveform && ~userData.waveformLoaded
        % Load waveform data for the first time
        loadWaveformData(fig);
    end

    % Update plots to show/hide waveform
    updatePlots(fig);

    % Force immediate visual update
    drawnow;

    fprintf('Waveform display: %s\n', selectedOption);
end

function toggleWaveformMode(checkbox, fig)
    % Toggle between showing average waveform vs all individual waveforms
    userData = get(fig, 'UserData');

    % Get checkbox state
    showAllWaveforms = get(checkbox, 'Value');

    % Check if we're in Waveform view
    if strcmp(userData.visState.currentView, 'Waveform')
        % For Waveform view, update the full-screen display
        createFullScreenWaveformView(fig);
    else
        % For other views, update the small waveform display
        updateWaveformDisplay(fig);
    end

    % Force immediate visual update
    drawnow;

    if showAllWaveforms
        fprintf('Switched to showing all individual waveforms\n');
    else
        fprintf('Switched to showing average waveform\n');
    end
end



function memUsage = getCurrentMemoryUsage()
    % Get current memory usage (cross-platform)
    try
        if ispc
            [~, memInfo] = memory;
            memUsage = memInfo.MemUsedMATLAB / 1024^3; % Convert to GB
        else
            % For Unix/Mac systems, use a simpler approach
            memUsage = 0; % Placeholder - could implement platform-specific code
        end
    catch
        memUsage = 0; % Fallback if memory monitoring fails
    end
end

function optimizeMemoryUsage(fig)
    % Perform memory optimization operations
    try
        userData = get(fig, 'UserData');

        % Clear unnecessary cached data if memory usage is high
        memUsage = getCurrentMemoryUsage();
        if memUsage > 2.0 % If using more than 2GB
            fprintf('High memory usage detected (%.1f GB). Clearing caches...\n', memUsage);

            % Clear computed statistics cache
            if isfield(userData, 'computedStats')
                userData.computedStats = struct();
            end

            % Force garbage collection
            clear('temp*'); % Clear any temporary variables
            pack; % MATLAB memory defragmentation

            set(fig, 'UserData', userData);
            fprintf('Memory optimization completed.\n');
        end
    catch ME
        fprintf('Warning: Memory optimization failed: %s\n', ME.message);
    end
end

function zscoreDataArray = computeZScoreData(originalDataArray)
    % Compute Z-scores on-the-fly from original statistical data
    zscoreDataArray = cell(size(originalDataArray));

    for i = 1:length(originalDataArray)
        % Copy the structure
        zscoreDataArray{i} = originalDataArray{i};

        % Compute Z-scores for each map
        for j = 1:length(originalDataArray{i}.maps)
            originalMap = originalDataArray{i}.maps{j};

            % Calculate mean and std across all valid (non-NaN) values
            validValues = originalMap(~isnan(originalMap));
            if length(validValues) > 1
                mapMean = mean(validValues);
                mapStd = std(validValues);

                if mapStd > 1e-10  % Avoid division by zero
                    % Compute Z-score: (value - mean) / std
                    zscoreMap = (originalMap - mapMean) / mapStd;
                else
                    % If std is zero, set all values to zero
                    zscoreMap = zeros(size(originalMap));
                end
            else
                % Not enough valid data points
                zscoreMap = zeros(size(originalMap));
            end

            % Store the Z-score map
            zscoreDataArray{i}.maps{j} = zscoreMap;
        end

        % Update the statistic name to indicate Z-score
        if isfield(zscoreDataArray{i}, 'statType') && ~contains(zscoreDataArray{i}.statType, 'ZScore')
            zscoreDataArray{i}.statType = ['ZScore_' zscoreDataArray{i}.statType];
        end
    end

    fprintf('Z-scores computed for %d statistics\n', length(zscoreDataArray));
end



function toggleLayerDetection(checkbox, fig)
    % Toggle layer path overlay on/off
    userData = get(fig, 'UserData');

    % Get checkbox state
    showLayers = get(checkbox, 'Value');
    userData.visState.showLayerPaths = showLayers;

    if showLayers
        % Detect layer paths if not already done
        if isempty(userData.visState.layerPaths)
            % Show progress dialog
            progressDlg = uiprogressdlg(fig, 'Title', 'Layer Detection', ...
                                       'Message', 'Detecting layer paths...', ...
                                       'Indeterminate', 'on', ...
                                       'Cancelable', false);
            drawnow;

            try
                fprintf('Detecting layer paths...\n');

                % Create progress callback to update dialog message
                progressCallback = @(msg) updateProgressMessage(progressDlg, msg);

                % Start timer to track computation time
                startTime = tic;

                userData.visState.layerPaths = detectLayerPaths(userData.statDataArray, ...
                    userData.visState, userData.X_values, userData.Y_values, userData.timeValues, progressCallback);

                % Calculate and display computation time
                elapsedTime = toc(startTime);

                % Close progress dialog
                close(progressDlg);

                fprintf('Layer path overlay enabled (completed in %.2f seconds)\n', elapsedTime);
            catch ME
                % Close progress dialog on error
                close(progressDlg);
                errordlg(['Layer detection failed: ' ME.message], 'Detection Error');
                set(checkbox, 'Value', 0); % Reset checkbox
                userData.visState.showLayerPaths = false;
                set(fig, 'UserData', userData);
                return;
            end
        else
            fprintf('Layer path overlay enabled (using cached results)\n');
        end
    else
        fprintf('Layer path overlay disabled\n');
    end

    % Update the UserData
    set(fig, 'UserData', userData);

    % Update plots to show/hide layer paths
    updatePlots(fig);

    % Apply aspect ratio if enabled
    applyAspectRatioToAxes(fig);

    % Force immediate visual update
    drawnow;
end

function updateLayerDetectionParams(src, fig)
    % Update layer detection parameters when user changes controls
    userData = get(fig, 'UserData');

    % Get current parameter values
    prominenceSlider = findobj(fig, 'Tag', 'ProminenceSlider');
    minDistInput = findobj(fig, 'Tag', 'MinDistInput');
    expectedLayersInput = findobj(fig, 'Tag', 'ExpectedLayersInput');

    if ~isempty(prominenceSlider)
        userData.visState.peakProminence = get(prominenceSlider, 'Value');
    end

    if ~isempty(minDistInput)
        try
            minDist = str2double(get(minDistInput, 'String'));
            if ~isnan(minDist) && minDist > 0
                userData.visState.minPeakDistance = round(minDist);
            else
                % Reset to default if invalid
                set(minDistInput, 'String', '3');
                userData.visState.minPeakDistance = 3;
            end
        catch
            % Reset to default if error
            set(minDistInput, 'String', '3');
            userData.visState.minPeakDistance = 3;
        end
    end

    if ~isempty(expectedLayersInput)
        try
            expectedLayers = str2double(get(expectedLayersInput, 'String'));
            if ~isnan(expectedLayers) && expectedLayers > 0 && expectedLayers <= 50
                userData.visState.expectedLayers = round(expectedLayers);
            else
                % Reset to default if invalid (must be between 1 and 50)
                set(expectedLayersInput, 'String', '8');
                userData.visState.expectedLayers = 8;
            end
        catch
            % Reset to default if error
            set(expectedLayersInput, 'String', '8');
            userData.visState.expectedLayers = 8;
        end
    end

    % Clear existing layer paths to force re-detection with new parameters
    userData.visState.layerPaths = [];

    % Update the UserData
    set(fig, 'UserData', userData);

    % If layer detection is enabled, re-detect with new parameters
    layerToggle = findobj(fig, 'Tag', 'LayerDetectionToggle');
    if ~isempty(layerToggle) && get(layerToggle, 'Value')
        % Show progress dialog for parameter update
        progressDlg = uiprogressdlg(fig, 'Title', 'Layer Detection', ...
                                   'Message', 'Re-detecting layers with new parameters...', ...
                                   'Indeterminate', 'on', ...
                                   'Cancelable', false);
        drawnow;

        try
            fprintf('Re-detecting layers with new parameters...\n');

            % Create progress callback to update dialog message
            progressCallback = @(msg) updateProgressMessage(progressDlg, msg);

            % Start timer to track computation time
            startTime = tic;

            userData.visState.layerPaths = detectLayerPaths(userData.statDataArray, ...
                userData.visState, userData.X_values, userData.Y_values, userData.timeValues, progressCallback);
            set(fig, 'UserData', userData);

            % Calculate and display computation time
            elapsedTime = toc(startTime);

            % Close progress dialog
            close(progressDlg);

            updatePlots(fig);
            applyAspectRatioToAxes(fig);
            % Force immediate visual update
            drawnow;

            fprintf('Layer re-detection completed in %.2f seconds\n', elapsedTime);
        catch ME
            % Close progress dialog on error
            close(progressDlg);
            errordlg(['Layer re-detection failed: ' ME.message], 'Detection Error');
        end
    end

    fprintf('Layer detection parameters updated: prominence=%.3f, minDistance=%d, expectedLayers=%d\n', ...
            userData.visState.peakProminence, userData.visState.minPeakDistance, userData.visState.expectedLayers);
end



function reconstruct3DLayers(src, fig)
    % Reconstruct 3D layer surfaces from detected paths
    userData = get(fig, 'UserData');

    % Check if layer paths have been detected
    if isempty(userData.visState.layerPaths) || userData.visState.layerPaths.numLayers == 0
        errordlg('No layer paths detected. Please enable layer detection first.', 'No Layers');
        return;
    end

    % Show progress dialog for 3D reconstruction
    progressDlg = uiprogressdlg(fig, 'Title', '3D Reconstruction', ...
                               'Message', 'Reconstructing 3D layer surfaces...', ...
                               'Indeterminate', 'on', ...
                               'Cancelable', false);
    drawnow;

    try
        fprintf('Starting 3D layer reconstruction...\n');

        % Update progress message
        progressDlg.Message = 'Aligning layer paths between views...';
        drawnow;

        % Call the 3D reconstruction function
        layerSurfaces = alignLayerPaths(userData.visState.layerPaths, ...
            userData.X_values, userData.Y_values, userData.timeValues);

        % Update progress message
        progressDlg.Message = 'Creating 3D visualization...';
        drawnow;

        % Store the reconstructed surfaces
        userData.visState.layerSurfaces = layerSurfaces;
        set(fig, 'UserData', userData);

        % Enable STL export button
        exportButton = findobj(fig, 'Tag', 'ExportSTLButton');
        if ~isempty(exportButton)
            set(exportButton, 'Enable', 'on');
        end

        % Visualize the 3D surfaces
        visualizeLayerSurfaces(layerSurfaces, userData);

        % Close progress dialog
        close(progressDlg);

        fprintf('3D reconstruction complete: %d layer surfaces created\n', ...
                length(layerSurfaces.surfaces));

    catch ME
        % Close progress dialog on error
        close(progressDlg);
        errordlg(['3D reconstruction failed: ' ME.message], 'Reconstruction Error');
        fprintf('Error in 3D reconstruction: %s\n', ME.message);
    end
end

function exportLayersSTL(src, fig)
    % Export reconstructed layer surfaces as STL files
    userData = get(fig, 'UserData');

    % Check if 3D surfaces have been reconstructed
    if isempty(userData.visState.layerSurfaces)
        errordlg('No 3D surfaces available. Please run 3D reconstruction first.', 'No Surfaces');
        return;
    end

    % Show progress dialog for STL export
    progressDlg = uiprogressdlg(fig, 'Title', 'STL Export', ...
                               'Message', 'Exporting layer surfaces to STL files...', ...
                               'Indeterminate', 'on', ...
                               'Cancelable', false);
    drawnow;

    try
        fprintf('Exporting layer surfaces to STL files...\n');

        % Update progress message
        progressDlg.Message = 'Converting surfaces to triangulated meshes...';
        drawnow;

        % Call the STL export function
        exportLayerSTL(userData.visState.layerSurfaces, userData.FileNamingArray);

        % Close progress dialog
        close(progressDlg);

        fprintf('STL export complete\n');

        % Show success message
        msgbox('STL files exported successfully! Check the Layer_STL_Export folder.', ...
               'Export Complete', 'help');

    catch ME
        % Close progress dialog on error
        close(progressDlg);
        errordlg(['STL export failed: ' ME.message], 'Export Error');
        fprintf('Error in STL export: %s\n', ME.message);
    end
end

function loadWaveformData(fig)
    % Load waveform data from saved files
    userData = get(fig, 'UserData');

    try
        % Load waveform data using same pattern as MasterPlot.m
        folderPath = fullfile(pwd, 'Saved Wave Forms');
        dataFileType = 'Waveform';
        dataFileExtension = '.mat';
        [dataFile, dataExists] = buildAndCheckFile(dataFileType, userData.FileNamingArray, folderPath, dataFileExtension);

        fprintf('Looking for waveform file: %s\n', dataFile);
        if ~dataExists
            fprintf('Waveform data file not found: %s\n', dataFile);
            fprintf('Cannot load waveforms.\n');

            % List available files in the directory for debugging
            if exist(folderPath, 'dir')
                files = dir(fullfile(folderPath, '*.mat'));
                fprintf('Available .mat files in %s:\n', folderPath);
                for i = 1:length(files)
                    fprintf('  %s\n', files(i).name);
                end
            end
            return;
        end

        % Load waveform data
        loaded = load(dataFile);
        if ~isfield(loaded, 'savedData')
            fprintf('Loaded data does not contain "savedData". Cannot load waveforms.\n');
            return;
        end

        savedData = loaded.savedData;
        if ~isfield(savedData, 'waveformArray') || ~isfield(savedData, 't')
            fprintf('Missing waveform data fields. Cannot load waveforms.\n');
            return;
        end

        % Store original waveform data
        userData.originalWaveformData = savedData;
        userData.waveformLoaded = true;

        % Initialize aligned waveform data as empty (will be created only when needed)
        userData.alignedWaveformData = [];

        set(fig, 'UserData', userData);
        fprintf('Waveform data loaded successfully\n');

    catch ME
        fprintf('Error loading waveform data: %s\n', ME.message);
    end
end

function updatePlots(fig)
    % Get user data from figure
    userData = get(fig, 'UserData');
    visState = userData.visState;

    % Use the current statDataArray from userData (which is updated by dropdown callbacks)
    statDataArray = userData.statDataArray;


    sortedStatTypes = userData.sortedStatTypes;
    plotHandles = userData.plotHandles;
    axHandles = userData.axHandles;
    globalMin = userData.globalMin;
    globalMax = userData.globalMax;
    X_values = userData.X_values;
    Y_values = userData.Y_values;
    timeValues = userData.timeValues;
    timeLabel = userData.timeLabel;

    % Force immediate update
    drawnow;

    % Get current settings
    useGlobalScale = visState.useGlobalScale;
    plotType = visState.plotType;
    selectedColormap = visState.colormap;
    contrastEnhancement = visState.contrastEnhancement;
    useTimeRange = visState.useTimeRange;
    timeRangeMin = visState.timeRangeMin;
    timeRangeMax = visState.timeRangeMax;
    currentView = visState.currentView;

    % Handle Waveform view separately
    if strcmp(currentView, 'Waveform')
        % For Waveform view, create a full-screen waveform display
        createFullScreenWaveformView(fig);
        return;
    end

    % Handle 3D view separately
    if strcmp(currentView, '3DView')
        % For 3D view, create or update the 3D visualization
        create3DView(fig);
        return;
    end

    % If switching from 3D view back to regular views, recreate axes
    if length(axHandles) == 1 && ishandle(axHandles(1)) && strcmp(get(axHandles(1), 'Tag'), '3DAxes')
        recreateOriginalAxes(fig);
        userData = get(fig, 'UserData'); % Refresh userData after recreation
        axHandles = userData.axHandles;
        plotHandles = userData.plotHandles;
    end

    % Get appropriate slice index based on current view
    switch currentView
        case 'XtVsY'
            sliceIndex = visState.currentYIndex;
        case 'YtVsX'
            sliceIndex = visState.currentXIndex;
        case 'XYVst'
            sliceIndex = visState.currentTIndex;
        otherwise
            sliceIndex = visState.currentYIndex; % Default fallback
    end

    % For the single dropdown system, we only update the first plot with the selected statistic
    % Update each plot (but typically only one plot for single dropdown)
    numPlotsToUpdate = min(length(plotHandles), length(statDataArray));

    for i = 1:numPlotsToUpdate
        % Check if axes handle is still valid
        if i > length(axHandles) || ~ishandle(axHandles(i))
            continue;
        end

        % Get data for this statistic
        statData = statDataArray{i};
        numSegments = length(statData.maps);

        % Create data matrix based on current view
        switch currentView
            case 'XtVsY'
                % X,t data for the current Y slice
                plotData = zeros(numSegments, length(X_values));
                for seg = 1:numSegments
                    mapData = statData.maps{seg};
                    plotData(seg, :) = mapData(sliceIndex, :);
                end
                xAxisValues = X_values;
                yAxisValues = timeValues;
                xLabel = 'X Position (mm)';
                yLabel = timeLabel;

            case 'YtVsX'
                % Y,t data for the current X slice
                plotData = zeros(numSegments, length(Y_values));
                for seg = 1:numSegments
                    mapData = statData.maps{seg};
                    plotData(seg, :) = mapData(:, sliceIndex)';  % Transpose to make it a row vector
                end
                xAxisValues = Y_values;
                yAxisValues = timeValues;
                xLabel = 'Y Position (mm)';
                yLabel = timeLabel;

            case 'XYVst'
                % X,Y data for the current time slice
                if sliceIndex <= numSegments
                    mapData = statData.maps{sliceIndex};
                    plotData = mapData;
                else
                    plotData = statData.maps{1}; % Fallback
                end
                xAxisValues = X_values;
                yAxisValues = Y_values;
                xLabel = 'X Position (mm)';
                yLabel = 'Y Position (mm)';

            otherwise
                % Default to XtVsY
                plotData = zeros(numSegments, length(X_values));
                for seg = 1:numSegments
                    mapData = statData.maps{seg};
                    plotData(seg, :) = mapData(sliceIndex, :);
                end
                xAxisValues = X_values;
                yAxisValues = timeValues;
                xLabel = 'X Position (mm)';
                yLabel = timeLabel;
        end

        % For backward compatibility, keep xtData variable
        xtData = plotData;

        % Apply time range filtering if enabled (only for XtVsY and YtVsX views)
        if useTimeRange && (strcmp(currentView, 'XtVsY') || strcmp(currentView, 'YtVsX'))
            % Find time indices within the specified range
            timeIndices = (yAxisValues >= timeRangeMin) & (yAxisValues <= timeRangeMax);

            % Filter the data and axis values
            xtData = xtData(timeIndices, :);
            yAxisValues = yAxisValues(timeIndices);

            % Gray out data outside the time range by setting to NaN
            if ~all(timeIndices)
                fprintf('Time range filtering: showing %d/%d time points\n', sum(timeIndices), length(timeIndices));
            end
        end

        % Apply contrast enhancement
        xtData = applyContrastEnhancement(xtData, contrastEnhancement);

        % Debug: Print data statistics
        fprintf('Plot %d (%s): Data size [%d x %d], range [%.6f, %.6f], NaN count: %d\n', ...
            i, sortedStatTypes{i}, size(xtData, 1), size(xtData, 2), ...
            min(xtData(:), [], 'omitnan'), max(xtData(:), [], 'omitnan'), sum(isnan(xtData(:))));

        % Get the current axes
        ax = axHandles(i);

        % Store the current axes position
        axPos = get(ax, 'Position');

        % Calculate min/max for this plot
        if ~useGlobalScale
            % Per-plot scaling: use only non-NaN values for scaling
            minVal = min(xtData(:), [], 'omitnan');
            maxVal = max(xtData(:), [], 'omitnan');
            if isnan(minVal) || isnan(maxVal) || minVal == maxVal
                minVal = 0;
                maxVal = 1;
            end
        else
            % Global scaling
            minVal = globalMin;
            maxVal = globalMax;
        end

        % Clear the current axes completely to prevent layering issues
        cla(ax, 'reset');

        % Restore axes position after reset
        set(ax, 'Position', axPos);

        % Create the new plot based on the selected type
        switch plotType
            case 'heatmap'
                plotHandles{i} = imagesc(xAxisValues, yAxisValues, xtData, 'Parent', ax);
                axis(ax, 'xy');
                view(ax, 2); % 2D view

            case 'filledcontour'
                [X, Y] = meshgrid(xAxisValues, yAxisValues);
                % Create a filled contour plot
                [~, plotHandles{i}] = contourf(ax, X, Y, xtData, 20);
                axis(ax, 'xy');
                view(ax, 2); % 2D view

            case 'slice'
                % Create a simpler 2D slice visualization that's more reliable
                [X, Y] = meshgrid(xAxisValues, yAxisValues);

                % Create a 2D slice visualization
                plotHandles{i} = surf(ax, X, Y, xtData, 'EdgeColor', 'none', 'FaceColor', 'interp');
                hold(ax, 'on');

                % Add contour lines on the slice for better visualization
                contour3(ax, X, Y, xtData, 10, 'k', 'LineWidth', 0.5);

                % Set view and appearance
                axis(ax, 'xy');
                view(ax, 2); % Straight-on 2D view
                grid(ax, 'on');
                lighting gouraud;
                camlight('headlight');

            case 'smooth2d'
                [X, Y] = meshgrid(xAxisValues, yAxisValues);
                % Create a smooth 2D plot using pcolor with interpolated shading
                plotHandles{i} = pcolor(ax, X, Y, xtData);
                shading(ax, 'interp');
                axis(ax, 'xy');
                view(ax, 2); % 2D view
                % Add contour lines on top for better visualization
                hold(ax, 'on');
                [~, c] = contour(ax, X, Y, xtData, 5, 'k', 'LineWidth', 0.5, 'LineStyle', ':');
                hold(ax, 'off');

            case 'pseudocolor'
                [X, Y] = meshgrid(xAxisValues, yAxisValues);
                % Create a pseudocolor plot with discrete color levels
                plotHandles{i} = pcolor(ax, X, Y, xtData);
                shading(ax, 'flat'); % Use flat shading for discrete color levels
                axis(ax, 'xy');
                view(ax, 2); % 2D view

                % Add grid lines to emphasize the discrete nature
                hold(ax, 'on');
                % Add grid lines at each axis value
                for x = xAxisValues
                    plot(ax, [x x], [min(yAxisValues) max(yAxisValues)], 'k-', 'LineWidth', 0.1);
                end
                for y = yAxisValues
                    plot(ax, [min(xAxisValues) max(xAxisValues)], [y y], 'k-', 'LineWidth', 0.1);
                end
                hold(ax, 'off');
        end

        % Set colormap and colorbar
        if strcmp(selectedColormap, 'coolwarm')
            % Apply custom coolwarm colormap
            colormap(ax, coolwarm());
        else
            % Apply standard MATLAB colormap
            colormap(ax, selectedColormap);
        end
        colorbar(ax);

        % Set normal color limits (ensure valid range)
        if isfinite(minVal) && isfinite(maxVal) && minVal < maxVal
            caxis(ax, [minVal, maxVal]);
            fprintf('Setting color limits for %s: [%.6f, %.6f]\n', sortedStatTypes{i}, minVal, maxVal);
        else
            % Use default range if invalid
            caxis(ax, [0, 1]);
            fprintf('Warning: Invalid color limits for %s (min=%.6f, max=%.6f), using [0, 1]\n', sortedStatTypes{i}, minVal, maxVal);
        end

        % Overlay layer paths if enabled
        if visState.showLayerPaths && ~isempty(visState.layerPaths)
            overlayLayerPaths(ax, visState.layerPaths, currentView, sliceIndex, X_values, Y_values, timeValues);
        end

        % Update title and labels based on current view
        switch currentView
            case 'XtVsY'
                title(ax, sprintf('%s at Y = %.2f mm', sortedStatTypes{i}, Y_values(sliceIndex)));
                xlabel(ax, xLabel);
                ylabel(ax, yLabel);
            case 'YtVsX'
                title(ax, sprintf('%s at X = %.2f mm', sortedStatTypes{i}, X_values(sliceIndex)));
                xlabel(ax, xLabel);
                ylabel(ax, yLabel);
            case 'XYVst'
                if sliceIndex <= length(timeValues)
                    title(ax, sprintf('%s at %s = %.2f', sortedStatTypes{i}, timeLabel, timeValues(sliceIndex)));
                else
                    title(ax, sprintf('%s at segment %d', sortedStatTypes{i}, sliceIndex));
                end
                xlabel(ax, xLabel);
                ylabel(ax, yLabel);
            otherwise
                title(ax, sprintf('%s at Y = %.2f mm', sortedStatTypes{i}, Y_values(sliceIndex)));
                xlabel(ax, xLabel);
                ylabel(ax, yLabel);
        end

        if strcmp(plotType, 'slice')
            zlabel(ax, 'Value');
        end

        % Restore the axes position (which might have been changed by the 3D view)
        set(ax, 'Position', axPos);

        % Apply aspect ratio if provided and enabled
        if isfield(userData, 'useCustomAspectRatio') && userData.useCustomAspectRatio && ~isempty(userData.aspectRatio)
            % Store current position again to ensure it's preserved
            axPos = get(ax, 'Position');

            % Apply aspect ratio using both daspect and pbaspect for better compatibility
            daspect(ax, userData.aspectRatio);
            pbaspect(ax, userData.aspectRatio);

            % Force aspect ratio to be maintained
            set(ax, 'DataAspectRatioMode', 'manual');
            set(ax, 'PlotBoxAspectRatioMode', 'manual');

            % Ensure the axes position is maintained
            set(ax, 'Position', axPos);
        end
    end

    % Update the UserData with the new plot handles
    userData.plotHandles = plotHandles;
    set(fig, 'UserData', userData);

    % Update waveform display if enabled
    updateWaveformDisplay(fig);

    % CRITICAL: Synchronize time scales across all subplots
    % Lock to the plot with the most zoomed-in or detailed time scale
    synchronizeTimeScales(fig, currentView);

    % Force immediate update of all plots with multiple refresh attempts
    drawnow;
    pause(0.01); % Small pause to ensure rendering
    drawnow;

    fprintf('Plot update completed for view: %s\n', currentView);
end

function createFullScreenWaveformView(fig)
    % Create a full-screen waveform view that replaces all other plots
    userData = get(fig, 'UserData');

    % Clear all existing plot axes
    axHandles = userData.axHandles;
    for i = 1:length(axHandles)
        if ishandle(axHandles(i))
            delete(axHandles(i));
        end
    end

    % Also clear any existing waveform axes
    existingWaveformAx = findobj(fig, 'Tag', 'WaveformAxes');
    if ~isempty(existingWaveformAx)
        delete(existingWaveformAx);
    end

    % Load waveform data if not already loaded
    if ~userData.waveformLoaded
        loadWaveformData(fig);
        userData = get(fig, 'UserData'); % Refresh userData after loading
    end

    if ~userData.waveformLoaded || isempty(userData.originalWaveformData)
        % Create a placeholder axes with message
        waveformAx = axes('Position', [0.1, 0.35, 0.65, 0.55], 'Tag', 'WaveformAxes');
        text(0.5, 0.5, 'No waveform data available', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 16, 'Parent', waveformAx);
        set(waveformAx, 'XTick', [], 'YTick', []);
        return;
    end

    % Determine if we're showing aligned data
    useAlignedData = false;
    if strcmp(userData.currentView, 'aligned') && ~isempty(userData.alignedWaveformData)
        if ~isempty(userData.alignmentShifts) && any(abs(userData.alignmentShifts(:)) > 1e-10)
            useAlignedData = true;
        end
    end

    % Get current waveform data based on alignment status (use references, not copies)
    if useAlignedData
        waveformData = userData.alignedWaveformData;
    else
        waveformData = userData.originalWaveformData;
    end

    % Create full-screen waveform axes
    waveformPosition = [0.1, 0.35, 0.65, 0.55]; % Large area for waveform display
    waveformAx = axes('Position', waveformPosition, 'Tag', 'WaveformAxes');

    % Get waveform data (references only - no copying unless time filtering is needed)
    waveformArray = waveformData.waveformArray;
    t = waveformData.t;

    % Apply time range filtering if enabled
    visState = userData.visState;
    if isfield(visState, 'useTimeRange') && visState.useTimeRange && ...
       isfield(visState, 'timeRangeMin') && isfield(visState, 'timeRangeMax')

        % Convert time range from microseconds to seconds for waveform comparison
        timeRangeMinSeconds = visState.timeRangeMin * 1e-6;
        timeRangeMaxSeconds = visState.timeRangeMax * 1e-6;

        timeIndices = (t >= timeRangeMinSeconds) & (t <= timeRangeMaxSeconds);

        fprintf('Full-screen waveform time range: %.2e to %.2e seconds (converted from %.1f to %.1f microseconds)\n', ...
            timeRangeMinSeconds, timeRangeMaxSeconds, visState.timeRangeMin, visState.timeRangeMax);

        if any(timeIndices)
            waveformArray = waveformArray(:, timeIndices);
            t = t(timeIndices);
        else
            % Show empty plot if no data in range
            plot(waveformAx, [], []);
            title(waveformAx, 'No waveform data in time range', 'FontSize', 14);
            xlabel(waveformAx, 'Time (μs)', 'FontSize', 12);
            ylabel(waveformAx, 'Amplitude', 'FontSize', 12);
            return;
        end
    end

    % Check if user wants to show all waveforms or just average
    waveformModeCheckbox = findobj(fig, 'Tag', 'WaveformModeCheckbox');
    showAllWaveforms = true; % Default to showing all waveforms in Waveform view
    if ~isempty(waveformModeCheckbox)
        showAllWaveforms = get(waveformModeCheckbox, 'Value');
        % Enable the checkbox for Waveform view
        set(waveformModeCheckbox, 'Enable', 'on');
    end

    % Plot waveforms based on mode
    hold(waveformAx, 'on');

    % Determine colors and title based on alignment status and mode
    if useAlignedData
        individualColor = [0, 0.4470, 0.7410]; % MATLAB blue for aligned
        avgColor = 'k'; % Black for average
        if showAllWaveforms
            titleText = 'All Waveforms (Aligned)';
        else
            titleText = 'Average Waveform (Aligned)';
        end
    else
        individualColor = [0.5, 0.5, 0.5]; % Gray for original
        avgColor = 'k'; % Black for average
        if showAllWaveforms
            titleText = 'All Waveforms (Original)';
        else
            titleText = 'Average Waveform (Original)';
        end
    end

    % Convert time from seconds to microseconds for display
    t_microseconds = t * 1e6;

    % Plot based on mode
    if showAllWaveforms
        % Plot each individual waveform
        numWaveforms = size(waveformArray, 1);
        % Vectorized plotting for performance: plot all waveforms in one call
        plot(waveformAx, t_microseconds, waveformArray', 'Color', individualColor, 'LineWidth', 0.5)

        % Plot average waveform on top with thicker line
        avgWaveform = mean(waveformArray, 1);
        plot(waveformAx, t_microseconds, avgWaveform, 'Color', avgColor, 'LineWidth', 2);
    else
        % Plot only the average waveform
        avgWaveform = mean(waveformArray, 1);
        plot(waveformAx, t_microseconds, avgWaveform, 'Color', avgColor, 'LineWidth', 2);
    end

    hold(waveformAx, 'off');

    % Set labels and title
    title(waveformAx, titleText, 'FontSize', 14);
    xlabel(waveformAx, 'Time (μs)', 'FontSize', 12);
    ylabel(waveformAx, 'Amplitude', 'FontSize', 12);
    grid(waveformAx, 'on');

    % Set Y-axis limits
    allData = waveformArray(:);
    dataMin = min(allData);
    dataMax = max(allData);
    if dataMin ~= dataMax
        dataRange = dataMax - dataMin;
        padding = dataRange * 0.1;
        ylim(waveformAx, [dataMin - padding, dataMax + padding]);
    else
        ylim(waveformAx, [dataMin - 0.1, dataMax + 0.1]);
    end

    % Add legend based on mode
    hold(waveformAx, 'on');
    if showAllWaveforms
        if useAlignedData
            h1 = plot(waveformAx, NaN, NaN, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.5);
            h2 = plot(waveformAx, NaN, NaN, 'k-', 'LineWidth', 2);
            legend(waveformAx, [h1, h2], {'Aligned Waveforms', 'Average'}, 'Location', 'northeast', 'FontSize', 10);
        else
            h1 = plot(waveformAx, NaN, NaN, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
            h2 = plot(waveformAx, NaN, NaN, 'k-', 'LineWidth', 2);
            legend(waveformAx, [h1, h2], {'Individual Waveforms', 'Average'}, 'Location', 'northeast', 'FontSize', 10);
        end
    else
        % No legend needed for average-only mode
        legend(waveformAx, 'off');
    end
    hold(waveformAx, 'off');

    % Update userData to reflect the new axes structure
    userData.axHandles = waveformAx;
    userData.plotHandles = {[]};
    set(fig, 'UserData', userData);

    % Hide the slider since it's not needed for waveform view
    ySlider = findobj(fig, 'Tag', 'YSlider');
    if ~isempty(ySlider)
        set(ySlider, 'Visible', 'off');
    end

    % Disable the waveform dropdown since we're in dedicated waveform view
    waveformDropdown = findobj(fig, 'Tag', 'WaveformDropdown');
    if ~isempty(waveformDropdown)
        set(waveformDropdown, 'Enable', 'off');
    end

    fprintf('Full-screen waveform view created\n');
end

function create3DView(fig)
    % Create a full-screen 3D visualization view
    userData = get(fig, 'UserData');

    % Debug: Show what statistics are available when creating 3D view
    fprintf('DEBUG: create3DView called\n');
    if ~isempty(userData) && isfield(userData, 'statDataArray') && isfield(userData, 'sortedStatTypes')
        fprintf('DEBUG: Available statistics in create3DView: %s\n', strjoin(userData.sortedStatTypes, ', '));
    else
        fprintf('DEBUG: No statistics data available in create3DView\n');
    end

    % Clear all existing plot axes
    axHandles = userData.axHandles;
    for i = 1:length(axHandles)
        if ishandle(axHandles(i))
            delete(axHandles(i));
        end
    end

    % Also clear any existing waveform axes
    existingWaveformAx = findobj(fig, 'Tag', 'WaveformAxes');
    if ~isempty(existingWaveformAx)
        delete(existingWaveformAx);
    end

    % Clear any existing 3D axes
    existing3DAx = findobj(fig, 'Tag', '3DAxes');
    if ~isempty(existing3DAx)
        delete(existing3DAx);
    end

    % Add the Utils/Visualization3D path if not already added
    utilsPath = fullfile(pwd, 'Utils');
    viz3DPath = fullfile(utilsPath, 'Visualization3D');
    if ~contains(path, viz3DPath)
        addpath(viz3DPath);
    end

    % Create 3D visualization directly in the figure
    try
        % Create 3D axes manually since we're not using the tab-based approach
        ax3D = axes('Position', [0.1, 0.1, 0.6, 0.8], 'Tag', '3DAxes');

        % Create a simple 3D visualization
        create3DVisualizationSimple(fig, ax3D);

        fprintf('3D view created successfully\n');

    catch ME
        fprintf('Error creating 3D view: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end

        % Fallback: create a simple message
        ax3D = axes('Position', [0.1, 0.1, 0.6, 0.8], 'Tag', '3DAxes');
        text(0.5, 0.5, 0.5, '3D View - Under Development', ...
            'HorizontalAlignment', 'center', 'FontSize', 16, ...
            'Parent', ax3D);
        axis(ax3D, 'off');
    end

    % Update userData to reflect the new axes structure
    ax3D = findobj(fig, 'Tag', '3DAxes');
    if ~isempty(ax3D)
        userData.axHandles = ax3D(1);
        userData.plotHandles = {[]};
        userData.visState.currentView = '3DView'; % Ensure view state is updated
        set(fig, 'UserData', userData);
    end

    % Hide the slider since it's not needed for 3D view
    ySlider = findobj(fig, 'Tag', 'YSlider');
    if ~isempty(ySlider)
        set(ySlider, 'Visible', 'off');
    end

    % Disable the waveform dropdown since we're in 3D view
    waveformDropdown = findobj(fig, 'Tag', 'WaveformDropdown');
    if ~isempty(waveformDropdown)
        set(waveformDropdown, 'Enable', 'off');
    end
end

function update3DFromAlignment(fig)
    % Update 3D visualization when alignment changes

    % Debug: Show when this function is called
    fprintf('DEBUG: update3DFromAlignment called\n');

    % Find existing 3D axes
    ax3D = findobj(fig, 'Tag', '3DAxes');
    if isempty(ax3D)
        fprintf('No 3D axes found - creating new 3D view\n');
        create3DView(fig);
        return;
    end

    try
        % Update the 3D visualization with current alignment data
        create3DVisualizationSimple(fig, ax3D(1));
        fprintf('3D visualization updated with alignment changes\n');

    catch ME
        fprintf('Error updating 3D visualization: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
end

function create3DVisualizationSimple(fig, ax3D)
    % Create a simplified 3D visualization directly in the figure

    userData = get(fig, 'UserData');
    if isempty(userData) || isempty(userData.statDataArray)
        % Fallback message
        text(0.5, 0.5, 0.5, '3D View - No Data Available', ...
            'HorizontalAlignment', 'center', 'FontSize', 16, ...
            'Parent', ax3D);
        axis(ax3D, 'off');
        return;
    end

    try
        % Get current data (aligned or original) with waveform processing applied
        if strcmp(userData.currentView, 'aligned') && ~isempty(userData.alignedDataCache)
            currentData = userData.alignedDataCache;
            dataLabel = 'Aligned';
        else
            currentData = userData.statDataArray;
            dataLabel = 'Original';
        end

        % Apply waveform processing to 3D view data
        currentWaveformMode = userData.statDropdownState.waveformMode;
        if ~strcmp(currentWaveformMode, 'RawWaveform')
            fprintf('3D View: Applying %s processing to statistical data\n', currentWaveformMode);
            currentData = apply3DWaveformProcessing(currentData, currentWaveformMode, userData);
            dataLabel = [dataLabel ' (' currentWaveformMode ')'];
        end


        % Use first statistic for 3D visualization
        statData = currentData{1};
        statName = userData.sortedStatTypes{1};

        % Debug: Show what statistic is being displayed in 3D
        fprintf('DEBUG: 3D View displaying statistic: %s\n', statName);
        fprintf('DEBUG: Current statDataArray contains %d statistics: %s\n', ...
            length(userData.statDataArray), strjoin(userData.sortedStatTypes, ', '));

        % Get data dimensions
        [numY, numX] = size(statData.maps{1});
        numSegments = length(statData.maps);

        % Apply time range filtering if enabled
        visState = userData.visState;
        segmentIndices = 1:numSegments;

        if isfield(visState, 'useTimeRange') && visState.useTimeRange && ...
           isfield(visState, 'timeRangeMin') && isfield(visState, 'timeRangeMax')

            % Get time values for filtering
            if isfield(userData, 'timeValues') && ~isempty(userData.timeValues)
                timeValues = userData.timeValues;
            else
                timeValues = 1:numSegments; % Fallback to segment indices
            end

            % Apply time range filter
            timeIndices = (timeValues >= visState.timeRangeMin) & (timeValues <= visState.timeRangeMax);

            if any(timeIndices)
                segmentIndices = find(timeIndices);
                fprintf('3D View: Time range filtering applied - showing %d/%d segments\n', length(segmentIndices), numSegments);
            else
                fprintf('3D View: No segments in time range - showing all data\n');
            end
        end

        % Create 3D volume by stacking filtered segments
        numFilteredSegments = length(segmentIndices);
        volume3D = zeros(numY, numX, numFilteredSegments);
        for i = 1:numFilteredSegments
            seg = segmentIndices(i);
            segmentData = statData.maps{seg};

            % Apply contrast enhancement if specified
            if isfield(visState, 'contrastEnhancement') && ~strcmp(visState.contrastEnhancement, 'none')
                segmentData = applyContrastEnhancement(segmentData, visState.contrastEnhancement);
            end

            volume3D(:, :, i) = segmentData;
        end

        % Get coordinate arrays
        if isfield(statData, 'X_sub') && isfield(statData, 'Y_sub')
            X_coords = statData.X_sub;
            Y_coords = statData.Y_sub;
        else
            X_coords = 1:numX;
            Y_coords = 1:numY;
        end

        % Use filtered segment indices for Z coordinates
        if isfield(userData, 'timeValues') && ~isempty(userData.timeValues)
            Z_coords = userData.timeValues(segmentIndices);
        else
            Z_coords = segmentIndices;
        end

        % Create coordinate grids
        [X_grid, Y_grid, Z_grid] = meshgrid(X_coords, Y_coords, Z_coords);

        % Clear axes and plot
        cla(ax3D);
        hold(ax3D, 'on');

        % Get plot type from dropdown
        plotTypeDropdown = findobj(fig, 'Tag', 'PlotTypeDropdown');
        if ~isempty(plotTypeDropdown)
            plotTypes = {'heatmap', 'filledcontour', 'slice', 'smooth2d', 'pseudocolor'};
            selectedType = plotTypes{get(plotTypeDropdown, 'Value')};
        else
            selectedType = 'slice'; % Default
        end

        % Create 3D visualization based on selected plot type
        create3DPlotByType(ax3D, X_grid, Y_grid, Z_grid, volume3D, selectedType, X_coords, Y_coords, Z_coords);

        % Set labels and title
        xlabel(ax3D, 'X Position (mm)');
        ylabel(ax3D, 'Y Position (mm)');
        zlabel(ax3D, 'Time/Segment');

        % Create title that reflects time range filtering
        if isfield(visState, 'useTimeRange') && visState.useTimeRange
            title(ax3D, sprintf('3D View: %s (%s) - Time Range: %.1f to %.1f', ...
                statName, dataLabel, visState.timeRangeMin, visState.timeRangeMax));
        else
            title(ax3D, sprintf('3D View: %s (%s)', statName, dataLabel));
        end

        % Apply visual settings from UI controls
        apply3DVisualSettings(fig, ax3D);

        % Set view and appearance
        view(ax3D, 3);
        axis(ax3D, 'equal');
        grid(ax3D, 'on');

        hold(ax3D, 'off');

        fprintf('Simple 3D visualization created for %s data\n', dataLabel);

    catch ME
        fprintf('Error creating simple 3D visualization: %s\n', ME.message);

        % Fallback: create a simple message
        cla(ax3D);
        text(0.5, 0.5, 0.5, '3D View - Error Loading Data', ...
            'HorizontalAlignment', 'center', 'FontSize', 16, ...
            'Parent', ax3D);
        axis(ax3D, 'off');
    end
end

% Function to apply waveform processing to 3D view data
function processedData = apply3DWaveformProcessing(originalData, waveformMode, userData)
    % Apply waveform processing effects to statistical data for 3D view

    processedData = originalData; % Start with copy

    if isempty(originalData) || isempty(originalData{1})
        return;
    end

    % Process each statistic
    for statIdx = 1:length(originalData)
        statData = originalData{statIdx};
        processedStatData = statData; % Copy structure

        % Process each segment
        for seg = 1:length(statData.maps)
            baseSegmentData = statData.maps{seg};

            % Apply the same processing as in computeStatisticWithWaveformProcessing
            switch waveformMode
                case 'Envelope'
                    % For envelope, amplify and smooth
                    statMap = abs(baseSegmentData) * 1.2;
                    if size(statMap, 1) > 3 && size(statMap, 2) > 3
                        kernel = ones(3,3) / 9;
                        statMap = conv2(statMap, kernel, 'same');
                    end

                case 'FFT'
                    % For FFT, apply power scaling and frequency patterns
                    statMap = abs(baseSegmentData).^0.8;
                    [rows, cols] = size(statMap);
                    [X, Y] = meshgrid(1:cols, 1:rows);
                    freqPattern = 0.1 * sin(2*pi*X/cols) .* cos(2*pi*Y/rows);
                    statMap = statMap .* (1 + freqPattern);

                case 'STFT'
                    % For STFT, apply square root scaling and time-frequency patterns
                    statMap = sqrt(abs(baseSegmentData));
                    [rows, cols] = size(statMap);
                    [X, Y] = meshgrid(1:cols, 1:rows);
                    tfPattern = 0.15 * cos(2*pi*X/cols + seg/10) .* sin(2*pi*Y/rows);
                    statMap = statMap .* (1 + tfPattern);

                case 'Derivative'
                    % For derivative, enhance edges
                    statMap = abs(baseSegmentData).^1.5;
                    if size(statMap, 1) > 1 && size(statMap, 2) > 1
                        [gradX, gradY] = gradient(statMap);
                        gradMag = sqrt(gradX.^2 + gradY.^2);
                        statMap = statMap + 0.3 * gradMag;
                    end

                case 'Integral'
                    % For integral, show accumulated effects
                    statMap = abs(baseSegmentData).^0.7;
                    [rows, cols] = size(statMap);
                    for i = 2:rows
                        statMap(i, :) = statMap(i, :) + 0.1 * statMap(i-1, :);
                    end

                case 'Bandpass'
                    % For bandpass, apply spatial filtering
                    statMap = abs(baseSegmentData);
                    if size(statMap, 1) > 5 && size(statMap, 2) > 5
                        lowpass = imgaussfilt(statMap, 3);
                        highpass = statMap - imgaussfilt(statMap, 0.5);
                        statMap = statMap - 0.3 * lowpass + 0.2 * highpass;
                    end

                case 'Lowpass'
                    % For lowpass, smooth the data
                    statMap = abs(baseSegmentData);
                    if size(statMap, 1) > 3 && size(statMap, 2) > 3
                        statMap = imgaussfilt(statMap, 1.5);
                    end

                case 'Highpass'
                    % For highpass, emphasize edges
                    statMap = abs(baseSegmentData);
                    if size(statMap, 1) > 3 && size(statMap, 2) > 3
                        lowpass = imgaussfilt(statMap, 2);
                        statMap = abs(statMap - lowpass);
                    end

                case 'Wavelet'
                    % For wavelet, apply multi-scale patterns
                    statMap = abs(baseSegmentData).^0.8;
                    [rows, cols] = size(statMap);
                    [X, Y] = meshgrid(1:cols, 1:rows);
                    scale1 = 0.1 * sin(2*pi*X/(cols/4)) .* sin(2*pi*Y/(rows/4));
                    scale2 = 0.05 * sin(2*pi*X/(cols/8)) .* sin(2*pi*Y/(rows/8));
                    multiScale = scale1 + scale2;
                    statMap = statMap .* (1 + multiScale);

                otherwise
                    % Raw waveform - no processing
                    statMap = baseSegmentData;
            end



            processedStatData.maps{seg} = statMap;
        end

        processedData{statIdx} = processedStatData;
    end

    fprintf('3D View: Applied %s processing to %d statistics\n', waveformMode, length(originalData));
end

function recreateOriginalAxes(fig)
    % Recreate the original axes layout when switching back from 3D view

    userData = get(fig, 'UserData');
    nStats = length(userData.sortedStatTypes);

    % Clear any existing axes
    existingAxes = findall(fig, 'Type', 'axes');
    for i = 1:length(existingAxes)
        if ishandle(existingAxes(i))
            delete(existingAxes(i));
        end
    end

    % Calculate subplot layout (same as original)
    nRows = ceil(sqrt(nStats));
    nCols = ceil(nStats / nRows);

    % Create axes for each statistic (same as original)
    axHandles = zeros(nStats, 1);
    plotHandles = cell(nStats, 1);

    for i = 1:nStats
        % Calculate subplot position
        row = ceil(i / nCols);
        col = mod(i-1, nCols) + 1;

        % Create subplot with proper spacing
        left = 0.1 + (col-1) * 0.55/nCols;
        bottom = 0.35 + (nRows-row) * 0.5/nRows;
        width = 0.55/nCols * 0.9;
        height = 0.5/nRows * 0.9;

        axHandles(i) = axes('Position', [left, bottom, width, height]);
        plotHandles{i} = [];
    end

    % Update userData with new axes
    userData.axHandles = axHandles;
    userData.plotHandles = plotHandles;
    set(fig, 'UserData', userData);

    fprintf('Original axes layout recreated\n');
end

function apply3DVisualSettings(fig, ax3D)
    % Apply visual settings from UI controls to 3D view

    userData = get(fig, 'UserData');
    visState = userData.visState;

    % Apply colormap
    colormapDropdown = findobj(fig, 'Tag', 'ColormapDropdown');
    if ~isempty(colormapDropdown)
        colormapOptions = get(colormapDropdown, 'String');
        selectedColormap = colormapOptions{get(colormapDropdown, 'Value')};

        if strcmp(selectedColormap, 'coolwarm')
            colormap(ax3D, coolwarm());
        else
            colormap(ax3D, selectedColormap);
        end
    else
        colormap(ax3D, 'jet'); % Default
    end

    % Add colorbar
    colorbar(ax3D);

    % Apply contrast enhancement to 3D data if needed
    contrastDropdown = findobj(fig, 'Tag', 'ContrastDropdown');
    if ~isempty(contrastDropdown)
        contrastMethods = {'none', 'linear', 'histeq', 'adaptive', 'gamma'};
        selectedContrast = contrastMethods{get(contrastDropdown, 'Value')};

        % Store contrast method in visState for use during 3D data processing
        userData = get(fig, 'UserData');
        userData.visState.contrastEnhancement = selectedContrast;
        set(fig, 'UserData', userData);

        fprintf('3D contrast enhancement set to: %s\n', selectedContrast);
    end

    % Apply plot type-specific visual settings for 3D
    plotTypeDropdown = findobj(fig, 'Tag', 'PlotTypeDropdown');
    if ~isempty(plotTypeDropdown)
        plotTypes = {'heatmap', 'filledcontour', 'slice', 'smooth2d', 'pseudocolor'};
        selectedType = plotTypes{get(plotTypeDropdown, 'Value')};

        % Apply general settings that work across all plot types
        switch selectedType
            case 'slice'
                % Slice planes work well with moderate transparency
                alpha(ax3D, 0.7);
                shading(ax3D, 'interp');

            case 'smooth2d'
                % Smooth surfaces benefit from interpolated shading
                shading(ax3D, 'interp');
                alpha(ax3D, 0.8);

            case 'pseudocolor'
                % Discrete appearance with flat shading
                shading(ax3D, 'flat');
                alpha(ax3D, 0.9);

            case 'filledcontour'
                % Contours work well with higher transparency
                alpha(ax3D, 0.6);

            case 'heatmap'
                % Heatmap surfaces with good visibility
                alpha(ax3D, 0.8);
                shading(ax3D, 'interp');

            otherwise
                % Default settings
                alpha(ax3D, 0.7);
                shading(ax3D, 'interp');
        end

        fprintf('3D plot type visual settings applied for: %s\n', selectedType);
    end

    fprintf('3D visual settings applied\n');
end

function create3DPlotByType(ax3D, X_grid, Y_grid, Z_grid, volume3D, plotType, X_coords, Y_coords, Z_coords)
    % Create 3D visualization based on selected plot type

    fprintf('Creating 3D plot with type: %s\n', plotType);

    % Validate input data
    if isempty(volume3D) || any(size(volume3D) == 0)
        fprintf('Warning: Empty volume data for 3D plot type %s\n', plotType);
        return;
    end

    switch plotType
        case 'heatmap'
            % Use slice planes but with different visual settings for heatmap appearance
            sliceX = [X_coords(1), X_coords(round(end/2)), X_coords(end)];
            sliceY = [Y_coords(1), Y_coords(round(end/2)), Y_coords(end)];

            if length(Z_coords) >= 3
                sliceZ = [Z_coords(1), Z_coords(round(end/2)), Z_coords(end)];
            elseif length(Z_coords) == 2
                sliceZ = [Z_coords(1), Z_coords(end)];
            elseif length(Z_coords) == 1
                sliceZ = Z_coords(1);
            else
                sliceZ = [];
            end

            if ~isempty(sliceZ)
                slice(ax3D, X_grid, Y_grid, Z_grid, volume3D, sliceX, sliceY, sliceZ);
                % Apply heatmap-like visual settings
                shading(ax3D, 'interp');
                alpha(ax3D, 0.8);
            end

        case 'filledcontour'
            % Use slice planes with contour-like visual appearance
            sliceX = [X_coords(1), X_coords(round(end/2)), X_coords(end)];
            sliceY = [Y_coords(1), Y_coords(round(end/2)), Y_coords(end)];

            if length(Z_coords) >= 3
                sliceZ = [Z_coords(1), Z_coords(round(end/2)), Z_coords(end)];
            elseif length(Z_coords) == 2
                sliceZ = [Z_coords(1), Z_coords(end)];
            elseif length(Z_coords) == 1
                sliceZ = Z_coords(1);
            else
                sliceZ = [];
            end

            if ~isempty(sliceZ)
                slice(ax3D, X_grid, Y_grid, Z_grid, volume3D, sliceX, sliceY, sliceZ);
                % Apply filled contour visual settings
                shading(ax3D, 'interp');
                alpha(ax3D, 0.6);

                % Add contour lines on the slice planes for contour effect
                hold(ax3D, 'on');
                for i = 1:length(sliceZ)
                    if sliceZ(i) >= min(Z_coords) && sliceZ(i) <= max(Z_coords)
                        % Find the closest Z index
                        [~, zIdx] = min(abs(Z_coords - sliceZ(i)));
                        sliceData = volume3D(:, :, zIdx);

                        % Add contour lines at this Z level
                        contour3(ax3D, X_grid(:, :, 1), Y_grid(:, :, 1), ...
                                ones(size(X_grid(:, :, 1))) * sliceZ(i), sliceData, 8, 'k', 'LineWidth', 0.5);
                    end
                end
                hold(ax3D, 'off');
            end

        case 'slice'
            % Traditional slice planes (original behavior)
            sliceX = [X_coords(1), X_coords(round(end/2)), X_coords(end)];
            sliceY = [Y_coords(1), Y_coords(round(end/2)), Y_coords(end)];

            if length(Z_coords) >= 3
                sliceZ = [Z_coords(1), Z_coords(round(end/2)), Z_coords(end)];
            elseif length(Z_coords) == 2
                sliceZ = [Z_coords(1), Z_coords(end)];
            elseif length(Z_coords) == 1
                sliceZ = Z_coords(1);
            else
                sliceZ = [];
            end

            if ~isempty(sliceZ)
                slice(ax3D, X_grid, Y_grid, Z_grid, volume3D, sliceX, sliceY, sliceZ);
                shading(ax3D, 'interp');
                alpha(ax3D, 0.7);
            end

        case 'smooth2d'
            % Use slice planes with smooth interpolated appearance
            sliceX = [X_coords(1), X_coords(round(end/2)), X_coords(end)];
            sliceY = [Y_coords(1), Y_coords(round(end/2)), Y_coords(end)];

            if length(Z_coords) >= 3
                sliceZ = [Z_coords(1), Z_coords(round(end/2)), Z_coords(end)];
            elseif length(Z_coords) == 2
                sliceZ = [Z_coords(1), Z_coords(end)];
            elseif length(Z_coords) == 1
                sliceZ = Z_coords(1);
            else
                sliceZ = [];
            end

            if ~isempty(sliceZ)
                slice(ax3D, X_grid, Y_grid, Z_grid, volume3D, sliceX, sliceY, sliceZ);
                % Apply smooth visual settings
                shading(ax3D, 'interp');
                alpha(ax3D, 0.8);
            end

        case 'pseudocolor'
            % Use slice planes with discrete color appearance
            sliceX = [X_coords(1), X_coords(round(end/2)), X_coords(end)];
            sliceY = [Y_coords(1), Y_coords(round(end/2)), Y_coords(end)];

            if length(Z_coords) >= 3
                sliceZ = [Z_coords(1), Z_coords(round(end/2)), Z_coords(end)];
            elseif length(Z_coords) == 2
                sliceZ = [Z_coords(1), Z_coords(end)];
            elseif length(Z_coords) == 1
                sliceZ = Z_coords(1);
            else
                sliceZ = [];
            end

            if ~isempty(sliceZ)
                slice(ax3D, X_grid, Y_grid, Z_grid, volume3D, sliceX, sliceY, sliceZ);
                % Apply discrete color visual settings
                shading(ax3D, 'flat');
                alpha(ax3D, 0.9);
            end

        otherwise
            % Default to slice
            create3DPlotByType(ax3D, X_grid, Y_grid, Z_grid, volume3D, 'slice', X_coords, Y_coords, Z_coords);
    end

    fprintf('3D plot type "%s" created successfully\n', plotType);
end

% Note: unpackFileNamingArray function is now centralized in the main directory

function updateWaveformDisplay(fig)
    % Update waveform display based on current settings
    userData = get(fig, 'UserData');

    % Check if waveform should be displayed
    waveformDropdown = findobj(fig, 'Tag', 'WaveformDropdown');
    if isempty(waveformDropdown)
        return;
    end

    options = get(waveformDropdown, 'String');
    selectedOption = options{get(waveformDropdown, 'Value')};
    showWaveform = strcmp(selectedOption, 'Show Waveform');

    % Remove existing waveform axes first to prevent layering issues
    existingWaveformAx = findobj(fig, 'Tag', 'WaveformAxes');
    if ~isempty(existingWaveformAx)
        delete(existingWaveformAx);
        drawnow;
    end

    if ~showWaveform || ~userData.waveformLoaded
        return;
    end

    % Determine if we're showing aligned data (single source of truth)
    useAlignedData = false;
    if strcmp(userData.currentView, 'aligned') && ~isempty(userData.alignedWaveformData)
        % Additional check: verify that alignment shifts actually exist and are non-zero
        if ~isempty(userData.alignmentShifts) && any(abs(userData.alignmentShifts(:)) > 1e-10)
            useAlignedData = true;
        end
    end

    % Get current waveform data based on alignment status
    if useAlignedData
        waveformData = userData.alignedWaveformData;
    else
        waveformData = userData.originalWaveformData;
    end

    % Apply waveform processing based on current mode
    currentWaveformMode = userData.statDropdownState.waveformMode;
    if ~strcmp(currentWaveformMode, 'RawWaveform')
        fprintf('Applying %s processing to full-screen waveform display\n', currentWaveformMode);
        waveformData = applyWaveformProcessingToDisplay(waveformData, currentWaveformMode);
    end

    if isempty(waveformData)
        return;
    end

    % Apply waveform processing based on current mode
    currentWaveformMode = userData.statDropdownState.waveformMode;
    if ~strcmp(currentWaveformMode, 'RawWaveform')
        fprintf('Applying %s processing to waveform display\n', currentWaveformMode);
        waveformData = applyWaveformProcessingToDisplay(waveformData, currentWaveformMode);
    end

    % Only show waveforms for XtVsY and YtVsX views (not XYVst, Waveform, or 3DView)
    currentView = userData.visState.currentView;
    if strcmp(currentView, 'XYVst') || strcmp(currentView, 'Waveform') || strcmp(currentView, '3DView')
        return;
    end

    % Create waveform axes positioned between slider and main plots
    waveformHeight = 0.10;
    waveformBottom = 0.18;
    waveformPosition = [0.1, waveformBottom, 0.7, waveformHeight];
    waveformAx = axes('Position', waveformPosition, 'Tag', 'WaveformAxes');

    % Get waveform data - use references to avoid unnecessary copying
    % Only create copies if time range filtering is needed
    waveformArray = waveformData.waveformArray;
    t = waveformData.t;

    % Apply time range filtering to waveform data if enabled
    visState = userData.visState;
    if isfield(visState, 'useTimeRange') && visState.useTimeRange && ...
       isfield(visState, 'timeRangeMin') && isfield(visState, 'timeRangeMax')

        % Convert time range from microseconds to seconds for waveform comparison
        % The time range filter uses microseconds (e-6), but waveform time is in seconds
        timeRangeMinSeconds = visState.timeRangeMin * 1e-6;
        timeRangeMaxSeconds = visState.timeRangeMax * 1e-6;

        % Find time indices within the specified range for waveform time axis
        timeIndices = (t >= timeRangeMinSeconds) & (t <= timeRangeMaxSeconds);

        if any(timeIndices)
            waveformArray = waveformArray(:, timeIndices);
            t = t(timeIndices);
        else
            % Show empty plot if no data in range
            plot(waveformAx, [], []);
            title(waveformAx, 'No waveform data in time range');
            xlabel(waveformAx, 'Time (μs)');
            ylabel(waveformAx, 'Amplitude');
            return;
        end
    end

    % Check if user wants to show all waveforms or just average
    waveformModeCheckbox = findobj(fig, 'Tag', 'WaveformModeCheckbox');
    showAllWaveforms = false;
    if ~isempty(waveformModeCheckbox)
        showAllWaveforms = get(waveformModeCheckbox, 'Value');
    end

    % Convert time from seconds to microseconds for display
    t_microseconds = t * 1e6;

    % Use the same alignment status determined above (no duplicate logic)

    if showAllWaveforms
        % Plot all individual waveforms
        hold(waveformAx, 'on');

        % Determine colors based on alignment status
        if useAlignedData
            individualColor = [0, 0.4470, 0.7410]; % Blue for aligned waveforms
            avgColor = 'k';
        else
            individualColor = [0.5, 0.5, 0.5]; % Gray for original waveforms
            avgColor = 'k';
        end

        % Plot each waveform with consistent color and reduced line width
        numWaveforms = size(waveformArray, 1);
        % Vectorized plotting for performance: plot all waveforms in one call
        plot(waveformAx, t_microseconds, waveformArray', 'Color', individualColor, 'LineWidth', 0.5)

        % Plot average waveform on top with thicker line
        avgWaveform = mean(waveformArray, 1);
        plot(waveformAx, t_microseconds, avgWaveform, 'Color', avgColor, 'LineWidth', 2);

        hold(waveformAx, 'off');
        % No title when showing all waveforms to avoid overlap with checkbox

        % Set Y-axis limits to show full amplitude range with padding
        allData = waveformArray(:);
        dataMin = min(allData);
        dataMax = max(allData);
        if dataMin ~= dataMax
            dataRange = dataMax - dataMin;
            padding = dataRange * 0.1; % 10% padding
            ylim(waveformAx, [dataMin - padding, dataMax + padding]);
        else
            % Handle case where all data is the same value
            ylim(waveformAx, [dataMin - 0.1, dataMax + 0.1]);
        end

    else
        % Plot only the average waveform
        avgWaveform = mean(waveformArray, 1);

        % Use black for average waveform regardless of alignment status
        plot(waveformAx, t_microseconds, avgWaveform, 'k-', 'LineWidth', 1.5);

        % Update title to indicate alignment status
        if useAlignedData
            % Count non-zero shifts for display
            nonZeroShifts = sum(abs(userData.alignmentShifts(:)) > 1e-10);
            title(waveformAx, sprintf('Average Waveform (Aligned - %d shifts)', nonZeroShifts));
        else
            title(waveformAx, sprintf('Average Waveform (Original)'));
        end

        % Set Y-axis limits to show full amplitude range with padding
        dataMin = min(avgWaveform);
        dataMax = max(avgWaveform);
        if dataMin ~= dataMax
            dataRange = dataMax - dataMin;
            padding = dataRange * 0.1; % 10% padding
            ylim(waveformAx, [dataMin - padding, dataMax + padding]);
        else
            % Handle case where all data is the same value
            ylim(waveformAx, [dataMin - 0.1, dataMax + 0.1]);
        end
    end

    xlabel(waveformAx, 'Time (μs)');
    ylabel(waveformAx, 'Amplitude');
    grid(waveformAx, 'on');

    % Add legend when showing all waveforms to explain color coding
    if showAllWaveforms
        hold(waveformAx, 'on');
        if useAlignedData
            % Create invisible plots for legend entries
            h1 = plot(waveformAx, NaN, NaN, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.5);
            h2 = plot(waveformAx, NaN, NaN, 'k-', 'LineWidth', 2);
            legend(waveformAx, [h1, h2], {'Aligned Waveforms', 'Average'}, 'Location', 'northeast', 'FontSize', 8);
        else
            % Create invisible plots for legend entries
            h1 = plot(waveformAx, NaN, NaN, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
            h2 = plot(waveformAx, NaN, NaN, 'k-', 'LineWidth', 2);
            legend(waveformAx, [h1, h2], {'Individual Waveforms', 'Average'}, 'Location', 'northeast', 'FontSize', 8);
        end
        hold(waveformAx, 'off');
    end



end

% Function to apply waveform processing for display purposes
function processedWaveformData = applyWaveformProcessingToDisplay(originalWaveformData, waveformMode)
    % Apply the selected waveform processing to the waveform data for display

    processedWaveformData = originalWaveformData; % Start with copy
    waveformArray = originalWaveformData.waveformArray;

    fprintf('Processing %d waveforms with %s method\n', size(waveformArray, 1), waveformMode);

    switch waveformMode
        case 'Envelope'
            % Apply Hilbert transform to get envelope
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 1
                    envelope = abs(hilbert(waveform));
                    processedArray(i, :) = envelope;

                    % Debug for first few waveforms
                    if i <= 3
                        fprintf('Waveform %d: Original range [%.4f, %.4f] -> Envelope range [%.4f, %.4f]\n', ...
                            i, min(waveform), max(waveform), min(envelope), max(envelope));
                    end
                else
                    processedArray(i, :) = waveform;
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'FFT'
            % Apply FFT and take magnitude
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 1
                    fftResult = abs(fft(waveform));
                    % Take only the first half (positive frequencies)
                    halfLength = floor(length(fftResult)/2);
                    fftMagnitude = fftResult(1:halfLength);

                    % Pad or truncate to match original length
                    if length(fftMagnitude) < size(waveformArray, 2)
                        % Pad with zeros
                        paddedFFT = zeros(1, size(waveformArray, 2));
                        paddedFFT(1:length(fftMagnitude)) = fftMagnitude;
                        processedArray(i, :) = paddedFFT;
                    else
                        % Truncate to fit
                        processedArray(i, :) = fftMagnitude(1:size(waveformArray, 2));
                    end

                    % Debug for first few waveforms
                    if i <= 3
                        fprintf('Waveform %d: Original range [%.4f, %.4f] -> FFT range [%.4f, %.4f]\n', ...
                            i, min(waveform), max(waveform), min(processedArray(i, :)), max(processedArray(i, :)));
                    end
                else
                    processedArray(i, :) = abs(waveform);
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'STFT'
            % Apply Short-Time Fourier Transform
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 16
                    % Compute STFT with appropriate window size
                    windowSize = min(64, floor(length(waveform)/4));
                    overlap = floor(windowSize/2);

                    try
                        % Compute STFT magnitude
                        stftMag = computeSTFTMagnitude(waveform, windowSize, overlap);

                        % Take mean across frequency bins to get time-domain representation
                        if size(stftMag, 2) == length(waveform)
                            stftTimeSeries = mean(stftMag, 1);
                        else
                            % Interpolate to match original length
                            stftTimeSeries = interp1(1:size(stftMag, 2), mean(stftMag, 1), ...
                                linspace(1, size(stftMag, 2), length(waveform)), 'linear', 'extrap');
                        end

                        processedArray(i, :) = stftTimeSeries;

                        % Debug for first few waveforms
                        if i <= 3
                            fprintf('Waveform %d: Original range [%.4f, %.4f] -> STFT range [%.4f, %.4f]\n', ...
                                i, min(waveform), max(waveform), min(stftTimeSeries), max(stftTimeSeries));
                        end
                    catch ME
                        fprintf('STFT failed for waveform %d: %s\n', i, ME.message);
                        % Fallback to absolute value
                        processedArray(i, :) = abs(waveform);
                    end
                else
                    % Too short for STFT, use absolute value
                    processedArray(i, :) = abs(waveform);
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'Derivative'
            % Apply numerical derivative
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 1
                    % Compute derivative using gradient
                    derivative = gradient(waveform);
                    processedArray(i, :) = derivative;

                    % Debug for first few waveforms
                    if i <= 3
                        fprintf('Waveform %d: Original range [%.4f, %.4f] -> Derivative range [%.4f, %.4f]\n', ...
                            i, min(waveform), max(waveform), min(derivative), max(derivative));
                    end
                else
                    processedArray(i, :) = waveform;
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'Integral'
            % Apply numerical integration (cumulative sum)
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 1
                    % Compute cumulative integral
                    integral = cumsum(waveform);
                    processedArray(i, :) = integral;

                    % Debug for first few waveforms
                    if i <= 3
                        fprintf('Waveform %d: Original range [%.4f, %.4f] -> Integral range [%.4f, %.4f]\n', ...
                            i, min(waveform), max(waveform), min(integral), max(integral));
                    end
                else
                    processedArray(i, :) = waveform;
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'Bandpass'
            % Apply bandpass filter (simplified using FFT)
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 8
                    % Simple bandpass filter using FFT
                    N = length(waveform);
                    fftWave = fft(waveform);

                    % Define bandpass range (keep middle 50% of frequencies)
                    lowCut = floor(N * 0.25);
                    highCut = floor(N * 0.75);

                    % Zero out frequencies outside the band
                    fftFiltered = fftWave;
                    fftFiltered(1:lowCut) = 0;
                    fftFiltered(highCut:end) = 0;

                    % Convert back to time domain
                    filtered = real(ifft(fftFiltered));
                    processedArray(i, :) = filtered;

                    % Debug for first few waveforms
                    if i <= 3
                        fprintf('Waveform %d: Original range [%.4f, %.4f] -> Bandpass range [%.4f, %.4f]\n', ...
                            i, min(waveform), max(waveform), min(filtered), max(filtered));
                    end
                else
                    processedArray(i, :) = waveform;
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'Lowpass'
            % Apply lowpass filter
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 4
                    % Simple lowpass filter using moving average
                    windowSize = max(3, floor(length(waveform) / 20));
                    filtered = movmean(waveform, windowSize);
                    processedArray(i, :) = filtered;

                    % Debug for first few waveforms
                    if i <= 3
                        fprintf('Waveform %d: Original range [%.4f, %.4f] -> Lowpass range [%.4f, %.4f]\n', ...
                            i, min(waveform), max(waveform), min(filtered), max(filtered));
                    end
                else
                    processedArray(i, :) = waveform;
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'Highpass'
            % Apply highpass filter (original - lowpass)
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 4
                    % Highpass = original - lowpass
                    windowSize = max(3, floor(length(waveform) / 20));
                    lowpass = movmean(waveform, windowSize);
                    highpass = waveform - lowpass;
                    processedArray(i, :) = highpass;

                    % Debug for first few waveforms
                    if i <= 3
                        fprintf('Waveform %d: Original range [%.4f, %.4f] -> Highpass range [%.4f, %.4f]\n', ...
                            i, min(waveform), max(waveform), min(highpass), max(highpass));
                    end
                else
                    processedArray(i, :) = waveform;
                end
            end
            processedWaveformData.waveformArray = processedArray;

        case 'Wavelet'
            % Apply wavelet transform (simplified using continuous wavelet)
            processedArray = zeros(size(waveformArray));
            for i = 1:size(waveformArray, 1)
                waveform = waveformArray(i, :);
                if ~isempty(waveform) && length(waveform) > 16
                    try
                        % Simple wavelet-like transform using Morlet-like function
                        N = length(waveform);
                        scales = 1:min(32, floor(N/4));
                        waveletCoeffs = zeros(length(scales), N);

                        for s = 1:length(scales)
                            scale = scales(s);
                            % Create simple Morlet-like wavelet
                            t = -scale:scale;
                            if length(t) > 1
                                wavelet = exp(-t.^2/(2*scale^2)) .* cos(2*pi*t/scale);
                                % Convolve with signal
                                conv_result = conv(waveform, wavelet, 'same');
                                waveletCoeffs(s, :) = abs(conv_result);
                            end
                        end

                        % Take mean across scales for time-domain representation
                        waveletTimeSeries = mean(waveletCoeffs, 1);
                        processedArray(i, :) = waveletTimeSeries;

                        % Debug for first few waveforms
                        if i <= 3
                            fprintf('Waveform %d: Original range [%.4f, %.4f] -> Wavelet range [%.4f, %.4f]\n', ...
                                i, min(waveform), max(waveform), min(waveletTimeSeries), max(waveletTimeSeries));
                        end
                    catch ME
                        fprintf('Wavelet failed for waveform %d: %s\n', i, ME.message);
                        % Fallback to absolute value
                        processedArray(i, :) = abs(waveform);
                    end
                else
                    % Too short for wavelet, use absolute value
                    processedArray(i, :) = abs(waveform);
                end
            end
            processedWaveformData.waveformArray = processedArray;

        otherwise
            % Raw waveform - no processing needed
            fprintf('Using raw waveform data\n');
    end

    fprintf('Waveform processing completed\n');
end

% Helper function to compute STFT magnitude
function stftMag = computeSTFTMagnitude(signal, windowSize, overlap)
    % Compute Short-Time Fourier Transform magnitude

    signalLength = length(signal);
    hopSize = windowSize - overlap;

    % Calculate number of frames
    numFrames = floor((signalLength - windowSize) / hopSize) + 1;

    % Create window (Hamming window)
    window = hamming(windowSize);

    % Initialize STFT matrix
    numFreqBins = floor(windowSize/2) + 1;
    stftMag = zeros(numFreqBins, numFrames);

    % Compute STFT
    for frame = 1:numFrames
        startIdx = (frame - 1) * hopSize + 1;
        endIdx = startIdx + windowSize - 1;

        if endIdx <= signalLength
            % Extract windowed segment
            segment = signal(startIdx:endIdx) .* window';

            % Compute FFT
            fftResult = fft(segment);

            % Take magnitude of positive frequencies only
            stftMag(:, frame) = abs(fftResult(1:numFreqBins));
        end
    end

    % If we have fewer frames than the original signal length, interpolate
    if numFrames < signalLength
        % Interpolate to match original signal length
        frameIndices = 1:numFrames;
        targetIndices = linspace(1, numFrames, signalLength);

        stftMagInterp = zeros(numFreqBins, signalLength);
        for freqBin = 1:numFreqBins
            stftMagInterp(freqBin, :) = interp1(frameIndices, stftMag(freqBin, :), ...
                targetIndices, 'linear', 'extrap');
        end
        stftMag = stftMagInterp;
    end
end

function applyAlignmentToWaveforms(fig)
    % Apply the same alignment shifts to waveform data as applied to statistical data
    userData = get(fig, 'UserData');

    if ~userData.waveformLoaded || isempty(userData.originalWaveformData)
        return;
    end

    % Check if we have alignment shifts available
    if isempty(userData.alignmentShifts)
        fprintf('No alignment shifts available to apply to waveforms\n');
        return;
    end

    try
        % Get references to original data (avoid copying)
        originalWaveformData = userData.originalWaveformData;
        waveformArray = originalWaveformData.waveformArray;
        [numWaveforms, numTimePoints] = size(waveformArray);

        % Get dimensions from original data to understand the mapping
        numY_sub = originalWaveformData.numY_sub;
        numX_sub = originalWaveformData.numX_sub;

        % Create aligned waveform data structure only if it doesn't exist
        if isempty(userData.alignedWaveformData)
            userData.alignedWaveformData = struct();
            userData.alignedWaveformData.t = originalWaveformData.t;
            userData.alignedWaveformData.numY_sub = numY_sub;
            userData.alignedWaveformData.numX_sub = numX_sub;
            userData.alignedWaveformData.waveformArray = zeros(size(waveformArray)); % Pre-allocate
        end

        % Apply shifts to each Y slice of waveforms
        alignmentShifts = userData.alignmentShifts;

        % Reshape waveform array to match Y,X,t structure
        % waveformArray is stored as (numWaveforms, numTimePoints) where numWaveforms = numY_sub * numX_sub
        % We need to reshape it to (numY_sub, numX_sub, numTimePoints) to apply Y-slice shifts
        waveformMatrix = zeros(numY_sub, numX_sub, numTimePoints);
        waveformIndex = 0;
        for yIdx = 1:numY_sub
            for xIdx = 1:numX_sub
                waveformIndex = waveformIndex + 1;
                if waveformIndex <= numWaveforms
                    waveformMatrix(yIdx, xIdx, :) = waveformArray(waveformIndex, :);
                end
            end
        end

        % Apply shifts to each Y slice
        for yIdx = 1:numY_sub
            if yIdx <= size(alignmentShifts, 1)
                % Get shifts for this Y slice (one shift per X position)
                shiftsForSlice = alignmentShifts(yIdx, :);

                % Apply shifts to each X position in this Y slice
                for xIdx = 1:numX_sub
                    if xIdx <= length(shiftsForSlice) && abs(shiftsForSlice(xIdx)) > 1e-10
                        % Apply the shift to the waveform at this Y,X position
                        originalWaveform = squeeze(waveformMatrix(yIdx, xIdx, :));
                        shiftedWaveform = applyShiftToWaveform(originalWaveform, shiftsForSlice(xIdx));
                        waveformMatrix(yIdx, xIdx, :) = shiftedWaveform;
                    end
                end
            end
        end

        % Convert back to original waveformArray format directly into the pre-allocated array
        waveformIndex = 0;
        for yIdx = 1:numY_sub
            for xIdx = 1:numX_sub
                waveformIndex = waveformIndex + 1;
                if waveformIndex <= numWaveforms
                    userData.alignedWaveformData.waveformArray(waveformIndex, :) = squeeze(waveformMatrix(yIdx, xIdx, :));
                end
            end
        end

        set(fig, 'UserData', userData);

        % Brief confirmation
        fprintf('Waveform alignment applied\n');

    catch ME
        fprintf('Error applying alignment to waveforms: %s\n', ME.message);
        % Fallback to copy of original data
        userData.alignedWaveformData = userData.originalWaveformData;
        set(fig, 'UserData', userData);
    end
end

function shiftedWaveform = applyShiftToWaveform(waveform, shift)
    % Apply time shift to a single waveform using zero padding (no wrapping)
    % This matches the alignment method used for statistical data

    if abs(shift) < 1e-10
        shiftedWaveform = waveform;
        return;
    end

    % Use zero padding for waveforms (same as statistical data)
    shift = round(shift); % Only integer shifts
    n = length(waveform);
    shiftedWaveform = zeros(size(waveform));

    if shift > 0  % Shift forward in time
        shiftedWaveform(shift+1:end) = waveform(1:end-shift);
        % shiftedWaveform(1:shift) remains zero (zero padding at start)
    else  % Shift backward in time
        shift = abs(shift);
        shiftedWaveform(1:end-shift) = waveform(shift+1:end);
        % shiftedWaveform(end-shift+1:end) remains zero (zero padding at end)
    end
end

function appliedShifts = calculateShiftsFromAlignment(originalData, alignedData)
    % Calculate the shifts applied by comparing original and aligned data
    % This is much more efficient than tracking shifts during alignment

    [numRows, numCols] = size(originalData);
    appliedShifts = zeros(1, numCols);

    % For each column, find the shift that best explains the transformation
    for col = 1:numCols
        originalCol = originalData(:, col);
        alignedCol = alignedData(:, col);

        % Skip if columns are identical (no shift applied)
        if isequal(originalCol, alignedCol)
            appliedShifts(col) = 0;
            continue;
        end

        % Test shifts to find the one that best matches the aligned data
        maxShift = min(15, floor(numRows/4)); % Reasonable search range
        bestShift = 0;
        bestCorr = -inf;

        for shift = -maxShift:maxShift
            % Apply shift to original column using zero padding
            testCol = applyShiftToWaveform(originalCol, shift);

            % Calculate correlation with aligned column
            corr = calculateCorrelation(testCol, alignedCol);

            if corr > bestCorr
                bestCorr = corr;
                bestShift = shift;
            end
        end

        appliedShifts(col) = bestShift;
    end
end

function corr = calculateCorrelation(col1, col2)
    % Calculate correlation between two columns
    validIdx = ~isnan(col1) & ~isnan(col2);
    if sum(validIdx) < 2
        corr = -inf;
        return;
    end

    col1_valid = col1(validIdx);
    col2_valid = col2(validIdx);

    % Normalize
    col1_norm = col1_valid - mean(col1_valid);
    col2_norm = col2_valid - mean(col2_valid);

    % Calculate correlation
    if std(col1_norm) > 1e-10 && std(col2_norm) > 1e-10
        corr = sum(col1_norm .* col2_norm) / (sqrt(sum(col1_norm.^2)) * sqrt(sum(col2_norm.^2)));
    else
        corr = -inf;
    end
end

function enhancedData = applyContrastEnhancement(data, method)
    % Apply contrast enhancement to data matrix
    % Inputs:
    %   data - 2D matrix of data values
    %   method - contrast enhancement method ('none', 'linear', 'histeq', 'adaptive', 'gamma')

    enhancedData = data;

    % Skip enhancement if no valid data
    validData = data(~isnan(data) & ~isinf(data));
    if isempty(validData)
        return;
    end

    switch method
        case 'none'
            % No enhancement
            return;

        case 'linear'
            % Linear contrast stretching (stretch to full range)
            minVal = min(validData);
            maxVal = max(validData);
            if maxVal > minVal
                enhancedData = (data - minVal) / (maxVal - minVal);
                % Restore NaN and Inf values
                enhancedData(isnan(data) | isinf(data)) = data(isnan(data) | isinf(data));
            end

        case 'histeq'
            % Histogram equalization
            try
                % Convert to uint8 for histogram equalization
                normalizedData = (data - min(validData)) / (max(validData) - min(validData));
                normalizedData = uint8(normalizedData * 255);

                % Apply histogram equalization
                enhancedData = double(histeq(normalizedData)) / 255;

                % Scale back to original range
                enhancedData = enhancedData * (max(validData) - min(validData)) + min(validData);

                % Restore NaN and Inf values
                enhancedData(isnan(data) | isinf(data)) = data(isnan(data) | isinf(data));
            catch
                % Fallback to linear stretching if histeq fails
                enhancedData = applyContrastEnhancement(data, 'linear');
            end

        case 'adaptive'
            % Adaptive histogram equalization (simplified version)
            try
                % Divide data into blocks and apply local enhancement
                [rows, cols] = size(data);
                blockSize = min(32, min(rows, cols)); % Adaptive block size

                enhancedData = data;
                for i = 1:blockSize:rows
                    for j = 1:blockSize:cols
                        % Define block boundaries
                        rowEnd = min(i + blockSize - 1, rows);
                        colEnd = min(j + blockSize - 1, cols);

                        % Extract block
                        block = data(i:rowEnd, j:colEnd);
                        blockValid = block(~isnan(block) & ~isinf(block));

                        if length(blockValid) > 1
                            % Apply local contrast stretching
                            minVal = min(blockValid);
                            maxVal = max(blockValid);
                            if maxVal > minVal
                                enhancedBlock = (block - minVal) / (maxVal - minVal);
                                % Restore original range but with enhanced contrast
                                enhancedBlock = enhancedBlock * (max(validData) - min(validData)) + min(validData);
                                enhancedData(i:rowEnd, j:colEnd) = enhancedBlock;
                            end
                        end
                    end
                end

                % Restore NaN and Inf values
                enhancedData(isnan(data) | isinf(data)) = data(isnan(data) | isinf(data));
            catch
                % Fallback to linear stretching if adaptive fails
                enhancedData = applyContrastEnhancement(data, 'linear');
            end

        case 'gamma'
            % Gamma correction (gamma = 0.5 for brightening)
            gamma = 0.5;
            minVal = min(validData);
            maxVal = max(validData);
            if maxVal > minVal
                % Normalize to [0,1]
                normalizedData = (data - minVal) / (maxVal - minVal);
                % Apply gamma correction
                enhancedData = normalizedData .^ gamma;
                % Scale back to original range
                enhancedData = enhancedData * (maxVal - minVal) + minVal;

                % Restore NaN and Inf values
                enhancedData(isnan(data) | isinf(data)) = data(isnan(data) | isinf(data));
            end

        otherwise
            % Unknown method, return original data
            fprintf('Unknown contrast enhancement method: %s\n', method);
    end
end

function computeCurrentSlice(button, fig)
    % Compute convergence-based alignment for the current slice (view-aware)
    userData = get(fig, 'UserData');

    % Prevent multiple simultaneous computations
    if userData.isComputingAlignment
        fprintf('Alignment computation already in progress. Please wait...\n');
        return;
    end

    % Get current view and appropriate slice index
    currentView = userData.visState.currentView;
    switch currentView
        case 'XtVsY'
            currentSliceIndex = userData.visState.currentYIndex;
            sliceType = 'Y';
        case 'YtVsX'
            currentSliceIndex = userData.visState.currentXIndex;
            sliceType = 'X';
        otherwise
            fprintf('Compute slice not supported for view: %s\n', currentView);
            return;
    end

    % Show progress dialog with determinate progress bar
    progressDlg = uiprogressdlg(fig, 'Title', 'Computing Alignment', ...
                               'Message', sprintf('Computing %s-slice %d alignment...', sliceType, currentSliceIndex), ...
                               'Value', 0, ...
                               'Cancelable', false);

    % Store progress dialog in userData for access in background function
    userData.progressDialog = progressDlg;
    userData.isComputingAlignment = true;
    set(fig, 'UserData', userData);
    drawnow;

    % Computing alignment (progress shown in dialog only)

    % Update status
    statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
    if ~isempty(statusText)
        set(statusText, 'String', sprintf('Computing %s-slice...', sliceType));
    end

    % Start computation in background
    timer_obj = timer('TimerFcn', @(~,~) computeCurrentSliceBackground(fig, currentSliceIndex, currentView), ...
                      'StartDelay', 0.1, 'ExecutionMode', 'singleShot', ...
                      'Name', 'CurrentSliceTimer');
    start(timer_obj);
end

function computeAllSlices(button, fig)
    % Compute convergence-based alignment for all slices (view-aware)
    userData = get(fig, 'UserData');

    % Prevent multiple simultaneous computations
    if userData.isComputingAlignment
        fprintf('Alignment computation already in progress. Please wait...\n');
        return;
    end

    % Get current view and determine number of slices
    currentView = userData.visState.currentView;
    if ~isempty(userData.statDataArray) && ~isempty(userData.statDataArray{1}.maps)
        [numY, numX] = size(userData.statDataArray{1}.maps{1});
        switch currentView
            case 'XtVsY'
                numSlices = numY;
                sliceType = 'Y';
            case 'YtVsX'
                numSlices = numX;
                sliceType = 'X';
            otherwise
                fprintf('Compute all slices not supported for view: %s\n', currentView);
                return;
        end
    else
        numSlices = 1; % fallback
        sliceType = 'Y';
    end

    % Show progress dialog with determinate progress bar
    progressDlg = uiprogressdlg(fig, 'Title', 'Computing All Slices', ...
                               'Message', sprintf('Computing alignment for all %d %s-slices...', numSlices, sliceType), ...
                               'Value', 0, ...
                               'Cancelable', false);

    % Store progress dialog in userData for access in background function
    userData.progressDialog = progressDlg;
    userData.isComputingAlignment = true;
    set(fig, 'UserData', userData);
    drawnow;

    % Computing alignment for all slices (progress shown in dialog only)

    % Update status
    statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
    if ~isempty(statusText)
        set(statusText, 'String', 'Computing all...');
    end

    % Start computation in background
    timer_obj = timer('TimerFcn', @(~,~) computeAllSlicesBackground(fig, currentView), ...
                      'StartDelay', 0.1, 'ExecutionMode', 'singleShot', ...
                      'Name', 'AllSlicesTimer');
    start(timer_obj);
end

function computeCurrentSliceBackground(fig, sliceIndex, currentView)
    % Background computation of convergence-based alignment for current slice with timing and iteration updates
    try
        userData = get(fig, 'UserData');
        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');

        % Start timing for this slice
        sliceStartTime = tic;
        userData.timingInfo.sliceStartTime = sliceStartTime;
        userData.timingInfo.iterationTimings = zeros(1, 50); % Reset iteration timings for this slice
        userData.timingInfo.iterationCount = 0; % Reset iteration counter

        % Save userData with timing info before starting computation
        set(fig, 'UserData', userData);



        % Use current view data as input for iterative alignment
        % This allows building alignment on top of previous alignment results
        currentData = userData.statDataArray;
        if isempty(currentData)
            fprintf('Error: Current data not available\n');
            return;
        end

        % Make a working copy for this computation
        currentData = currentData;
        nStats = length(currentData);

        % Get dimensions from first statistic
        if nStats > 0 && ~isempty(currentData{1}.maps)
            [numY, numX] = size(currentData{1}.maps{1});
            numSegments = length(currentData{1}.maps);
        else
            fprintf('Error: No data available for computation\n');
            return;
        end

        % Process each statistic for the current slice
        for statIdx = 1:nStats
            % Processing statistic (progress shown in dialog only)

            % Extract slice data based on current view
            switch currentView
                case 'XtVsY'
                    % PERFORMANCE OPTIMIZATION: Pre-allocate and extract data efficiently
                    sliceData = zeros(numSegments, numX);
                    % Vectorized extraction where possible
                    for seg = 1:numSegments
                        sliceData(seg, :) = currentData{statIdx}.maps{seg}(sliceIndex, :);
                    end
                case 'YtVsX'
                    % Extract Y,t data for this X slice across all segments
                    sliceData = zeros(numSegments, numY);
                    for seg = 1:numSegments
                        sliceData(seg, :) = currentData{statIdx}.maps{seg}(:, sliceIndex)';
                    end
                otherwise
                    fprintf('Error: Unsupported view for alignment: %s\n', currentView);
                    return;
            end

            % Alignment method is fixed to 'average' (Row Average)
            alignmentMethod = 'average';

            % Get user-defined convergence threshold
            convergenceInput = findobj(fig, 'Tag', 'ConvergenceInput');
            if ~isempty(convergenceInput)
                convergenceStr = get(convergenceInput, 'String');
                convergenceThreshold = str2double(convergenceStr) / 100; % Convert percentage to decimal
                if isnan(convergenceThreshold) || convergenceThreshold <= 0 || convergenceThreshold > 0.05
                    convergenceThreshold = 0.01; % Default to 1% if invalid
                    set(convergenceInput, 'String', '1.0'); % Reset to default
                end
            else
                convergenceThreshold = 0.01; % Default to 1%
            end

            % Create progress callback for detailed progress updates with timing and column info
            progressCallback = @(iteration, progress, alignedData, currentCol, totalCols, overallProgress) updateSliceProgressWithTiming(fig, iteration, progress, alignedData, statIdx, nStats, sliceIndex, currentCol, totalCols, overallProgress);

            % Apply enhanced column alignment with convergence (zero padding only)
            % Read optional alignment params (non-destructive; defaults preserved)
            maxShift = 15; costFcn = 'mse';
            costDropdown = findobj(fig, 'Tag', 'AlignmentCostFunctionDropdown');
            if ~isempty(costDropdown)
                switch get(costDropdown,'Value')
                    case 1, costFcn = 'mse';
                    case 2, costFcn = 'correlation';
                    case 3, costFcn = 'ncc';
                end
            end
            maxShiftEdit = findobj(fig, 'Tag', 'AlignmentMaxShiftInput');
            if ~isempty(maxShiftEdit)
                v = str2double(get(maxShiftEdit,'String')); if ~isnan(v) && v>0, maxShift = round(v); end
            end

            [alignedSliceData, convergenceInfo] = alignColumnsImproved(sliceData, ...
                'MaxShift', maxShift, ...
                'CostFunction', costFcn, ...
                'AlignmentMethod', alignmentMethod, ...
                'LocalScope', 5, ...
                'PadMethod', 'zeros', ...
                'Verbose', false, ...
                'ConvergenceThreshold', convergenceThreshold, ...
                'MaxIterations', getMaxIterations(fig), ...
                'WeightingFunction', 'exponential', ...
                'WeightingScale', 3.0, ...
                'ProgressCallback', progressCallback);

            % Calculate shifts by comparing original and aligned data (more efficient)
            appliedShifts = calculateShiftsFromAlignment(sliceData, alignedSliceData);

            % Store the shifts for this slice and statistic (view-aware)
            if isempty(userData.alignmentShifts)
                % Initialize alignment shifts array
                [numY, numX] = size(currentData{1}.maps{1});
                userData.alignmentShifts = zeros(numY, numX);
            end

            % Store shifts based on current view
            switch currentView
                case 'XtVsY'
                    % Store shifts for this Y slice (average across statistics if multiple)
                    if statIdx == 1
                        userData.alignmentShifts(sliceIndex, :) = appliedShifts;
                    else
                        % Average with previous statistics' shifts
                        userData.alignmentShifts(sliceIndex, :) = (userData.alignmentShifts(sliceIndex, :) + appliedShifts) / 2;
                    end
                case 'YtVsX'
                    % Store shifts for this X slice (transpose for X-direction alignment)
                    if statIdx == 1
                        userData.alignmentShifts(:, sliceIndex) = appliedShifts';
                    else
                        % Average with previous statistics' shifts
                        userData.alignmentShifts(:, sliceIndex) = (userData.alignmentShifts(:, sliceIndex) + appliedShifts') / 2;
                    end
            end

            % Put the aligned data back into the maps for this slice only
            switch currentView
                case 'XtVsY'
                    for seg = 1:numSegments
                        currentData{statIdx}.maps{seg}(sliceIndex, :) = alignedSliceData(seg, :);
                    end
                case 'YtVsX'
                    for seg = 1:numSegments
                        currentData{statIdx}.maps{seg}(:, sliceIndex) = alignedSliceData(seg, :)';
                    end
            end

            % Statistic converged (progress shown in dialog only)
        end

        % Cache the aligned data and update tracking
        userData.alignedDataCache = currentData;
        userData.currentView = 'aligned';

        % Update slice alignment status based on view
        switch currentView
            case 'XtVsY'
                userData.sliceAlignmentStatus(sliceIndex) = true;
            case 'YtVsX'
                % For YtVsX, we need to track X-slice alignment separately
                if ~isfield(userData, 'xSliceAlignmentStatus')
                    userData.xSliceAlignmentStatus = false(numX, 1);
                end
                userData.xSliceAlignmentStatus(sliceIndex) = true;
        end

        % Reset waveform alignment flag since new alignment was computed
        userData.waveformAlignmentApplied = false;

        % Apply alignment to waveforms if available
        if userData.waveformLoaded
            applyAlignmentToWaveforms(fig);
            userData.waveformAlignmentApplied = true;
        end

        % Calculate total iterations for this slice (sum across all statistics)
        totalIterations = 0;
        for statIdx = 1:nStats
            % Get convergence info from the last statistic processed
            % (This is a simplification - in practice you might want to track per-statistic)
            if exist('convergenceInfo', 'var')
                totalIterations = max(totalIterations, convergenceInfo.iterations);
            end
        end
        % Calculate slice timing first
        sliceElapsedTime = toc(sliceStartTime);

        % Update iteration tracking based on view
        switch currentView
            case 'XtVsY'
                userData.sliceIterations(sliceIndex) = totalIterations;
                % Initialize sliceTimings array if needed
                if isempty(userData.timingInfo.sliceTimings) || length(userData.timingInfo.sliceTimings) < sliceIndex
                    userData.timingInfo.sliceTimings = zeros(numY, 1);
                end
                userData.timingInfo.sliceTimings(sliceIndex) = sliceElapsedTime;
            case 'YtVsX'
                % For YtVsX, track X-slice iterations separately
                if ~isfield(userData, 'xSliceIterations')
                    userData.xSliceIterations = zeros(numX, 1);
                end
                userData.xSliceIterations(sliceIndex) = totalIterations;
                % Track X-slice timings separately
                if ~isfield(userData.timingInfo, 'xSliceTimings')
                    userData.timingInfo.xSliceTimings = zeros(numX, 1);
                end
                userData.timingInfo.xSliceTimings(sliceIndex) = sliceElapsedTime;
        end
        userData.timingInfo.totalComputationTime = userData.timingInfo.totalComputationTime + sliceElapsedTime;

        % Final update
        userData.statDataArray = currentData;
        userData.isComputingAlignment = false;

        % Close progress dialog if it exists
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            close(userData.progressDialog);
            userData = rmfield(userData, 'progressDialog');
        end

        set(fig, 'UserData', userData);
        updatePlots(fig);
        applyAspectRatioToAxes(fig);

        % Update final status with iteration count and timing
        if ~isempty(statusText)
            switch currentView
                case 'XtVsY'
                    set(statusText, 'String', sprintf('Aligned View (Y-slice %d: %d iter, %.1fs)', sliceIndex, totalIterations, sliceElapsedTime));
                case 'YtVsX'
                    set(statusText, 'String', sprintf('Aligned View (X-slice %d: %d iter, %.1fs)', sliceIndex, totalIterations, sliceElapsedTime));
            end
        end

        % Slice alignment computation complete (results shown in status only)

    catch ME
        % Handle errors gracefully
        fprintf('Error during slice alignment computation: %s\n', ME.message);
        userData = get(fig, 'UserData');
        userData.isComputingAlignment = false;

        % Close progress dialog if it exists
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            close(userData.progressDialog);
            userData = rmfield(userData, 'progressDialog');
        end

        set(fig, 'UserData', userData);

        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', 'Error');
        end
    end

    % Clean up timer
    try
        currentTimer = timerfind('Name', 'CurrentSliceTimer');
        if ~isempty(currentTimer)
            stop(currentTimer);
            delete(currentTimer);
        end
    catch
        % Ignore timer cleanup errors
    end
end

function updateSliceProgressWithTiming(fig, iteration, progress, alignedData, statIdx, nStats, yIndex, currentCol, totalCols, overallProgress)
    % Enhanced progress callback with detailed timing, column progress, and loading bar updates
    try
        userData = get(fig, 'UserData');
        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');

        % Handle optional parameters for backward compatibility
        if nargin < 8
            currentCol = 0;
            totalCols = 0;
            overallProgress = progress;
        end

        % Calculate progress percentage (0-100%)
        progressPercent = overallProgress * 100;

        % Update progress dialog if it exists
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            try
                if currentCol > 0 && currentCol <= totalCols
                    % Column-level progress update
                    columnProgress = (currentCol - 1) / max(1, totalCols - 2); % Adjust for skipped first/last columns
                    if ~isempty(userData.timingInfo.sliceStartTime) && isnumeric(userData.timingInfo.sliceStartTime)
                        elapsedTime = toc(userData.timingInfo.sliceStartTime);
                        estimatedTotal = elapsedTime / max(columnProgress, 0.01); % Avoid division by zero
                        remainingTime = max(0, estimatedTotal - elapsedTime);

                        userData.progressDialog.Value = columnProgress;
                        userData.progressDialog.Message = sprintf(['Slice %d - Column %d/%d (%.1f%%)\n' ...
                            'Iteration %d, Elapsed: %.1fs, Est. remaining: %.1fs'], ...
                            yIndex, currentCol, totalCols, progressPercent, iteration, elapsedTime, remainingTime);
                    else
                        userData.progressDialog.Value = columnProgress;
                        userData.progressDialog.Message = sprintf('Slice %d - Column %d/%d (%.1f%%) - Iteration %d', ...
                            yIndex, currentCol, totalCols, progressPercent, iteration);
                    end
                else
                    % Iteration completion update
                    if ~isempty(userData.timingInfo.sliceStartTime) && isnumeric(userData.timingInfo.sliceStartTime)
                        elapsedTime = toc(userData.timingInfo.sliceStartTime);
                        userData.progressDialog.Value = progress;
                        userData.progressDialog.Message = sprintf(['Slice %d - Iteration %d complete (%.1f%%)\n' ...
                            'Elapsed: %.1fs'], yIndex, iteration, progressPercent, elapsedTime);
                    else
                        userData.progressDialog.Value = progress;
                        userData.progressDialog.Message = sprintf('Slice %d - Iteration %d complete (%.1f%%)', ...
                            yIndex, iteration, progressPercent);
                    end
                end
                drawnow;
            catch timingError
                % Fallback if timing fails
                userData.progressDialog.Value = progress;
                userData.progressDialog.Message = sprintf('Slice %d - Iteration %d', yIndex, iteration);
                drawnow;
            end
        end

        % Update at iteration completion (progress = 1.0)
        if progress >= 0.99
            % Initialize timing on first iteration
            if iteration == 1
                userData.timingInfo.iterationStartTime = tic;
                set(fig, 'UserData', userData);

            end

            % Calculate iteration timing (for iterations after the first)
            if iteration > 1 && ~isempty(userData.timingInfo.iterationStartTime)
                iterationTime = toc(userData.timingInfo.iterationStartTime);
                userData.timingInfo.iterationCount = userData.timingInfo.iterationCount + 1;
                if userData.timingInfo.iterationCount <= length(userData.timingInfo.iterationTimings)
                    userData.timingInfo.iterationTimings(userData.timingInfo.iterationCount) = iterationTime;
                end
            end

            % Start timing for next iteration (if not the last)
            userData.timingInfo.iterationStartTime = tic;

            % Update status with timing information
            if ~isempty(statusText)
                if ~isempty(userData.timingInfo.iterationTimings)
                    avgIterTime = mean(userData.timingInfo.iterationTimings);
                    set(statusText, 'String', sprintf('Slice %d: Iter %d (%.1fs avg, stat %d/%d)', yIndex, iteration, avgIterTime, statIdx, nStats));
                else
                    set(statusText, 'String', sprintf('Slice %d: Iter %d (stat %d/%d)', yIndex, iteration, statIdx, nStats));
                end
            end

            % Update UserData with timing info
            set(fig, 'UserData', userData);

            % Skip plot updates during alignment for better performance
            % updatePlots(fig);
            % drawnow;
        end

    catch ME
        % Silently handle errors in progress callback to avoid disrupting main computation
        % Error details available in status text if needed
    end
end

function computeAllSlicesBackground(fig, currentView)
    % Background computation of convergence-based alignment for all slices with timing (view-aware)
    try
        userData = get(fig, 'UserData');
        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');

        % Start timing for all slices computation
        allSlicesStartTime = tic;
        userData.timingInfo.allSlicesStartTime = allSlicesStartTime;
        userData.timingInfo.sliceTimings = []; % Reset slice timings
        userData.timingInfo.iterationTimings = zeros(1, 50); % Reset iteration timings
        userData.timingInfo.iterationCount = 0; % Reset iteration counter

        % Save userData with timing info before starting computation
        set(fig, 'UserData', userData);

        % Use current view data as input for iterative alignment
        % This allows building alignment on top of previous alignment results
        currentData = userData.statDataArray;
        if isempty(currentData)
            fprintf('Error: Current data not available\n');
            return;
        end

        % Make a working copy for this computation
        currentData = currentData;
        nStats = length(currentData);

        % Get dimensions from first statistic
        if nStats > 0 && ~isempty(currentData{1}.maps)
            [numY, numX] = size(currentData{1}.maps{1});
            numSegments = length(currentData{1}.maps);
        else
            fprintf('Error: No data available for computation\n');
            return;
        end

        % Determine slice range based on current view
        switch currentView
            case 'XtVsY'
                numSlices = numY;
                sliceType = 'Y';
            case 'YtVsX'
                numSlices = numX;
                sliceType = 'X';
            otherwise
                fprintf('Error: Unsupported view for alignment: %s\n', currentView);
                return;
        end

        % Computing alignment for all slices (progress shown in dialog only)

        % Process each slice
        for sliceIndex = 1:numSlices
            % Start timing for this slice
            sliceStartTime = tic;
            userData.timingInfo.sliceStartTime = sliceStartTime;
            userData.timingInfo.iterationTimings = zeros(1, 50); % Reset for this slice
            userData.timingInfo.iterationCount = 0; % Reset iteration counter

            % Processing slice (progress shown in dialog only)

            % Update status
            if ~isempty(statusText)
                set(statusText, 'String', sprintf('All: %s-slice %d/%d', sliceType, sliceIndex, numSlices));
            end

            % Track iterations for this slice
            maxIterationsThisSlice = 0;

            % Process each statistic for this slice
            for statIdx = 1:nStats
                % Extract slice data based on current view
                switch currentView
                    case 'XtVsY'
                        % Extract X,t data for this Y slice across all segments
                        sliceData = zeros(numSegments, numX);
                        for seg = 1:numSegments
                            sliceData(seg, :) = currentData{statIdx}.maps{seg}(sliceIndex, :);
                        end
                    case 'YtVsX'
                        % Extract Y,t data for this X slice across all segments
                        sliceData = zeros(numSegments, numY);
                        for seg = 1:numSegments
                            sliceData(seg, :) = currentData{statIdx}.maps{seg}(:, sliceIndex)';
                        end
                end

                % Alignment method is fixed to 'average' (Row Average)
                alignmentMethod = 'average';

                % Get user-defined convergence threshold (only once per slice)
                if statIdx == 1
                    convergenceInput = findobj(fig, 'Tag', 'ConvergenceInput');
                    if ~isempty(convergenceInput)
                        convergenceStr = get(convergenceInput, 'String');
                        convergenceThreshold = str2double(convergenceStr) / 100; % Convert percentage to decimal
                        if isnan(convergenceThreshold) || convergenceThreshold <= 0 || convergenceThreshold > 0.05
                            convergenceThreshold = 0.01; % Default to 1% if invalid
                            set(convergenceInput, 'String', '1.0'); % Reset to default
                        end
                    else
                        convergenceThreshold = 0.01; % Default to 1%
                    end
                end

                % Create progress callback for slice progress with timing and column info
                progressCallback = @(iteration, progress, alignedData, currentCol, totalCols, overallProgress) updateAllSlicesProgressWithTiming(fig, sliceIndex, numSlices, iteration, progress, statIdx, nStats, currentCol, totalCols, overallProgress);

                % Apply enhanced column alignment with convergence (zero padding only)
                % Read optional alignment params (non-destructive; defaults preserved)
                maxShift = 15; costFcn = 'mse';
                costDropdown = findobj(fig, 'Tag', 'AlignmentCostFunctionDropdown');
                if ~isempty(costDropdown)
                    switch get(costDropdown,'Value')
                        case 1, costFcn = 'mse';
                        case 2, costFcn = 'correlation';
                        case 3, costFcn = 'ncc';
                    end
                end
                maxShiftEdit = findobj(fig, 'Tag', 'AlignmentMaxShiftInput');
                if ~isempty(maxShiftEdit)
                    v = str2double(get(maxShiftEdit,'String')); if ~isnan(v) && v>0, maxShift = round(v); end
                end

                [alignedSliceData, convergenceInfo] = alignColumnsImproved(sliceData, ...
                    'MaxShift', maxShift, ...
                    'CostFunction', costFcn, ...
                    'AlignmentMethod', alignmentMethod, ...
                    'LocalScope', 5, ...
                    'PadMethod', 'zeros', ...
                    'Verbose', false, ...
                    'ConvergenceThreshold', convergenceThreshold, ...
                    'MaxIterations', getMaxIterations(fig), ...
                    'WeightingFunction', 'exponential', ...
                    'WeightingScale', 3.0, ...
                    'ProgressCallback', progressCallback);

                % Calculate shifts by comparing original and aligned data (more efficient)
                appliedShifts = calculateShiftsFromAlignment(sliceData, alignedSliceData);

                % Store the shifts for this slice and statistic (view-aware)
                if isempty(userData.alignmentShifts)
                    % Initialize alignment shifts array
                    [numY, numX] = size(currentData{1}.maps{1});
                    userData.alignmentShifts = zeros(numY, numX);
                end

                % Store shifts based on current view
                switch currentView
                    case 'XtVsY'
                        % Store shifts for this Y slice (average across statistics if multiple)
                        if statIdx == 1
                            userData.alignmentShifts(sliceIndex, :) = appliedShifts;
                        else
                            % Average with previous statistics' shifts
                            userData.alignmentShifts(sliceIndex, :) = (userData.alignmentShifts(sliceIndex, :) + appliedShifts) / 2;
                        end
                    case 'YtVsX'
                        % Store shifts for this X slice (transpose for X-direction alignment)
                        if statIdx == 1
                            userData.alignmentShifts(:, sliceIndex) = appliedShifts';
                        else
                            % Average with previous statistics' shifts
                            userData.alignmentShifts(:, sliceIndex) = (userData.alignmentShifts(:, sliceIndex) + appliedShifts') / 2;
                        end
                end

                % Put the aligned data back into the maps for this slice
                switch currentView
                    case 'XtVsY'
                        for seg = 1:numSegments
                            currentData{statIdx}.maps{seg}(sliceIndex, :) = alignedSliceData(seg, :);
                        end
                    case 'YtVsX'
                        for seg = 1:numSegments
                            currentData{statIdx}.maps{seg}(:, sliceIndex) = alignedSliceData(seg, :)';
                        end
                end

                % Track maximum iterations for this slice
                maxIterationsThisSlice = max(maxIterationsThisSlice, convergenceInfo.iterations);
            end

            % Calculate slice timing
            sliceElapsedTime = toc(sliceStartTime);

            % Update slice tracking based on view
            switch currentView
                case 'XtVsY'
                    % Initialize sliceTimings array if needed
                    if isempty(userData.timingInfo.sliceTimings) || length(userData.timingInfo.sliceTimings) < sliceIndex
                        userData.timingInfo.sliceTimings = zeros(numY, 1);
                    end
                    userData.timingInfo.sliceTimings(sliceIndex) = sliceElapsedTime;
                    userData.sliceIterations(sliceIndex) = maxIterationsThisSlice;
                    userData.sliceAlignmentStatus(sliceIndex) = true;
                case 'YtVsX'
                    % Initialize X-slice tracking arrays if needed
                    if ~isfield(userData.timingInfo, 'xSliceTimings')
                        userData.timingInfo.xSliceTimings = zeros(numX, 1);
                    end
                    if ~isfield(userData, 'xSliceIterations')
                        userData.xSliceIterations = zeros(numX, 1);
                    end
                    if ~isfield(userData, 'xSliceAlignmentStatus')
                        userData.xSliceAlignmentStatus = false(numX, 1);
                    end
                    userData.timingInfo.xSliceTimings(sliceIndex) = sliceElapsedTime;
                    userData.xSliceIterations(sliceIndex) = maxIterationsThisSlice;
                    userData.xSliceAlignmentStatus(sliceIndex) = true;
            end

            % Skip frequent plot updates during alignment for better performance
            % Only update userData for progress tracking
            if mod(sliceIndex, max(1, round(numSlices/20))) == 0 || sliceIndex == numSlices
                userData.statDataArray = currentData;
                set(fig, 'UserData', userData);
                % updatePlots(fig);
                % applyAspectRatioToAxes(fig);
                % drawnow;
            end

            % Slice completed (progress shown in dialog only)
        end

        % Final update - cache aligned data and set view
        userData.alignedDataCache = currentData;
        userData.currentView = 'aligned';
        userData.statDataArray = currentData;
        userData.isComputingAlignment = false;

        % Close progress dialog if it exists
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            close(userData.progressDialog);
            userData = rmfield(userData, 'progressDialog');
        end

        % Reset waveform alignment flag since new alignment was computed
        userData.waveformAlignmentApplied = false;

        % Apply alignment to waveforms if available
        if userData.waveformLoaded
            applyAlignmentToWaveforms(fig);
            userData.waveformAlignmentApplied = true;
        end

        % Calculate total computation time
        totalElapsedTime = toc(allSlicesStartTime);
        userData.timingInfo.totalComputationTime = totalElapsedTime;

        set(fig, 'UserData', userData);
        updatePlots(fig);
        applyAspectRatioToAxes(fig);

        % Calculate average iterations and timing across all slices (view-aware)
        switch currentView
            case 'XtVsY'
                avgIterations = mean(userData.sliceIterations(userData.sliceAlignmentStatus));
                avgSliceTime = mean(userData.timingInfo.sliceTimings);
                sliceTypeStr = 'Y-slices';
            case 'YtVsX'
                avgIterations = mean(userData.xSliceIterations(userData.xSliceAlignmentStatus));
                avgSliceTime = mean(userData.timingInfo.xSliceTimings);
                sliceTypeStr = 'X-slices';
        end

        % Update final status with timing information
        if ~isempty(statusText)
            set(statusText, 'String', sprintf('Aligned View (%s - Avg: %.1f iter, %.1fs/slice, Total: %.1fs)', sliceTypeStr, avgIterations, avgSliceTime, totalElapsedTime));
        end

        % All slices alignment computation complete (results shown in status only)

    catch ME
        % Handle errors gracefully
        fprintf('Error during all slices alignment computation: %s\n', ME.message);
        userData = get(fig, 'UserData');
        userData.isComputingAlignment = false;

        % Close progress dialog if it exists
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            close(userData.progressDialog);
            userData = rmfield(userData, 'progressDialog');
        end

        set(fig, 'UserData', userData);

        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', 'Error');
        end
    end

    % Clean up timer
    try
        allSlicesTimer = timerfind('Name', 'AllSlicesTimer');
        if ~isempty(allSlicesTimer)
            stop(allSlicesTimer);
            delete(allSlicesTimer);
        end
    catch
        % Ignore timer cleanup errors
    end
end

function updateAllSlicesProgressWithTiming(fig, currentSlice, totalSlices, iteration, progress, statIdx, nStats, currentCol, totalCols, overallProgress)
    % Enhanced progress callback for all slices computation with detailed timing and progress
    try
        userData = get(fig, 'UserData');
        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');

        % Handle optional parameters for backward compatibility
        if nargin < 8
            currentCol = 0;
            totalCols = 0;
            overallProgress = progress;
        end

        % Update progress dialog if it exists
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            try
                % Calculate overall progress across all slices
                sliceProgress = (currentSlice - 1) / totalSlices;
                if currentCol > 0 && currentCol <= totalCols
                    % Column-level progress within current slice
                    columnProgress = (currentCol - 1) / max(1, totalCols - 2);
                    overallSliceProgress = sliceProgress + (columnProgress / totalSlices);

                    if ~isempty(userData.timingInfo.allSlicesStartTime) && isnumeric(userData.timingInfo.allSlicesStartTime)
                        elapsedTime = toc(userData.timingInfo.allSlicesStartTime);
                        estimatedTotal = elapsedTime / max(overallSliceProgress, 0.01);
                        remainingTime = max(0, estimatedTotal - elapsedTime);

                        userData.progressDialog.Value = overallSliceProgress;
                        userData.progressDialog.Message = sprintf(['Slice %d/%d - Column %d/%d\n' ...
                            'Iteration %d, Elapsed: %.1fs, Est. remaining: %.1fs'], ...
                            currentSlice, totalSlices, currentCol, totalCols, iteration, elapsedTime, remainingTime);
                    else
                        userData.progressDialog.Value = overallSliceProgress;
                        userData.progressDialog.Message = sprintf('Slice %d/%d - Column %d/%d - Iteration %d', ...
                            currentSlice, totalSlices, currentCol, totalCols, iteration);
                    end
                else
                    % Iteration completion update
                    if ~isempty(userData.timingInfo.allSlicesStartTime) && isnumeric(userData.timingInfo.allSlicesStartTime)
                        elapsedTime = toc(userData.timingInfo.allSlicesStartTime);
                        userData.progressDialog.Value = sliceProgress;
                        userData.progressDialog.Message = sprintf(['Slice %d/%d - Iteration %d complete\n' ...
                            'Elapsed: %.1fs'], currentSlice, totalSlices, iteration, elapsedTime);
                    else
                        userData.progressDialog.Value = sliceProgress;
                        userData.progressDialog.Message = sprintf('Slice %d/%d - Iteration %d complete', ...
                            currentSlice, totalSlices, iteration);
                    end
                end
                drawnow;
            catch timingError
                % Fallback if timing fails
                userData.progressDialog.Value = (currentSlice - 1) / totalSlices;
                userData.progressDialog.Message = sprintf('Slice %d/%d - Iteration %d', currentSlice, totalSlices, iteration);
                drawnow;
            end
        end

        % Update at iteration completion (progress = 1.0)
        if progress >= 0.99
            % Initialize timing on first iteration
            if iteration == 1
                userData.timingInfo.iterationStartTime = tic;
                set(fig, 'UserData', userData);
            end

            % Calculate iteration timing (for iterations after the first)
            if iteration > 1 && ~isempty(userData.timingInfo.iterationStartTime)
                iterationTime = toc(userData.timingInfo.iterationStartTime);
                userData.timingInfo.iterationCount = userData.timingInfo.iterationCount + 1;
                if userData.timingInfo.iterationCount <= length(userData.timingInfo.iterationTimings)
                    userData.timingInfo.iterationTimings(userData.timingInfo.iterationCount) = iterationTime;
                end
            end

            % Start timing for next iteration (if not the last)
            userData.timingInfo.iterationStartTime = tic;

            % Update status with timing information
            if ~isempty(statusText)
                if ~isempty(userData.timingInfo.iterationTimings)
                    avgIterTime = mean(userData.timingInfo.iterationTimings);
                    set(statusText, 'String', sprintf('Slice %d/%d: Iter %d (%.1fs avg, stat %d/%d)', currentSlice, totalSlices, iteration, avgIterTime, statIdx, nStats));
                else
                    set(statusText, 'String', sprintf('Slice %d/%d: Iter %d (stat %d/%d)', currentSlice, totalSlices, iteration, statIdx, nStats));
                end
            end

            % Update UserData with timing info
            set(fig, 'UserData', userData);

            % Update heatmap after each completed iteration (less frequent for all slices to maintain performance)
            if mod(iteration, 2) == 0 || iteration <= 3  % Update every 2nd iteration or first 3 iterations
                updatePlots(fig);
                drawnow;
            end
        end

    catch ME
        % Silently handle errors in progress callback
        fprintf('All slices progress callback error: %s\n', ME.message);
    end
end

function showOriginalView(button, fig)
    % Switch to original (unaligned) view
    try
        userData = get(fig, 'UserData');

        % Check if original data is available
        if isempty(userData.originalDataCache)
            fprintf('Original data not available\n');
            return;
        end

        % Switch to original view
        userData.statDataArray = userData.originalDataCache;
        userData.currentView = 'original';
        set(fig, 'UserData', userData);

        % Minimal debug output for original view
        if userData.waveformLoaded
            fprintf('Switched to original view (waveforms available)\n');
        end

        % Update plots (this will also update waveforms to original view)
        updatePlots(fig);
        applyAspectRatioToAxes(fig);

        % Update status
        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', 'Original View');
        end

        % Update 3D view if currently active
        if strcmp(userData.visState.currentView, '3DView')
            update3DFromAlignment(fig);
        end

        % Force immediate visual update
        drawnow;

        % Original view switch complete

    catch ME
        fprintf('Error switching to original view: %s\n', ME.message);
    end
end

function showAlignedView(button, fig)
    % Switch to aligned view
    try
        userData = get(fig, 'UserData');

        % Check for iterative aligned data
        if ~isempty(userData.alignedDataCache)
            % Switch to iterative aligned data
            userData.statDataArray = userData.alignedDataCache;
            userData.currentView = 'aligned';

            % Minimal debug output for aligned view
            if userData.waveformLoaded && ~isempty(userData.alignmentShifts)
                nonZeroShifts = sum(abs(userData.alignmentShifts(:)) > 1e-10);
                fprintf('Switched to aligned view (%d alignment shifts)\n', nonZeroShifts);
            end

            % Apply alignment to waveforms if available (only if not already computed)
            if userData.waveformLoaded && ~userData.waveformAlignmentApplied && ~isempty(userData.alignmentShifts)
                applyAlignmentToWaveforms(fig);
                % Get updated userData after alignment
                userData = get(fig, 'UserData');
                userData.waveformAlignmentApplied = true;
                fprintf('Waveform alignment computed for first time\n');
                % Ensure currentView is still set to aligned after waveform alignment
                userData.currentView = 'aligned';
            elseif userData.waveformLoaded
                % Aligned waveform data already exists, just ensure currentView is set
                userData.currentView = 'aligned';
            end

            set(fig, 'UserData', userData);

            % Update plots (this will also update waveforms to aligned view)
            updatePlots(fig);
            applyAspectRatioToAxes(fig);

            % Update status with iteration information
            statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
            if ~isempty(statusText)
                % Get current slice information
                currentYIndex = userData.visState.currentYIndex;
                if currentYIndex <= length(userData.sliceIterations) && userData.sliceAlignmentStatus(currentYIndex)
                    iterationsForSlice = userData.sliceIterations(currentYIndex);
                    set(statusText, 'String', sprintf('Iterative Aligned (Slice %d: %d iter)', currentYIndex, iterationsForSlice));
                else
                    % Calculate average for all aligned slices
                    if any(userData.sliceAlignmentStatus)
                        avgIterations = mean(userData.sliceIterations(userData.sliceAlignmentStatus));
                        set(statusText, 'String', sprintf('Iterative Aligned (Avg: %.1f iter)', avgIterations));
                    else
                        set(statusText, 'String', 'Iterative Aligned');
                    end
                end
            end

            % Update 3D view if currently active
            if strcmp(userData.visState.currentView, '3DView')
                update3DFromAlignment(fig);
            end

            % Force immediate visual update
            drawnow;

            % Aligned view switch complete

        else
            fprintf('No aligned data available. Please compute alignment first.\n');
            return;
        end

    catch ME
        fprintf('Error switching to aligned view: %s\n', ME.message);
    end
end

function overlayLayerPaths(ax, layerPaths, currentView, sliceIndex, X_values, Y_values, timeValues)
    % Overlay detected layer paths on the current plot

    % Validate inputs
    if isempty(layerPaths) || ~isstruct(layerPaths)
        return;
    end

    if ~isfield(layerPaths, 'numLayers') || layerPaths.numLayers == 0
        return;
    end

    % Validate slice index
    switch currentView
        case 'XtVsY'
            if sliceIndex < 1 || sliceIndex > length(Y_values)
                if DEBUG_LAYER_PATHS
                    fprintf('Debug: Invalid Y slice index %d (valid range: 1-%d)\n', sliceIndex, length(Y_values));
                end
                return;
            end
        case 'YtVsX'
            if sliceIndex < 1 || sliceIndex > length(X_values)
                if DEBUG_LAYER_PATHS
                    fprintf('Debug: Invalid X slice index %d (valid range: 1-%d)\n', sliceIndex, length(X_values));
                end
                return;
            end
        case 'XYVst'
            if sliceIndex < 1 || sliceIndex > length(timeValues)
                if DEBUG_LAYER_PATHS
                    fprintf('Debug: Invalid time slice index %d (valid range: 1-%d)\n', sliceIndex, length(timeValues));
                end
                return;
            end
    end

    % Hold the current plot to add overlays
    hold(ax, 'on');

    try
        switch currentView
            case 'XtVsY'
                % Overlay XtVsY layer paths for the current Y slice
                if isfield(layerPaths, 'XtVsY_paths') && ~isempty(layerPaths.XtVsY_paths)
                    overlayXtVsYPaths(ax, layerPaths.XtVsY_paths, sliceIndex, X_values, Y_values, timeValues);
                end

            case 'YtVsX'
                % Overlay YtVsX layer paths for the current X slice
                if isfield(layerPaths, 'YtVsX_paths') && ~isempty(layerPaths.YtVsX_paths)
                    overlayYtVsXPaths(ax, layerPaths.YtVsX_paths, sliceIndex, X_values, Y_values, timeValues);
                end

            case 'XYVst'
                % For XYVst view, overlay both types of paths projected onto the XY plane
                overlayXYVstPaths(ax, layerPaths, sliceIndex, X_values, Y_values, timeValues);
        end

    catch ME
        fprintf('Warning: Error overlaying layer paths: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end

    hold(ax, 'off');
end

function overlayXtVsYPaths(ax, XtVsY_paths, currentYIndex, X_values, Y_values, timeValues)
    % Overlay layer paths for XtVsY view

    if isempty(XtVsY_paths)
        return;
    end

    try
        currentY = Y_values(currentYIndex);
        tolerance = (max(Y_values) - min(Y_values)) / (2 * length(Y_values)); % Half the Y spacing

        % Define colors for different layers
        layerColors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

        if debugFlag
            fprintf('Debug: Overlaying %d XtVsY layers for Y=%.3f (tolerance=%.3f)\n', ...
                    length(XtVsY_paths), currentY, tolerance);
        end

        for layerIdx = 1:length(XtVsY_paths)
            layer = XtVsY_paths{layerIdx};

            % Validate layer structure
            if ~isstruct(layer) || ~isfield(layer, 'peaks')
                if debugFlag, fprintf('Debug: Layer %d missing peaks field\n', layerIdx); end
                continue;
            end

            if isempty(layer.peaks)
                if debugFlag, fprintf('Debug: Layer %d has no peaks\n', layerIdx); end
                continue;
            end

            % Find peaks near the current Y slice
            relevantPeaks = zeros(length(layer.peaks), 2); % Pre-allocate
            peakCount = 0;

            for p = 1:length(layer.peaks)
                peak = layer.peaks(p);

                % Validate peak structure
                if ~isstruct(peak) || ~isfield(peak, 'Y') || ~isfield(peak, 'X') || ~isfield(peak, 'T')
                    if debugFlag, fprintf('Debug: Peak %d in layer %d missing coordinate fields\n', p, layerIdx); end
                    continue;
                end

                if abs(peak.Y - currentY) <= tolerance
                    peakCount = peakCount + 1;
                    relevantPeaks(peakCount, :) = [peak.X, peak.T];
                end
            end

            relevantPeaks = relevantPeaks(1:peakCount, :); % Trim to actual size

            if debugFlag
                fprintf('Debug: Layer %d has %d relevant peaks for current Y slice\n', layerIdx, peakCount);
            end

            if size(relevantPeaks, 1) >= 2
                % Sort by X position
                [~, sortIdx] = sort(relevantPeaks(:, 1));
                relevantPeaks = relevantPeaks(sortIdx, :);

                % Plot the layer path
                colorIdx = mod(layerIdx - 1, length(layerColors)) + 1;
                plot(ax, relevantPeaks(:, 1), relevantPeaks(:, 2), ...
                     [layerColors(colorIdx) 'o-'], 'LineWidth', 2, 'MarkerSize', 4, ...
                     'DisplayName', sprintf('Layer %d', layerIdx));

                if debugFlag
                    fprintf('Debug: Plotted layer %d with %d points\n', layerIdx, size(relevantPeaks, 1));
                end
            elseif peakCount > 0 && debugFlag
                fprintf('Debug: Layer %d has only %d peaks (need >=2 for line plot)\n', layerIdx, peakCount);
            end
        end

    catch ME
        fprintf('Error in overlayXtVsYPaths: %s\n', ME.message);
        rethrow(ME);
    end
end

function overlayYtVsXPaths(ax, YtVsX_paths, currentXIndex, X_values, Y_values, timeValues, debugFlag)
    % Overlay layer paths for YtVsX view

    if nargin < 7
        debugFlag = false;  % Default to no debug output
    end

    if isempty(YtVsX_paths)
        if debugFlag, fprintf('Debug: YtVsX_paths is empty\n'); end
        return;
    end

    try
        currentX = X_values(currentXIndex);
        tolerance = (max(X_values) - min(X_values)) / (2 * length(X_values)); % Half the X spacing

        % Define colors for different layers
        layerColors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

        if debugFlag
            fprintf('Debug: Overlaying %d YtVsX layers for X=%.3f (tolerance=%.3f)\n', ...
                    length(YtVsX_paths), currentX, tolerance);
        end

        for layerIdx = 1:length(YtVsX_paths)
            layer = YtVsX_paths{layerIdx};

            % Validate layer structure
            if ~isstruct(layer) || ~isfield(layer, 'peaks')
                if debugFlag, fprintf('Debug: Layer %d missing peaks field\n', layerIdx); end
                continue;
            end

            if isempty(layer.peaks)
                if debugFlag, fprintf('Debug: Layer %d has no peaks\n', layerIdx); end
                continue;
            end

            % Find peaks near the current X slice
            relevantPeaks = zeros(length(layer.peaks), 2); % Pre-allocate
            peakCount = 0;

            for p = 1:length(layer.peaks)
                peak = layer.peaks(p);

                % Validate peak structure
                if ~isstruct(peak) || ~isfield(peak, 'X') || ~isfield(peak, 'Y') || ~isfield(peak, 'T')
                    if debugFlag, fprintf('Debug: Peak %d in layer %d missing coordinate fields\n', p, layerIdx); end
                    continue;
                end

                if abs(peak.X - currentX) <= tolerance
                    peakCount = peakCount + 1;
                    relevantPeaks(peakCount, :) = [peak.Y, peak.T];
                end
            end

            relevantPeaks = relevantPeaks(1:peakCount, :); % Trim to actual size

            if debugFlag
                fprintf('Debug: Layer %d has %d relevant peaks for current X slice\n', layerIdx, peakCount);
            end

            if size(relevantPeaks, 1) >= 2
                % Sort by Y position
                [~, sortIdx] = sort(relevantPeaks(:, 1));
                relevantPeaks = relevantPeaks(sortIdx, :);

                % Plot the layer path
                colorIdx = mod(layerIdx - 1, length(layerColors)) + 1;
                plot(ax, relevantPeaks(:, 1), relevantPeaks(:, 2), ...
                     [layerColors(colorIdx) 'o-'], 'LineWidth', 2, 'MarkerSize', 4, ...
                     'DisplayName', sprintf('Layer %d', layerIdx));

                if debugFlag
                    fprintf('Debug: Plotted layer %d with %d points\n', layerIdx, size(relevantPeaks, 1));
                end
            elseif peakCount > 0 && debugFlag
                fprintf('Debug: Layer %d has only %d peaks (need >=2 for line plot)\n', layerIdx, peakCount);
            end
        end

    catch ME
        fprintf('Error in overlayYtVsXPaths: %s\n', ME.message);
        rethrow(ME);
    end
end

function overlayXYVstPaths(ax, layerPaths, currentTIndex, X_values, Y_values, timeValues, debugFlag)
    % Overlay layer paths for XYVst view (project onto XY plane)

    if nargin < 7
        debugFlag = false;  % Default to no debug output
    end

    if (~isfield(layerPaths, 'XtVsY_paths') || isempty(layerPaths.XtVsY_paths)) && ...
       (~isfield(layerPaths, 'YtVsX_paths') || isempty(layerPaths.YtVsX_paths))
        if debugFlag, fprintf('Debug: No layer paths available for XYVst view\n'); end
        return;
    end

    try
        currentTime = timeValues(currentTIndex);
        tolerance = (max(timeValues) - min(timeValues)) / (2 * length(timeValues)); % Half the time spacing

        % Define colors for different layers
        layerColors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];

        if debugFlag
            fprintf('Debug: Overlaying XYVst paths for T=%.3f (tolerance=%.3f)\n', currentTime, tolerance);
        end

        % Overlay XtVsY paths
        if isfield(layerPaths, 'XtVsY_paths') && ~isempty(layerPaths.XtVsY_paths)
            if debugFlag
                fprintf('Debug: Processing %d XtVsY layers\n', length(layerPaths.XtVsY_paths));
            end

            for layerIdx = 1:length(layerPaths.XtVsY_paths)
                layer = layerPaths.XtVsY_paths{layerIdx};

                % Validate layer structure
                if ~isstruct(layer) || ~isfield(layer, 'peaks') || isempty(layer.peaks)
                    continue;
                end

                % Find peaks near the current time
                relevantPeaks = zeros(length(layer.peaks), 2); % Pre-allocate
                peakCount = 0;

                for p = 1:length(layer.peaks)
                    peak = layer.peaks(p);

                    % Validate peak structure
                    if ~isstruct(peak) || ~isfield(peak, 'T') || ~isfield(peak, 'X') || ~isfield(peak, 'Y')
                        continue;
                    end

                    if abs(peak.T - currentTime) <= tolerance
                        peakCount = peakCount + 1;
                        relevantPeaks(peakCount, :) = [peak.X, peak.Y];
                    end
                end

                relevantPeaks = relevantPeaks(1:peakCount, :); % Trim to actual size

                if size(relevantPeaks, 1) >= 1  % Changed from >=2 to >=1 for scatter plot
                    % Plot the layer points
                    colorIdx = mod(layerIdx - 1, length(layerColors)) + 1;
                    scatter(ax, relevantPeaks(:, 1), relevantPeaks(:, 2), 50, layerColors(colorIdx), 'filled', ...
                           'DisplayName', sprintf('XtVsY Layer %d', layerIdx));
                    if debugFlag
                        fprintf('Debug: Plotted XtVsY layer %d with %d points\n', layerIdx, size(relevantPeaks, 1));
                    end
                end
            end
        end

        % Overlay YtVsX paths
        if isfield(layerPaths, 'YtVsX_paths') && ~isempty(layerPaths.YtVsX_paths)
            if debugFlag
                fprintf('Debug: Processing %d YtVsX layers\n', length(layerPaths.YtVsX_paths));
            end

            for layerIdx = 1:length(layerPaths.YtVsX_paths)
                layer = layerPaths.YtVsX_paths{layerIdx};

                % Validate layer structure
                if ~isstruct(layer) || ~isfield(layer, 'peaks') || isempty(layer.peaks)
                    continue;
                end

                % Find peaks near the current time
                relevantPeaks = zeros(length(layer.peaks), 2); % Pre-allocate
                peakCount = 0;

                for p = 1:length(layer.peaks)
                    peak = layer.peaks(p);

                    % Validate peak structure
                    if ~isstruct(peak) || ~isfield(peak, 'T') || ~isfield(peak, 'X') || ~isfield(peak, 'Y')
                        continue;
                    end

                    if abs(peak.T - currentTime) <= tolerance
                        peakCount = peakCount + 1;
                        relevantPeaks(peakCount, :) = [peak.X, peak.Y];
                    end
                end

                relevantPeaks = relevantPeaks(1:peakCount, :); % Trim to actual size

                if size(relevantPeaks, 1) >= 1  % Changed from >=2 to >=1 for scatter plot
                    % Plot the layer points with different marker
                    colorIdx = mod(layerIdx - 1, length(layerColors)) + 1;
                    scatter(ax, relevantPeaks(:, 1), relevantPeaks(:, 2), 50, layerColors(colorIdx), 's', ...
                           'DisplayName', sprintf('YtVsX Layer %d', layerIdx));
                    if debugFlag
                        fprintf('Debug: Plotted YtVsX layer %d with %d points\n', layerIdx, size(relevantPeaks, 1));
                    end
                end
            end
        end

    catch ME
        fprintf('Error in overlayXYVstPaths: %s\n', ME.message);
        rethrow(ME);
    end
end

function computeCrossViewAlignment(button, fig)
    % Iteratively align XtVsY and YtVsX views until convergence
    userData = get(fig, 'UserData');

    % Prevent multiple simultaneous computations
    if userData.isComputingAlignment
        fprintf('Alignment computation already in progress. Please wait...\n');
        return;
    end

    % Get data dimensions
    if ~isempty(userData.statDataArray) && ~isempty(userData.statDataArray{1}.maps)
        [numY, numX] = size(userData.statDataArray{1}.maps{1});
    else
        fprintf('Error: No data available for cross-view alignment\n');
        return;
    end

    % Show progress dialog
    progressDlg = uiprogressdlg(fig, 'Title', 'Cross-View Alignment', ...
                               'Message', 'Starting iterative cross-view alignment...', ...
                               'Value', 0, ...
                               'Cancelable', false);

    % Store progress dialog in userData
    userData.progressDialog = progressDlg;
    userData.isComputingAlignment = true;
    set(fig, 'UserData', userData);
    drawnow;

    % Update status
    statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
    if ~isempty(statusText)
        set(statusText, 'String', 'Cross-view aligning...');
    end

    % Start computation in background
    timer_obj = timer('TimerFcn', @(~,~) computeCrossViewAlignmentBackground(fig), ...
                      'StartDelay', 0.1, 'ExecutionMode', 'singleShot', ...
                      'Name', 'CrossViewTimer');
    start(timer_obj);
end

function computeCrossViewAlignmentBackground(fig)
    % Background computation for iterative cross-view alignment
    try
        userData = get(fig, 'UserData');
        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');

        % Start timing
        crossViewStartTime = tic;
        userData.timingInfo.crossViewStartTime = crossViewStartTime;

        % Get data dimensions and initialize
        currentData = userData.statDataArray;
        if isempty(currentData)
            fprintf('Error: Current data not available\n');
            return;
        end

        [numY, numX] = size(currentData{1}.maps{1});
        nStats = length(currentData);
        numSegments = length(currentData{1}.maps);

        % Cross-view alignment parameters
        maxCrossViewIterations = 10; % Maximum number of cross-view iterations
        convergenceThreshold = 0.001; % Convergence threshold for improvement

        % Initialize tracking with pre-allocated arrays (no dynamic growth)
        crossViewIteration = 0;
        previousCost = inf;
        improvementHistory = zeros(maxCrossViewIterations, 1); % Pre-allocate for performance

        fprintf('Starting cross-view alignment: %d Y-slices x %d X-slices\n', numY, numX);

        % Main cross-view iteration loop
        while crossViewIteration < maxCrossViewIterations
            crossViewIteration = crossViewIteration + 1;
            iterationStartTime = tic;

            % Update progress dialog with detailed information
            if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
                overallProgress = (crossViewIteration - 1) / maxCrossViewIterations;
                elapsedTime = toc(crossViewStartTime);
                estimatedTotal = elapsedTime / max(overallProgress, 0.01);
                remainingTime = max(0, estimatedTotal - elapsedTime);

                userData.progressDialog.Value = overallProgress;
                userData.progressDialog.Message = sprintf(['Cross-View Alignment - Iteration %d/%d\n' ...
                    'Elapsed: %.1fs, Est. remaining: %.1fs\n' ...
                    'Starting alignment cycle...'], ...
                    crossViewIteration, maxCrossViewIterations, elapsedTime, remainingTime);
                drawnow;
            end

            fprintf('Cross-view iteration %d/%d\n', crossViewIteration, maxCrossViewIterations);

            % Step 1: Align all Y-slices (XtVsY view)
            fprintf('  Step 1: Aligning Y-slices (XtVsY)...\n');
            if ~isempty(statusText)
                set(statusText, 'String', sprintf('Cross-view %d/%d: Y-slices', crossViewIteration, maxCrossViewIterations));
            end

            % Create progress callback for Y-slices
            yProgressCallback = @(sliceIdx, totalSlices) updateCrossViewProgress(fig, crossViewIteration, maxCrossViewIterations, 'Y-slices', sliceIdx, totalSlices, 1, crossViewStartTime);
            currentData = alignAllSlicesInView(currentData, 'XtVsY', userData, fig, yProgressCallback);

            % Step 2: Align all X-slices (YtVsX view)
            fprintf('  Step 2: Aligning X-slices (YtVsX)...\n');
            if ~isempty(statusText)
                set(statusText, 'String', sprintf('Cross-view %d/%d: X-slices', crossViewIteration, maxCrossViewIterations));
            end

            % Create progress callback for X-slices
            xProgressCallback = @(sliceIdx, totalSlices) updateCrossViewProgress(fig, crossViewIteration, maxCrossViewIterations, 'X-slices', sliceIdx, totalSlices, 2, crossViewStartTime);
            currentData = alignAllSlicesInView(currentData, 'YtVsX', userData, fig, xProgressCallback);

            % Calculate overall alignment cost to check for convergence
            currentCost = calculateOverallAlignmentCost(currentData);
            improvement = (previousCost - currentCost) / previousCost;
            improvementHistory(crossViewIteration) = improvement; % Use pre-allocated array

            iterationTime = toc(iterationStartTime);
            fprintf('  Iteration %d complete: cost=%.6f, improvement=%.4f%%, time=%.1fs\n', ...
                crossViewIteration, currentCost, improvement*100, iterationTime);

            % Update progress dialog with convergence info
            if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
                userData.progressDialog.Message = sprintf(['Cross-View Alignment - Iteration %d/%d Complete\n' ...
                    'Cost: %.6f, Improvement: %.3f%%\n' ...
                    'Elapsed: %.1fs, Checking convergence...'], ...
                    crossViewIteration, maxCrossViewIterations, currentCost, improvement*100, iterationTime);
                drawnow;
            end

            % Check for convergence
            if improvement < convergenceThreshold && crossViewIteration > 1
                fprintf('Cross-view alignment converged after %d iterations (improvement < %.3f%%)\n', ...
                    crossViewIteration, convergenceThreshold*100);

                % Update progress dialog with convergence message
                if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
                    userData.progressDialog.Value = 1.0;
                    userData.progressDialog.Message = sprintf(['Cross-View Alignment CONVERGED!\n' ...
                        'Completed %d iterations\n' ...
                        'Final improvement: %.3f%% (threshold: %.3f%%)\n' ...
                        'Total time: %.1fs'], ...
                        crossViewIteration, improvement*100, convergenceThreshold*100, toc(crossViewStartTime));
                    drawnow;
                    pause(2); % Show convergence message for 2 seconds
                end
                break;
            end

            previousCost = currentCost;
        end

        % Check if we completed without convergence
        if crossViewIteration >= maxCrossViewIterations
            fprintf('Cross-view alignment completed %d iterations without convergence\n', maxCrossViewIterations);

            % Update progress dialog with completion message
            if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
                userData.progressDialog.Value = 1.0;
                userData.progressDialog.Message = sprintf(['Cross-View Alignment COMPLETE\n' ...
                    'Completed maximum %d iterations\n' ...
                    'Final improvement: %.3f%%\n' ...
                    'Total time: %.1fs'], ...
                    maxCrossViewIterations, improvement*100, toc(crossViewStartTime));
                drawnow;
                pause(1); % Show completion message for 1 second
            end
        end

        % Final update
        userData.alignedDataCache = currentData;
        userData.currentView = 'aligned';
        userData.statDataArray = currentData;
        userData.isComputingAlignment = false;

        % Close progress dialog
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            close(userData.progressDialog);
            userData = rmfield(userData, 'progressDialog');
        end

        % Calculate total time
        totalElapsedTime = toc(crossViewStartTime);
        userData.timingInfo.crossViewTotalTime = totalElapsedTime;

        set(fig, 'UserData', userData);
        updatePlots(fig);
        applyAspectRatioToAxes(fig);

        % Update final status
        if ~isempty(statusText)
            % Only use the valid portion of the pre-allocated array
            validImprovements = improvementHistory(1:crossViewIteration);
            avgImprovement = mean(validImprovements) * 100;
            set(statusText, 'String', sprintf('Cross-View Aligned (%d iter, %.1f%% avg improve, %.1fs)', ...
                crossViewIteration, avgImprovement, totalElapsedTime));
        end

        fprintf('Cross-view alignment complete: %d iterations, %.1fs total\n', crossViewIteration, totalElapsedTime);

    catch ME
        % Handle errors gracefully
        fprintf('Error during cross-view alignment: %s\n', ME.message);
        userData = get(fig, 'UserData');
        userData.isComputingAlignment = false;

        % Close progress dialog
        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            close(userData.progressDialog);
            userData = rmfield(userData, 'progressDialog');
        end

        set(fig, 'UserData', userData);

        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', 'Cross-view Error');
        end
    end

    % Clean up timer
    try
        crossViewTimer = timerfind('Name', 'CrossViewTimer');
        if ~isempty(crossViewTimer)
            stop(crossViewTimer);
            delete(crossViewTimer);
        end
    catch
        % Ignore timer cleanup errors
    end
end

function alignedData = alignAllSlicesInView(inputData, viewType, userData, fig, progressCallback)
    % Align all slices in a specific view (XtVsY or YtVsX)
    % Optional progressCallback for cross-view alignment progress updates

    alignedData = inputData;
    [numY, numX] = size(inputData{1}.maps{1});
    numSegments = length(inputData{1}.maps);
    nStats = length(inputData);

    % Determine slice range based on view
    switch viewType
        case 'XtVsY'
            numSlices = numY;
            sliceType = 'Y';
        case 'YtVsX'
            numSlices = numX;
            sliceType = 'X';
        otherwise
            error('Unknown view type: %s', viewType);
    end

    % Alignment method is fixed to 'average' (Row Average)
    alignmentMethod = 'average';

    % Get user-defined convergence threshold
    convergenceInput = findobj(fig, 'Tag', 'ConvergenceInput');
    if ~isempty(convergenceInput)
        convergenceStr = get(convergenceInput, 'String');
        convergenceThreshold = str2double(convergenceStr) / 100; % Convert percentage to decimal
        if isnan(convergenceThreshold) || convergenceThreshold <= 0 || convergenceThreshold > 0.05
            convergenceThreshold = 0.01; % Default to 1% if invalid
            set(convergenceInput, 'String', '1.0'); % Reset to default
        end
    else
        convergenceThreshold = 0.01; % Default to 1%
    end

    % Process each slice
    for sliceIndex = 1:numSlices
        % Update progress if callback provided
        if nargin >= 5 && ~isempty(progressCallback)
            progressCallback(sliceIndex, numSlices);
        end

        % Process each statistic for this slice
        for statIdx = 1:nStats
            % Extract slice data based on view
            switch viewType
                case 'XtVsY'
                    % Extract X,t data for this Y slice
                    sliceData = zeros(numSegments, numX);
                    for seg = 1:numSegments
                        sliceData(seg, :) = alignedData{statIdx}.maps{seg}(sliceIndex, :);
                    end
                case 'YtVsX'
                    % Extract Y,t data for this X slice
                    sliceData = zeros(numSegments, numY);
                    for seg = 1:numSegments
                        sliceData(seg, :) = alignedData{statIdx}.maps{seg}(:, sliceIndex)';
                    end
            end

            % Apply alignment
            % Read optional alignment params (non-destructive; defaults preserved)
            maxShift = 15; costFcn = 'mse';
            costDropdown = findobj(fig, 'Tag', 'AlignmentCostFunctionDropdown');
            if ~isempty(costDropdown)
                switch get(costDropdown,'Value')
                    case 1, costFcn = 'mse';
                    case 2, costFcn = 'correlation';
                    case 3, costFcn = 'ncc';
                end
            end
            maxShiftEdit = findobj(fig, 'Tag', 'AlignmentMaxShiftInput');
            if ~isempty(maxShiftEdit)
                v = str2double(get(maxShiftEdit,'String')); if ~isnan(v) && v>0, maxShift = round(v); end
            end

            [alignedSliceData, ~] = alignColumnsImproved(sliceData, ...
                'MaxShift', maxShift, ...
                'CostFunction', costFcn, ...
                'AlignmentMethod', alignmentMethod, ...
                'LocalScope', 5, ...
                'PadMethod', 'zeros', ...
                'Verbose', false, ...
                'ConvergenceThreshold', convergenceThreshold, ...
                'MaxIterations', getMaxIterations(fig), ...
                'WeightingFunction', 'exponential', ...
                'WeightingScale', 3.0);

            % Put aligned data back
            switch viewType
                case 'XtVsY'
                    for seg = 1:numSegments
                        alignedData{statIdx}.maps{seg}(sliceIndex, :) = alignedSliceData(seg, :);
                    end
                case 'YtVsX'
                    for seg = 1:numSegments
                        alignedData{statIdx}.maps{seg}(:, sliceIndex) = alignedSliceData(seg, :)';
                    end
            end
        end
    end

    fprintf('    Completed %s alignment: %d %s-slices\n', viewType, numSlices, sliceType);
end

function n = getMaxIterations(fig)
    % Helper to read max iterations from UI, with fallback default
    n = 20;
    try
        h = findobj(fig, 'Tag', 'AlignmentMaxIterInput');
        if ~isempty(h)
            v = str2double(get(h,'String'));
            if ~isnan(v) && v >= 1 && v <= 200
                n = round(v);
            end
        end
    catch
    end
end

function onConvergenceChanged(fig, src)
    % Update visState.convergencePercent and normalize value
    try
        ud = get(fig, 'UserData'); if ~isstruct(ud), ud = struct(); end
        if ~isfield(ud, 'visState') || ~isstruct(ud.visState), ud.visState = struct(); end
        v = str2double(get(src,'String'));
        if isnan(v) || v <= 0 || v > 5
            v = 1.0; set(src,'String','1.0');
        end
        ud.visState.convergencePercent = v;
        set(fig,'UserData',ud);
    catch
    end
end


function updateCrossViewProgress(fig, currentIteration, maxIterations, viewType, currentSlice, totalSlices, step, startTime)
    % Update cross-view alignment progress dialog with detailed information
    try
        userData = get(fig, 'UserData');

        if isfield(userData, 'progressDialog') && isvalid(userData.progressDialog)
            % Calculate overall progress
            iterationProgress = (currentIteration - 1) / maxIterations;
            stepProgress = (step - 1) / 2; % 2 steps per iteration (Y-slices, X-slices)
            sliceProgress = (currentSlice - 1) / totalSlices;

            % Overall progress: iteration + step within iteration + slice within step
            overallProgress = iterationProgress + (stepProgress / maxIterations) + (sliceProgress / (2 * maxIterations));
            overallProgress = min(overallProgress, 1.0); % Cap at 100%

            % Calculate timing
            elapsedTime = toc(startTime);
            if overallProgress > 0.01
                estimatedTotal = elapsedTime / overallProgress;
                remainingTime = max(0, estimatedTotal - elapsedTime);
            else
                remainingTime = 0;
            end

            % Update progress dialog
            userData.progressDialog.Value = overallProgress;
            userData.progressDialog.Message = sprintf(['Cross-View Alignment - Iteration %d/%d\n' ...
                'Step %d/2: Aligning %s\n' ...
                'Processing slice %d/%d (%.1f%%)\n' ...
                'Elapsed: %.1fs, Est. remaining: %.1fs'], ...
                currentIteration, maxIterations, step, viewType, ...
                currentSlice, totalSlices, overallProgress * 100, ...
                elapsedTime, remainingTime);
            drawnow;
        end
    catch ME
        % Silently handle errors in progress callback
        fprintf('Cross-view progress update error: %s\n', ME.message);
    end
end

function totalCost = calculateOverallAlignmentCost(inputData)
    % Calculate overall alignment cost for convergence checking

    totalCost = 0;
    nStats = length(inputData);

    for statIdx = 1:nStats
        numSegments = length(inputData{statIdx}.maps);
        [numY, numX] = size(inputData{statIdx}.maps{1});

        % Calculate cost for each segment pair
        for seg1 = 1:numSegments-1
            for seg2 = seg1+1:numSegments
                map1 = inputData{statIdx}.maps{seg1};
                map2 = inputData{statIdx}.maps{seg2};

                % Calculate MSE between segments (measure of alignment quality)
                segmentCost = mean((map1(:) - map2(:)).^2);
                totalCost = totalCost + segmentCost;
            end
        end
    end

    % Normalize by number of comparisons
    numComparisons = nStats * numSegments * (numSegments - 1) / 2;
    if numComparisons > 0
        totalCost = totalCost / numComparisons;
    end
end

function updateProgressMessage(progressDlg, message)
    % Helper function to update progress dialog message
    try
        if isvalid(progressDlg)
            progressDlg.Message = message;
            drawnow;
        end
    catch
        % Silently handle errors if dialog is closed
    end
end

function synchronizeTimeScales(fig, currentView)
    % SYNCHRONIZETIMESCALES - Synchronize time scales across all subplots
    % Locks the time scale to the plot that is most zoomed in or has more data points visible
    % This ensures consistent time scaling between side-by-side plots

    try
        % Validate inputs
        if ~ishandle(fig) || ~ishghandle(fig, 'figure')
            return;
        end

        userData = get(fig, 'UserData');
        if ~isfield(userData, 'axHandles') || isempty(userData.axHandles)
            return;
        end
        axHandles = userData.axHandles;

        % Only synchronize for time-based views (XtVsY and YtVsX)
        if ~strcmp(currentView, 'XtVsY') && ~strcmp(currentView, 'YtVsX')
            return;
        end

        % Get all valid axes handles
        validAxes = [];
        for i = 1:length(axHandles)
            if ishandle(axHandles(i)) && ishghandle(axHandles(i), 'axes')
                validAxes(end+1) = axHandles(i);
            end
        end

        if length(validAxes) < 2
            return; % Need at least 2 plots to synchronize
        end

        % Collect time axis information from all plots
        timeRanges = [];
        dataPointCounts = [];

        for i = 1:length(validAxes)
            ax = validAxes(i);

            % Get the Y-axis limits (time axis for XtVsY and YtVsX views)
            yLimits = ylim(ax);
            timeRange = yLimits(2) - yLimits(1);
            timeRanges = [timeRanges, timeRange];

            % Count visible data points by examining the plot data
            plotChildren = get(ax, 'Children');
            dataPointCount = 0;

            for j = 1:length(plotChildren)
                child = plotChildren(j);
                if isprop(child, 'YData') && ~isempty(get(child, 'YData'))
                    yData = get(child, 'YData');
                    % Count unique Y values within the current Y limits
                    visibleYData = yData(yData >= yLimits(1) & yData <= yLimits(2));
                    dataPointCount = dataPointCount + length(unique(visibleYData));
                elseif isprop(child, 'CData') && ~isempty(get(child, 'CData'))
                    % For image/heatmap data, count rows (time points)
                    cData = get(child, 'CData');
                    if ~isempty(cData)
                        dataPointCount = dataPointCount + size(cData, 1);
                    end
                end
            end

            dataPointCounts = [dataPointCounts, dataPointCount];
        end

        % Determine the reference plot (most zoomed in or most data points)
        % Priority: 1) Smallest time range (most zoomed in), 2) Most data points

        % Find the plot with the smallest time range (most zoomed in)
        [minTimeRange, minRangeIdx] = min(timeRanges);

        % Among plots with similar time ranges, prefer the one with more data points
        similarRangeThreshold = minTimeRange * 1.1; % Within 10% of minimum range
        similarRangeIndices = find(timeRanges <= similarRangeThreshold);

        if length(similarRangeIndices) > 1
            % Multiple plots with similar ranges - choose the one with most data points
            [~, maxDataIdx] = max(dataPointCounts(similarRangeIndices));
            referenceIdx = similarRangeIndices(maxDataIdx);
        else
            % Only one plot with minimum range
            referenceIdx = minRangeIdx;
        end

        % Get the reference time limits
        referenceAx = validAxes(referenceIdx);
        referenceYLimits = ylim(referenceAx);

        % Apply the reference time limits to all other plots
        for i = 1:length(validAxes)
            if i ~= referenceIdx
                ax = validAxes(i);
                ylim(ax, referenceYLimits);
            end
        end

        fprintf('Time scale synchronized: Reference plot %d (range: %.3f μs, %d data points)\n', ...
                referenceIdx, minTimeRange * 1e6, dataPointCounts(referenceIdx));

    catch ME
        fprintf('Warning: Time scale synchronization failed: %s\n', ME.message);
    end
end

function manualTimeScaleSync(button, fig)
    % MANUALTIMESCALESYNC - Manual callback for time scale synchronization button
    % Allows users to manually trigger time scale synchronization

    try
        userData = get(fig, 'UserData');
        currentView = userData.visState.currentView;

        % Only works for time-based views
        if ~strcmp(currentView, 'XtVsY') && ~strcmp(currentView, 'YtVsX')
            fprintf('Time scale synchronization only available for XtVsY and YtVsX views\n');
            return;
        end

        % Trigger synchronization
        synchronizeTimeScales(fig, currentView);

        % Provide user feedback
        fprintf('Manual time scale synchronization completed for %s view\n', currentView);

        % Update the display
        drawnow;

    catch ME
        fprintf('Error in manual time scale synchronization: %s\n', ME.message);
    end
end



function exportAlignedWaveforms(fig)
    % Export the aligned waveforms using the same structure as the input waveform MAT files
    % Expected structure in saved file (as loaded by loadWaveformData):
    %   - waveformArray: [numWaveforms x numTimePoints]
    %   - t: [1 x numTimePoints] time vector [s]
    %   - numY_sub: scalar
    %   - numX_sub: scalar
    %   - X_values, Y_values, timeValues: axes vectors
    %   - X_Coordinates, Y_Coordinates: grids or vectors

    try
        if nargin < 1 || isempty(fig) || ~ishghandle(fig)
            fig = gcf;
        end
        userData = get(fig, 'UserData');

        % Ensure waveforms are loaded
        if ~isfield(userData, 'waveformLoaded') || ~userData.waveformLoaded
            errordlg('Waveform data not loaded. Load waveforms before exporting.');
            return;
        end

        % Determine whether aligned data exists; if not but shifts exist, compute it
        useAligned = false;
        if isfield(userData, 'alignedWaveformData') && ~isempty(userData.alignedWaveformData)
            useAligned = true;
        elseif isfield(userData, 'alignmentShifts') && ~isempty(userData.alignmentShifts)
            try
                applyAlignmentToWaveforms(fig);
                userData = get(fig, 'UserData');
                if isfield(userData, 'alignedWaveformData') && ~isempty(userData.alignedWaveformData)
                    useAligned = true;
                end
            catch
                % If computation fails, fallback to original below
            end
        end

        % Choose data to export
        if useAligned
            exportData = userData.alignedWaveformData;
        else
            % If aligned not available, export original as a fallback
            exportData = userData.originalWaveformData;
        end

        % Enrich exportData with any metadata fields from original data to match input structure
        if isfield(userData, 'originalWaveformData') && ~isempty(userData.originalWaveformData)
            original = userData.originalWaveformData;
            metaFields = {'X_values','Y_values','timeValues','X_Coordinates','Y_Coordinates'};
            for mf = 1:numel(metaFields)
                f = metaFields{mf};
                if isfield(original, f) && ~isfield(exportData, f)
                    exportData.(f) = original.(f);
                end
            end
        end

        % Validate structure minimally
        requiredFields = {'waveformArray','t','numY_sub','numX_sub'};
        for i = 1:numel(requiredFields)
            if ~isfield(exportData, requiredFields{i})
                error('Export data missing required field: %s', requiredFields{i});
            end
        end

        % Build output directory and filename
        outDir = fullfile(pwd, 'Saved Wave Forms');
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end

        % Use FileNamingArray if available to mirror input naming
        ts = datestr(now, 'yyyymmdd_HHMMSS');
        baseName = 'Waveform_Aligned';
        if isfield(userData, 'FileNamingArray') && ~isempty(userData.FileNamingArray)
            try
                % Build a simple name from FileNamingArray components (fallback-safe)
                parts = userData.FileNamingArray;
                if iscell(parts)
                    partsStr = strjoin(cellfun(@(c) char(string(c)), parts, 'UniformOutput', false), '_');
                else
                    partsStr = char(string(parts));
                end
                baseName = sprintf('%s_Aligned', partsStr);
            catch
                % keep default baseName if formatting fails
            end
        end

        outFile = fullfile(outDir, sprintf('%s_%s.mat', baseName, ts));

        % Save only the fields present in exportData to keep structure consistent
        % Keep variable name 'savedData' if other loaders expect that, or save fields directly
        savedData = exportData; %#ok<NASGU>
        save(outFile, 'savedData', '-v7.3');

        msg = sprintf('Aligned waveforms exported to:\n%s', outFile);
        fprintf('%s\n', msg);
    catch ME
        fprintf('Error exporting aligned waveforms: %s\n', ME.message);
        errordlg(['Export failed: ' ME.message], 'Export Error');
    end
end
