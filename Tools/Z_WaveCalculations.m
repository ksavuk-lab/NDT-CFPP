%% Z_WaveCalculations - Extract and visualize waveforms from specific locations
%
% This function extracts waveforms from user-defined locations in the dataset
% and visualizes them with interactive toggle controls.

%% --- User-Defined Settings ---

clc, clear, %close all;

% Add directories to path for utility functions
addpath(fullfile(pwd, 'Utils'));
addpath(fullfile(pwd, 'Data and Computation Scripts'));
if exist('addUtilsPath', 'file')
    addUtilsPath();
end

% I. Dataset Selection
DATASET = 4; % Choose from 1 to 5
% 1: L0P5S25_10Mhz  (No Defect, Sample 25)
% 2: L8P12S5_10Mhz  (8 Degree Defect, Sample 5)
% 3: L16P1S9_10MHz  (16 Degree Defect, Sample 9)
% 4: L16P1S10_10MHz (16 Degree Defect, Sample 10)
% 5: L16P34S3_10MHz (16 Degree Defect, Sample 3)

% II. Data Processing Options
TrimWaveData  = 1; % 1 = trim spatial data, 0 = keep all
TrimTimeRange = 1; % 1 = trim time data, 0 = keep all

% II.a Data Alignment
AlignWaveFormsAtFirstPeak = 1; % 1 = align at first peak, 0 = keep original
Align_Between_Time_Range = [0.09e-6 0.18e-6]; % Find first peak between these time ranges.
AlignPeakType = 'positive';  % Options: 'positive', 'negative', or 'both'
AlignDiagnostics = false;  % Enable diagnostic logging for alignment
AlignDiagnosticSamples = 50;  % Number of random waveforms to sample for diagnostics

% IV. Peak Detection Options
AnalyzePeaks = 0; % 1 = analyze peaks, 0 = don't analyze peaks
PeakOptions = struct();
PeakOptions.minPeakHeight = 0.02; % Minimum peak height as fraction of max amplitude (lower value to detect more peaks)
PeakOptions.minPeakDistance = 2; % Minimum distance between peaks in samples (smaller value to detect closer peaks)
PeakOptions.minPeakProminence = 0.02; % Minimum peak prominence as fraction of max amplitude
PeakOptions.useAbsoluteValue = true; % Use absolute value of waveform for peak detection
PeakOptions.findValleys = true; % Also find valleys (negative peaks)

% III. Trim Ranges
xRange = [20 70]; % X range in mm
yRange = [5 20]; % Y range in mm
TimeRange = [0.0e-6 2e-6]; % Time range in seconds

%% --- Define Waveform Extraction Locations ---

% Extract 5 waveforms from the left side
Left_wave_1 = [10, 20];
Left_wave_2 = [15, 20];
Left_wave_3 = [20, 20];
Left_wave_4 = [25, 20];
Left_wave_5 = [30, 20];

% Extract 5 waveforms from the right side
Right_wave_1 = [50, 20];
Right_wave_2 = [55, 20];
Right_wave_3 = [60, 20];
Right_wave_4 = [65, 20];
Right_wave_5 = [70, 20];

% Extract 2 waveforms from the Center left side
Center_Left_wave_1 = [37.5, 20];
Center_Left_wave_2 = [39.0, 20];

% Extract 2 waveforms from the Center Right side
Center_Right_wave_1 = [41.0, 20];
Center_Right_wave_2 = [42.5, 20];

% Extract Center Wave (center of trimmed ranges)
Center_wave = [40, 20]; % Center of xRange and yRange

%% --- Data Loading and Preprocessing ---

disp('Loading dataset...');
[DataStructure] = loadData(DATASET);

% Get raw data directly
waveform3DMatrix = DataStructure.waveform3DMatrix;
waveformFull = waveform3DMatrix;

% Generate time vector for alignment
SampRate = DataStructure.sampRate * 1e6;  % Convert to Hz
dt = 1 / SampRate;
tSize = DataStructure.tSize;
t_full = dt:dt:(tSize * dt);

% Align waveforms if enabled (before trimming)
if AlignWaveFormsAtFirstPeak == 1
    % Check if time range for alignment is specified
    if exist('Align_Between_Time_Range', 'var') && ~isempty(Align_Between_Time_Range) && length(Align_Between_Time_Range) == 2
        % Check if peak type is specified
        if exist('AlignPeakType', 'var') && ~isempty(AlignPeakType)
            if exist('AlignDiagnostics', 'var') && ~isempty(AlignDiagnostics)
                if exist('AlignDiagnosticSamples', 'var') && ~isempty(AlignDiagnosticSamples)
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, Align_Between_Time_Range, AlignPeakType, AlignDiagnostics, AlignDiagnosticSamples);
                else
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, Align_Between_Time_Range, AlignPeakType, AlignDiagnostics);
                end
            else
                [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, Align_Between_Time_Range, AlignPeakType);
            end
        else
            % Use default peak type (positive)
            if exist('AlignDiagnostics', 'var') && ~isempty(AlignDiagnostics)
                if exist('AlignDiagnosticSamples', 'var') && ~isempty(AlignDiagnosticSamples)
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, Align_Between_Time_Range, 'positive', AlignDiagnostics, AlignDiagnosticSamples);
                else
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, Align_Between_Time_Range, 'positive', AlignDiagnostics);
                end
            else
                [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, Align_Between_Time_Range);
            end
        end
    else
        % No time range specified
        if exist('AlignPeakType', 'var') && ~isempty(AlignPeakType)
            if exist('AlignDiagnostics', 'var') && ~isempty(AlignDiagnostics)
                if exist('AlignDiagnosticSamples', 'var') && ~isempty(AlignDiagnosticSamples)
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, [], AlignPeakType, AlignDiagnostics, AlignDiagnosticSamples);
                else
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, [], AlignPeakType, AlignDiagnostics);
                end
            else
                [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, [], AlignPeakType);
            end
        else
            % Use default peak type (positive)
            if exist('AlignDiagnostics', 'var') && ~isempty(AlignDiagnostics)
                if exist('AlignDiagnosticSamples', 'var') && ~isempty(AlignDiagnosticSamples)
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, [], 'positive', AlignDiagnostics, AlignDiagnosticSamples);
                else
                    [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, [], 'positive', AlignDiagnostics);
                end
            else
                [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0);
            end
        end
    end

    % Update the original matrix with aligned data
    waveform3DMatrix = waveformFull;
    DataStructure.waveform3DMatrix = waveformFull;
end

disp('Extracting and trimming data...');
[SampRate, ScanRes, indexRes, xSize, ySize, tSize, waveform3DMatrix, ~, t, X_Coordinates, Y_Coordinates] = ...
    ExtractAndTrimData(DataStructure, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange);

% Update waveformFull with trimmed data
waveformFull = waveform3DMatrix;

%% --- Extract Waveforms from Specified Locations ---

% Function to find the nearest index for a given coordinate
findNearestIndex = @(coord, coordArray) find(abs(coordArray - coord) == min(abs(coordArray - coord)), 1);

% Store all waveform locations for easier processing
waveLocations = {
    'Left 1', Left_wave_1;
    'Left 2', Left_wave_2;
    'Left 3', Left_wave_3;
    'Left 4', Left_wave_4;
    'Left 5', Left_wave_5;
    'Right 1', Right_wave_1;
    'Right 2', Right_wave_2;
    'Right 3', Right_wave_3;
    'Right 4', Right_wave_4;
    'Right 5', Right_wave_5;
    'Center 1 (Left)', Center_Left_wave_1;
    'Center 2 (Left)', Center_Left_wave_2;
    'Center 3 (Right)', Center_Right_wave_1;
    'Center 4 (Right)', Center_Right_wave_2;
    'Center 5 (Middle)', Center_wave
};

% Extract waveforms
numWaveforms = size(waveLocations, 1);
waveforms = cell(numWaveforms, 1);
actualCoordinates = zeros(numWaveforms, 2);

disp('Extracting waveforms from specified locations...');
for i = 1:numWaveforms
    % Get the location
    location = waveLocations{i, 2};

    % Find nearest indices in X and Y coordinates
    xIdx = findNearestIndex(location(1), X_Coordinates);
    yIdx = findNearestIndex(location(2), Y_Coordinates);

    % Store actual coordinates (might be slightly different from requested due to discretization)
    actualCoordinates(i, :) = [X_Coordinates(xIdx), Y_Coordinates(yIdx)];

    % Extract the waveform
    waveforms{i} = squeeze(waveformFull(yIdx, xIdx, :));

    % Display the extraction information
    disp(['Extracted ' waveLocations{i, 1} ' waveform from coordinates [' ...
          num2str(actualCoordinates(i, 1)) ', ' num2str(actualCoordinates(i, 2)) ']']);
end

%% --- Create Interactive Visualization ---

% Create figure with enough space for plots and controls
fig = figure('Name', 'Waveform Visualization', 'Position', [100, 100, 1200, 800]);

% Create axes for the plot
ax = axes('Position', [0.05, 0.1, 0.7, 0.85]);
hold(ax, 'on');

% Define colors for different groups
leftColors = winter(5);
rightColors = autumn(5);
centerLeftColors = [0, 0.7, 0; 0, 0.5, 0];
centerRightColors = [0.7, 0, 0.7; 0.5, 0, 0.5];
centerColor = [0, 0, 0];

% Plot all waveforms initially
plotHandles = gobjects(numWaveforms, 1);
for i = 1:numWaveforms
    % Determine color based on waveform group
    if i <= 5 % Left waveforms
        color = leftColors(i, :);
    elseif i <= 10 % Right waveforms
        color = rightColors(i-5, :);
    elseif i <= 12 % Center Left waveforms
        color = centerLeftColors(i-10, :);
    elseif i <= 14 % Center Right waveforms
        color = centerRightColors(i-12, :);
    else % Center waveform
        color = centerColor;
    end

    % Plot the waveform
    plotHandles(i) = plot(ax, t * 1e6, waveforms{i}, 'Color', color, 'LineWidth', 1.5);
end

% Set up the plot
xlabel(ax, 'Time (μs)');
ylabel(ax, 'Amplitude');
title(ax, 'Extracted Waveforms');
grid(ax, 'on');
box(ax, 'on');

% Create toggle buttons for each waveform
buttonWidth = 0.15;
buttonHeight = 0.03;
startX = 0.8;
startY = 0.9;
spacing = 0.035;

% Group labels
uicontrol('Style', 'text', 'String', 'Left Waveforms', ...
    'Position', [startX*fig.Position(3), (startY+spacing)*fig.Position(4), buttonWidth*fig.Position(3), buttonHeight*fig.Position(4)], ...
    'FontWeight', 'bold');

uicontrol('Style', 'text', 'String', 'Right Waveforms', ...
    'Position', [startX*fig.Position(3), (startY-5*spacing)*fig.Position(4), buttonWidth*fig.Position(3), buttonHeight*fig.Position(4)], ...
    'FontWeight', 'bold');

uicontrol('Style', 'text', 'String', 'Center Waveforms', ...
    'Position', [startX*fig.Position(3), (startY-10*spacing)*fig.Position(4), buttonWidth*fig.Position(3), buttonHeight*fig.Position(4)], ...
    'FontWeight', 'bold');

% Create toggle buttons
toggleButtons = gobjects(numWaveforms, 1);
for i = 1:numWaveforms
    % Calculate position based on waveform group
    if i <= 5 % Left waveforms
        posY = startY - (i-1)*spacing;
    elseif i <= 10 % Right waveforms
        posY = startY - 6*spacing - (i-6)*spacing;
    else % Center waveforms (including Center Left, Center Right, and Center)
        posY = startY - 11*spacing - (i-11)*spacing;
    end

    % Create the toggle button
    toggleButtons(i) = uicontrol('Style', 'checkbox', ...
        'String', waveLocations{i, 1}, ...
        'Value', 1, ... % Initially checked
        'Position', [startX*fig.Position(3), posY*fig.Position(4), buttonWidth*fig.Position(3), buttonHeight*fig.Position(4)], ...
        'Callback', {@toggleWaveform, i, plotHandles});
end

% Add a legend
legend(ax, waveLocations(:, 1), 'Location', 'eastoutside');

% Display a small map showing the waveform locations in the bottom right corner
mapAx = axes('Position', [0.8, 0.1, 0.15, 0.2]);
hold(mapAx, 'on');

% Plot the boundaries of the trimmed region
plot(mapAx, [xRange(1), xRange(2), xRange(2), xRange(1), xRange(1)], ...
         [yRange(1), yRange(1), yRange(2), yRange(2), yRange(1)], 'k-', 'LineWidth', 1.5);

% Plot the waveform locations
for i = 1:5 % Left waveforms
    plot(mapAx, actualCoordinates(i, 1), actualCoordinates(i, 2), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', leftColors(i, :), 'MarkerEdgeColor', 'k');
end
for i = 6:10 % Right waveforms
    plot(mapAx, actualCoordinates(i, 1), actualCoordinates(i, 2), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', rightColors(i-5, :), 'MarkerEdgeColor', 'k');
end
for i = 11:12 % Center Left waveforms
    plot(mapAx, actualCoordinates(i, 1), actualCoordinates(i, 2), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', centerLeftColors(i-10, :), 'MarkerEdgeColor', 'k');
end
for i = 13:14 % Center Right waveforms
    plot(mapAx, actualCoordinates(i, 1), actualCoordinates(i, 2), 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', centerRightColors(i-12, :), 'MarkerEdgeColor', 'k');
end
% Center waveform
plot(mapAx, actualCoordinates(15, 1), actualCoordinates(15, 2), 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', centerColor, 'MarkerEdgeColor', 'k');

% Set up the map
title(mapAx, 'Waveform Locations');
xlabel(mapAx, 'X (mm)');
ylabel(mapAx, 'Y (mm)');
axis(mapAx, [xRange(1)-1, xRange(2)+1, yRange(1)-1, yRange(2)+1]);
grid(mapAx, 'on');
box(mapAx, 'on');

% Add buttons in a single row below the waveform locations map
% Calculate button widths for three buttons in a row
singleButtonWidth = buttonWidth/3 - 0.01;

% Add "Select All" button
uicontrol('Style', 'pushbutton', ...
    'String', 'Select All', ...
    'Position', [startX*fig.Position(3), 0.05*fig.Position(4), singleButtonWidth*fig.Position(3), buttonHeight*fig.Position(4)], ...
    'Callback', {@selectAll, toggleButtons, plotHandles, 1});

% Add "Deselect All" button
uicontrol('Style', 'pushbutton', ...
    'String', 'Deselect All', ...
    'Position', [(startX+singleButtonWidth+0.005)*fig.Position(3), 0.05*fig.Position(4), singleButtonWidth*fig.Position(3), buttonHeight*fig.Position(4)], ...
    'Callback', {@selectAll, toggleButtons, plotHandles, 0});

% Add "Save Figure" button
uicontrol('Style', 'pushbutton', ...
    'String', 'Save Figure', ...
    'Position', [(startX+2*singleButtonWidth+0.01)*fig.Position(3), 0.05*fig.Position(4), singleButtonWidth*fig.Position(3), buttonHeight*fig.Position(4)], ...
    'Callback', @saveFigureCallback);



% Display completion message
disp('Waveform extraction and visualization complete!');

%% --- Peak Analysis ---

if AnalyzePeaks
    disp('Analyzing peaks in waveforms...');

    % Extract peak data for each waveform
    peakDataArray = cell(numWaveforms, 1);

    for i = 1:numWaveforms
        % Adjust minimum peak height based on the maximum amplitude of this waveform
        currentOptions = PeakOptions;
        currentOptions.minPeakHeight = PeakOptions.minPeakHeight * max(abs(waveforms{i}));

        % Extract peak data
        peakDataArray{i} = extractPeakData(waveforms{i}, t, currentOptions);

        % Display number of peaks found
        disp([waveLocations{i, 1} ': Found ' num2str(height(peakDataArray{i})) ' peaks']);
    end

    % Create figure for peak analysis comparison
    peakFig = figure('Name', 'Peak Analysis Comparison', 'Position', [100, 100, 1200, 800]);

    % Define colors for different locations
    peakColors = lines(numWaveforms);

    % 1. Time Differences Between Consecutive Peaks
    subplot(2, 2, 1);
    hold on;

    % Check if any waveform has more than one peak
    hasPeakPairs = false;
    maxPeaks = 0;

    for i = 1:numWaveforms
        peakData = peakDataArray{i};
        if height(peakData) > 1
            hasPeakPairs = true;
            plot(2:height(peakData), peakData.TimeDifference(2:end) * 1e6, 'o-', ...
                'Color', peakColors(i,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', peakColors(i,:));
            maxPeaks = max(maxPeaks, height(peakData));
        end
    end

    title('Time Differences Between Consecutive Peaks');
    xlabel('Peak Pair');
    ylabel('Time Difference (\mus)');

    if hasPeakPairs && maxPeaks > 1
        xticks(2:maxPeaks);
        legend(waveLocations(:, 1), 'Location', 'best');
    else
        text(0.5, 0.5, 'No peak pairs found', 'HorizontalAlignment', 'center', ...
            'Units', 'normalized', 'FontSize', 14);
    end

    grid on;

    % 2. Attenuation Between Consecutive Peaks
    subplot(2, 2, 2);
    hold on;

    % Check if any waveform has more than one peak
    hasPeakPairs = false;
    maxPeaks = 0;

    for i = 1:numWaveforms
        peakData = peakDataArray{i};
        if height(peakData) > 1
            hasPeakPairs = true;
            plot(2:height(peakData), peakData.Attenuation(2:end), 'o-', ...
                'Color', peakColors(i,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', peakColors(i,:));
            maxPeaks = max(maxPeaks, height(peakData));
        end
    end

    title('Attenuation Between Consecutive Peaks');
    xlabel('Peak Pair');
    ylabel('Attenuation (dB)');

    if hasPeakPairs && maxPeaks > 1
        xticks(2:maxPeaks);
        legend(waveLocations(:, 1), 'Location', 'best');
    else
        text(0.5, 0.5, 'No peak pairs found', 'HorizontalAlignment', 'center', ...
            'Units', 'normalized', 'FontSize', 14);
    end

    grid on;

    % 3. Peak Amplitudes
    subplot(2, 2, 3);
    hold on;

    % Check if any peaks were found
    hasPeaks = false;
    maxPeaks = 0;

    % Add a horizontal line at zero
    plot([0, numWaveforms+1], [0, 0], 'k--', 'LineWidth', 1);

    % Create legend entries
    legendEntries = {};

    for i = 1:numWaveforms
        peakData = peakDataArray{i};
        if height(peakData) > 0
            hasPeaks = true;

            % Separate peaks and valleys for plotting
            isPeak = peakData.PeakType == 1;
            isValley = peakData.PeakType == -1;

            % Plot peaks with upward triangles
            if any(isPeak)
                peakIndices = find(isPeak);
                plot(peakData.PeakNumber(isPeak), peakData.PeakAmplitude(isPeak), '^-', ...
                    'Color', peakColors(i,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', peakColors(i,:));
            end

            % Plot valleys with downward triangles
            if any(isValley)
                valleyIndices = find(isValley);
                plot(peakData.PeakNumber(isValley), peakData.PeakAmplitude(isValley), 'v-', ...
                    'Color', peakColors(i,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', peakColors(i,:));
            end

            maxPeaks = max(maxPeaks, height(peakData));
            legendEntries{end+1} = waveLocations{i, 1};
        end
    end

    title('Peak and Valley Amplitudes');
    xlabel('Peak/Valley Number');
    ylabel('Amplitude');

    if hasPeaks && maxPeaks > 0
        xticks(1:maxPeaks);
        legend(legendEntries, 'Location', 'best');
    else
        text(0.5, 0.5, 'No peaks found', 'HorizontalAlignment', 'center', ...
            'Units', 'normalized', 'FontSize', 14);
    end

    grid on;

    % 4. Peak Times
    subplot(2, 2, 4);
    hold on;

    % Check if any peaks were found
    hasPeaks = false;
    maxPeaks = 0;

    % Create legend entries
    legendEntries = {};

    for i = 1:numWaveforms
        peakData = peakDataArray{i};
        if height(peakData) > 0
            hasPeaks = true;

            % Separate peaks and valleys for plotting
            isPeak = peakData.PeakType == 1;
            isValley = peakData.PeakType == -1;

            % Plot peaks with upward triangles
            if any(isPeak)
                peakIndices = find(isPeak);
                plot(peakData.PeakNumber(isPeak), peakData.PeakTime(isPeak) * 1e6, '^-', ...
                    'Color', peakColors(i,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', peakColors(i,:));
            end

            % Plot valleys with downward triangles
            if any(isValley)
                valleyIndices = find(isValley);
                plot(peakData.PeakNumber(isValley), peakData.PeakTime(isValley) * 1e6, 'v-', ...
                    'Color', peakColors(i,:), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', peakColors(i,:));
            end

            maxPeaks = max(maxPeaks, height(peakData));
            legendEntries{end+1} = waveLocations{i, 1};
        end
    end

    title('Peak and Valley Times');
    xlabel('Peak/Valley Number');
    ylabel('Time (\mus)');

    if hasPeaks && maxPeaks > 0
        xticks(1:maxPeaks);
        legend(legendEntries, 'Location', 'best');
    else
        text(0.5, 0.5, 'No peaks found', 'HorizontalAlignment', 'center', ...
            'Units', 'normalized', 'FontSize', 14);
    end

    grid on;

    % Create summary table
    summaryTable = table();

    % Check if any peaks were found
    totalPeaks = sum(cellfun(@height, peakDataArray));

    if totalPeaks > 0
        for i = 1:numWaveforms
            % Get peak data for this location
            peakData = peakDataArray{i};

            % Skip if no peaks were found for this waveform
            if height(peakData) == 0
                continue;
            end

            % Add location information
            locationColumn = repmat(string(waveLocations{i, 1}), height(peakData), 1);
            xCoordColumn = repmat(actualCoordinates(i, 1), height(peakData), 1);
            yCoordColumn = repmat(actualCoordinates(i, 2), height(peakData), 1);

            % Create a table with location information
            locationTable = table(locationColumn, xCoordColumn, yCoordColumn, ...
                'VariableNames', {'Location', 'X_Coordinate', 'Y_Coordinate'});

            % Combine with peak data
            combinedTable = [locationTable, peakData];

            % Append to summary table
            summaryTable = [summaryTable; combinedTable];
        end

        % Display summary table
        disp('Summary of Peak Data:');
        disp(summaryTable);
    else
        disp('No peaks were found in any waveform with the current settings.');

        % Create an empty table with the correct structure
        summaryTable = table([], [], [], [], [], [], [], [], ...
            'VariableNames', {'Location', 'X_Coordinate', 'Y_Coordinate', 'PeakNumber', 'PeakTime', 'PeakAmplitude', 'TimeDifference', 'Attenuation', 'PeakType'});
    end

    % Save the peak data to a MAT file
    save('PeakAnalysisResults.mat', 'summaryTable', 'peakDataArray', 'waveforms', ...
        't', 'actualCoordinates', 'waveLocations', 'PeakOptions');

    % Export the summary table to a CSV file
    writetable(summaryTable, 'PeakAnalysisResults.csv');

    disp('Peak analysis complete! Results saved to PeakAnalysisResults.mat and PeakAnalysisResults.csv');

    % Create figure for waveform visualization with peaks
    waveformPeaksFig = figure('Name', 'Waveforms and Peaks', 'Position', [100, 100, 1200, 800]);

    % Create a figure-wide legend
    legendHandles = [];
    legendLabels = {};

    % Plot all waveforms with peaks
    for i = 1:numWaveforms
        subplot(numWaveforms, 1, i);

        % Plot the waveform
        waveHandle = plot(t * 1e6, waveforms{i}, 'Color', peakColors(i,:), 'LineWidth', 1.5);
        hold on;

        % Save waveform handle for legend (only for the first subplot)
        if i == 1
            legendHandles = [legendHandles, waveHandle];
            legendLabels = [legendLabels, 'Waveform'];
        end

        % Plot the detected peaks
        peakData = peakDataArray{i};

        % Set up the plot
        title([waveLocations{i, 1} ' Waveform and Peaks/Valleys']);
        xlabel('Time (\mus)');
        ylabel('Amplitude');
        grid on;

        if height(peakData) > 0
            % Separate peaks and valleys
            isPeak = peakData.PeakType == 1;
            isValley = peakData.PeakType == -1;

            % Plot peaks (upward triangles)
            if any(isPeak)
                peakHandle = scatter(peakData.PeakTime(isPeak) * 1e6, peakData.PeakAmplitude(isPeak), 100, 'r', '^', 'filled', 'MarkerEdgeColor', 'k');

                % Save peak handle for legend (only for the first subplot)
                if i == 1
                    legendHandles = [legendHandles, peakHandle];
                    legendLabels = [legendLabels, 'Peaks'];
                end
            end

            % Plot valleys (downward triangles)
            if any(isValley)
                valleyHandle = scatter(peakData.PeakTime(isValley) * 1e6, peakData.PeakAmplitude(isValley), 100, 'b', 'v', 'filled', 'MarkerEdgeColor', 'k');

                % Save valley handle for legend (only for the first subplot)
                if i == 1
                    legendHandles = [legendHandles, valleyHandle];
                    legendLabels = [legendLabels, 'Valleys'];
                end
            end

            % Add peak/valley numbers
            for j = 1:height(peakData)
                text(peakData.PeakTime(j) * 1e6, peakData.PeakAmplitude(j), [' ' num2str(j)], ...
                    'FontSize', 10, 'FontWeight', 'bold');
            end
        else
            % Add text indicating no peaks found
            text(mean(xlim), mean(ylim), 'No peaks found', 'HorizontalAlignment', 'center', ...
                'FontSize', 12, 'Color', 'r');
        end
    end

    % Add a single legend to the figure
    if ~isempty(legendHandles)
        figLegend = legend(legendHandles, legendLabels, 'Orientation', 'horizontal');
        figLegend.Position = [0.3, 0.01, 0.4, 0.03]; % Position at the bottom of the figure
    end
end

%% --- Callback Functions ---

% Toggle waveform visibility
function toggleWaveform(src, ~, idx, plotHandles)
    % Set visibility based on checkbox state
    if src.Value
        set(plotHandles(idx), 'Visible', 'on');
    else
        set(plotHandles(idx), 'Visible', 'off');
    end
end

% Select or deselect all waveforms
function selectAll(~, ~, toggleButtons, plotHandles, value)
    for i = 1:length(toggleButtons)
        % Update checkbox
        toggleButtons(i).Value = value;

        % Update plot visibility
        if value
            set(plotHandles(i), 'Visible', 'on');
        else
            set(plotHandles(i), 'Visible', 'off');
        end
    end
end

% Save figure callback
function saveFigureCallback(~, ~)
    [filename, pathname] = uiputfile({'*.fig', 'MATLAB Figure (*.fig)'; ...
                                     '*.png', 'PNG Image (*.png)'; ...
                                     '*.jpg', 'JPEG Image (*.jpg)'}, ...
                                     'Save Figure As');
    if isequal(filename, 0) || isequal(pathname, 0)
        return; % User canceled
    end

    fullpath = fullfile(pathname, filename);
    [~, ~, ext] = fileparts(fullpath);

    if strcmpi(ext, '.fig')
        savefig(gcf, fullpath);
    else
        print(gcf, fullpath, '-dpng', '-r300');
    end

    disp(['Figure saved to: ' fullpath]);
end

function peakData = extractPeakData(waveform, t, options)
    % EXTRACTPEAKDATA - Extract peak data from a waveform
    %
    % This function extracts peaks from a waveform and calculates various
    % metrics such as peak time, time differences between consecutive peaks,
    % and attenuation between consecutive peaks.
    %
    % Inputs:
    %   waveform - The waveform signal (1D array)
    %   t        - Time vector corresponding to the waveform
    %   options  - (Optional) Structure with the following fields:
    %              .minPeakHeight - Minimum peak height (default: 0.1*max(abs(waveform)))
    %              .minPeakDistance - Minimum distance between peaks in samples (default: 10)
    %              .useAbsoluteValue - Whether to use absolute value of waveform (default: true)
    %
    % Outputs:
    %   peakData - Table with the following columns:
    %              .PeakNumber - Index of the peak
    %              .PeakTime - Time location of the peak (seconds)
    %              .PeakAmplitude - Amplitude of the peak
    %              .TimeDifference - Time difference from previous peak (seconds)
    %              .Attenuation - Attenuation from previous peak (dB)

    % Set default options if not provided
    if nargin < 3
        options = struct();
    end

    if ~isfield(options, 'minPeakHeight')
        options.minPeakHeight = 0.1 * max(abs(waveform));
    end

    if ~isfield(options, 'minPeakDistance')
        options.minPeakDistance = 10; % samples
    end

    if ~isfield(options, 'useAbsoluteValue')
        options.useAbsoluteValue = true;
    end

    % Prepare the waveform for peak detection
    if options.useAbsoluteValue
        waveformForPeaks = abs(waveform);
    else
        waveformForPeaks = waveform;
    end

    % Initialize arrays for peaks and valleys
    allPeakIndices = [];
    allPeakAmplitudes = [];
    peakTypes = []; % 1 for peak, -1 for valley

    % Find positive peaks without limiting the number
    minHeight = options.minPeakHeight * max(abs(waveformForPeaks));

    % Add prominence parameter if specified
    if isfield(options, 'minPeakProminence')
        minProminence = options.minPeakProminence * max(abs(waveformForPeaks));
        [peakAmplitudes, peakIndices] = findpeaks(waveformForPeaks, ...
            'MinPeakHeight', minHeight, ...
            'MinPeakDistance', options.minPeakDistance, ...
            'MinPeakProminence', minProminence);
    else
        [peakAmplitudes, peakIndices] = findpeaks(waveformForPeaks, ...
            'MinPeakHeight', minHeight, ...
            'MinPeakDistance', options.minPeakDistance);
    end

    % Add positive peaks to the arrays
    allPeakIndices = [allPeakIndices; peakIndices];
    allPeakAmplitudes = [allPeakAmplitudes; waveform(peakIndices)];
    peakTypes = [peakTypes; ones(length(peakIndices), 1)]; % Mark as peaks

    % Find valleys (negative peaks) if requested
    if isfield(options, 'findValleys') && options.findValleys
        % Invert the waveform to find valleys
        minHeight = options.minPeakHeight * max(abs(waveform));

        % Add prominence parameter if specified
        if isfield(options, 'minPeakProminence')
            minProminence = options.minPeakProminence * max(abs(waveform));
            [valleyAmplitudes, valleyIndices] = findpeaks(-waveform, ...
                'MinPeakHeight', minHeight, ...
                'MinPeakDistance', options.minPeakDistance, ...
                'MinPeakProminence', minProminence);
        else
            [valleyAmplitudes, valleyIndices] = findpeaks(-waveform, ...
                'MinPeakHeight', minHeight, ...
                'MinPeakDistance', options.minPeakDistance);
        end

        % Add valleys to the arrays
        allPeakIndices = [allPeakIndices; valleyIndices];
        allPeakAmplitudes = [allPeakAmplitudes; waveform(valleyIndices)];
        peakTypes = [peakTypes; -ones(length(valleyIndices), 1)]; % Mark as valleys
    end

    % Check if any peaks or valleys were found
    if isempty(allPeakIndices)
        % Create an empty table with the correct variable names
        peakData = table([], [], [], [], [], [], ...
            'VariableNames', {'PeakNumber', 'PeakTime', 'PeakAmplitude', 'TimeDifference', 'Attenuation', 'PeakType'});

        % Add units to variable descriptions
        peakData.Properties.VariableDescriptions{2} = 'seconds';
        peakData.Properties.VariableDescriptions{3} = 'amplitude units';
        peakData.Properties.VariableDescriptions{4} = 'seconds';
        peakData.Properties.VariableDescriptions{5} = 'dB';
        peakData.Properties.VariableDescriptions{6} = '1=peak, -1=valley';

        return;
    end

    % Sort all peaks and valleys by time (index)
    [sortedIndices, sortIdx] = sort(allPeakIndices);
    sortedAmplitudes = allPeakAmplitudes(sortIdx);
    sortedTypes = peakTypes(sortIdx);

    % Calculate peak times
    peakTimes = t(sortedIndices);

    % Initialize arrays for time differences and attenuation
    numPeaks = length(sortedIndices);
    timeDifferences = zeros(numPeaks, 1);
    attenuations = zeros(numPeaks, 1);

    % Calculate time differences and attenuation between consecutive peaks
    for i = 2:numPeaks
        % Time difference from previous peak
        timeDifferences(i) = peakTimes(i) - peakTimes(i-1);

        % Attenuation from previous peak in dB
        % Avoid division by zero or log of zero/negative
        if sortedAmplitudes(i-1) ~= 0 && sortedAmplitudes(i) ~= 0
            attenuations(i) = 20 * log10(abs(sortedAmplitudes(i)) / abs(sortedAmplitudes(i-1)));
        else
            attenuations(i) = NaN; % Not calculable
        end
    end

    % Create a table with the results
    peakData = table((1:numPeaks)', peakTimes(:), sortedAmplitudes(:), timeDifferences(:), attenuations(:), sortedTypes(:), ...
        'VariableNames', {'PeakNumber', 'PeakTime', 'PeakAmplitude', 'TimeDifference', 'Attenuation', 'PeakType'});

    % Add units to variable descriptions
    peakData.Properties.VariableDescriptions{2} = 'seconds';
    peakData.Properties.VariableDescriptions{3} = 'amplitude units';
    peakData.Properties.VariableDescriptions{4} = 'seconds';
    peakData.Properties.VariableDescriptions{5} = 'dB';
    peakData.Properties.VariableDescriptions{6} = '1=peak, -1=valley';
end
