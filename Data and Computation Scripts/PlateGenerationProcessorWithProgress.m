function plateData = PlateGenerationProcessorWithProgress(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions, fig, progressDlg)
% PlateGenerationProcessorWithProgress - Enhanced plate generation with live progress updates
%
% This function generates plates with real-time visual updates and progress feedback
%
% INPUTS:
%   peakData        - Cell array of peak/valley data from PeakExtractionProcessor
%   X_Coordinates   - X spatial coordinates array
%   Y_Coordinates   - Y spatial coordinates array  
%   numY_sub        - Number of Y subdivisions
%   numX_sub        - Number of X subdivisions
%   plateOptions    - Structure with plate generation settings
%   fig             - Figure handle for live updates
%   progressDlg     - Progress dialog handle
%
% OUTPUTS:
%   plateData       - Structure containing generated plate information

fprintf('\n=== Enhanced Plate Generation with Live Updates ===\n');

% Start timing
tic;

% Initialize output structure
plateData = struct();
plateData.plates = [];
plateData.numPlates = 0;
plateData.processingTime = 0;
plateData.statistics = struct();

% Check if plate generation is enabled
if ~plateOptions.enablePlateGeneration
    if isvalid(progressDlg)
        progressDlg.Message = 'Plate generation is disabled.';
        pause(0.5);
    end
    fprintf('Plate generation is disabled. Skipping processing.\n');
    plateData.processingTime = toc;
    return;
end

% Validate inputs
if isempty(peakData)
    if isvalid(progressDlg)
        progressDlg.Message = 'Error: No peak data provided.';
        pause(0.5);
    end
    fprintf('Error: No peak data provided for plate generation.\n');
    plateData.processingTime = toc;
    return;
end

% Update progress
if isvalid(progressDlg)
    progressDlg.Message = 'Extracting points from waveforms...';
    progressDlg.Indeterminate = 'off';
    progressDlg.Value = 0.1;
    drawnow;
end

fprintf('Starting plate generation with settings:\n');
fprintf('  Amplitude tolerance: ±%.3f\n', plateOptions.amplitudeTolerance);
fprintf('  Time tolerance: ±%.1f μs (±%.2e sec)\n', plateOptions.timeTolerance, plateOptions.timeTolerance * 1e-6);
fprintf('  Minimum points per plate: %d\n', plateOptions.minPointsPerPlate);
fprintf('  Plate type: %s\n', plateOptions.plateType);

% Extract all points from peak data
fprintf('Extracting points from %d waveforms...\n', length(peakData));
allPoints = extractAllPointsWithProgress(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions.plateType, progressDlg);

if isempty(allPoints)
    if isvalid(progressDlg)
        progressDlg.Message = 'Warning: No points extracted.';
        pause(0.5);
    end
    fprintf('Warning: No points extracted for plate generation.\n');
    plateData.processingTime = toc;
    return;
end

fprintf('Extracted %d total points for plate generation.\n', size(allPoints, 1));

% Update progress
if isvalid(progressDlg)
    progressDlg.Message = sprintf('Grouping %d points into plates...', size(allPoints, 1));
    progressDlg.Value = 0.3;
    drawnow;
end

% Group points into plates with live updates
plates = groupPointsIntoPlatesWithProgress(allPoints, plateOptions, fig, progressDlg);

% Update progress
if isvalid(progressDlg)
    progressDlg.Message = 'Filtering plates by size requirements...';
    progressDlg.Value = 0.8;
    drawnow;
end

% Filter plates by minimum point requirement
validPlates = filterPlatesBySize(plates, plateOptions.minPointsPerPlate);

% Update progress
if isvalid(progressDlg)
    progressDlg.Message = 'Finalizing results...';
    progressDlg.Value = 0.9;
    drawnow;
end

% Store results
plateData.plates = validPlates;
plateData.numPlates = length(validPlates);
plateData.processingTime = toc;

% Generate statistics
plateData.statistics = generatePlateStatistics(validPlates, allPoints);

% Final progress update
if isvalid(progressDlg)
    progressDlg.Message = sprintf('Completed! Generated %d plates.', plateData.numPlates);
    progressDlg.Value = 1.0;
    drawnow;
    pause(0.5);
end

% Display results
fprintf('\n=== Plate Generation Results ===\n');
fprintf('Total plates generated: %d\n', plateData.numPlates);
fprintf('Processing time: %.2f seconds\n', plateData.processingTime);
fprintf('Average points per plate: %.1f\n', plateData.statistics.avgPointsPerPlate);
fprintf('Largest plate: %d points\n', plateData.statistics.maxPointsPerPlate);
fprintf('Smallest plate: %d points\n', plateData.statistics.minPointsPerPlate);

end

function allPoints = extractAllPointsWithProgress(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateType, progressDlg)
% Extract all relevant points with progress updates

allPoints = [];
pointCount = 0;
numWaveforms = length(peakData);

% Estimate total points for preallocation
estimatedPoints = numWaveforms * 20; % Rough estimate
tempPoints = zeros(estimatedPoints, 6); % [X, Y, Time, Amplitude, Type, WaveformIndex]

% Process each waveform with progress updates
for waveformIdx = 1:numWaveforms
    % Update progress every 1000 waveforms
    if mod(waveformIdx, 1000) == 0 && isvalid(progressDlg)
        progress = 0.1 + 0.2 * (waveformIdx / numWaveforms);
        progressDlg.Value = progress;
        progressDlg.Message = sprintf('Extracting points: %d/%d waveforms...', waveformIdx, numWaveforms);
        drawnow;
    end
    
    waveformData = peakData{waveformIdx};
    
    if isempty(waveformData)
        continue;
    end
    
    % Calculate spatial coordinates for this waveform
    [yIdx, xIdx] = ind2sub([numY_sub, numX_sub], waveformIdx);
    
    % Validate indices
    if xIdx < 1 || xIdx > numX_sub || yIdx < 1 || yIdx > numY_sub
        continue;
    end
    
    spatialX = X_Coordinates(xIdx);
    spatialY = Y_Coordinates(yIdx);
    
    % Filter by plate type
    if strcmp(lower(plateType), 'peaks')
        validIndices = waveformData.TransitionType == 1;
    elseif strcmp(lower(plateType), 'valleys')
        validIndices = waveformData.TransitionType == -1;
    else % 'both'
        validIndices = true(height(waveformData), 1);
    end
    
    validData = waveformData(validIndices, :);
    numValidPoints = height(validData);
    
    if numValidPoints > 0
        % Store points: [X, Y, Time, Amplitude, Type, WaveformIndex]
        newPoints = [repmat(spatialX, numValidPoints, 1), ...
                    repmat(spatialY, numValidPoints, 1), ...
                    validData.TransitionTime, ...
                    validData.TransitionAmplitude, ...
                    validData.TransitionType, ...
                    repmat(waveformIdx, numValidPoints, 1)];
        
        tempPoints(pointCount+1:pointCount+numValidPoints, :) = newPoints;
        pointCount = pointCount + numValidPoints;
    end
end

% Trim to actual size
allPoints = tempPoints(1:pointCount, :);

end

function plates = groupPointsIntoPlatesWithProgress(allPoints, plateOptions, fig, progressDlg)
% Group points into plates with live visual updates

% Pre-allocate plates array to avoid dynamic growth
maxPossiblePlates = min(10000, floor(size(allPoints, 1) / plateOptions.minPointsPerPlate));
plates = cell(maxPossiblePlates, 1);
usedPoints = false(size(allPoints, 1), 1);
plateCount = 0;
totalPoints = size(allPoints, 1);

% Convert time tolerance to seconds
timeToleranceSeconds = plateOptions.timeTolerance * 1e-6;

fprintf('Grouping points with tolerances: Amp=±%.3f, Time=±%.2e sec\n', ...
        plateOptions.amplitudeTolerance, timeToleranceSeconds);

% Get figure data for live updates
figData = get(fig, 'UserData');

% Process each unused point as a potential plate seed
for i = 1:totalPoints
    % Update progress every 1000 points
    if mod(i, 1000) == 0 && isvalid(progressDlg)
        progress = 0.3 + 0.5 * (i / totalPoints);
        progressDlg.Value = progress;
        progressDlg.Message = sprintf('Processing point %d/%d (found %d plates)...', i, totalPoints, plateCount);
        drawnow;
        
        % Live visual update every 2000 points or every 10 plates
        if (mod(i, 2000) == 0 || mod(plateCount, 10) == 0) && plateCount > 0
            % Check if figure is still valid before updating
            if isvalid(fig) && ishghandle(fig)
                updateLivePlateVisualization(fig, plates(1:plateCount), figData);
            end
        end
    end
    
    if usedPoints(i)
        continue;
    end
    
    % Current point as seed
    seedPoint = allPoints(i, :);
    seedAmp = seedPoint(4);
    seedTime = seedPoint(3);
    seedType = seedPoint(5);
    
    % Find all points within tolerance of this seed
    ampDiff = abs(allPoints(:, 4) - seedAmp);
    timeDiff = abs(allPoints(:, 3) - seedTime);
    typeSame = allPoints(:, 5) == seedType;
    
    % Points that match criteria and aren't already used
    matchingPoints = ~usedPoints & ...
                    ampDiff <= plateOptions.amplitudeTolerance & ...
                    timeDiff <= timeToleranceSeconds & ...
                    typeSame;
    
    if sum(matchingPoints) >= plateOptions.minPointsPerPlate && plateCount < maxPossiblePlates
        plateCount = plateCount + 1;

        % Create plate structure
        platePoints = allPoints(matchingPoints, :);

        plates{plateCount} = struct();
        plates{plateCount}.id = plateCount;
        plates{plateCount}.points = platePoints;
        plates{plateCount}.numPoints = size(platePoints, 1);
        plates{plateCount}.centerAmplitude = mean(platePoints(:, 4));
        plates{plateCount}.centerTime = mean(platePoints(:, 3));
        plates{plateCount}.amplitudeRange = [min(platePoints(:, 4)), max(platePoints(:, 4))];
        plates{plateCount}.timeRange = [min(platePoints(:, 3)), max(platePoints(:, 3))];
        plates{plateCount}.spatialBounds = [min(platePoints(:, 1)), max(platePoints(:, 1)), ...
                                           min(platePoints(:, 2)), max(platePoints(:, 2))];
        plates{plateCount}.plateType = seedType; % 1 for peaks, -1 for valleys

        % Mark these points as used
        usedPoints(matchingPoints) = true;
    end
end

% Trim plates array to actual size
plates = plates(1:plateCount);

% Final live update
if isvalid(fig) && ishghandle(fig)
    updateLivePlateVisualization(fig, plates, figData);
end

fprintf('Completed grouping. Generated %d candidate plates.\n', plateCount);

end

function updateLivePlateVisualization(fig, plates, figData)
% Update the 3D plot to show plates as they are generated

% Validate figure handle first
if ~isvalid(fig) || ~ishghandle(fig)
    return;
end

if isempty(plates) || ~isfield(figData, 'mainAx') || ~isvalid(figData.mainAx)
    return;
end

try
    % Clear previous plate visualization and redraw
    cla(figData.mainAx);
    hold(figData.mainAx, 'on');

    % First plot original peak/valley data (faded)
    if isfield(figData, 'allPeakX') && ~isempty(figData.allPeakX)
        scatter3(figData.mainAx, figData.allPeakX, figData.allPeakY, figData.allPeakZ, ...
                10, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.3);
    end

    % Plot each plate with a different bright color
    colors = lines(length(plates));
    for i = 1:length(plates)
        if ~isempty(plates{i}) && isfield(plates{i}, 'points')
            platePoints = plates{i}.points;
            if ~isempty(platePoints)
                % Plot plate points with bright colors and larger markers
                scatter3(figData.mainAx, platePoints(:, 1), platePoints(:, 2), platePoints(:, 3) * 1e6, ...
                        50, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'k');
            end
        end
    end

    % Restore axis labels and formatting
    xlabel(figData.mainAx, 'X Position (mm)');
    ylabel(figData.mainAx, 'Y Position (mm)');
    zlabel(figData.mainAx, 'Aligned Time (μs)');
    title(figData.mainAx, sprintf('Live Plate Generation - %d Plates Found', length(plates)));
    grid(figData.mainAx, 'on');

    hold(figData.mainAx, 'off');
    drawnow;

catch ME
    % Silently handle visualization errors (figure may have been closed)
    % Don't print error to avoid cluttering output during normal operation
end

end

% Include the original helper functions
function validPlates = filterPlatesBySize(plates, minPointsPerPlate)
% Filter plates that meet minimum point requirements

validPlates = {};
validCount = 0;

for i = 1:length(plates)
    if plates{i}.numPoints >= minPointsPerPlate
        validCount = validCount + 1;
        validPlates{validCount} = plates{i};
        validPlates{validCount}.id = validCount; % Renumber
    end
end

fprintf('Filtered plates: %d valid plates (minimum %d points each)\n', ...
        validCount, minPointsPerPlate);

end

function stats = generatePlateStatistics(plates, allPoints)
% Generate statistics about the plate generation process

stats = struct();

if isempty(plates)
    stats.avgPointsPerPlate = 0;
    stats.maxPointsPerPlate = 0;
    stats.minPointsPerPlate = 0;
    stats.totalPointsInPlates = 0;
    stats.pointUtilization = 0;
    return;
end

% Calculate statistics
pointCounts = cellfun(@(p) p.numPoints, plates);
stats.avgPointsPerPlate = mean(pointCounts);
stats.maxPointsPerPlate = max(pointCounts);
stats.minPointsPerPlate = min(pointCounts);
stats.totalPointsInPlates = sum(pointCounts);
stats.pointUtilization = stats.totalPointsInPlates / size(allPoints, 1) * 100;

% Type distribution
peakPlates = sum(cellfun(@(p) p.plateType == 1, plates));
valleyPlates = sum(cellfun(@(p) p.plateType == -1, plates));
stats.peakPlates = peakPlates;
stats.valleyPlates = valleyPlates;

end
