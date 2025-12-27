function plateData = PlateGenerationProcessor(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions, varargin)
% PlateGenerationProcessor - Support script for Main.m to generate flat plates from peak/valley data
%
% This function groups similar peak/valley points into flat plate structures based on
% amplitude and time tolerances, creating a layered representation of NDT data.
%
% INPUTS:
%   peakData        - Cell array of peak/valley data from PeakExtractionProcessor
%   X_Coordinates   - X spatial coordinates array
%   Y_Coordinates   - Y spatial coordinates array
%   numY_sub        - Number of Y subdivisions
%   numX_sub        - Number of X subdivisions
%   plateOptions    - Structure with plate generation settings:
%                     .amplitudeTolerance - Amplitude grouping tolerance (e.g., 0.1)
%                     .timeTolerance      - Time grouping tolerance in units of 1e-6 sec (e.g., 1 = 1e-6 sec)
%                     .enablePlateGeneration - Boolean to enable/disable plate generation
%                     .minPointsPerPlate  - Minimum points required to form a plate
%                     .plateType          - 'peaks', 'valleys', or 'both'
%   varargin        - Optional aligned coordinate data:
%                     varargin{1} = allPeakX (aligned X coordinates)
%                     varargin{2} = allPeakY (aligned Y coordinates)
%                     varargin{3} = allPeakZ (aligned Z/time coordinates)
%                     varargin{4} = allPeakAmps (aligned amplitudes)
%                     varargin{5} = allPeakTypes (peak/valley types)
%
% OUTPUTS:
%   plateData       - Structure containing generated plate information:
%                     .plates            - Array of plate structures
%                     .numPlates         - Total number of plates generated
%                     .processingTime    - Time taken for plate generation
%                     .statistics        - Processing statistics

fprintf('\n=== Plate Generation Processor ===\n');

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
    fprintf('Plate generation is disabled. Skipping processing.\n');
    plateData.processingTime = toc;
    return;
end

% Validate inputs
if isempty(peakData)
    fprintf('Error: No peak data provided for plate generation.\n');
    plateData.processingTime = toc;
    return;
end

fprintf('Starting plate generation with settings:\n');
fprintf('  Amplitude tolerance: ±%.3f\n', plateOptions.amplitudeTolerance);
fprintf('  Time tolerance: ±%.1f μs (±%.2e sec)\n', plateOptions.timeTolerance, plateOptions.timeTolerance * 1e-6);
fprintf('  Minimum points per plate: %d\n', plateOptions.minPointsPerPlate);
fprintf('  Plate type: %s\n', plateOptions.plateType);

% Check if aligned coordinate data is provided
if length(varargin) >= 5
    % Use aligned coordinate data directly
    allPeakX = varargin{1};
    allPeakY = varargin{2};
    allPeakZ = varargin{3};
    allPeakAmps = varargin{4};
    allPeakTypes = varargin{5};

    fprintf('Using aligned coordinate data (%d points)...\n', length(allPeakX));
    allPoints = createAllPointsFromAligned(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, plateOptions.plateType);
else
    % Extract all points from peak data (original method)
    fprintf('Extracting points from %d waveforms...\n', length(peakData));
    allPoints = extractAllPoints(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions.plateType);
end

if isempty(allPoints)
    fprintf('Warning: No points extracted for plate generation.\n');
    plateData.processingTime = toc;
    return;
end

fprintf('Extracted %d total points for plate generation.\n', size(allPoints, 1));

% Group points into plates based on amplitude and time tolerances
fprintf('Grouping points into plates...\n');

% DISABLED: Progress dialog causing issues - use console output only
progressHandle = [];
fprintf('Using console-only progress reporting (popup dialog disabled)\n');

plates = groupPointsIntoPlatesWithProgress(allPoints, plateOptions, progressHandle, X_Coordinates, Y_Coordinates);

% Close progress dialog
if ~isempty(progressHandle)
    try
        closeProgress(progressHandle);
    catch ME
        fprintf('Warning: Failed to close progress dialog: %s\n', ME.message);
    end
end

% Filter plates by minimum point requirement
validPlates = filterPlatesBySize(plates, plateOptions.minPointsPerPlate);

% Store results
plateData.plates = validPlates;
plateData.numPlates = length(validPlates);
plateData.processingTime = toc;

% Generate statistics
plateData.statistics = generatePlateStatistics(validPlates, allPoints);

% Display results
fprintf('\n=== Plate Generation Results ===\n');
fprintf('Total plates generated: %d\n', plateData.numPlates);
fprintf('Processing time: %.2f seconds\n', plateData.processingTime);
fprintf('Average points per plate: %.1f\n', plateData.statistics.avgPointsPerPlate);
fprintf('Largest plate: %d points\n', plateData.statistics.maxPointsPerPlate);
fprintf('Smallest plate: %d points\n', plateData.statistics.minPointsPerPlate);

end

function allPoints = createAllPointsFromAligned(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, plateType)
% Create allPoints array from aligned coordinate data
% This ensures plate generation uses the same aligned data as visualization

fprintf('Creating points array from aligned coordinate data...\n');

% Convert Z from microseconds back to seconds for consistency with original data
allPeakZ_seconds = allPeakZ * 1e-6;

% Filter by plate type
switch plateType
    case 'peaks'
        validIndices = allPeakTypes == 1;
    case 'valleys'
        validIndices = allPeakTypes == -1;
    case 'both'
        validIndices = true(size(allPeakTypes));
    otherwise
        validIndices = true(size(allPeakTypes));
end

% Create allPoints array: [X, Y, Time, Amplitude, Type, WaveformIndex]
% Note: WaveformIndex is set to 0 since we don't have that info from aligned data
allPoints = [allPeakX(validIndices), ...
             allPeakY(validIndices), ...
             allPeakZ_seconds(validIndices), ...
             allPeakAmps(validIndices), ...
             allPeakTypes(validIndices), ...
             zeros(sum(validIndices), 1)]; % WaveformIndex set to 0

fprintf('✓ Created %d points from aligned data (filtered for %s)\n', size(allPoints, 1), plateType);
end

function allPoints = extractAllPoints(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateType)
% Extract all relevant points from peak data based on plate type

allPoints = [];
pointCount = 0;

% Estimate total points for preallocation
estimatedPoints = length(peakData) * 20; % Rough estimate
tempPoints = zeros(estimatedPoints, 6); % [X, Y, Time, Amplitude, Type, WaveformIndex]

for waveformIdx = 1:length(peakData)
    if isempty(peakData{waveformIdx})
        continue;
    end
    
    % Calculate spatial coordinates
    [yIdx, xIdx] = ind2sub([numY_sub, numX_sub], waveformIdx);
    if xIdx > length(X_Coordinates) || yIdx > length(Y_Coordinates)
        continue;
    end
    
    spatialX = X_Coordinates(xIdx);
    spatialY = Y_Coordinates(yIdx);
    
    % Get peak/valley data for this waveform
    waveformData = peakData{waveformIdx};
    
    % Filter by plate type
    switch plateType
        case 'peaks'
            validIndices = waveformData.TransitionType == 1;
        case 'valleys'
            validIndices = waveformData.TransitionType == -1;
        case 'both'
            validIndices = true(height(waveformData), 1);
        otherwise
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

function plates = groupPointsIntoPlatesWithProgress(allPoints, plateOptions, progressHandle, X_Coordinates, Y_Coordinates)
% Group points into plates based on amplitude and time tolerances with STRICT spatial constraints
% ALGORITHM #1: Amplitude + Time + Strict Grid Connectivity with Progress Tracking
% - Points must be within amplitude ± tolerance AND time ± tolerance
% - NO INTERSECTING GROUPS (critical constraint)
% - Links only allowed within 1 step size (0.1mm) in any direction (3x3x3 cube)
% - Each point can only connect to immediate neighbors in grid

plates = {};
usedPoints = false(size(allPoints, 1), 1);
plateCount = 0;

% Convert time tolerance to seconds
timeToleranceSeconds = plateOptions.timeTolerance * 1e-6;

% CORRECTED: Calculate actual step size from coordinate data
if length(X_Coordinates) > 1
    stepSizeX = abs(X_Coordinates(2) - X_Coordinates(1));
else
    stepSizeX = 0.1; % Fallback
end

if length(Y_Coordinates) > 1
    stepSizeY = abs(Y_Coordinates(2) - Y_Coordinates(1));
else
    stepSizeY = 0.1; % Fallback
end

% Use the actual step size (should be consistent for X and Y)
stepSize = max(stepSizeX, stepSizeY); % Use the larger step size for safety

fprintf('Detected step size from data: X=%.3f mm, Y=%.3f mm, Using=%.3f mm\n', ...
    stepSizeX, stepSizeY, stepSize);

totalPoints = size(allPoints, 1);

fprintf('Algorithm #1: Amplitude + Time + Strict Grid Connectivity\n');
fprintf('Grouping points with tolerances: Amp=±%.3f, Time=±%.2e sec, Step=%.1f mm (3x3x3 cube)\n', ...
        plateOptions.amplitudeTolerance, timeToleranceSeconds, stepSize);
fprintf('CONSTRAINT: No intersecting groups allowed - strict grid connectivity only\n');

% Process each unused point as a potential plate seed
progressInterval = max(100, floor(totalPoints / 500)); % Show progress ~500 times (much more frequent)
fprintf('Processing %d points to find grid-connected plates...\n', totalPoints);

% Add safety timeout to prevent infinite hanging
startTime = tic;
maxProcessingTime = 1800; % 30 minutes maximum
timeoutWarningShown = false;
lastTimeUpdate = tic; % For time-based progress updates

for i = 1:totalPoints
    % Check for timeout every 1000 iterations
    if mod(i, 1000) == 0
        elapsedTime = toc(startTime);
        if elapsedTime > maxProcessingTime
            fprintf('⚠️  Processing timeout reached (%.1f minutes). Stopping early.\n', elapsedTime/60);
            fprintf('   Processed %d/%d points, found %d plates\n', i, totalPoints, plateCount);
            break;
        elseif elapsedTime > maxProcessingTime * 0.8 && ~timeoutWarningShown
            fprintf('⚠️  Processing is taking a long time (%.1f minutes elapsed)\n', elapsedTime/60);
            fprintf('   Consider using a smaller dataset or different algorithm\n');
            timeoutWarningShown = true;
        end
    end

    % Show progress every progressInterval points OR every 10 seconds
    timeBasedUpdate = (toc(lastTimeUpdate) > 10); % Every 10 seconds
    pointBasedUpdate = (mod(i, progressInterval) == 0);

    if pointBasedUpdate || timeBasedUpdate
        percentComplete = (i / totalPoints) * 100;
        elapsedTime = toc(startTime);
        statusMessage = sprintf('Processing point %d/%d (%d plates found)', i, totalPoints, plateCount);
        currentFunction = 'groupPointsIntoPlatesWithProgress';

        % Progress dialog disabled - using console output only

        fprintf('Progress: %.1f%% (%d/%d points processed, %d plates found, %.1f min elapsed)\n', ...
                percentComplete, i, totalPoints, plateCount, elapsedTime/60);

        if timeBasedUpdate
            lastTimeUpdate = tic; % Reset time-based update timer
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

    % CORRECTED: Only check immediate grid neighbors (not entire dataset)
    % Start region growing from seed point using immediate neighbors only
    plateIndices = regionGrowImmediateNeighbors(allPoints, i, usedPoints, plateOptions, stepSize);

        if length(plateIndices) >= plateOptions.minPointsPerPlate
            plateCount = plateCount + 1;

            % Create plate structure
            platePoints = allPoints(plateIndices, :);
        
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
            plates{plateCount}.isGridConnected = true; % Algorithm #1 ensures strict grid connectivity
            plates{plateCount}.stepSize = stepSize; % Record the step size used

            % Mark these points as used
            usedPoints(plateIndices) = true;

            % Show plate generation progress more frequently
            if mod(plateCount, 5) == 0  % Every 5 plates instead of 10
                fprintf('Generated %d plates... (%.1f%% points processed)\n', plateCount, (i/totalPoints)*100);
            end
        end
end

fprintf('✓ Completed grouping. Generated %d candidate plates.\n', plateCount);
if plateCount > 0
    fprintf('  → Plates will be filtered by minimum size requirement (%d points)\n', plateOptions.minPointsPerPlate);
end

% Final progress update
if ~isempty(progressHandle)
    try
        updateProgress(progressHandle, totalPoints, sprintf('Complete: Generated %d plates with strict grid connectivity', plateCount), 'groupPointsIntoPlatesWithProgress');
    catch ME
        fprintf('Warning: Final progress update failed: %s\n', ME.message);
    end
end

end % End of groupPointsIntoPlatesWithProgress function

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

fprintf('✓ Filtered plates: %d valid plates (minimum %d points each)\n', ...
        validCount, minPointsPerPlate);
if validCount > 0
    fprintf('  → Plates are ready for visualization\n');
else
    fprintf('  → No plates met the minimum size requirement\n');
end

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

function connectedIndices = findGridConnectedGroup(allPoints, candidateIndices, seedIdx, stepSize)
% Find grid-connected group starting from seed point using STRICT 3x3x3 cube connectivity
% CONSTRAINTS:
% - Only points within 1 step size (0.1mm) in ANY direction can connect
% - Uses 3x3x3 cube neighborhood (max distance = stepSize in X, Y, Z)
% - Prevents intersecting groups by enforcing strict grid connectivity
% - No diagonal connections beyond step size

numCandidates = length(candidateIndices);
visited = false(numCandidates, 1);
connectedIndices = [];

% Find seed in candidate list
seedInCandidates = find(candidateIndices == seedIdx);
if isempty(seedInCandidates)
    return;
end

% Start region growing from seed with STRICT grid constraints
toVisit = [seedInCandidates];
visited(seedInCandidates) = true;

while ~isempty(toVisit)
    currentCandidateIdx = toVisit(1);
    toVisit(1) = [];

    currentPointIdx = candidateIndices(currentCandidateIdx);
    connectedIndices = [connectedIndices, currentPointIdx];

    currentPoint = allPoints(currentPointIdx, :);
    currentX = currentPoint(1);
    currentY = currentPoint(2);
    currentZ = currentPoint(3); % Include Z (time) in grid connectivity

    % Find unvisited neighbors within STRICT 3x3x3 cube (step size in each direction)
    for i = 1:numCandidates
        if visited(i)
            continue;
        end

        neighborPointIdx = candidateIndices(i);
        neighborPoint = allPoints(neighborPointIdx, :);
        neighborX = neighborPoint(1);
        neighborY = neighborPoint(2);
        neighborZ = neighborPoint(3);

        % Check STRICT 3x3x3 cube connectivity - must be within step size in ALL directions
        deltaX = abs(neighborX - currentX);
        deltaY = abs(neighborY - currentY);
        deltaZ = abs(neighborZ - currentZ);

        % STRICT constraint: within step size in X, Y, AND Z directions
        if deltaX <= stepSize && deltaY <= stepSize && deltaZ <= (stepSize * 1e-6) % Convert Z to seconds
            visited(i) = true;
            toVisit = [toVisit, i];
        end
    end
end

end

function plateIndices = regionGrowImmediateNeighbors(allPoints, seedIdx, usedPoints, plateOptions, stepSize)
% CORRECTED: Region growing that only checks immediate grid neighbors
% This is the correct implementation of amplitude + time method

plateIndices = [];
toVisit = [seedIdx];
visited = false(size(allPoints, 1), 1);
visited(seedIdx) = true;

seedPoint = allPoints(seedIdx, :);
seedAmp = seedPoint(4);
seedTime = seedPoint(3);
seedType = seedPoint(5);

timeToleranceSeconds = plateOptions.timeTolerance * 1e-6;

while ~isempty(toVisit)
    currentIdx = toVisit(1);
    toVisit(1) = [];

    if usedPoints(currentIdx)
        continue;
    end

    plateIndices = [plateIndices, currentIdx];
    currentPoint = allPoints(currentIdx, :);

    % Find immediate grid neighbors only (not entire dataset)
    immediateNeighbors = findImmediateGridNeighbors(allPoints, currentIdx, stepSize);

    for neighborIdx = immediateNeighbors
        if visited(neighborIdx) || usedPoints(neighborIdx)
            continue;
        end

        neighborPoint = allPoints(neighborIdx, :);

        % Check amplitude, time, and type tolerances
        ampDiff = abs(neighborPoint(4) - seedAmp);
        timeDiff = abs(neighborPoint(3) - seedTime);
        typeSame = neighborPoint(5) == seedType;

        if ampDiff <= plateOptions.amplitudeTolerance && ...
           timeDiff <= timeToleranceSeconds && ...
           typeSame
            visited(neighborIdx) = true;
            toVisit = [toVisit, neighborIdx];
        end
    end
end
end

function immediateNeighbors = findImmediateGridNeighbors(allPoints, centerIdx, stepSize)
% Find only immediate grid neighbors (adjacent points within one step size)

centerPoint = allPoints(centerIdx, :);
centerX = centerPoint(1);
centerY = centerPoint(2);

immediateNeighbors = [];

% Check all points for immediate spatial neighbors
for i = 1:size(allPoints, 1)
    if i == centerIdx
        continue;
    end

    point = allPoints(i, :);
    deltaX = abs(point(1) - centerX);
    deltaY = abs(point(2) - centerY);

    % Only immediate neighbors (within one step size in X and Y)
    if deltaX <= stepSize && deltaY <= stepSize
        immediateNeighbors = [immediateNeighbors, i];
    end
end
end
