function plateData = PlateGenerationProcessorWithSpatial(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions)
% PlateGenerationProcessorWithSpatial - Amplitude + Time tolerance with spatial connectivity
%
% This function combines the simple amplitude + time tolerance approach with
% spatial connectivity enforcement to ensure plates are spatially continuous.
%
% ALGORITHM:
% 1. Use amplitude + time tolerance to find candidate points (like simple method)
% 2. Apply spatial connectivity check to ensure plates are continuous
% 3. Split disconnected regions into separate plates
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
%                     .minPointsPerPlate  - Minimum points required to form a plate
%                     .plateType          - 'peaks', 'valleys', or 'both'
%
% OUTPUTS:
%   plateData       - Structure containing generated plate information

fprintf('\n=== Amplitude + Time + Spatial Plate Generation ===\n');

% Start timing
tic;

% Initialize output structure
plateData = struct();
plateData.plates = [];
plateData.numPlates = 0;
plateData.processingTime = 0;
plateData.statistics = struct();

fprintf('Starting plate generation with spatial connectivity:\n');
fprintf('  Amplitude tolerance: ±%.3f\n', plateOptions.amplitudeTolerance);
fprintf('  Time tolerance: ±%.1f μs (±%.2e sec)\n', plateOptions.timeTolerance, plateOptions.timeTolerance * 1e-6);
fprintf('  Minimum points per plate: %d\n', plateOptions.minPointsPerPlate);
fprintf('  Plate type: %s\n', plateOptions.plateType);
fprintf('  Spatial connectivity: ENABLED\n');

% Extract all points from peak data
fprintf('Extracting points from %d waveforms...\n', length(peakData));
allPoints = extractAllPoints(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions.plateType);

if isempty(allPoints)
    fprintf('Warning: No points extracted for plate generation.\n');
    plateData.processingTime = toc;
    return;
end

fprintf('Extracted %d total points for plate generation.\n', size(allPoints, 1));

% Group points using amplitude + time tolerance with spatial connectivity
fprintf('Grouping points with spatial connectivity...\n');
plates = groupPointsWithSpatialConnectivity(allPoints, plateOptions);

% Filter plates by minimum point requirement
validPlates = filterPlatesBySize(plates, plateOptions.minPointsPerPlate);

% Store results
plateData.plates = validPlates;
plateData.numPlates = length(validPlates);
plateData.processingTime = toc;

% Generate statistics
plateData.statistics = generatePlateStatistics(validPlates, allPoints);

% Display results
fprintf('\n=== Amplitude + Time + Spatial Results ===\n');
fprintf('Total plates generated: %d\n', plateData.numPlates);
fprintf('Processing time: %.2f seconds\n', plateData.processingTime);
if plateData.numPlates > 0
    fprintf('Average points per plate: %.1f\n', plateData.statistics.avgPointsPerPlate);
    fprintf('Largest plate: %d points\n', plateData.statistics.maxPointsPerPlate);
    fprintf('Smallest plate: %d points\n', plateData.statistics.minPointsPerPlate);
    fprintf('Point utilization: %.1f%%\n', plateData.statistics.pointUtilization);
    fprintf('✓ All plates are spatially continuous\n');
end

end

function allPoints = extractAllPoints(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateType)
% Extract all relevant points from peak data based on plate type
% Same as original PlateGenerationProcessor

allPoints = [];
pointCount = 0;

% Pre-allocate for efficiency
maxPoints = length(peakData) * 10; % Estimate
tempPoints = zeros(maxPoints, 6); % [X, Y, Time, Amplitude, Type, WaveformIndex]

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

function plates = groupPointsWithSpatialConnectivity(allPoints, plateOptions)
% Group points using amplitude + time tolerance, then enforce spatial connectivity

% Pre-allocate for better performance
maxPlates = 2000;
plates = cell(maxPlates, 1);
plateCount = 0;
usedPoints = false(size(allPoints, 1), 1);

% Convert time tolerance to seconds
timeToleranceSeconds = plateOptions.timeTolerance * 1e-6;
spatialTolerance = 2.0; % mm - maximum distance for spatial connectivity

fprintf('Grouping with tolerances: Amp=±%.3f, Time=±%.2e sec, Spatial=%.1f mm\n', ...
        plateOptions.amplitudeTolerance, timeToleranceSeconds, spatialTolerance);

% Process each unused point as a potential plate seed
totalPoints = size(allPoints, 1);
fprintf('Processing %d points...\n', totalPoints);

for i = 1:totalPoints
    if usedPoints(i)
        continue;
    end
    
    % Current point as seed
    seedPoint = allPoints(i, :);
    seedAmp = seedPoint(4);
    seedTime = seedPoint(3);
    seedType = seedPoint(5);
    
    % Step 1: Find all points within amplitude + time tolerance (like simple method)
    ampDiff = abs(allPoints(:, 4) - seedAmp);
    timeDiff = abs(allPoints(:, 3) - seedTime);
    typeSame = allPoints(:, 5) == seedType;
    
    candidatePoints = ~usedPoints & ...
                     ampDiff <= plateOptions.amplitudeTolerance & ...
                     timeDiff <= timeToleranceSeconds & ...
                     typeSame;
    
    if sum(candidatePoints) >= plateOptions.minPointsPerPlate
        % Step 2: Apply spatial connectivity to split into connected components
        candidateIndices = find(candidatePoints);
        candidateData = allPoints(candidateIndices, :);
        
        % Find spatially connected components
        connectedComponents = findConnectedComponents(candidateData, spatialTolerance);
        
        % Create plates for each connected component
        for compIdx = 1:length(connectedComponents)
            component = connectedComponents{compIdx};
            if length(component) >= plateOptions.minPointsPerPlate
                plateCount = plateCount + 1;
                
                % Expand if needed
                if plateCount > maxPlates
                    plates = [plates; cell(maxPlates, 1)];
                    maxPlates = maxPlates * 2;
                end
                
                % Get actual point indices
                actualIndices = candidateIndices(component);
                platePoints = allPoints(actualIndices, :);
                
                % Create plate structure
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
                plates{plateCount}.plateType = seedType;
                plates{plateCount}.isSpatiallyConnected = true;
                
                % Mark points as used
                usedPoints(actualIndices) = true;
            end
        end
    end
    
    % Progress reporting
    if mod(i, 5000) == 0
        fprintf('Processed %d/%d points, found %d plates...\n', i, totalPoints, plateCount);
    end
end

% Trim to actual size
plates = plates(1:plateCount);

fprintf('✓ Completed grouping. Generated %d spatially connected plates.\n', plateCount);

end

function connectedComponents = findConnectedComponents(points, spatialTolerance)
% Find spatially connected components using simple region growing
% This is much simpler than the complex topological method

numPoints = size(points, 1);
visited = false(numPoints, 1);
connectedComponents = {};
componentCount = 0;

for i = 1:numPoints
    if visited(i)
        continue;
    end

    % Start new connected component
    componentCount = componentCount + 1;
    component = [];
    toVisit = [i];

    while ~isempty(toVisit)
        currentIdx = toVisit(1);
        toVisit(1) = [];

        if visited(currentIdx)
            continue;
        end

        visited(currentIdx) = true;
        component = [component, currentIdx];

        % Find neighbors within spatial tolerance
        currentPoint = points(currentIdx, :);
        currentX = currentPoint(1);
        currentY = currentPoint(2);

        for j = 1:numPoints
            if visited(j)
                continue;
            end

            neighborPoint = points(j, :);
            neighborX = neighborPoint(1);
            neighborY = neighborPoint(2);

            % Check spatial distance
            spatialDist = sqrt((neighborX - currentX)^2 + (neighborY - currentY)^2);

            if spatialDist <= spatialTolerance
                toVisit = [toVisit, j];
            end
        end
    end

    connectedComponents{componentCount} = component;
end

end

function validPlates = filterPlatesBySize(plates, minPointsPerPlate)
% Filter plates that meet minimum point requirements

maxPlates = length(plates);
validPlates = cell(maxPlates, 1);
validCount = 0;

for i = 1:length(plates)
    if plates{i}.numPoints >= minPointsPerPlate
        validCount = validCount + 1;
        validPlates{validCount} = plates{i};
        validPlates{validCount}.id = validCount; % Renumber
    end
end

% Trim to actual size
validPlates = validPlates(1:validCount);

fprintf('  ✓ Filtered plates: %d valid plates (minimum %d points each)\n', ...
        validCount, minPointsPerPlate);
end

function stats = generatePlateStatistics(plates, allPoints)
% Generate statistics for the generated plates

stats = struct();

if isempty(plates)
    stats.avgPointsPerPlate = 0;
    stats.maxPointsPerPlate = 0;
    stats.minPointsPerPlate = 0;
    stats.totalPointsInPlates = 0;
    stats.pointUtilization = 0;
    return;
end

% Basic statistics
pointCounts = cellfun(@(p) p.numPoints, plates);
stats.avgPointsPerPlate = mean(pointCounts);
stats.maxPointsPerPlate = max(pointCounts);
stats.minPointsPerPlate = min(pointCounts);
stats.totalPointsInPlates = sum(pointCounts);
stats.pointUtilization = stats.totalPointsInPlates / size(allPoints, 1) * 100;

end
