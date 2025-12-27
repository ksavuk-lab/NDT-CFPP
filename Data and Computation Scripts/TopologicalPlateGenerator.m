function plateData = TopologicalPlateGenerator(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions, varargin)
% TOPOLOGICALPLATEGENERATOR - Generate non-intersecting layered plates
%
% This function creates plates that follow proper topological rules:
% 1. Groups cannot intersect each other
% 2. Groups can be layered like an onion (depth-based layers)
% 3. Groups must be spatially continuous
% 4. No "tunneling" through other groups without proper connections
%
% INPUTS:
%   peakData        - Cell array of peak/valley data from PeakExtractionProcessor
%   X_Coordinates   - X spatial coordinates array
%   Y_Coordinates   - Y spatial coordinates array  
%   numY_sub        - Number of Y subdivisions
%   numX_sub        - Number of X subdivisions
%   plateOptions    - Structure with plate generation settings

tic;
fprintf('\n=== Topological Plate Generation ===\n');
fprintf('Generating non-intersecting layered plates...\n');

% Initialize output structure
plateData = struct();
plateData.plates = {};
plateData.numPlates = 0;
plateData.processingTime = 0;

fprintf('Starting topological plate generation with settings:\n');
fprintf('  Amplitude tolerance: ±%.3f\n', plateOptions.amplitudeTolerance);
fprintf('  Time tolerance: ±%.1f μs (±%.2e sec)\n', plateOptions.timeTolerance, plateOptions.timeTolerance * 1e-6);
fprintf('  Minimum points per plate: %d\n', plateOptions.minPointsPerPlate);
fprintf('  Plate type: %s\n', plateOptions.plateType);
fprintf('  Spatial connectivity: ENABLED\n');
fprintf('  Layer topology: NON-INTERSECTING\n');

% Check if aligned coordinate data is provided
if length(varargin) >= 5
    % Use aligned coordinate data directly
    allPeakX = varargin{1};
    allPeakY = varargin{2};
    allPeakZ = varargin{3};
    allPeakAmps = varargin{4};
    allPeakTypes = varargin{5};

    fprintf('Using aligned coordinate data (%d points)...\n', length(allPeakX));
    allPoints = createAllPointsFromAlignedTopo(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, plateOptions.plateType);
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

% Step 1: Create depth-based layers (onion layers)
fprintf('Step 1: Creating depth-based layers...\n');
step1_start = tic;
depthLayers = createDepthLayers(allPoints, plateOptions);
step1_time = toc(step1_start);
fprintf('  → Step 1 completed in %.2f seconds\n', step1_time);

% Step 2: Within each layer, create spatially continuous groups
fprintf('Step 2: Creating spatially continuous groups within layers...\n');
step2_start = tic;
layeredPlates = createSpatialGroupsWithinLayers(depthLayers, plateOptions);
step2_time = toc(step2_start);
fprintf('  → Step 2 completed in %.2f seconds\n', step2_time);

% Step 3: Validate non-intersection constraint
fprintf('Step 3: Validating non-intersection constraints...\n');
step3_start = tic;
validatedPlates = validateNonIntersection(layeredPlates, plateOptions);
step3_time = toc(step3_start);
fprintf('  → Step 3 completed in %.2f seconds\n', step3_time);

% Step 4: Filter by minimum size and finalize
validPlates = filterPlatesBySize(validatedPlates, plateOptions.minPointsPerPlate);

% Store results
plateData.plates = validPlates;
plateData.numPlates = length(validPlates);
plateData.processingTime = toc;

% Generate statistics
plateData.statistics = generatePlateStatistics(validPlates, allPoints);

% Display results
fprintf('\n=== Topological Plate Generation Results ===\n');
fprintf('Total plates generated: %d\n', plateData.numPlates);
fprintf('Processing time: %.2f seconds\n', plateData.processingTime);
if plateData.numPlates > 0
    fprintf('Average points per plate: %.1f\n', plateData.statistics.avgPointsPerPlate);
    fprintf('Largest plate: %d points\n', plateData.statistics.maxPointsPerPlate);
    fprintf('Smallest plate: %d points\n', plateData.statistics.minPointsPerPlate);
    fprintf('Depth layers: %d\n', plateData.statistics.numDepthLayers);
    fprintf('✓ All plates are non-intersecting and topologically valid\n');
end

end

function allPoints = createAllPointsFromAlignedTopo(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, plateType)
% Create allPoints array from aligned coordinate data for topological processing
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
% Same as original but with enhanced spatial indexing

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

function depthLayers = createDepthLayers(allPoints, plateOptions)
% Create depth-based layers (onion layers) based on time/depth
% Points are grouped into layers that don't intersect in the depth dimension

fprintf('  Creating depth layers based on time values...\n');

% Sort points by time (depth)
[sortedTimes, sortIdx] = sort(allPoints(:, 3)); % Column 3 is time
sortedPoints = allPoints(sortIdx, :);

% Create layers with time tolerance
timeToleranceSeconds = plateOptions.timeTolerance * 1e-6;

% Pre-allocate cell array for better performance
maxLayers = 200; % Reasonable estimate
depthLayers = cell(maxLayers, 1);
layerCount = 0;
usedPoints = false(size(sortedPoints, 1), 1);

for i = 1:size(sortedPoints, 1)
    if usedPoints(i)
        continue;
    end

    currentTime = sortedPoints(i, 3);
    currentType = sortedPoints(i, 5);

    % Find all points within time tolerance and same type (vectorized)
    timeDiff = abs(sortedPoints(:, 3) - currentTime);
    typeSame = sortedPoints(:, 5) == currentType;
    layerPoints = ~usedPoints & timeDiff <= timeToleranceSeconds & typeSame;

    if sum(layerPoints) >= plateOptions.minPointsPerPlate
        layerCount = layerCount + 1;
        depthLayers{layerCount} = sortedPoints(layerPoints, :);
        usedPoints(layerPoints) = true;

        fprintf('    Layer %d: %.2e sec, %d points (type %d)\n', ...
                layerCount, currentTime, sum(layerPoints), currentType);
    end
end

% Trim to actual size
depthLayers = depthLayers(1:layerCount);

fprintf('  ✓ Created %d depth layers\n', layerCount);
end

function layeredPlates = createSpatialGroupsWithinLayers(depthLayers, plateOptions)
% Within each depth layer, create spatially continuous groups
% This ensures groups don't "tunnel" through each other

fprintf('  Creating spatially continuous groups within each layer...\n');

% Pre-allocate for better performance
maxPlates = 2000; % Reasonable estimate
layeredPlates = cell(maxPlates, 1);
plateCount = 0;

for layerIdx = 1:length(depthLayers)
    layerPoints = depthLayers{layerIdx};

    if size(layerPoints, 1) < plateOptions.minPointsPerPlate
        continue;
    end

    fprintf('    Processing layer %d (%d points)...\n', layerIdx, size(layerPoints, 1));

    % Create spatially connected groups within this layer
    layer_start = tic;
    spatialGroups = createSpatiallyConnectedGroupsOptimized(layerPoints, plateOptions);
    layer_time = toc(layer_start);
    fprintf('      → Layer %d processed in %.2f seconds\n', layerIdx, layer_time);

    % Add layer information to each group
    for groupIdx = 1:length(spatialGroups)
        plateCount = plateCount + 1;
        if plateCount > maxPlates
            % Expand if needed
            layeredPlates = [layeredPlates; cell(maxPlates, 1)];
            maxPlates = maxPlates * 2;
        end
        layeredPlates{plateCount} = spatialGroups{groupIdx};
        layeredPlates{plateCount}.layerIndex = layerIdx;
        layeredPlates{plateCount}.layerDepth = mean(layerPoints(:, 3)); % Average time for layer
    end

    fprintf('      → Created %d spatial groups in layer %d\n', length(spatialGroups), layerIdx);
end

% Trim to actual size
layeredPlates = layeredPlates(1:plateCount);

fprintf('  ✓ Created %d total spatial groups across all layers\n', plateCount);
end

function spatialGroups = createSpatiallyConnectedGroupsOptimized(layerPoints, plateOptions)
% Create spatially connected groups using optimized region growing
% Ensures groups are continuous and don't have spatial gaps

% Pre-allocate for better performance
maxGroups = 500; % Reasonable estimate
spatialGroups = cell(maxGroups, 1);
groupCount = 0;
usedPoints = false(size(layerPoints, 1), 1);

% Spatial connectivity parameters
spatialTolerance = 2.0; % mm - maximum distance for spatial connectivity
ampTolerance = plateOptions.amplitudeTolerance;

% Create spatial index for faster neighbor lookup
spatialIndex = createSpatialIndex(layerPoints, spatialTolerance);

for i = 1:size(layerPoints, 1)
    if usedPoints(i)
        continue;
    end

    % Start a new spatial group with optimized region growing
    seedPoint = layerPoints(i, :);
    groupPoints = regionGrowSpatialGroupOptimized(layerPoints, i, usedPoints, spatialTolerance, ampTolerance, spatialIndex);

    if length(groupPoints) >= plateOptions.minPointsPerPlate
        groupCount = groupCount + 1;

        % Expand if needed
        if groupCount > maxGroups
            spatialGroups = [spatialGroups; cell(maxGroups, 1)];
            maxGroups = maxGroups * 2;
        end

        % Create group structure
        groupData = layerPoints(groupPoints, :);
        spatialGroups{groupCount} = struct();
        spatialGroups{groupCount}.id = groupCount;
        spatialGroups{groupCount}.points = groupData;
        spatialGroups{groupCount}.numPoints = size(groupData, 1);
        spatialGroups{groupCount}.centerAmplitude = mean(groupData(:, 4));
        spatialGroups{groupCount}.centerTime = mean(groupData(:, 3));
        spatialGroups{groupCount}.amplitudeRange = [min(groupData(:, 4)), max(groupData(:, 4))];
        spatialGroups{groupCount}.timeRange = [min(groupData(:, 3)), max(groupData(:, 3))];
        spatialGroups{groupCount}.spatialBounds = [min(groupData(:, 1)), max(groupData(:, 1)), ...
                                                  min(groupData(:, 2)), max(groupData(:, 2))];
        spatialGroups{groupCount}.plateType = seedPoint(5);
        spatialGroups{groupCount}.isSpatiallyConnected = true;

        % Mark points as used
        usedPoints(groupPoints) = true;
    end
end

% Trim to actual size
spatialGroups = spatialGroups(1:groupCount);
end

function spatialIndex = createSpatialIndex(layerPoints, spatialTolerance)
% Create spatial index for faster neighbor lookup using grid-based hashing

% Calculate grid parameters
xMin = min(layerPoints(:, 1));
yMin = min(layerPoints(:, 2));
xMax = max(layerPoints(:, 1));
yMax = max(layerPoints(:, 2));

gridSize = spatialTolerance; % Use spatial tolerance as grid size
nGridX = ceil((xMax - xMin) / gridSize) + 1;
nGridY = ceil((yMax - yMin) / gridSize) + 1;

% Initialize spatial index
spatialIndex = struct();
spatialIndex.xMin = xMin;
spatialIndex.yMin = yMin;
spatialIndex.gridSize = gridSize;
spatialIndex.nGridX = nGridX;
spatialIndex.nGridY = nGridY;
spatialIndex.grid = cell(nGridX, nGridY);

% Populate grid with point indices
for i = 1:size(layerPoints, 1)
    x = layerPoints(i, 1);
    y = layerPoints(i, 2);

    gridX = floor((x - xMin) / gridSize) + 1;
    gridY = floor((y - yMin) / gridSize) + 1;

    % Clamp to grid bounds
    gridX = max(1, min(nGridX, gridX));
    gridY = max(1, min(nGridY, gridY));

    spatialIndex.grid{gridX, gridY} = [spatialIndex.grid{gridX, gridY}, i];
end
end

function groupPoints = regionGrowSpatialGroupOptimized(layerPoints, seedIdx, usedPoints, spatialTolerance, ampTolerance, spatialIndex)
% Optimized region growing algorithm using spatial indexing and pre-allocated arrays

% Pre-allocate arrays to avoid dynamic growth
maxGroupSize = size(layerPoints, 1);
groupPointsArray = zeros(maxGroupSize, 1);
toProcessArray = zeros(maxGroupSize, 1);

groupCount = 1;
processCount = 1;
processIdx = 1;

groupPointsArray(1) = seedIdx;
toProcessArray(1) = seedIdx;

processed = false(size(layerPoints, 1), 1);
processed(seedIdx) = true;

seedAmp = layerPoints(seedIdx, 4);
seedType = layerPoints(seedIdx, 5);

while processIdx <= processCount
    currentIdx = toProcessArray(processIdx);
    processIdx = processIdx + 1;

    currentPoint = layerPoints(currentIdx, :);
    currentX = currentPoint(1);
    currentY = currentPoint(2);

    % Get candidate points from spatial index (much faster than checking all points)
    candidates = getSpatialNeighbors(currentX, currentY, spatialIndex);

    % Check each candidate
    for j = 1:length(candidates)
        i = candidates(j);

        if processed(i) || usedPoints(i)
            continue;
        end

        candidatePoint = layerPoints(i, :);
        candidateX = candidatePoint(1);
        candidateY = candidatePoint(2);
        candidateAmp = candidatePoint(4);
        candidateType = candidatePoint(5);

        % Check spatial distance (only for nearby candidates)
        spatialDist = sqrt((candidateX - currentX)^2 + (candidateY - currentY)^2);

        % Check amplitude similarity and type match
        ampDiff = abs(candidateAmp - seedAmp);

        if spatialDist <= spatialTolerance && ampDiff <= ampTolerance && candidateType == seedType
            % Add to group
            groupCount = groupCount + 1;
            processCount = processCount + 1;
            groupPointsArray(groupCount) = i;
            toProcessArray(processCount) = i;
            processed(i) = true;
        end
    end
end

% Return trimmed array
groupPoints = groupPointsArray(1:groupCount);
end

function neighbors = getSpatialNeighbors(x, y, spatialIndex)
% Get all points in neighboring grid cells

neighbors = [];

% Calculate grid position
gridX = floor((x - spatialIndex.xMin) / spatialIndex.gridSize) + 1;
gridY = floor((y - spatialIndex.yMin) / spatialIndex.gridSize) + 1;

% Check 3x3 neighborhood of grid cells
for dx = -1:1
    for dy = -1:1
        checkX = gridX + dx;
        checkY = gridY + dy;

        % Check bounds
        if checkX >= 1 && checkX <= spatialIndex.nGridX && ...
           checkY >= 1 && checkY <= spatialIndex.nGridY
            cellPoints = spatialIndex.grid{checkX, checkY};
            neighbors = [neighbors, cellPoints];
        end
    end
end
end

function validatedPlates = validateNonIntersection(layeredPlates, plateOptions)
% Validate that plates don't intersect each other
% Remove or modify plates that violate non-intersection constraint

fprintf('  Validating non-intersection constraints...\n');

% Pre-allocate for better performance
maxPlates = length(layeredPlates);
validatedPlates = cell(maxPlates, 1);
validCount = 0;

% Sort plates by layer depth (deepest first)
if ~isempty(layeredPlates)
    layerDepths = cellfun(@(p) p.layerDepth, layeredPlates);
    [~, sortIdx] = sort(layerDepths);
    sortedPlates = layeredPlates(sortIdx);
else
    sortedPlates = {};
end

% Check each plate against previously validated plates
for i = 1:length(sortedPlates)
    currentPlate = sortedPlates{i};
    isValid = true;

    % Check against all previously validated plates
    for j = 1:validCount
        if platesIntersect(currentPlate, validatedPlates{j})
            fprintf('    Plate %d intersects with plate %d - removing\n', i, j);
            isValid = false;
            break;
        end
    end

    if isValid
        validCount = validCount + 1;
        validatedPlates{validCount} = currentPlate;
        validatedPlates{validCount}.id = validCount; % Renumber
    end
end

% Trim to actual size
validatedPlates = validatedPlates(1:validCount);

fprintf('  ✓ Validated %d non-intersecting plates (removed %d intersecting)\n', ...
        validCount, length(layeredPlates) - validCount);
end

function intersects = platesIntersect(plate1, plate2)
% Check if two plates intersect in 3D space
% Uses spatial bounds and depth overlap to determine intersection

% Get spatial bounds [minX, maxX, minY, maxY]
bounds1 = plate1.spatialBounds;
bounds2 = plate2.spatialBounds;

% Check XY spatial overlap
xyOverlap = ~(bounds1(2) < bounds2(1) || bounds2(2) < bounds1(1) || ...
              bounds1(4) < bounds2(3) || bounds2(4) < bounds1(3));

if ~xyOverlap
    intersects = false;
    return;
end

% Check depth (time) overlap
timeRange1 = plate1.timeRange;
timeRange2 = plate2.timeRange;

timeOverlap = ~(timeRange1(2) < timeRange2(1) || timeRange2(2) < timeRange1(1));

% Plates intersect if they overlap in both XY and time
intersects = xyOverlap && timeOverlap;
end

function validPlates = filterPlatesBySize(plates, minPointsPerPlate)
% Filter plates that meet minimum point requirements (optimized)

% Pre-allocate for better performance
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
% Generate statistics about the plate generation process (enhanced)

stats = struct();

if isempty(plates)
    stats.avgPointsPerPlate = 0;
    stats.maxPointsPerPlate = 0;
    stats.minPointsPerPlate = 0;
    stats.totalPointsInPlates = 0;
    stats.pointUtilization = 0;
    stats.numDepthLayers = 0;
    return;
end

% Basic statistics
pointCounts = cellfun(@(p) p.numPoints, plates);
stats.avgPointsPerPlate = mean(pointCounts);
stats.maxPointsPerPlate = max(pointCounts);
stats.minPointsPerPlate = min(pointCounts);
stats.totalPointsInPlates = sum(pointCounts);
stats.pointUtilization = stats.totalPointsInPlates / size(allPoints, 1) * 100;

% Layer statistics
if isfield(plates{1}, 'layerIndex')
    layerIndices = cellfun(@(p) p.layerIndex, plates);
    stats.numDepthLayers = length(unique(layerIndices));
else
    stats.numDepthLayers = 1;
end
end
