function plateData = FastPlateGenerator(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions, varargin)
% FASTPLATEGENERATOR - High-performance plate generation using spatial indexing
% This replaces the slow O(n²) algorithm with an O(n log n) spatial index approach
%
% PERFORMANCE IMPROVEMENTS:
% - Spatial indexing for O(log n) neighbor lookup instead of O(n)
% - Grid-based hashing for instant spatial queries
% - Vectorized operations where possible
% - Early termination and pruning

fprintf('=== FAST PLATE GENERATOR ===\n');
fprintf('Using optimized spatial indexing for high-performance plate generation\n');

tic;

% Check if aligned coordinate data is provided
if length(varargin) >= 5
    % Use aligned coordinate data directly
    allPeakX = varargin{1};
    allPeakY = varargin{2};
    allPeakZ = varargin{3};
    allPeakAmps = varargin{4};
    allPeakTypes = varargin{5};

    fprintf('Using aligned coordinate data (%d points)...\n', length(allPeakX));
    allPoints = createAllPointsFromAlignedFast(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, plateOptions.plateType);
else
    % Extract all points (original method)
    allPoints = extractAllPointsFast(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions.plateType);
end
fprintf('Extracted %d total points for fast plate generation.\n', size(allPoints, 1));

% Create spatial index for O(log n) neighbor lookup
fprintf('Building spatial index...\n');
spatialIndex = buildSpatialIndex(allPoints, plateOptions, X_Coordinates, Y_Coordinates);
fprintf('✓ Spatial index built in %.2f seconds\n', toc);

% Generate plates using fast algorithm
fprintf('Generating plates with fast neighbor lookup...\n');
startTime = tic;
plates = generatePlatesFast(allPoints, spatialIndex, plateOptions);
plateGenTime = toc(startTime);

% Filter and finalize
validPlates = filterPlatesBySize(plates, plateOptions.minPointsPerPlate);

% Store results
plateData.plates = validPlates;
plateData.numPlates = length(validPlates);
plateData.processingTime = toc;
plateData.plateGenTime = plateGenTime;

% Display results
fprintf('\n=== FAST PLATE GENERATION RESULTS ===\n');
fprintf('Total plates generated: %d\n', plateData.numPlates);
fprintf('Total processing time: %.2f seconds\n', plateData.processingTime);
fprintf('Plate generation time: %.2f seconds\n', plateGenTime);
if plateData.numPlates > 0
    avgPoints = mean(cellfun(@(p) p.numPoints, validPlates));
    fprintf('Average points per plate: %.1f\n', avgPoints);
end
fprintf('Performance: %.0f points/second\n', size(allPoints,1)/plateData.processingTime);

end

function allPoints = createAllPointsFromAlignedFast(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, plateType)
% Create allPoints array from aligned coordinate data for fast processing
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

function allPoints = extractAllPointsFast(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateType)
% Fast point extraction with pre-allocation

% Pre-allocate based on estimated points per waveform
estimatedPointsPerWaveform = 25; % Based on your data
maxPoints = length(peakData) * estimatedPointsPerWaveform;
tempPoints = zeros(maxPoints, 5); % [X, Y, Time, Amplitude, Type]
pointCount = 0;

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
    
    % Get peak/valley data
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
    
    % Extract valid transitions
    validTimes = waveformData.TransitionTime(validIndices);
    validAmps = waveformData.TransitionAmplitude(validIndices);
    validTypes = waveformData.TransitionType(validIndices);
    
    numValidPoints = length(validTimes);
    if numValidPoints == 0
        continue;
    end
    
    % Add to points array
    endIdx = pointCount + numValidPoints;
    tempPoints(pointCount+1:endIdx, 1) = spatialX;
    tempPoints(pointCount+1:endIdx, 2) = spatialY;
    tempPoints(pointCount+1:endIdx, 3) = validTimes;
    tempPoints(pointCount+1:endIdx, 4) = validAmps;
    tempPoints(pointCount+1:endIdx, 5) = validTypes;
    
    pointCount = endIdx;
end

% Trim to actual size
allPoints = tempPoints(1:pointCount, :);
end

function spatialIndex = buildSpatialIndex(allPoints, plateOptions, X_Coordinates, Y_Coordinates)
% Build spatial index using grid-based hashing for O(1) neighbor lookup

% Calculate grid parameters
xMin = min(allPoints(:, 1));
yMin = min(allPoints(:, 2));
xMax = max(allPoints(:, 1));
yMax = max(allPoints(:, 2));

% Calculate actual step size from coordinate data
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

actualStepSize = max(stepSizeX, stepSizeY);

% Grid size should be larger than search radius for efficiency
searchRadius = actualStepSize * 2; % mm - 2x actual step size for neighbor search
gridSize = searchRadius * 2; % Each grid cell covers 2x search radius

fprintf('FastPlateGenerator: Detected step size X=%.3f mm, Y=%.3f mm, Using=%.3f mm\n', ...
    stepSizeX, stepSizeY, actualStepSize);

nGridX = ceil((xMax - xMin) / gridSize) + 1;
nGridY = ceil((yMax - yMin) / gridSize) + 1;

% Initialize grid
spatialIndex = struct();
spatialIndex.grid = cell(nGridY, nGridX);
spatialIndex.xMin = xMin;
spatialIndex.yMin = yMin;
spatialIndex.gridSize = gridSize;
spatialIndex.nGridX = nGridX;
spatialIndex.nGridY = nGridY;
spatialIndex.searchRadius = searchRadius;

% Populate grid with point indices
for i = 1:size(allPoints, 1)
    x = allPoints(i, 1);
    y = allPoints(i, 2);
    
    % Calculate grid coordinates
    gridX = floor((x - xMin) / gridSize) + 1;
    gridY = floor((y - yMin) / gridSize) + 1;
    
    % Clamp to valid range
    gridX = max(1, min(nGridX, gridX));
    gridY = max(1, min(nGridY, gridY));
    
    % Add point index to grid cell
    spatialIndex.grid{gridY, gridX} = [spatialIndex.grid{gridY, gridX}, i];
end

fprintf('✓ Spatial index: %dx%d grid, %.1f mm cells\n', nGridY, nGridX, gridSize);
end

function neighborIndices = findNeighborsFast(allPoints, spatialIndex, pointIdx, plateOptions)
% Fast neighbor lookup using spatial index - O(1) average case

point = allPoints(pointIdx, :);
x = point(1);
y = point(2);
amp = point(4);
time = point(3);
type = point(5);

% Calculate search bounds in grid coordinates
searchRadius = spatialIndex.searchRadius;
xMin = spatialIndex.xMin;
yMin = spatialIndex.yMin;
gridSize = spatialIndex.gridSize;

gridXMin = max(1, floor((x - searchRadius - xMin) / gridSize) + 1);
gridXMax = min(spatialIndex.nGridX, floor((x + searchRadius - xMin) / gridSize) + 1);
gridYMin = max(1, floor((y - searchRadius - yMin) / gridSize) + 1);
gridYMax = min(spatialIndex.nGridY, floor((y + searchRadius - yMin) / gridSize) + 1);

% Collect candidate points from nearby grid cells
candidateIndices = [];
for gx = gridXMin:gridXMax
    for gy = gridYMin:gridYMax
        candidateIndices = [candidateIndices, spatialIndex.grid{gy, gx}];
    end
end

if isempty(candidateIndices)
    neighborIndices = [];
    return;
end

% Filter candidates by actual distance and tolerances
candidatePoints = allPoints(candidateIndices, :);

% Spatial distance check
spatialDist = sqrt((candidatePoints(:, 1) - x).^2 + (candidatePoints(:, 2) - y).^2);
spatialValid = spatialDist <= searchRadius;

% Amplitude tolerance check
ampDiff = abs(candidatePoints(:, 4) - amp);
ampValid = ampDiff <= plateOptions.amplitudeTolerance;

% Time tolerance check
timeDiff = abs(candidatePoints(:, 3) - time);
timeValid = timeDiff <= (plateOptions.timeTolerance * 1e-6);

% Type check
typeValid = candidatePoints(:, 5) == type;

% Combine all criteria
validMask = spatialValid & ampValid & timeValid & typeValid;
neighborIndices = candidateIndices(validMask);

% Remove self
neighborIndices = neighborIndices(neighborIndices ~= pointIdx);
end

function plates = generatePlatesFast(allPoints, spatialIndex, plateOptions)
% Fast plate generation using region growing with spatial index

plates = {};
usedPoints = false(size(allPoints, 1), 1);
plateCount = 0;
totalPoints = size(allPoints, 1);

fprintf('Processing %d points with fast neighbor lookup...\n', totalPoints);
progressInterval = max(1000, floor(totalPoints / 100));
startTime = tic;

for i = 1:totalPoints
    if usedPoints(i)
        continue;
    end
    
    % Start region growing from this seed point
    platePoints = regionGrowFast(allPoints, spatialIndex, i, usedPoints, plateOptions);
    
    if length(platePoints) >= plateOptions.minPointsPerPlate
        plateCount = plateCount + 1;
        
        % Create plate structure
        platePointData = allPoints(platePoints, :);
        
        plates{plateCount} = struct();
        plates{plateCount}.id = plateCount;
        plates{plateCount}.points = platePointData;
        plates{plateCount}.numPoints = length(platePoints);
        plates{plateCount}.centerAmplitude = mean(platePointData(:, 4));
        plates{plateCount}.centerTime = mean(platePointData(:, 3));
        plates{plateCount}.amplitudeRange = [min(platePointData(:, 4)), max(platePointData(:, 4))];
        plates{plateCount}.timeRange = [min(platePointData(:, 3)), max(platePointData(:, 3))];
        plates{plateCount}.spatialBounds = [min(platePointData(:, 1)), max(platePointData(:, 1)), ...
                                           min(platePointData(:, 2)), max(platePointData(:, 2))];
        plates{plateCount}.plateType = platePointData(1, 5);
        plates{plateCount}.isFastGenerated = true;
        
        % Mark points as used
        usedPoints(platePoints) = true;
        
        % Progress reporting
        if mod(plateCount, 10) == 0
            fprintf('Generated %d plates... (%.1f%% points processed, %.1f sec)\n', ...
                plateCount, (i/totalPoints)*100, toc(startTime));
        end
    end
    
    % Overall progress
    if mod(i, progressInterval) == 0
        fprintf('Progress: %.1f%% (%d/%d points, %d plates, %.1f sec)\n', ...
            (i/totalPoints)*100, i, totalPoints, plateCount, toc(startTime));
    end
end

fprintf('✓ Fast generation complete: %d plates from %d points\n', plateCount, totalPoints);
end

function platePoints = regionGrowFast(allPoints, spatialIndex, seedIdx, usedPoints, plateOptions)
% Fast region growing using spatial index for neighbor lookup

platePoints = [];
toVisit = [seedIdx];
visited = false(size(allPoints, 1), 1);
visited(seedIdx) = true;

while ~isempty(toVisit)
    currentIdx = toVisit(1);
    toVisit(1) = [];
    
    if usedPoints(currentIdx)
        continue;
    end
    
    platePoints = [platePoints, currentIdx];
    
    % Find neighbors using fast spatial index
    neighbors = findNeighborsFast(allPoints, spatialIndex, currentIdx, plateOptions);
    
    % Add unvisited neighbors to queue
    for neighborIdx = neighbors
        if ~visited(neighborIdx) && ~usedPoints(neighborIdx)
            visited(neighborIdx) = true;
            toVisit = [toVisit, neighborIdx];
        end
    end
end
end

function validPlates = filterPlatesBySize(plates, minPointsPerPlate)
% Filter plates by minimum size requirement

validPlates = {};
validCount = 0;

for i = 1:length(plates)
    if plates{i}.numPoints >= minPointsPerPlate
        validCount = validCount + 1;
        validPlates{validCount} = plates{i};
        validPlates{validCount}.id = validCount;
    end
end

fprintf('✓ Filtered plates: %d valid plates (minimum %d points each)\n', ...
    validCount, minPointsPerPlate);
end
