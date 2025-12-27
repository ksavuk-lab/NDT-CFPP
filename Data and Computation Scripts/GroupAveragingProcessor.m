function plateData = GroupAveragingProcessor(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions)
% GroupAveragingProcessor - Algorithm #3: Group Averaging
%
% This function implements group averaging for binding:
% 1. Start with initial groups based on amplitude + time tolerance
% 2. Calculate group averages (position, amplitude, time)
% 3. Use group average as basis for binding new points
% 4. Points must be within tolerance of group averages (not individual points)
% 5. Creates more stable, continuous groups for large defect areas
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

fprintf('\n=== Algorithm #3: Group Averaging ===\n');

% Start timing
tic;

% Initialize output structure
plateData = struct();
plateData.plates = [];
plateData.numPlates = 0;
plateData.processingTime = 0;
plateData.statistics = struct();

% Group averaging parameters
maxIterations = 5; % Maximum refinement iterations
convergenceThreshold = 0.01; % Convergence threshold for group averages
spatialTolerance = 3.0; % mm - spatial tolerance for group membership

fprintf('Starting group averaging with settings:\n');
fprintf('  Amplitude tolerance: ±%.3f\n', plateOptions.amplitudeTolerance);
fprintf('  Time tolerance: ±%.1f μs (±%.2e sec)\n', plateOptions.timeTolerance, plateOptions.timeTolerance * 1e-6);
fprintf('  Minimum points per plate: %d\n', plateOptions.minPointsPerPlate);
fprintf('  Plate type: %s\n', plateOptions.plateType);
fprintf('  Max iterations: %d\n', maxIterations);
fprintf('  Convergence threshold: %.3f\n', convergenceThreshold);
fprintf('  Spatial tolerance: %.1f mm\n', spatialTolerance);

% Extract all points from peak data
fprintf('Extracting points from %d waveforms...\n', length(peakData));
allPoints = extractAllPoints(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateOptions.plateType);

if isempty(allPoints)
    fprintf('Warning: No points extracted for plate generation.\n');
    plateData.processingTime = toc;
    return;
end

fprintf('Extracted %d total points for group averaging.\n', size(allPoints, 1));

% Apply group averaging algorithm
fprintf('Applying group averaging algorithm...\n');

% Create progress dialog
progressHandle = ProgressDialog('Algorithm #3: Group Averaging', ...
                               'Starting group averaging...', ...
                               size(allPoints, 1));

plates = groupAveragingBindingWithProgress(allPoints, plateOptions, maxIterations, convergenceThreshold, spatialTolerance, progressHandle);

% Close progress dialog
closeProgress(progressHandle);

% Filter plates by minimum point requirement
validPlates = filterPlatesBySize(plates, plateOptions.minPointsPerPlate);

% Store results
plateData.plates = validPlates;
plateData.numPlates = length(validPlates);
plateData.processingTime = toc;

% Generate statistics
plateData.statistics = generatePlateStatistics(validPlates, allPoints);

% Display results
fprintf('\n=== Group Averaging Results ===\n');
fprintf('Total plates generated: %d\n', plateData.numPlates);
fprintf('Processing time: %.2f seconds\n', plateData.processingTime);
if plateData.numPlates > 0
    fprintf('Average points per plate: %.1f\n', plateData.statistics.avgPointsPerPlate);
    fprintf('Largest plate: %d points\n', plateData.statistics.maxPointsPerPlate);
    fprintf('Smallest plate: %d points\n', plateData.statistics.minPointsPerPlate);
    fprintf('Point utilization: %.1f%%\n', plateData.statistics.pointUtilization);
    fprintf('Average iterations: %.1f\n', plateData.statistics.avgIterations);
    fprintf('✓ All plates use group averaging for stable binding\n');
end

end

function allPoints = extractAllPoints(peakData, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, plateType)
% Extract all relevant points from peak data based on plate type
% Same as other algorithms for consistency

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

function plates = groupAveragingBindingWithProgress(allPoints, plateOptions, maxIterations, convergenceThreshold, spatialTolerance, progressHandle)
% ALGORITHM #3: Group Averaging with Progress Tracking
% Start with initial groups, refine using group averages as binding criteria

% Pre-allocate for better performance
maxPlates = 1000;
plates = cell(maxPlates, 1);
plateCount = 0;
usedPoints = false(size(allPoints, 1), 1);

% Convert time tolerance to seconds
timeToleranceSeconds = plateOptions.timeTolerance * 1e-6;
totalPoints = size(allPoints, 1);

fprintf('Processing %d points with group averaging...\n', totalPoints);

% Step 1: Create initial groups using simple amplitude + time tolerance
fprintf('Step 1: Creating initial groups...\n');
updateProgress(progressHandle, totalPoints * 0.1, 'Step 1: Starting initial group creation...', 'createInitialGroupsWithProgress');
initialGroups = createInitialGroupsWithProgress(allPoints, plateOptions.amplitudeTolerance, timeToleranceSeconds, spatialTolerance, progressHandle, totalPoints * 0.1, totalPoints * 0.5);
fprintf('  → Created %d initial groups\n', length(initialGroups));
updateProgress(progressHandle, totalPoints * 0.5, sprintf('Step 1 complete: Created %d initial groups', length(initialGroups)), 'createInitialGroupsWithProgress');

% Step 2: Refine groups using iterative averaging
fprintf('Step 2: Refining groups using iterative averaging...\n');
updateProgress(progressHandle, totalPoints * 0.5, 'Step 2: Starting iterative group refinement...', 'refineGroupsWithAveragingWithProgress');
refinedGroups = refineGroupsWithAveragingWithProgress(allPoints, initialGroups, plateOptions.amplitudeTolerance, timeToleranceSeconds, spatialTolerance, maxIterations, convergenceThreshold, progressHandle, totalPoints * 0.5, totalPoints * 0.8);
fprintf('  → Refined to %d stable groups\n', length(refinedGroups));
updateProgress(progressHandle, totalPoints * 0.8, sprintf('Step 2 complete: Refined to %d stable groups', length(refinedGroups)), 'refineGroupsWithAveragingWithProgress');

% Step 3: Convert groups to plates
fprintf('Step 3: Converting groups to plates...\n');
for groupIdx = 1:length(refinedGroups)
    group = refinedGroups{groupIdx};
    
    if length(group.members) >= plateOptions.minPointsPerPlate
        plateCount = plateCount + 1;
        
        % Expand if needed
        if plateCount > maxPlates
            plates = [plates; cell(maxPlates, 1)];
            maxPlates = maxPlates * 2;
        end
        
        % Create plate structure
        platePoints = allPoints(group.members, :);
        
        plates{plateCount} = struct();
        plates{plateCount}.id = plateCount;
        plates{plateCount}.points = platePoints;
        plates{plateCount}.numPoints = size(platePoints, 1);
        plates{plateCount}.centerAmplitude = group.avgAmplitude;
        plates{plateCount}.centerTime = group.avgTime;
        plates{plateCount}.centerPosition = [group.avgX, group.avgY];
        plates{plateCount}.amplitudeRange = [min(platePoints(:, 4)), max(platePoints(:, 4))];
        plates{plateCount}.timeRange = [min(platePoints(:, 3)), max(platePoints(:, 3))];
        plates{plateCount}.spatialBounds = [min(platePoints(:, 1)), max(platePoints(:, 1)), ...
                                           min(platePoints(:, 2)), max(platePoints(:, 2))];
        plates{plateCount}.plateType = platePoints(1, 5); % All same type in group
        plates{plateCount}.algorithm = 'group_averaging';
        plates{plateCount}.iterations = group.iterations;
        plates{plateCount}.convergence = group.converged;
        
        % Mark points as used
        usedPoints(group.members) = true;
    end
end

% Trim to actual size
plates = plates(1:plateCount);

fprintf('✓ Completed group averaging. Generated %d plates.\n', plateCount);

% Final progress update
updateProgress(progressHandle, totalPoints, sprintf('Complete: Generated %d plates with group averaging', plateCount), 'groupAveragingBindingWithProgress');

end

function initialGroups = createInitialGroups(allPoints, ampTolerance, timeTolerance, spatialTolerance)
% Create initial groups using simple amplitude + time + spatial tolerance

groups = {};
groupCount = 0;
usedPoints = false(size(allPoints, 1), 1);

for i = 1:size(allPoints, 1)
    if usedPoints(i)
        continue;
    end

    seedPoint = allPoints(i, :);
    seedAmp = seedPoint(4);
    seedTime = seedPoint(3);
    seedType = seedPoint(5);

    % Find all points within tolerance
    ampDiff = abs(allPoints(:, 4) - seedAmp);
    timeDiff = abs(allPoints(:, 3) - seedTime);
    typeSame = allPoints(:, 5) == seedType;

    candidates = ~usedPoints & ...
                ampDiff <= ampTolerance & ...
                timeDiff <= timeTolerance & ...
                typeSame;

    if sum(candidates) >= 3 % Minimum for initial group
        % Apply spatial connectivity
        candidateIndices = find(candidates);
        connectedIndices = findSpatiallyConnectedGroup(allPoints, candidateIndices, i, spatialTolerance);

        if length(connectedIndices) >= 3
            groupCount = groupCount + 1;

            % Calculate initial group averages
            groupPoints = allPoints(connectedIndices, :);
            groups{groupCount} = struct();
            groups{groupCount}.members = connectedIndices;
            groups{groupCount}.avgX = mean(groupPoints(:, 1));
            groups{groupCount}.avgY = mean(groupPoints(:, 2));
            groups{groupCount}.avgTime = mean(groupPoints(:, 3));
            groups{groupCount}.avgAmplitude = mean(groupPoints(:, 4));
            groups{groupCount}.type = seedType;
            groups{groupCount}.iterations = 0;
            groups{groupCount}.converged = false;

            usedPoints(connectedIndices) = true;
        end
    end
end

initialGroups = groups;

end

function refinedGroups = refineGroupsWithAveraging(allPoints, initialGroups, ampTolerance, timeTolerance, spatialTolerance, maxIterations, convergenceThreshold)
% Refine groups using iterative averaging - key part of Algorithm #3

refinedGroups = initialGroups;

for groupIdx = 1:length(refinedGroups)
    group = refinedGroups{groupIdx};

    for iteration = 1:maxIterations
        oldAvgX = group.avgX;
        oldAvgY = group.avgY;
        oldAvgTime = group.avgTime;
        oldAvgAmp = group.avgAmplitude;

        % Find all points within tolerance of current group averages
        spatialDist = sqrt((allPoints(:, 1) - group.avgX).^2 + (allPoints(:, 2) - group.avgY).^2);
        timeDiff = abs(allPoints(:, 3) - group.avgTime);
        ampDiff = abs(allPoints(:, 4) - group.avgAmplitude);
        typeSame = allPoints(:, 5) == group.type;

        newMembers = find(spatialDist <= spatialTolerance & ...
                         timeDiff <= timeTolerance & ...
                         ampDiff <= ampTolerance & ...
                         typeSame);

        if length(newMembers) >= 3
            % Update group averages
            newGroupPoints = allPoints(newMembers, :);
            group.avgX = mean(newGroupPoints(:, 1));
            group.avgY = mean(newGroupPoints(:, 2));
            group.avgTime = mean(newGroupPoints(:, 3));
            group.avgAmplitude = mean(newGroupPoints(:, 4));
            group.members = newMembers;

            % Check convergence
            positionChange = sqrt((group.avgX - oldAvgX)^2 + (group.avgY - oldAvgY)^2);
            timeChange = abs(group.avgTime - oldAvgTime);
            ampChange = abs(group.avgAmplitude - oldAvgAmp);

            if positionChange < convergenceThreshold && ...
               timeChange < (timeTolerance * convergenceThreshold) && ...
               ampChange < (ampTolerance * convergenceThreshold)
                group.converged = true;
                group.iterations = iteration;
                break;
            end
        else
            % Group became too small, revert to previous iteration
            break;
        end

        group.iterations = iteration;
    end

    refinedGroups{groupIdx} = group;
end

end

function connectedIndices = findSpatiallyConnectedGroup(allPoints, candidateIndices, seedIdx, spatialTolerance)
% Find spatially connected group starting from seed point
% Same as Algorithm #1 for consistency - OPTIMIZED to avoid dynamic arrays

numCandidates = length(candidateIndices);
visited = false(numCandidates, 1);

% Pre-allocate arrays to avoid dynamic growth
maxConnected = numCandidates;
connectedArray = zeros(maxConnected, 1);
toVisitArray = zeros(maxConnected, 1);
connectedCount = 0;
toVisitCount = 0;
toVisitIdx = 1;

% Find seed in candidate list
seedInCandidates = find(candidateIndices == seedIdx);
if isempty(seedInCandidates)
    connectedIndices = [];
    return;
end

% Start region growing from seed
toVisitCount = 1;
toVisitArray(1) = seedInCandidates;
visited(seedInCandidates) = true;

while toVisitIdx <= toVisitCount
    currentCandidateIdx = toVisitArray(toVisitIdx);
    toVisitIdx = toVisitIdx + 1;

    currentPointIdx = candidateIndices(currentCandidateIdx);
    connectedCount = connectedCount + 1;
    connectedArray(connectedCount) = currentPointIdx;

    currentPoint = allPoints(currentPointIdx, :);
    currentX = currentPoint(1);
    currentY = currentPoint(2);

    % Find unvisited neighbors within spatial tolerance
    for i = 1:numCandidates
        if visited(i)
            continue;
        end

        neighborPointIdx = candidateIndices(i);
        neighborPoint = allPoints(neighborPointIdx, :);
        neighborX = neighborPoint(1);
        neighborY = neighborPoint(2);

        % Check spatial distance
        spatialDist = sqrt((neighborX - currentX)^2 + (neighborY - currentY)^2);

        if spatialDist <= spatialTolerance
            visited(i) = true;
            toVisitCount = toVisitCount + 1;
            toVisitArray(toVisitCount) = i;
        end
    end
end

% Return trimmed array
connectedIndices = connectedArray(1:connectedCount);

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
    stats.avgIterations = 0;
    return;
end

% Basic statistics
pointCounts = cellfun(@(p) p.numPoints, plates);
stats.avgPointsPerPlate = mean(pointCounts);
stats.maxPointsPerPlate = max(pointCounts);
stats.minPointsPerPlate = min(pointCounts);
stats.totalPointsInPlates = sum(pointCounts);
stats.pointUtilization = stats.totalPointsInPlates / size(allPoints, 1) * 100;

% Group averaging specific statistics
if isfield(plates{1}, 'iterations')
    iterations = cellfun(@(p) p.iterations, plates);
    stats.avgIterations = mean(iterations);
else
    stats.avgIterations = 0;
end

end

function initialGroups = createInitialGroupsWithProgress(allPoints, ampTolerance, timeTolerance, spatialTolerance, progressHandle, startProgress, endProgress)
% Create initial groups using simple amplitude + time + spatial tolerance with progress tracking

groups = {};
groupCount = 0;
usedPoints = false(size(allPoints, 1), 1);
totalPoints = size(allPoints, 1);

for i = 1:totalPoints
    if usedPoints(i)
        continue;
    end

    % Update progress every 2000 points
    if mod(i, 2000) == 0
        progressValue = startProgress + (i / totalPoints) * (endProgress - startProgress);
        statusMessage = sprintf('Scanning point %d/%d for initial groups (%d groups found)', i, totalPoints, groupCount);
        updateProgress(progressHandle, progressValue, statusMessage, 'createInitialGroupsWithProgress');
    end

    seedPoint = allPoints(i, :);
    seedAmp = seedPoint(4);
    seedTime = seedPoint(3);
    seedType = seedPoint(5);

    % Find all points within tolerance
    ampDiff = abs(allPoints(:, 4) - seedAmp);
    timeDiff = abs(allPoints(:, 3) - seedTime);
    typeSame = allPoints(:, 5) == seedType;

    candidates = ~usedPoints & ...
                ampDiff <= ampTolerance & ...
                timeDiff <= timeTolerance & ...
                typeSame;

    if sum(candidates) >= 3 % Minimum for initial group
        % Apply spatial connectivity
        candidateIndices = find(candidates);
        connectedIndices = findSpatiallyConnectedGroup(allPoints, candidateIndices, i, spatialTolerance);

        if length(connectedIndices) >= 3
            groupCount = groupCount + 1;

            % Calculate initial group averages
            groupPoints = allPoints(connectedIndices, :);
            groups{groupCount} = struct();
            groups{groupCount}.members = connectedIndices;
            groups{groupCount}.avgX = mean(groupPoints(:, 1));
            groups{groupCount}.avgY = mean(groupPoints(:, 2));
            groups{groupCount}.avgTime = mean(groupPoints(:, 3));
            groups{groupCount}.avgAmplitude = mean(groupPoints(:, 4));
            groups{groupCount}.type = seedType;
            groups{groupCount}.iterations = 0;
            groups{groupCount}.converged = false;

            usedPoints(connectedIndices) = true;
        end
    end
end

initialGroups = groups;

end

function refinedGroups = refineGroupsWithAveragingWithProgress(allPoints, initialGroups, ampTolerance, timeTolerance, spatialTolerance, maxIterations, convergenceThreshold, progressHandle, startProgress, endProgress)
% Refine groups using iterative averaging with progress tracking

refinedGroups = initialGroups;
totalGroups = length(refinedGroups);

for groupIdx = 1:totalGroups
    % Update progress for each group
    progressValue = startProgress + (groupIdx / totalGroups) * (endProgress - startProgress);
    statusMessage = sprintf('Refining group %d/%d (iterative averaging)', groupIdx, totalGroups);
    updateProgress(progressHandle, progressValue, statusMessage, 'refineGroupsWithAveragingWithProgress');

    group = refinedGroups{groupIdx};

    for iteration = 1:maxIterations
        oldAvgX = group.avgX;
        oldAvgY = group.avgY;
        oldAvgTime = group.avgTime;
        oldAvgAmp = group.avgAmplitude;

        % Find all points within tolerance of current group averages
        spatialDist = sqrt((allPoints(:, 1) - group.avgX).^2 + (allPoints(:, 2) - group.avgY).^2);
        timeDiff = abs(allPoints(:, 3) - group.avgTime);
        ampDiff = abs(allPoints(:, 4) - group.avgAmplitude);
        typeSame = allPoints(:, 5) == group.type;

        newMembers = find(spatialDist <= spatialTolerance & ...
                         timeDiff <= timeTolerance & ...
                         ampDiff <= ampTolerance & ...
                         typeSame);

        if length(newMembers) >= 3
            % Update group averages
            newGroupPoints = allPoints(newMembers, :);
            group.avgX = mean(newGroupPoints(:, 1));
            group.avgY = mean(newGroupPoints(:, 2));
            group.avgTime = mean(newGroupPoints(:, 3));
            group.avgAmplitude = mean(newGroupPoints(:, 4));
            group.members = newMembers;

            % Check convergence
            positionChange = sqrt((group.avgX - oldAvgX)^2 + (group.avgY - oldAvgY)^2);
            timeChange = abs(group.avgTime - oldAvgTime);
            ampChange = abs(group.avgAmplitude - oldAvgAmp);

            if positionChange < convergenceThreshold && ...
               timeChange < (timeTolerance * convergenceThreshold) && ...
               ampChange < (ampTolerance * convergenceThreshold)
                group.converged = true;
                group.iterations = iteration;
                break;
            end
        else
            % Group became too small, revert to previous iteration
            break;
        end

        group.iterations = iteration;
    end

    refinedGroups{groupIdx} = group;
end

end
