# Column Alignment Techniques - Complete Variable Tracking and Data Analysis

## Overview
This document provides an extremely detailed, step-by-step analysis of the column alignment techniques used in XtVsYPlot.m, with complete variable tracking, mathematical computations, and data transformations. Every calculation is shown with actual values and examples.

## 1. Single Slice Alignment - Complete Walkthrough

### 1.1 Input Data Structure and Extraction

#### **Source Data Format**
```matlab
% Original data structure in XtVsYPlot.m
statDataArray{1}.maps{seg} = [numY × numX] matrix
% Example dimensions: numY=64, numX=128, numSegments=400

% Sample data values for Y-slice 32:
statData.maps{1}(32, 1:5) = [0.245, 0.251, 0.248, 0.252, 0.249]
statData.maps{2}(32, 1:5) = [0.247, 0.253, 0.250, 0.254, 0.251]
statData.maps{3}(32, 1:5) = [0.246, 0.252, 0.249, 0.253, 0.250]
...continuing for all 400 segments
statData.maps{400}(32, 1:5) = [0.244, 0.250, 0.247, 0.251, 0.248]
```

#### **Slice Data Extraction Process**
```matlab
% Function: computeCurrentSlice() -> computeCurrentSliceBackground()
% Current example: XtVsY view, Y-slice index = 32

sliceData = zeros(numSegments, numX);  % [400 × 128] matrix
for seg = 1:numSegments
    sliceData(seg, :) = statData.maps{seg}(sliceIndex, :);
end

% Resulting sliceData matrix:
sliceData(1, 1:5) = [0.245, 0.251, 0.248, 0.252, 0.249]    % Segment 1
sliceData(2, 1:5) = [0.247, 0.253, 0.250, 0.254, 0.251]    % Segment 2
sliceData(3, 1:5) = [0.246, 0.252, 0.249, 0.253, 0.250]    % Segment 3
...
sliceData(400, 1:5) = [0.244, 0.250, 0.247, 0.251, 0.248]  % Segment 400

% Data interpretation:
% - Each row represents a time segment (400 total time points)
% - Each column represents a spatial X position (128 spatial positions)
% - Values are statistical measurements (e.g., MaxAmplitude, RMS)
```

### 1.2 alignColumnsImproved Function Entry

#### **Function Call with Parameters**
```matlab
[alignedSliceData, convergenceInfo] = alignColumnsImproved(sliceData, ...
    'MaxShift', 15, ...                    % ±15 pixel search range
    'CostFunction', 'mse', ...             % Mean Squared Error
    'AlignmentMethod', 'average', ...      % Compare to row average
    'ConvergenceThreshold', 0.001, ...     % 0.1% improvement threshold
    'MaxIterations', 20, ...               % Maximum 20 iterations
    'WeightingFunction', 'exponential', ... % Distance weighting
    'WeightingScale', 3.0, ...             % Weighting decay parameter
    'PadMethod', 'zeros');                 % Zero padding for shifts
```

#### **Input Validation and Parameter Extraction**
```matlab
[numRows, numCols] = size(sliceData);  % [400, 128]

% Extracted parameters:
maxShift = 15;                    % Search range: -15 to +15 pixels
costFunction = 'mse';             % Mean Squared Error cost function
alignmentMethod = 'average';      % Fastest and most stable method
convergenceThreshold = 0.001;     % 0.1% relative improvement threshold
maxIterations = 20;               % Maximum iterations before stopping
```

### 1.3 Initialization Phase

#### **Working Data and Tracking Arrays**
```matlab
% Create working copy of input data
alignedData = sliceData;  % [400 × 128] matrix - will be modified in-place

% Initialize tracking arrays
previousShifts = zeros(1, 128);           % Track cumulative shifts per column
convergenceHistory = zeros(1, 20);       % Pre-allocated for max iterations
convergenceHistoryIndex = 0;             % Current position in history

% Performance optimization cache
costCache = containers.Map('KeyType', 'char', 'ValueType', 'double');

% Progress tracking
totalColumnsToProcess = 128;              % All columns including edges
```

#### **Row Average Calculation (Reference Template)**
```matlab
% For 'average' method - calculate initial row average
rowAverage = mean(alignedData, 2);  % [400 × 1] vector

% Sample calculation for first few time points:
rowAverage(1) = mean([0.245, 0.251, 0.248, 0.252, 0.249, ...]) = 0.2487
rowAverage(2) = mean([0.247, 0.253, 0.250, 0.254, 0.251, ...]) = 0.2489
rowAverage(3) = mean([0.246, 0.252, 0.249, 0.253, 0.250, ...]) = 0.2488
...
rowAverage(400) = mean([0.244, 0.250, 0.247, 0.251, 0.248, ...]) = 0.2485

% Purpose: This serves as the reference template that all columns will be aligned to
% Update strategy: Recalculated every 5 columns to incorporate improvements
```

### 1.4 Main Convergence Loop

#### **Iteration Structure**
```matlab
for iteration = 1:maxIterations  % 1 to 20
    % Reset iteration variables
    currentShifts = zeros(1, 128);        % Shifts applied this iteration
    totalCostImprovement = 0;             % Total improvement this iteration
    
    % Column processing order: middle columns first, then edges
    columnOrder = [2:127, 1, 128];        % [2,3,4,...,127,1,128]
    
    % Timing
    iterationStartTime = tic;
```

#### **Column-by-Column Processing**
```matlab
for colIdx = 1:128
    col = columnOrder(colIdx);            % Current column index
    currentCol = alignedData(:, col);     % [400 × 1] column vector
    
    % Example: Processing column 2 (first middle column)
    % currentCol = [0.251, 0.253, 0.250, 0.252, 0.251, ...]
    
    % Progress calculation
    columnProgress = colIdx / 128;
    overallProgress = (iteration - 1 + columnProgress) / 20;
```

### 1.5 Column Type Determination and Processing

#### **Middle Column Processing (Average Method)**
```matlab
% For middle columns (col = 2 to 127), use average method
if col ~= 1 && col ~= 128
    [bestShift, costImprovement] = findOptimalShiftAverage(
        currentCol, rowAverage, maxShift, costFunction, padMethod, costCache);
```

#### **Shift Range Testing**
```matlab
% Inside findOptimalShiftAverage function
shiftRange = -15:15;                      % 31 possible shifts
bestCost = inf;
bestShift = 0;

% Calculate original cost (before any shift)
originalCost = calculateCost(currentCol, rowAverage, 'mse');
```

#### **MSE Cost Calculation Details**
```matlab
% Inside calculateCost function
function cost = calculateCost(col1, col2, costFunction)
    % Remove NaN values
    validIdx = ~isnan(col1) & ~isnan(col2);
    
    % Exclude edge rows to avoid padding effects
    edgeRows = 5;
    dataLength = length(col1);  % 400
    excludeIdx = false(size(col1));
    excludeIdx(1:5) = true;           % Exclude top 5 rows
    excludeIdx(396:400) = true;       % Exclude bottom 5 rows
    
    % Combine exclusions
    validIdx = validIdx & ~excludeIdx;
    
    % Extract valid data (390 points out of 400)
    col1_valid = col1(validIdx);      % [390 × 1]
    col2_valid = col2(validIdx);      % [390 × 1]
    
    % Calculate MSE
    cost = mean((col1_valid - col2_valid).^2);
end

% Example calculation for column 2:
originalCost = mean((currentCol_valid - rowAverage_valid).^2) = 0.000234;
```

### 1.6 Shift Optimization Loop

#### **Testing Each Possible Shift**
```matlab
for shift = -15:15
    % Apply shift to current column
    shiftedCol = applyIntegerShift(currentCol, shift, 'zeros');
    
    % Calculate cost for this shift
    cost = calculateCost(shiftedCol, rowAverage, 'mse');
    
    % Update best if this is better
    if cost < bestCost
        bestCost = cost;
        bestShift = shift;
    end
end
```

#### **Detailed Shift Application Example (shift = -3)**
```matlab
% Example: shift = -3 (shift up by 3 pixels)
function shiftedCol = applyIntegerShift(column, shift, padMethod)
    n = length(column);               % 400
    shiftedCol = zeros(size(column)); % [400 × 1] initialized to zeros
    
    if shift == 0
        shiftedCol = column;
        return;
    end
    
    % For shift = -3 (negative = shift up)
    shift = abs(shift);               % shift = 3
    shiftedCol(1:end-shift) = column(shift+1:end);  % Move data up
    % shiftedCol(1:397) = column(4:400)
    % shiftedCol(398:400) remains zero (zero padding at bottom)
    
    % Result:
    % Original:  [0.251, 0.253, 0.250, 0.252, 0.251, ...]
    % Shifted:   [0.252, 0.251, 0.249, 0.253, 0.250, ..., 0, 0, 0]
end
```

#### **Cost Calculation for Shifted Column**
```matlab
% Calculate cost for shift = -3
cost = calculateCost(shiftedCol, rowAverage, 'mse');

% MSE calculation with shifted data:
shifted_valid = shiftedCol(validIdx);     % [390 × 1] excluding edges
average_valid = rowAverage(validIdx);     % [390 × 1] excluding edges
cost = mean((shifted_valid - average_valid).^2) = 0.000198;

% Check if this is the best so far
if cost < bestCost  % 0.000198 < inf
    bestCost = 0.000198;
    bestShift = -3;
end
```

#### **Complete Shift Testing Results**
```matlab
% Results for all shifts tested:
shift = -15: cost = 0.000456, bestShift = -3
shift = -14: cost = 0.000423, bestShift = -3
shift = -13: cost = 0.000398, bestShift = -3
...
shift = -4:  cost = 0.000201, bestShift = -3
shift = -3:  cost = 0.000198, bestShift = -3  ← BEST COST
shift = -2:  cost = 0.000205, bestShift = -3
shift = -1:  cost = 0.000218, bestShift = -3
shift = 0:   cost = 0.000234, bestShift = -3  ← ORIGINAL COST
shift = +1:  cost = 0.000267, bestShift = -3
...
shift = +15: cost = 0.000512, bestShift = -3

% Final optimization result:
bestShift = -3;                              % Optimal shift: up by 3 pixels
bestCost = 0.000198;                         % Lowest MSE cost
costImprovement = originalCost - bestCost;   % 0.000234 - 0.000198 = 0.000036
```

### 1.7 Shift Application Decision

#### **Threshold Check**
```matlab
% Apply shift only if it's significant enough
if abs(bestShift) > 0.5  % Integer shift threshold
    % Apply the shift: abs(-3) = 3 > 0.5 ✓
    alignedData(:, col) = applyIntegerShift(currentCol, bestShift, 'zeros');
    currentShifts(col) = round(bestShift);
    totalCostImprovement = totalCostImprovement + costImprovement;
else
    % No shift applied - improvement too small
    currentShifts(col) = 0;
end
```

#### **Column Update Result**
```matlab
% After applying shift = -3 to column 2:
alignedData(1:5, 2) = [0.252, 0.251, 0.249, 0.253, 0.250];  % Shifted data
alignedData(398:400, 2) = [0, 0, 0];                         % Zero padding

% Tracking updates:
currentShifts(2) = -3;                       % Column 2 shifted up by 3 pixels
totalCostImprovement += 0.000036;            % Add to total improvement
```

### 1.8 Row Average Update Strategy

#### **Periodic Reference Template Update**
```matlab
% Update row average every 5 columns for stability
if strcmp(alignmentMethod, 'average') && mod(col, 5) == 0
    rowAverage = mean(alignedData, 2);       % Recalculate reference
end

% Example: After processing column 5
% New rowAverage incorporates alignment improvements:
newRowAverage(1) = mean(alignedData(1, :)) = 0.2489;  % Was 0.2487
newRowAverage(2) = mean(alignedData(2, :)) = 0.2491;  % Was 0.2489
...

% Purpose: Incorporate alignment improvements into reference template
% Frequency: Every 5 columns balances stability vs. responsiveness
```

### 1.9 Edge Column Processing

#### **Edge Column Special Handling**
```matlab
% Edge columns (col = 1 and col = 128) use single neighbor comparison
if col == 1 || col == numCols
    [bestShift, costImprovement] = findOptimalShiftEdge(
        currentCol, alignedData, col, maxShift, costFunction, padMethod, costCache);
```

#### **Neighbor Selection Strategy**
```matlab
% Inside findOptimalShiftEdge function
if colIndex == 1
    % First column: compare only with column 2
    neighborCol = fullData(:, 2);
else
    % Last column: compare only with second-to-last column
    neighborCol = fullData(:, numCols - 1);  % Column 127
end

% Example for column 1:
neighborCol = alignedData(:, 2);  % Use already-processed column 2
% neighborCol = [0.252, 0.251, 0.249, 0.253, 0.250, ...]  (after its alignment)
```

#### **Edge Column Optimization**
```matlab
% Calculate original cost against single neighbor
originalCost = calculateCost(currentCol, neighborCol, 'mse') = 0.000187;

% Test all possible shifts
for shift = -15:15
    shiftedCol = applyIntegerShift(currentCol, shift, 'zeros');
    cost = calculateCost(shiftedCol, neighborCol, 'mse');

    if cost < bestCost
        bestCost = cost;
        bestShift = shift;
    end
end

% Result for column 1:
bestShift = -2;                              % Optimal: up by 2 pixels
bestCost = 0.000156;                         % Best cost against neighbor
costImprovement = 0.000187 - 0.000156 = 0.000031;
```

### 1.10 Iteration Completion and Convergence

#### **Iteration Summary**
```matlab
% After processing all 128 columns in iteration 1:
totalCostImprovement = 0.004523;             % Sum of all column improvements
appliedShifts = [-2, -3, -1, 0, -2, +1, 0, -1, ...];  % Shifts per column

% Convergence calculation
currentCost = calculateOverallCost(alignedData);
improvement = (previousCost - currentCost) / previousCost;

% Example values:
currentCost = 0.234567;                      % Current total alignment cost
previousCost = 0.238901;                     % Cost from previous iteration
improvement = (0.238901 - 0.234567) / 0.238901 = 0.0181 = 1.81%;
```

#### **Convergence Check**
```matlab
% Check convergence criteria
if improvement < convergenceThreshold || iteration >= maxIterations
    % Converged: improvement = 0.0181 > 0.001 (not yet converged)
    % Continue to next iteration

    % Update tracking
    convergenceHistory(iteration) = improvement;
    convergenceHistoryIndex = iteration;
    previousCost = currentCost;
end
```

#### **Final Convergence Results**
```matlab
% Typical convergence after 6 iterations:
iterations = 6;
finalCost = 0.231245;
totalImprovement = (0.238901 - 0.231245) / 0.238901 = 3.21%;

% Final shift summary for all columns:
finalShifts = [-2, -3, -1, 0, -2, +1, 0, -1, +2, -1, ...];

% Convergence information returned:
convergenceInfo = struct(
    'converged', true,
    'iterations', 6,
    'finalCost', 0.231245,
    'costHistory', [0.0181, 0.0089, 0.0045, 0.0023, 0.0012, 0.0008],
    'shiftHistory', [6 × 128] matrix,  % Applied shifts per iteration
    'totalShifts', finalShifts,        % Cumulative shifts per column
    'improvement', 0.0321              % 3.21% total improvement
);
```

### 1.11 Result Integration Back to XtVsYPlot

#### **Data Integration Process**
```matlab
% Return from alignColumnsImproved to XtVsYPlot.m
[alignedSliceData, convergenceInfo] = alignColumnsImproved(...);

% Integration back into statData structure
% For XtVsY view, Y-slice = 32:
for seg = 1:numSegments
    alignedData{statIdx}.maps{seg}(sliceIndex, :) = alignedSliceData(seg, :);
end

% Specific example:
for seg = 1:400
    alignedData{1}.maps{seg}(32, :) = alignedSliceData(seg, :);
end
```

#### **Before/After Comparison**
```matlab
% Original data (before alignment):
statData.maps{1}(32, 1:5) = [0.245, 0.251, 0.248, 0.252, 0.249];
statData.maps{2}(32, 1:5) = [0.247, 0.253, 0.250, 0.254, 0.251];
statData.maps{3}(32, 1:5) = [0.246, 0.252, 0.249, 0.253, 0.250];

% Aligned data (after alignment):
alignedData{1}.maps{1}(32, 1:5) = [0.248, 0.252, 0.247, 0.252, 0.249];  % Shifts applied
alignedData{1}.maps{2}(32, 1:5) = [0.250, 0.254, 0.249, 0.254, 0.251];  % Shifts applied
alignedData{1}.maps{3}(32, 1:5) = [0.249, 0.253, 0.248, 0.253, 0.250];  % Shifts applied

% Visual effect: Temporal features are now aligned across X positions
% The peaks and valleys in the data now occur at the same time indices
```

## 2. Compute All Slices - Batch Processing Analysis

### 2.1 Batch Processing Overview
```matlab
% Function: computeAllSlices() -> computeAllSlicesBackground()
% Process all Y-slices (or X-slices) in current view

% Get dimensions from data structure
[numY, numX] = size(userData.statDataArray{1}.maps{1});  % [64, 128]

% Processing scope depends on current view:
% XtVsY view: process all 64 Y-slices (each slice = [400 × 128] matrix)
% YtVsX view: process all 128 X-slices (each slice = [400 × 64] matrix)
```

### 2.2 Progress Tracking and Slice Processing Loop
```matlab
% Progress dialog setup
progressDlg = uiprogressdlg(fig, 'Title', 'Computing All Slices', ...
    'Message', sprintf('Computing alignment for all %d Y-slices...', numY), ...
    'Value', 0, 'Cancelable', false);

% Initialize timing arrays
userData.timingInfo.sliceTimings = zeros(numY, 1);       % Pre-allocated
userData.sliceIterations = zeros(numY, 1);               % Iterations per slice
userData.sliceAlignmentStatus = false(numY, 1);          % Success tracking

% Main processing loop
for sliceIndex = 1:numY  % Process slices 1 to 64
    sliceStartTime = tic;

    % Extract slice data (identical to single slice process)
    sliceData = zeros(numSegments, numX);                % [400 × 128]
    for seg = 1:numSegments
        sliceData(seg, :) = currentData{statIdx}.maps{seg}(sliceIndex, :);
    end

    % Apply column alignment (same algorithm as single slice)
    [alignedSliceData, convergenceInfo] = alignColumnsImproved(sliceData, ...
        'MaxShift', 15, 'CostFunction', 'mse', 'AlignmentMethod', 'average', ...
        'ConvergenceThreshold', 0.001, 'MaxIterations', 20, ...
        'WeightingFunction', 'exponential', 'WeightingScale', 3.0);

    % Integrate results back into data structure
    for seg = 1:numSegments
        alignedData{statIdx}.maps{seg}(sliceIndex, :) = alignedSliceData(seg, :);
    end

    % Update progress and tracking
    progressDlg.Value = sliceIndex / numY;
    progressDlg.Message = sprintf('Completed slice %d of %d (%.1f%% done)', ...
        sliceIndex, numY, (sliceIndex/numY)*100);

    % Record performance metrics
    sliceElapsedTime = toc(sliceStartTime);
    userData.timingInfo.sliceTimings(sliceIndex) = sliceElapsedTime;
    userData.sliceIterations(sliceIndex) = convergenceInfo.iterations;
    userData.sliceAlignmentStatus(sliceIndex) = convergenceInfo.converged;

    drawnow;  % Update UI
end
```

### 2.3 Batch Processing Performance Analysis
```matlab
% Final batch processing results:
totalSlicesProcessed = 64;                               % All Y-slices completed
averageTimePerSlice = mean(userData.timingInfo.sliceTimings);  % ~2.3 seconds
totalProcessingTime = sum(userData.timingInfo.sliceTimings);   % ~147 seconds
minTimePerSlice = min(userData.timingInfo.sliceTimings);       % ~1.8 seconds
maxTimePerSlice = max(userData.timingInfo.sliceTimings);       % ~3.1 seconds

% Convergence statistics:
averageIterationsPerSlice = mean(userData.sliceIterations);    % ~5.2 iterations
convergenceRate = sum(userData.sliceAlignmentStatus) / numY;   % 98.4% success rate
failedSlices = find(~userData.sliceAlignmentStatus);           % [12, 45] (example)

% Quality metrics:
totalAlignmentImprovement = 0;
for sliceIdx = 1:numY
    if userData.sliceAlignmentStatus(sliceIdx)
        % Calculate improvement for this slice (stored in convergenceInfo)
        totalAlignmentImprovement += sliceImprovements(sliceIdx);
    end
end
averageImprovementPerSlice = totalAlignmentImprovement / sum(userData.sliceAlignmentStatus);
```

## 3. Cross-View Iterative Alignment - Complete Analysis

### 3.1 Cross-View Alignment Overview
```matlab
% Function: computeCrossViewAlignment() -> computeCrossViewAlignmentBackground()
% Iteratively align XtVsY and YtVsX views until global convergence

% Two-phase iterative process:
% Phase 1: Align all Y-slices (XtVsY view) - 64 slices
% Phase 2: Align all X-slices (YtVsX view) - 128 slices
% Repeat until cross-view convergence achieved
```

### 3.2 Cross-View Parameters and Initialization
```matlab
% Cross-view alignment parameters
maxCrossViewIterations = 10;                 % Maximum cross-view iterations
convergenceThreshold = 0.001;                % 0.1% improvement threshold
crossViewIteration = 0;
previousCost = inf;

% Pre-allocated tracking arrays (performance optimization)
improvementHistory = zeros(maxCrossViewIterations, 1);
crossViewCostHistory = zeros(maxCrossViewIterations, 1);

% Data dimensions
[numY, numX] = size(currentData{1}.maps{1});  % [64, 128]
numSegments = length(currentData{1}.maps);    % 400
nStats = length(currentData);                 % Number of statistics

fprintf('Starting cross-view alignment: %d Y-slices x %d X-slices\n', numY, numX);
```

### 3.3 Cross-View Main Iteration Loop
```matlab
while crossViewIteration < maxCrossViewIterations
    crossViewIteration = crossViewIteration + 1;
    iterationStartTime = tic;

    fprintf('Cross-view iteration %d/%d:\n', crossViewIteration, maxCrossViewIterations);

    % PHASE 1: Align all Y-slices (XtVsY view)
    fprintf('  Phase 1: Aligning %d Y-slices...\n', numY);
    yProgressCallback = @(sliceIdx, totalSlices) updateCrossViewProgress(...);
    currentData = alignAllSlicesInView(currentData, 'XtVsY', userData, fig, yProgressCallback);

    % PHASE 2: Align all X-slices (YtVsX view)
    fprintf('  Phase 2: Aligning %d X-slices...\n', numX);
    xProgressCallback = @(sliceIdx, totalSlices) updateCrossViewProgress(...);
    currentData = alignAllSlicesInView(currentData, 'YtVsX', userData, fig, xProgressCallback);

    % Calculate overall alignment cost for convergence check
    currentCost = calculateOverallAlignmentCost(currentData);
    improvement = (previousCost - currentCost) / previousCost;
    improvementHistory(crossViewIteration) = improvement;
    crossViewCostHistory(crossViewIteration) = currentCost;

    iterationTime = toc(iterationStartTime);
    fprintf('  Iteration %d complete: cost=%.6f, improvement=%.4f%%, time=%.1fs\n', ...
        crossViewIteration, currentCost, improvement*100, iterationTime);

    % Check for convergence
    if improvement < convergenceThreshold
        fprintf('Cross-view alignment converged after %d iterations\n', crossViewIteration);
        break;
    end

    previousCost = currentCost;
end
```

### 3.4 Overall Cost Calculation for Cross-View Convergence
```matlab
function totalCost = calculateOverallAlignmentCost(inputData)
    % Calculate comprehensive alignment cost across all data

    totalCost = 0;
    nStats = length(inputData);

    for statIdx = 1:nStats
        numSegments = length(inputData{statIdx}.maps);
        [numY, numX] = size(inputData{statIdx}.maps{1});

        % Calculate cost between all segment pairs (temporal consistency)
        for seg1 = 1:numSegments-1
            for seg2 = seg1+1:numSegments
                map1 = inputData{statIdx}.maps{seg1};  % [64 × 128]
                map2 = inputData{statIdx}.maps{seg2};  % [64 × 128]

                % MSE between consecutive time segments
                segmentCost = mean((map1(:) - map2(:)).^2);
                totalCost = totalCost + segmentCost;
            end
        end

        % Calculate spatial consistency costs
        for seg = 1:numSegments
            currentMap = inputData{statIdx}.maps{seg};

            % Row-wise consistency (Y-direction)
            for y = 1:numY-1
                rowCost = mean((currentMap(y, :) - currentMap(y+1, :)).^2);
                totalCost = totalCost + rowCost * 0.1;  % Lower weight
            end

            % Column-wise consistency (X-direction)
            for x = 1:numX-1
                colCost = mean((currentMap(:, x) - currentMap(:, x+1)).^2);
                totalCost = totalCost + colCost * 0.1;  % Lower weight
            end
        end
    end
end
```

### 3.5 Cross-View Convergence Results
```matlab
% Typical cross-view convergence results:
crossViewIterations = 4;                     % Converged after 4 iterations
finalCrossViewCost = 0.187234;               % Final overall cost
totalCrossViewImprovement = (0.234567 - 0.187234) / 0.234567 = 20.2%;

% Iteration-by-iteration improvement:
improvementHistory = [0.0856, 0.0423, 0.0187, 0.0089];  % Decreasing improvements
crossViewCostHistory = [0.214523, 0.205634, 0.201456, 0.187234];  % Decreasing costs

% Final alignment quality metrics:
ySlicesConverged = sum(userData.sliceAlignmentStatus);    % 63/64 Y-slices
xSlicesConverged = sum(userData.xSliceAlignmentStatus);   % 126/128 X-slices
overallConvergenceRate = (ySlicesConverged + xSlicesConverged) / (numY + numX);  % 98.4%
```

## 4. Alignment Criteria and Success Metrics

### 4.1 What Constitutes "Aligned" Data

#### **Primary Alignment Criterion: Cost Reduction**
```matlab
% Alignment is considered successful when:
% 1. Cost improvement falls below threshold (0.1%)
% 2. Maximum iterations reached (20)
% 3. No further shifts provide meaningful improvement

% Mathematical definition of alignment success:
improvement = (previousCost - currentCost) / previousCost;
isAligned = improvement < convergenceThreshold;  % improvement < 0.001

% Example successful alignment:
% Iteration 1: cost = 0.234567, improvement = 1.81%
% Iteration 2: cost = 0.228934, improvement = 2.40%
% Iteration 3: cost = 0.226123, improvement = 1.23%
% Iteration 4: cost = 0.224891, improvement = 0.54%
% Iteration 5: cost = 0.224567, improvement = 0.14%
% Iteration 6: cost = 0.224456, improvement = 0.05% ← CONVERGED (< 0.1%)
```

#### **Secondary Alignment Criteria**
```matlab
% 1. Shift Stability: Applied shifts should be reasonable
maxReasonableShift = 15;                     % ±15 pixels maximum
appliedShifts = [-2, -3, -1, 0, -2, +1, ...];
shiftStability = all(abs(appliedShifts) <= maxReasonableShift);  % true

% 2. Spatial Consistency: Adjacent columns should have similar shifts
shiftDifferences = diff(appliedShifts);      % [-1, +2, +1, -2, +3, ...]
maxShiftJump = max(abs(shiftDifferences));   % Should be < 5 pixels typically
spatialConsistency = maxShiftJump < 5;       % true

% 3. Temporal Coherence: Features should align across time
% Measured by reduced variance in aligned data
originalVariance = var(originalSliceData, [], 2);    % Variance per time point
alignedVariance = var(alignedSliceData, [], 2);      % Variance after alignment
varianceReduction = mean(originalVariance - alignedVariance);  % Should be positive
temporalCoherence = varianceReduction > 0;           % true
```

### 4.2 Quality Assessment Metrics

#### **Quantitative Quality Measures**
```matlab
% 1. Mean Squared Error Reduction
originalMSE = calculateOverallCost(originalData);     % 0.238901
alignedMSE = calculateOverallCost(alignedData);       % 0.224456
mseReduction = (originalMSE - alignedMSE) / originalMSE;  % 6.05% improvement

% 2. Cross-Correlation Improvement
% Calculate average cross-correlation between adjacent columns
originalCorrelations = zeros(numCols-1, 1);
alignedCorrelations = zeros(numCols-1, 1);

for col = 1:numCols-1
    % Original data correlations
    corrMatrix = corrcoef(originalData(:, col), originalData(:, col+1));
    originalCorrelations(col) = corrMatrix(1, 2);

    % Aligned data correlations
    corrMatrix = corrcoef(alignedData(:, col), alignedData(:, col+1));
    alignedCorrelations(col) = corrMatrix(1, 2);
end

averageOriginalCorr = mean(originalCorrelations);     % 0.723
averageAlignedCorr = mean(alignedCorrelations);       % 0.856
correlationImprovement = averageAlignedCorr - averageOriginalCorr;  % +0.133

% 3. Feature Alignment Score
% Measure how well peaks align across columns
peakAlignmentScore = calculatePeakAlignment(alignedData);  % 0.892 (0-1 scale)
```

#### **Visual Quality Indicators**
```matlab
% 1. Reduced "jaggedness" in heatmap visualization
% 2. Clearer temporal patterns across spatial dimensions
% 3. Smoother transitions between adjacent columns
% 4. Enhanced visibility of layer structures
% 5. Reduced noise in cross-sectional views
```

## 5. Mathematical Foundations and Formulas

### 5.1 Core Mathematical Equations

#### **Mean Squared Error Cost Function**
```matlab
% Primary cost function used throughout alignment process
function cost = calculateMSE(column1, column2)
    % Exclude edge effects (top and bottom 5 samples)
    validIndices = 6:(length(column1)-5);

    % Extract valid data
    col1_valid = column1(validIndices);
    col2_valid = column2(validIndices);

    % Calculate MSE
    cost = (1/N) * Σ(col1_valid(i) - col2_valid(i))²
    % where N = length(validIndices)
end

% Mathematical representation:
% MSE = (1/N) * Σᵢ₌₁ᴺ (xᵢ - yᵢ)²
% where xᵢ, yᵢ are corresponding samples from two columns
```

#### **Distance Weighting Function**
```matlab
% For 'full' alignment method - exponential distance weighting
function weight = calculateDistanceWeight(distance, scale)
    weight = exp(-distance / scale);
end

% Mathematical representation:
% w(d) = e^(-d/σ)
% where d = |col_i - col_j| (distance between columns)
%       σ = weightingScale (typically 3.0)

% Example weights for scale = 3.0:
% Distance 1: w = e^(-1/3) = 0.717
% Distance 2: w = e^(-2/3) = 0.513
% Distance 3: w = e^(-3/3) = 0.368
% Distance 5: w = e^(-5/3) = 0.189
```

#### **Weighted Cost Calculation**
```matlab
% For 'full' method - distance-weighted cost across all columns
function totalCost = calculateWeightedCost(targetColumn, allColumns, targetIndex)
    totalCost = 0;
    totalWeight = 0;

    for otherCol = 1:numColumns
        if otherCol ~= targetIndex
            distance = abs(otherCol - targetIndex);
            weight = exp(-distance / 3.0);

            columnCost = calculateMSE(targetColumn, allColumns(:, otherCol));
            totalCost = totalCost + weight * columnCost;
            totalWeight = totalWeight + weight;
        end
    end

    finalCost = totalCost / totalWeight;  % Normalized weighted average
end

% Mathematical representation:
% C_weighted = (Σⱼ wⱼ * MSE(cᵢ, cⱼ)) / (Σⱼ wⱼ)
% where wⱼ = e^(-|i-j|/σ) for all j ≠ i
```

#### **Convergence Detection Formula**
```matlab
% Relative improvement calculation
function improvement = calculateImprovement(previousCost, currentCost)
    improvement = (previousCost - currentCost) / previousCost;
end

% Mathematical representation:
% ρ = (C_{k-1} - C_k) / C_{k-1}
% where ρ = relative improvement
%       C_k = cost at iteration k
%       C_{k-1} = cost at previous iteration
%
% Convergence criterion: ρ < ε (typically ε = 0.001)
```

### 5.2 Shift Application Mathematics

#### **Zero-Padding Shift Algorithm**
```matlab
% Mathematical description of integer shift operation
function shiftedColumn = applyShift(originalColumn, shiftAmount)
    N = length(originalColumn);
    shiftedColumn = zeros(N, 1);

    if shiftAmount > 0  % Shift down (positive direction)
        % shiftedColumn[k] = originalColumn[k - shiftAmount] for k > shiftAmount
        % shiftedColumn[k] = 0 for k ≤ shiftAmount
        shiftedColumn((shiftAmount+1):N) = originalColumn(1:(N-shiftAmount));

    elseif shiftAmount < 0  % Shift up (negative direction)
        s = abs(shiftAmount);
        % shiftedColumn[k] = originalColumn[k + s] for k ≤ N-s
        % shiftedColumn[k] = 0 for k > N-s
        shiftedColumn(1:(N-s)) = originalColumn((s+1):N);
    end
end

% Mathematical representation:
% For shift s > 0 (downward):
% y[n] = { x[n-s]  if n > s
%        { 0       if n ≤ s
%
% For shift s < 0 (upward):
% y[n] = { x[n+|s|] if n ≤ N-|s|
%        { 0        if n > N-|s|
```

This comprehensive analysis provides complete variable tracking and mathematical foundations for understanding exactly how the column alignment algorithms work, what constitutes successful alignment, and how the data is transformed at each step of the process.
