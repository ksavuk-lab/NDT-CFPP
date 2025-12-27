# AlignWaveformsAtFirstPeak.m - Complete Variable Tracking and Data Analysis

## Overview
This document provides an extremely detailed, step-by-step analysis of the AlignWaveformsAtFirstPeak.m function, with complete variable tracking, mathematical computations, and data transformations. Every calculation is shown with actual values and examples.

## 1. Function Entry and Input Data Structure

### 1.1 Function Signature and Parameters
```matlab
function [alignedWaveform, shiftIndices] = AlignWaveformsAtFirstPeak(waveform3DMatrix, t, isEnvelope, timeRange, peakType, enableDiagnostics, numDiagnosticSamples)
```

#### **Input Data Structure**
```matlab
% Primary input: 3D waveform matrix
waveform3DMatrix: [numY × numX × numT] = [64 × 128 × 2000]
% - numY = 64: Number of Y spatial positions (rows)
% - numX = 128: Number of X spatial positions (columns)  
% - numT = 2000: Number of time samples per waveform

% Time vector
t: [1 × 2000] = [0, 5e-9, 1e-8, 1.5e-8, ..., 9.995e-6]
% - Sample interval: dt = 5e-9 seconds = 5 nanoseconds
% - Total time span: 0 to 9.995 microseconds
% - Sampling frequency: 200 MHz

% Sample data values for position (32, 64):
waveform3DMatrix(32, 64, 1:10) = [0.012, 0.015, 0.018, 0.021, 0.019, 0.017, 0.020, 0.023, 0.025, 0.022]
waveform3DMatrix(32, 64, 996:1000) = [0.089, 0.156, 0.234, 0.187, 0.123]  % Peak region
waveform3DMatrix(32, 64, 1996:2000) = [0.008, 0.006, 0.004, 0.003, 0.002]  % End region
```

#### **Parameter Processing and Defaults**
```matlab
% Force envelope flag to false (envelope processing removed)
isEnvelope = false;  % Always process as raw waveform data

% Input validation
[numY, numX, numT] = size(waveform3DMatrix);  % [64, 128, 2000]

% Fixed threshold parameter
threshold = 1.0;  % Amplitude threshold multiplier

% Time range processing
if nargin < 4 || isempty(timeRange)
    timeRangeIndices = [1, length(t)];  % [1, 2000] - search entire waveform
    searchRangeMsg = 'searching entire waveform';
else
    % Convert time range to indices
    timeRangeIndices = [find(t >= timeRange(1), 1, 'first'), 
                       find(t <= timeRange(2), 1, 'last')];
    searchRangeMsg = sprintf('searching between %.3f and %.3f μs', 
                            timeRange(1)*1e6, timeRange(2)*1e6);
end

% Peak type parameter
if nargin < 5 || isempty(peakType)
    peakType = 'positive';  % Default to positive peak detection
end
% Options: 'positive', 'negative', 'both'

% Diagnostic parameters
if nargin < 6 || isempty(enableDiagnostics)
    enableDiagnostics = false;
end
if nargin < 7 || isempty(numDiagnosticSamples)
    numDiagnosticSamples = 100;
end
```

### 1.2 Output Initialization and Performance Tracking
```matlab
% Pre-allocate output arrays for performance
alignedWaveform = zeros(size(waveform3DMatrix));  % [64 × 128 × 2000]
shiftIndices = zeros(numY, numX);                 % [64 × 128]

% Performance tracking variables
startTime = tic;                                  % Start timing
processedCount = 0;                               % Waveforms processed counter
totalWaveforms = numY * numX;                     % 64 × 128 = 8192 total waveforms

% Diagnostic setup
if enableDiagnostics
    diagnosticSamples = min(numDiagnosticSamples, totalWaveforms);
    diagnosticIndices = randperm(totalWaveforms, diagnosticSamples);
    % Example: diagnosticIndices = [1247, 3891, 5632, 7234, ...]
end
```

## 2. Main Processing Loop - Waveform-by-Waveform Analysis

### 2.1 Nested Loop Structure
```matlab
% Process each spatial position
for y = 1:numY  % Y positions: 1 to 64
    for x = 1:numX  % X positions: 1 to 128
        processedCount = processedCount + 1;
        
        % Extract current waveform
        currentWaveform = squeeze(waveform3DMatrix(y, x, :));  % [2000 × 1]
        
        % Progress calculation
        progress = processedCount / totalWaveforms * 100;
        
        % Example for first waveform (y=1, x=1):
        % processedCount = 1
        % progress = 1/8192 * 100 = 0.012%
```

### 2.2 Current Waveform Extraction Example
```matlab
% Example: Processing position (32, 64)
y = 32; x = 64;
currentWaveform = squeeze(waveform3DMatrix(32, 64, :));  % [2000 × 1]

% Sample waveform data:
currentWaveform(1:10) = [0.012, 0.015, 0.018, 0.021, 0.019, 0.017, 0.020, 0.023, 0.025, 0.022];
currentWaveform(995:1005) = [0.067, 0.089, 0.156, 0.234, 0.187, 0.123, 0.089, 0.067, 0.045, 0.034];
currentWaveform(1996:2000) = [0.008, 0.006, 0.004, 0.003, 0.002];

% Key characteristics:
% - Early time: Low amplitude noise (~0.01-0.03)
% - Peak region: Strong signal around index 998 (amplitude 0.234)
% - Late time: Decay to near zero
```

### 2.3 NaN Validation Check
```matlab
% Check for invalid data
if any(isnan(currentWaveform))
    % Skip alignment - preserve original waveform
    alignedWaveform(y, x, :) = currentWaveform;
    shiftIndices(y, x) = 0;  % No shift applied
    continue;  % Move to next waveform
end

% For our example:
hasNaN = any(isnan(currentWaveform));  % false - proceed with alignment
```

## 3. Peak Detection Algorithm - Raw Waveform Processing

### 3.1 Range Extraction
```matlab
% Extract search range from waveform
rangeWaveform = currentWaveform(timeRangeIndices(1):timeRangeIndices(2));

% For entire waveform search (most common case):
rangeStart = timeRangeIndices(1);  % 1
rangeEnd = timeRangeIndices(2);    % 2000
rangeWaveform = currentWaveform(1:2000);  % Full waveform
numPoints = length(rangeWaveform);  % 2000

% Sample range data (same as currentWaveform for full search):
rangeWaveform(1:10) = [0.012, 0.015, 0.018, 0.021, 0.019, 0.017, 0.020, 0.023, 0.025, 0.022];
rangeWaveform(995:1005) = [0.067, 0.089, 0.156, 0.234, 0.187, 0.123, 0.089, 0.067, 0.045, 0.034];
```

### 3.2 Raw Waveform Peak Detection Function Call
```matlab
% Call peak detection algorithm
peakIdx = findFirstPeakInRawWaveform(currentWaveform, threshold, timeRangeIndices, peakType);

% Input parameters:
% - waveform: currentWaveform [2000 × 1]
% - threshold: 1.0
% - rangeIndices: [1, 2000]
% - peakType: 'positive'
```

## 4. Positive Peak Detection Algorithm (Detailed)

### 4.1 Threshold Calculation
```matlab
% Calculate adaptive threshold based on signal amplitude
adjustedThreshold = threshold × max(rangeWaveform);

% Sample calculation:
maxValue = max(rangeWaveform);  % 0.234 (at index 998)
adjustedThreshold = 1.0 × 0.234;  % 0.234

% This means we're looking for the first point that reaches the maximum amplitude
```

### 4.2 Threshold Crossing Search
```matlab
% Find first point above threshold
thresholdCrossings = find(rangeWaveform >= adjustedThreshold, 1, 'first');

% Sample calculation:
% rangeWaveform >= 0.234 returns logical array
% Only rangeWaveform(998) = 0.234 meets this criterion
% Result: thresholdCrossings = 998
```

### 4.3 Peak Verification Algorithm
```matlab
% Search for actual peak after threshold crossing
if ~isempty(thresholdCrossings)
    startIdx = thresholdCrossings;  % 998
    verificationWindow = min(5, numPoints - startIdx);  % min(5, 2000-998) = 5
    
    % Peak search loop
    for i = startIdx:(numPoints-1)  % 998 to 1999
        % Check if signal is decreasing (potential peak)
        if rangeWaveform(i) > rangeWaveform(i+1)
            % Verify peak stability over verification window
            windowEnd = min(i + verificationWindow, numPoints);
            windowValues = rangeWaveform((i+1):windowEnd);
            
            % Stability condition: current value >= all following values in window
            if all(rangeWaveform(i) >= windowValues)
                peakIdx = i + rangeStart - 1;  % Convert to global index
                break;  % Peak found and verified
            end
        end
    end
end
```

### 4.4 Peak Verification Example
```matlab
% Sample verification for i = 998:
currentValue = rangeWaveform(998);  % 0.234
nextValue = rangeWaveform(999);     % 0.187

% Decreasing signal check:
isDecreasing = currentValue > nextValue;  % 0.234 > 0.187 = true ✓

% Verification window:
windowEnd = min(998 + 5, 2000);  % min(1003, 2000) = 1003
windowValues = rangeWaveform(999:1003);  % [0.187, 0.123, 0.089, 0.067, 0.045]

% Stability check:
stabilityCheck = all(0.234 >= [0.187, 0.123, 0.089, 0.067, 0.045]);  % true ✓

% Result:
peakIdx = 998 + 1 - 1;  % 998 (global index)
peakTime = t(998);      % 4.985e-6 seconds = 4.985 μs
peakValue = currentWaveform(998);  % 0.234
```

## 5. Alternative Peak Detection Scenarios

### 5.1 Negative Peak Detection (Valley Detection)
```matlab
% For peakType = 'negative'
adjustedThreshold = threshold × abs(min(rangeWaveform));

% Sample calculation:
minValue = min(rangeWaveform);  % -0.156 (hypothetical valley at index 1205)
adjustedThreshold = 1.0 × abs(-0.156);  % 0.156

% Threshold crossing search for valleys:
thresholdCrossings = find(rangeWaveform <= -adjustedThreshold, 1, 'first');

% Valley verification (similar to peak but looking for signal increase):
for i = startIdx:(numPoints-1)
    if rangeWaveform(i) < rangeWaveform(i+1)  % Signal increasing (end of valley)
        windowValues = rangeWaveform((i+1):windowEnd);
        if all(rangeWaveform(i) <= windowValues)  % Valley stability
            peakIdx = i + rangeStart - 1;
            break;
        end
    end
end
```

### 5.2 Both Peak Types Detection
```matlab
% For peakType = 'both' - find both positive and negative peaks
posPeakIdx = findFirstPeakInRawWaveform(currentWaveform, threshold, timeRangeIndices, 'positive');
negPeakIdx = findFirstPeakInRawWaveform(currentWaveform, threshold, timeRangeIndices, 'negative');

% Magnitude comparison:
if ~isempty(posPeakIdx) && ~isempty(negPeakIdx)
    posValue = abs(currentWaveform(posPeakIdx));  % abs(0.234) = 0.234
    negValue = abs(currentWaveform(negPeakIdx));  % abs(-0.156) = 0.156

    % Selection criterion: choose larger absolute magnitude
    if posValue >= negValue  % 0.234 >= 0.156 ✓
        peakIdx = posPeakIdx;  % 998
        selectedType = 'positive';
    else
        peakIdx = negPeakIdx;
        selectedType = 'negative';
    end
elseif ~isempty(posPeakIdx)
    peakIdx = posPeakIdx;
    selectedType = 'positive';
elseif ~isempty(negPeakIdx)
    peakIdx = negPeakIdx;
    selectedType = 'negative';
else
    peakIdx = [];  % No peak found
end
```

### 5.3 Fallback to Global Maximum/Minimum
```matlab
% If no threshold crossing found, use global extremum
if isempty(thresholdCrossings)
    if strcmp(peakType, 'positive')
        [~, localPeakIdx] = max(rangeWaveform);
        peakIdx = localPeakIdx + rangeStart - 1;
    elseif strcmp(peakType, 'negative')
        [~, localPeakIdx] = min(rangeWaveform);
        peakIdx = localPeakIdx + rangeStart - 1;
    end
end

% Mathematical guarantee: max/min always exists for non-empty arrays
% Example: localPeakIdx = 998, peakIdx = 998 + 1 - 1 = 998
```

## 6. Peak Index Validation and Alignment Execution

### 6.1 Peak Index Validation
```matlab
% Validation conditions for successful peak detection:
isValidPeak = ~isempty(peakIdx) && peakIdx > 0 && peakIdx <= numT;

% Current example validation:
condition1 = ~isempty(peakIdx);  % ~isempty(998) = true ✓
condition2 = peakIdx > 0;        % 998 > 0 = true ✓
condition3 = peakIdx <= numT;    % 998 <= 2000 = true ✓

% Overall result:
isValidPeak = true && true && true;  % true ✓

% Additional information:
peakTime = t(peakIdx);           % t(998) = 4.985e-6 seconds = 4.985 μs
peakValue = currentWaveform(peakIdx);  % currentWaveform(998) = 0.234
```

### 6.2 Waveform Alignment Execution
```matlab
if isValidPeak
    % Record peak location
    shiftIndices(y, x) = peakIdx;  % shiftIndices(32, 64) = 998

    % Create aligned waveform
    shiftedWaveform = zeros(size(currentWaveform));  % [2000 × 1] zeros

    % Calculate copy parameters
    copyLength = numT - peakIdx + 1;  % 2000 - 998 + 1 = 1003

    % Copy data from peak onwards to beginning
    shiftedWaveform(1:copyLength) = currentWaveform(peakIdx:end);
    % shiftedWaveform(1:1003) = currentWaveform(998:2000)

    % Store aligned waveform
    alignedWaveform(y, x, :) = shiftedWaveform;
end
```

### 6.3 Alignment Transformation Details
```matlab
% Before alignment:
% currentWaveform(995:1005) = [0.067, 0.089, 0.156, 0.234, 0.187, 0.123, 0.089, 0.067, 0.045, 0.034]
%                              idx:995  996   997   998   999   1000  1001  1002  1003  1004
%                                                   ↑ Peak at index 998

% After alignment:
% shiftedWaveform(1:11) = [0.234, 0.187, 0.123, 0.089, 0.067, 0.045, 0.034, ..., 0.002, 0, 0]
%                          idx:1    2     3     4     5     6     7          1003  1004 1005
%                              ↑ Peak now at index 1 (t = 0)

% Mathematical transformation:
% shiftedWaveform(i) = currentWaveform(peakIdx + i - 1) for i = 1 to copyLength
% shiftedWaveform(i) = 0 for i = copyLength+1 to numT

% Specific examples:
shiftedWaveform(1) = currentWaveform(998);   % 0.234 (peak at t=0)
shiftedWaveform(2) = currentWaveform(999);   % 0.187
shiftedWaveform(3) = currentWaveform(1000);  % 0.123
...
shiftedWaveform(1003) = currentWaveform(2000);  % 0.002 (last data point)
shiftedWaveform(1004) = 0;  % Zero padding begins
shiftedWaveform(1005) = 0;  % Zero padding continues
...
shiftedWaveform(2000) = 0;  % Zero padding to end
```

## 7. Failure Handling and Edge Cases

### 7.1 Alignment Failure Scenarios
```matlab
% If peak validation fails:
if ~isValidPeak
    % Preserve original waveform
    alignedWaveform(y, x, :) = currentWaveform;
    shiftIndices(y, x) = 0;  % No shift recorded

    % Possible failure reasons:
    % 1. No peak found in specified range
    % 2. Peak index out of bounds (peakIdx <= 0 or peakIdx > numT)
    % 3. All values below threshold
    % 4. Peak verification failed (no stable peak found)
    % 5. Empty waveform or all NaN values
end
```

### 7.2 Diagnostic Logging (Optional)
```matlab
% If diagnostics enabled and current waveform selected for sampling:
if enableDiagnostics
    linearIndex = (y-1)*numX + x;  % Convert 2D position to linear index
    % Example: linearIndex = (32-1)*128 + 64 = 31*128 + 64 = 4032

    if ismember(linearIndex, diagnosticIndices)
        % Extract context around peak
        contextStart = max(1, peakIdx - 5);      % max(1, 998-5) = 993
        contextEnd = min(numT, peakIdx + 5);     % min(2000, 998+5) = 1003
        contextValues = currentWaveform(contextStart:contextEnd);

        % Log detailed information
        fprintf('Position [%d, %d]: Peak at index %d (%.3f μs), value %.6f\n', ...
                y, x, peakIdx, t(peakIdx)*1e6, currentWaveform(peakIdx));
        fprintf('  Context values: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', ...
                contextValues);

        % Example output:
        % Position [32, 64]: Peak at index 998 (4.985 μs), value 0.234000
        % Context values: [0.078, 0.089, 0.123, 0.156, 0.189, 0.234, 0.187, 0.123, 0.089, 0.067, 0.045]
    end
end
```

## 8. Processing Completion and Results

### 8.1 Final Processing Statistics
```matlab
% After processing all 8192 waveforms:
totalProcessingTime = toc(startTime);  % e.g., 12.34 seconds
successfulAlignments = sum(shiftIndices(:) > 0);  % e.g., 7834 successful alignments
failedAlignments = totalWaveforms - successfulAlignments;  % e.g., 358 failures
successRate = successfulAlignments / totalWaveforms * 100;  % 95.6%

% Performance metrics:
averageTimePerWaveform = totalProcessingTime / totalWaveforms;  % 1.5 milliseconds
waveformsPerSecond = totalWaveforms / totalProcessingTime;      % 664 waveforms/second

% Peak timing statistics:
validShifts = shiftIndices(shiftIndices > 0);
averagePeakTime = mean(t(validShifts));     % 4.2 μs (average peak arrival time)
stdPeakTime = std(t(validShifts));          % 1.8 μs (standard deviation)
minPeakTime = min(t(validShifts));          % 1.2 μs (earliest peak)
maxPeakTime = max(t(validShifts));          % 8.7 μs (latest peak)
```

### 8.2 Output Data Structure
```matlab
% Function returns:
% 1. alignedWaveform: [64 × 128 × 2000] - aligned waveforms
% 2. shiftIndices: [64 × 128] - peak locations for each position

% alignedWaveform characteristics:
% - All detected peaks aligned to t = 0 (first time sample)
% - Original amplitude characteristics preserved
% - Zero-padding applied after original data ends
% - Failed alignments preserve original waveform

% shiftIndices characteristics:
% - Contains peak index for each spatial position
% - Zero values indicate failed alignment
% - Can be used to calculate peak arrival times: t(shiftIndices)
% - Provides spatial map of peak timing variations
```

## 9. Alignment Criteria and Success Metrics

### 9.1 What Constitutes "Aligned" Data
```matlab
% Alignment success criteria:
% 1. Valid peak detected within search range
% 2. Peak passes verification stability test
% 3. Peak index within valid bounds [1, numT]
% 4. Waveform contains no NaN values

% Mathematical definition of successful alignment:
isAligned = ~isempty(peakIdx) && peakIdx > 0 && peakIdx <= numT && ...
            all(~isnan(currentWaveform));

% Quality indicators:
% - Peak at t = 0 in aligned waveform
% - Consistent temporal reference across all spatial positions
% - Preserved signal characteristics after peak
% - Appropriate zero-padding for data beyond original length
```

### 9.2 Alignment Quality Assessment
```matlab
% Quantitative quality measures:
% 1. Success rate: percentage of waveforms successfully aligned
successRate = sum(shiftIndices(:) > 0) / numel(shiftIndices) * 100;

% 2. Peak timing consistency: standard deviation of peak times
peakTimes = t(shiftIndices(shiftIndices > 0));
timingConsistency = std(peakTimes);  % Lower values indicate better consistency

% 3. Spatial coherence: similarity of peak times in neighboring positions
spatialCoherence = calculateSpatialCoherence(shiftIndices, t);

% 4. Signal preservation: correlation between original and aligned signals
signalPreservation = calculateSignalPreservation(waveform3DMatrix, alignedWaveform);
```

This comprehensive analysis provides complete variable tracking and mathematical foundations for understanding exactly how the AlignWaveformsAtFirstPeak.m function works, what constitutes successful alignment, and how the data is transformed at each step of the process.
