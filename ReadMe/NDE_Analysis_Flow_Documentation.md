# NDE Analysis Data Flow Documentation
## Complete Variable and Function Tracking from Main.m to XtVsYPlot.m

### Overview
This document provides an extremely detailed tracking of all variables, functions, equations, and mathematical operations in the NDE analysis pipeline from Main.m through XtVsYPlot.m.

---

## Chunk 1: Main.m Initialization

### Path Setup
```matlab
addpath(fullfile(pwd, 'Utils'))
addpath(fullfile(pwd, 'Plotting Scripts'))
addpath(fullfile(pwd, 'Plotting Scripts', 'Data'))
addpath(fullfile(pwd, 'Plotting Scripts', 'UI'))
addpath(fullfile(pwd, 'Plotting Scripts', 'Utils'))
addpath(fullfile(pwd, 'Plotting Scripts', 'Visualization'))
addpath(fullfile(pwd, 'Plotting Scripts', 'Visualization3D'))
addpath(fullfile(pwd, 'Data and Computation Scripts'))
addUtilsPath()
```

### User Configuration Variables
| Variable | Type | Value | Purpose |
|----------|------|-------|---------|
| `DATASET` | integer | 1 | Dataset selection (1-5) |
| `TrimWaveData` | logical | 1 | Enable spatial trimming |
| `TrimTimeRange` | logical | 1 | Enable temporal trimming |
| `xRange` | 1×2 double | [20 60] | X spatial range (mm) |
| `yRange` | 1×2 double | [5 10] | Y spatial range (mm) |
| `TimeRange` | 1×2 double | [3.0e-6 7.0e-6] | Time range (seconds) |
| `sampleThickness` | double | 10 | Z-height for 3D (mm) |
| `TotalSlices` | integer | 400 | Number of time segments |

### Dataset Options
1. `L0P5S25_10Mhz` - No Defect, Sample 25
2. `L8P12S5_10Mhz` - 8 Degree Defect, Sample 5
3. `L16P1S9_10MHz` - 16 Degree Defect, Sample 9
4. `L16P1S10_10MHz` - 16 Degree Defect, Sample 10
5. `L16P34S3_10MHz` - 16 Degree Defect, Sample 3

### Alignment Settings
| Variable | Type | Value | Purpose |
|----------|------|-------|---------|
| `AlignWaveFormsAtFirstPeak` | logical | 0 | Disable waveform alignment |
| `Align_Between_Time_Range` | 1×2 double | [3.0e-6 7.0e-6] | Peak search range |
| `AlignPeakType` | string | 'negative' | Peak type to align |
| `AlignDiagnostics` | logical | 0 | Disable diagnostics |
| `AlignDiagnosticSamples` | integer | 100 | Diagnostic sample count |

### Statistical Analysis Flags
| Variable | Type | Value | Purpose |
|----------|------|-------|---------|
| `ComputeMaxAmplitude` | logical | 1 | Enable MaxAmplitude |
| `ComputeRMS` | logical | 0 | Disable RMS |
| `ComputeVariance` | logical | 0 | Disable Variance |
| `ComputeSkewness` | logical | 0 | Disable Skewness |
| `ComputeKurtosis` | logical | 0 | Disable Kurtosis |

---

## Chunk 2: Data Loading Process

### Function Call: loadData(DATASET)
**File:** `Data and Computation Scripts/loadData.m`
**Input:** `DATASET = 1`
**Output:** `DataStructure`

### Data Loading Algorithm
```matlab
DataStructure = struct();
baseDirs = {
    '/Users/ksavuk/Documents/.../L0/',
    '/Users/kostia/Documents/.../L0/'
};

% Dataset 1 file selection
fileNames{1} = 'L0P5S25_10Mhz_TimingCorrected.mat';

% Load attempt
for i = 1:length(baseDirs)
    fullPath = [baseDirs{i}, fileNames{1}];
    temp = load(fullPath);
    DataStructure = temp.datastructure;
    break; % Success
end
```

### DataStructure Fields Extracted
| Field | Type | Purpose |
|-------|------|---------|
| `sampRate` | double | Sampling rate (MHz) |
| `scanRes` | double | Scan resolution |
| `indexRes` | double | Index resolution |
| `xSize` | integer | X dimension size |
| `ySize` | integer | Y dimension size |
| `tSize` | integer | Time dimension size |
| `waveform3DMatrix` | 3D double | Raw waveform data [ySize×xSize×tSize] |

### Time Vector Generation
```matlab
SampRate = DataStructure.sampRate × 1e6;  % Convert MHz to Hz
dt = 1 / SampRate;                        % Time step (seconds)
tSize = DataStructure.tSize;              % Number of time samples
t_full = dt:dt:(tSize × dt);              % Time vector
```

### Variable Assignments
```matlab
waveform3DMatrix = DataStructure.waveform3DMatrix;
waveformFull = waveform3DMatrix;  % Working copy
```

---

## Chunk 3: Data Extraction and Trimming

### Function Call: ExtractAndTrimData()
**File:** `Data and Computation Scripts/ExtractAndTrimData.m`
**Inputs:** `DataStructure, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange`
**Outputs:** `SampRate, ScanRes, indexRes, xSize, ySize, tSize, waveform3DMatrix, waveformEnvelope3DMatrix, t, X_Coordinates, Y_Coordinates`

### Field Extraction
```matlab
SampRate = DataStructure.sampRate × 1e6;  % MHz to Hz conversion
ScanRes = DataStructure.scanRes;
indexRes = DataStructure.indexRes;
xSize = DataStructure.xSize;
ySize = DataStructure.ySize;
tSize = DataStructure.tSize;
waveform3DMatrix = DataStructure.waveform3DMatrix;
```

### Coordinate Generation
```matlab
X_Coordinates = linspace(based on xSize and ScanRes);  % mm
Y_Coordinates = linspace(based on ySize and ScanRes);  % mm
t = time_vector(based on tSize and SampRate);          % seconds
```

### Spatial Trimming (TrimWaveData = 1)
```matlab
% X Range Trimming
xIndices = find(X_Coordinates >= xRange(1) & X_Coordinates <= xRange(2));
X_Coordinates = X_Coordinates(xIndices);

% Y Range Trimming  
yIndices = find(Y_Coordinates >= yRange(1) & Y_Coordinates <= yRange(2));
Y_Coordinates = Y_Coordinates(yIndices);

% Apply spatial trimming to waveform matrix
waveform3DMatrix = waveform3DMatrix(yIndices, xIndices, :);
```

### Temporal Trimming (TrimTimeRange = 1)
```matlab
% Time Range Trimming
tIndices = find(t >= TimeRange(1) & t <= TimeRange(2));
t = t(tIndices);
waveform3DMatrix = waveform3DMatrix(:, :, tIndices);
```

### Size Updates
```matlab
xSize = length(X_Coordinates);  % New X dimension after trimming
ySize = length(Y_Coordinates);  % New Y dimension after trimming
tSize = length(t);              % New time dimension after trimming
```

### Aspect Ratio Calculation
```matlab
xlength = max(X_Coordinates) - min(X_Coordinates);  % X spatial extent (mm)
ylength = max(Y_Coordinates) - min(Y_Coordinates);  % Y spatial extent (mm)
aspectRatio = [xlength ylength 1];                  % Aspect ratio vector
```

---

## Chunk 4: Waveform Preprocessing and Segmentation

### FileNamingArray Construction
```matlab
alignFlag = 0;  % No alignment applied
caseNumber = DATASET;  % Case identifier = 1
FileNamingArray = [DATASET, caseNumber, TrimWaveData, TrimTimeRange, ...
                   xRange, yRange, TimeRange, alignFlag];
% Result: [1, 1, 1, 1, 20, 60, 5, 10, 3.0e-6, 7.0e-6, 0]
```

### Function Call: PreComputeWaveformData()
**File:** `Data and Computation Scripts/PreComputeWaveformData.m`
**Purpose:** Extract and save waveform data

### Waveform Extraction Algorithm
```matlab
numY = length(Y_Coordinates);
numX = length(X_Coordinates);
totalWaveforms = numY × numX;
waveformArray = zeros(totalWaveforms, size(dataToPlot, 3));

currentWaveform = 0;
for i = 1:numY
    for j = 1:numX
        currentWaveform = currentWaveform + 1;
        waveformArray(currentWaveform, :) = squeeze(dataToPlot(i, j, :))';
    end
end
```

### Data Structure Creation
```matlab
savedData.waveformArray = waveformArray;
savedData.t = t;
savedData.X_Coordinates = X_Coordinates;
savedData.Y_Coordinates = Y_Coordinates;
savedData.numY_sub = numY;
savedData.numX_sub = numX;
savedData.dataType = 'Raw';
savedData.FileNamingArray = FileNamingArray;
```

### File Save Operation
**File:** `Saved Wave Forms/WaveformCase1_Raw_X20.0-60.0_Y5.0-10.0_T3.0e-06-7.0e-06.mat`

### Function Call: SegmentPrecomputedWaveform()
**File:** `Data and Computation Scripts/SegmentPrecomputedWaveform.m`
**Inputs:** `dataFile, 'totalSlices', [], TotalSlices=400`
**Output:** `segmentedDataTotalSlices`

### Segmentation Algorithm
```matlab
M = size(waveformArray, 2);  % Number of time points
dt = t(2) - t(1);            % Time step
numSlices = 400;             % Number of segments

% Create segment boundaries
edges = linspace(1, M+1, numSlices+1);  % 401 elements for 400 segments

% Segmentation loop
for i = 1:numSlices
    time_indices = edges(i):edges(i+1)-1;
    segment_waveform = waveformArray(:, time_indices);
    segmentedData{i}.waveform = reshape(segment_waveform, ySize, xSize, []);
    segmentedData{i}.time = t(time_indices);
end
```

### Result
- **segmentedData:** Cell array with 400 elements
- **Each cell contains:**
  - `.waveform`: [ySize × xSize × segmentTimePoints]
  - `.time`: Time vector for segment

**File Save:** `Saved Wave Forms/SegmentedWaveformCase1_Raw_X20-60_Y5-10_T3.0e-06-7.0e-06_TotalSlices_400.mat`

---

## Chunk 5: Statistical Computations

### Statistics Selection Logic
```matlab
StatsToCompute = {};
if ComputeMaxAmplitude == 1
    StatsToCompute{end+1} = 'MaxAmplitude';  % Added
end
if ComputeRMS == 1
    StatsToCompute{end+1} = 'RMS';          % Skipped (RMS = 0)
end
% ... other statistics skipped
% Result: StatsToCompute = {'MaxAmplitude'}
```

### Function Call: ComputeAndTransformStats()
**File:** `Data and Computation Scripts/ComputeAndTransformStats.m`
**Inputs:** `waveformFull, t, X_Coordinates, Y_Coordinates, FileNamingArray, StatsToCompute, TotalSlices, 'raw'`

### Constants and Setup
```matlab
method = 'totalSlices';           % Hardcoded
Speed_of_Wave = 2968.67;         % m/s (hardcoded)
folderPath = 'Statistical Analysis';
transformType = 'raw';
```

### TOF and Depth Computation
```matlab
% Initialize maps
tofMap = nan(numY, numX);
depthMap = nan(numY, numX);

% TOF computation loop
for y = 1:numY
    for x = 1:numX
        waveform = squeeze(waveformFull(y, x, :));
        [pks, locs] = findpeaks(abs(waveform), 'SortStr', 'descend', 'NPeaks', 1);
        if ~isempty(locs)
            tof = t(locs(1));                              % First peak time
            tofMap(y, x) = tof;
            depthMap(y, x) = Speed_of_Wave * tof / 2 * 1000;  % Depth equation
        end
    end
end
```

### Depth Calculation Equation
```
depth = (Speed_of_Wave × tof / 2) × 1000
Where:
• Speed_of_Wave = 2968.67 m/s
• tof = time of flight (seconds)
• Factor of 2: round-trip time
• Factor of 1000: convert m to mm
Result: depth in mm
```

### RelativeTOF Computation
```matlab
relativeTofMaps = cell(numSegments, 1);

for seg = 1:numSegments
    segmentTofMap = nan(numY, numX);
    segmentTime = segmentedData{seg}.time;
    segmentWaveform = segmentedData{seg}.waveform;

    for y = 1:numY
        for x = 1:numX
            waveform = squeeze(segmentWaveform(y, x, :));
            [pks, locs] = findpeaks(abs(waveform), 'SortStr', 'descend', 'NPeaks', 1);
            if ~isempty(locs)
                segmentTofMap(y, x) = segmentTime(locs(1));
            end
        end
    end

    meanTof = mean(segmentTofMap(:), 'omitnan');
    relativeTofMaps{seg} = segmentTofMap - meanTof;
end
```

### MaxAmplitude Computation Algorithm
**Function:** `computeSignedMaxAmplitude(data)`
**File:** `Data and Computation Scripts/computeSignedMaxAmplitude.m`

```matlab
% Preserve sign of maximum amplitude
maxVals = max(data, [], 3, 'omitnan');
minVals = min(data, [], 3, 'omitnan');
signedMax = zeros(rows, cols);

for i = 1:rows
    for j = 1:cols
        if abs(maxVals(i,j)) >= abs(minVals(i,j))
            signedMax(i,j) = maxVals(i,j);
        else
            signedMax(i,j) = minVals(i,j);
        end
    end
end
```

### Statistical Function Definitions
```matlab
computeStat = @(data, statType, seg) struct(...
    'RMS', sqrt(mean(data.^2, 3, 'omitnan')), ...
    'MaxAmplitude', computeSignedMaxAmplitude(data), ...
    'Variance', var(data, 0, 3, 'omitnan'), ...
    'Skewness', skewness(data, 1, 3), ...
    'Kurtosis', kurtosis(data, 1, 3), ...
    'TOF', tofMap, ...
    'Depth', depthMap, ...
    'RelativeTOF', relativeTofMaps{seg}).(statType);
```

### Segment-wise Statistical Processing
```matlab
statMaps = cell(numSegments, 1);
for seg = 1:numSegments  % 400 segments
    waveformSegment = segmentedData{seg}.waveform;
    statVal = computeStat(waveformSegment, 'MaxAmplitude', seg);
    transformedStat = statVal;  % No transform for 'raw'
    statMaps{seg} = transformedStat;
end
```

### Statistical Data Structure Assembly
```matlab
statData.maps = statMaps;                    % 400 segments
statData.X_sub = X_Coordinates;
statData.Y_sub = Y_Coordinates;
statData.segmentTimes = cellfun(@(x) x.time, segmentedData);
statData.statType = 'MaxAmplitude';
statData.method = 'totalSlices';
statData.param = 400;
statData.transform_type = 'raw';
statData.FileNamingArray = FileNamingArray;
```

**File Save:** `Statistical Analysis/Raw_TotalSlices(400)_MaxAmplitude_GlobalScaling_No_Z-Score_Case1_X20-60_Y5-10_T3.0e-06-7.0e-06.mat`

---

## Chunk 6: XtVsYPlot Entry and Initialization

### Main Plotting Application Call
```matlab
if EnableMainPlottingApplication == 1  % true
    Main_Plotting_Application(FileNamingArray, StatsToCompute, ...
                             'totalSlices', TotalSlices, aspectRatio, savePlot);
end
```

### XtVsYPlot Function Entry
**File:** `XtVsYPlot.m`
**Function:** `XtVsYPlot(FileNamingArray, statTypes, method, param, aspectRatio, savePlot)`

**Inputs:**
- `FileNamingArray`: [1,1,1,1,20,60,5,10,3.0e-6,7.0e-6,0]
- `statTypes`: {'MaxAmplitude'}
- `method`: 'totalSlices'
- `param`: 400
- `aspectRatio`: [xlength, ylength, 1]
- `savePlot`: 0

### Input Validation
```matlab
if ~iscell(statTypes)
    statTypes = {statTypes};
end
% Result: statTypes = {'MaxAmplitude'} (already cell)
```

### Statistical Data Loading Setup
```matlab
nStats = length(statTypes);        % = 1
statDataArray = cell(nStats, 1);
sortedStatTypes = cell(nStats, 1);
```

### Function Call: loadStatData()
**Internal function in XtVsYPlot.m**
```matlab
for i = 1:nStats  % i = 1
    statType = statTypes{i};  % = 'MaxAmplitude'
    statData = loadStatData(FileNamingArray, statType, method, param);
    statDataArray{i} = statData;
    sortedStatTypes{i} = statType;
end
```

### loadStatData Function Execution
```matlab
folderPath = 'Statistical Analysis';
fileExtension = '.mat';

% Unpack FileNamingArray
caseNumber = 1;
xRange = [20, 60];
yRange = [5, 10];
TimeRange = [3.0e-6, 7.0e-6];

% File name construction
dataType = 'Raw';
xStr = 'X20-60';
yStr = 'Y5-10';
tStr = 'T3.0e-06-7.0e-06';
fileName = 'Raw_TotalSlices(400)_MaxAmplitude_GlobalScaling_No_Z-Score_Case1_X20-60_Y5-10_T3.0e-06-7.0e-06';
filePath = fullfile('Statistical Analysis', [fileName, '.mat']);

% Load data
loadedData = load(filePath);
statData = loadedData.statData;
```

### Loaded Statistical Data Structure
```matlab
statData.maps           % 400 MaxAmplitude maps
statData.X_sub          % X coordinates
statData.Y_sub          % Y coordinates
statData.segmentTimes   % Time vectors for each segment
statData.statType       % 'MaxAmplitude'
statData.method         % 'totalSlices'
statData.param          % 400

### Coordinate Extraction and Time Unit Conversion
```matlab
X_values = statDataArray{1}.X_sub;
Y_values = statDataArray{1}.Y_sub;
timeValues = statDataArray{1}.segmentTimes{1};

% Time unit conversion for display
if max(timeValues) > 1e-3
    timeLabel = 'Time (ms)';
    timeValues = timeValues * 1000;
else
    timeLabel = 'Time (μs)';
    timeValues = timeValues * 1e6;
end
```

### Visualization State Initialization
```matlab
visState.currentView = 'XtVsY';
visState.currentYIndex = round(length(Y_values) / 2);
visState.currentXIndex = round(length(X_values) / 2);
visState.currentSliceIndex = round(length(timeValues) / 2);
visState.showAligned = false;
visState.peakProminence = 0.1;
visState.minPeakDistance = 5;
visState.expectedLayers = 3;
visState.layerPaths = struct();
visState.layerSurfaces = [];
```

### Statistical Dropdown State Initialization
```matlab
statDropdownState.waveformMode = 'RawWaveform';
statDropdownState.linearProcessingMode = 'No Change';
statDropdownState.computationType = 'Amplitude';
statDropdownState.alignmentMethod = 'iterative';
statDropdownState.convergenceThreshold = 1e-6;
```

### Figure and Subplot Creation
```matlab
% Main figure
fig = figure('Name', 'X,t vs Y Visualization', 'NumberTitle', 'off', ...
             'Position', [100, 100, 1200, 800]);
set(fig, 'ResizeFcn', @(src, event) maintainUIPositions(src));

% Subplot layout calculation
nStats = 1;
nRows = ceil(sqrt(nStats));     % = 1
nCols = ceil(nStats / nRows);   % = 1

% Subplot margins and positioning
leftMargin = 0.1;
rightMargin = 0.76;
topMargin = 0.9;
bottomMargin = 0.35;
width = (rightMargin - leftMargin) / nCols;   % = 0.66
height = (topMargin - bottomMargin) / nRows;  % = 0.55

% Axes creation
axHandles = zeros(nStats, 1);
plotHandles = cell(nStats, 1);
for i = 1:nStats  % i = 1
    row = ceil(i / nCols);                    % = 1
    col = mod(i-1, nCols) + 1;               % = 1
    left = leftMargin + (col-1) * width;     % = 0.1
    bottom = bottomMargin + (nRows-row) * height;  % = 0.35
    axHandles(i) = axes('Position', [0.1, 0.35, 0.66, 0.55]);
end
```

### Initial Plot Creation (XtVsY View)
```matlab
currentYIndex = visState.currentYIndex;
statData = statDataArray{1};
numSegments = length(statData.maps);  % = 400
xtData = zeros(numSegments, length(X_values));

% Extract data for XtVsY view
for seg = 1:numSegments
    mapData = statData.maps{seg};
    xtData(seg, :) = mapData(currentYIndex, :);
end

% Create heatmap
plotHandles{1} = imagesc(X_values, timeValues, xtData, 'Parent', axHandles(1));
axis(axHandles(1), 'xy');
colormap(axHandles(1), 'jet');
colorbar(axHandles(1));
title(axHandles(1), sprintf('MaxAmplitude at Y = %.2f mm', Y_values(currentYIndex)));
xlabel(axHandles(1), 'X (mm)');
ylabel(axHandles(1), timeLabel);

% Apply aspect ratio
if ~isempty(aspectRatio)
    pbaspect(axHandles(1), aspectRatio);
    set(axHandles(1), 'DataAspectRatioMode', 'manual');
    set(axHandles(1), 'PlotBoxAspectRatioMode', 'manual');
end
```

---

## Chunk 7: XtVsYPlot View Options and Mathematical Processing

### View Options Available

#### 1. XtVsY View Processing
```matlab
plotData = zeros(numSegments, length(X_values));
for seg = 1:numSegments
    mapData = statData.maps{seg};
    plotData(seg, :) = mapData(currentYIndex, :);
end
xAxisValues = X_values;
yAxisValues = timeValues;
% Display: X Position (mm) vs Time/Depth
```

#### 2. YtVsX View Processing
```matlab
plotData = zeros(numSegments, length(Y_values));
for seg = 1:numSegments
    mapData = statData.maps{seg};
    plotData(seg, :) = mapData(:, currentXIndex);
end
xAxisValues = Y_values;
yAxisValues = timeValues;
% Display: Y Position (mm) vs Time/Depth
```

#### 3. XYVst View Processing
```matlab
if sliceIndex <= numSegments
    mapData = statData.maps{sliceIndex};
    plotData = mapData;
else
    plotData = statData.maps{1};
end
xAxisValues = X_values;
yAxisValues = Y_values;
% Display: X vs Y Position at time slice
```

#### 4. Waveform View Processing
```matlab
% Load original waveform data
loadWaveformData(fig);
% Apply waveform processing modes
% Display individual waveforms
```

#### 5. 3D View Processing
```matlab
% Create 3D scatter/surface visualization
create3DView(fig);
% X,Y,Time coordinates with amplitude coloring
```

### Statistical Analysis Options

#### Computation Types with Equations
1. **Amplitude:** `statMap = abs(baseSegmentData)` → `|data|`
2. **RMS:** `statMap = sqrt(baseSegmentData²)` → `√(data²)`
3. **Variance:** `statMap = computeLocalVariance(baseSegmentData)` → 3×3 sliding window variance
4. **Skewness:** `statMap = computeLocalSkewness(baseSegmentData)` → 3×3 sliding window skewness
5. **Kurtosis:** `statMap = computeLocalKurtosis(baseSegmentData)` → 3×3 sliding window kurtosis

### Waveform Processing Modes

#### 1. RawWaveform
```matlab
% No processing applied
```

#### 2. Envelope
```matlab
envelope = abs(hilbert(waveform));
% Amplify by 1.2, apply 3×3 smoothing kernel
kernel = ones(3,3) / 9;
statMap = conv2(statMap, kernel, 'same');
```

#### 3. FFT
```matlab
% Apply power scaling
statMap = abs(data).^0.8;
% Add frequency patterns
freqPattern = 0.1 * sin(2*pi*X/cols) .* cos(2*pi*Y/rows);
statMap = statMap * 1 + freqPattern;
```

#### 4. STFT
```matlab
% Apply sqrt scaling
statMap = sqrt(abs(data));
% Add time-frequency patterns
tfPattern = 0.15 * cos(2*pi*X/cols + seg/10) .* sin(2*pi*Y/rows);
statMap = statMap * 1 + tfPattern;
```

#### 5. Derivative
```matlab
% Enhance edges
statMap = abs(data).^1.5;
[gradX, gradY] = gradient(statMap);
gradMag = sqrt(gradX.^2 + gradY.^2);
statMap = statMap + 0.3 * gradMag;
```

#### 6. Integral
```matlab
% Accumulate effects
statMap = abs(data).^0.7;
for i = 2:rows
    statMap(i,:) = statMap(i,:) + 0.1 * statMap(i-1,:);
end
```

#### 7. Wavelet
```matlab
% Multi-scale analysis
statMap = abs(data).^0.8;
scale1 = 0.1 * sin(2*pi*X/cols/4) .* sin(2*pi*Y/rows/4);
scale2 = 0.05 * sin(2*pi*X/cols/8) .* sin(2*pi*Y/rows/8);
statMap = statMap * 1 + scale1 + scale2;
```
```
