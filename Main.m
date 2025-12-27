%% Main Script - NDE Analysis
% This script serves as the central control for Non-Destructive Evaluation (NDE) analysis of C-scan data.
% It allows users to load a specified dataset, process it according to user-defined settings, and visualize
% the results through 2D waveform plots, FFT plots, and statistical heatmaps with interactive sliders for
% segmented data.

%% --- User-Defined Settings ---

clc, clear, close all;

% Add directories to path for utility functions (robust bootstrap)
% Ensure Utils is on path based on this file's location (independent of pwd)
thisDir  = fileparts(mfilename('fullpath'));
utilsDir = fullfile(thisDir, 'Utils');
if exist(utilsDir, 'dir') && ~contains(path, utilsDir)
    addpath(utilsDir);
end
% Centralized project path setup
addUtilsPath();
%% 
% I. Dataset Selection
DATASET = 4;

% Dataset number 1 : L0P5S25_10Mhz  ( No Defect, Sample 25)
% Dataset number 2 : L8P12S5_10Mhz  ( 8 Degree Defect, Sample 5)
% Dataset number 3 : L16P1S9_10MHz  ( 16 Degree Defect, Sample 9)
% Dataset number 4 : L16P1S10_10MHz ( 16 Degree Defect, Sample 10)
% Dataset number 5 : L16P34S3_10MHz ( 16 Degree Defect, Sample 3)

%% II. Data Processing Options
% Trim Ranges:
TrimWaveData        = 1;    % 1 = trim spatial data, 0 = keep all
TrimTimeRange       = 1;    % 1 = trim time data, 0 = keep all

xRange          = [30 60];   % X range in mm
yRange          = [10 20];    % Y range in mm
TimeRange       = [0.0e-6 2.00e-6]; % Time range in seconds
sampleThickness = 10;  % (mm) Replace 10 with your actual Z-height

%% III. Waveform Segmentation Settings
TotalSlices                     = 200; % Number of slices

%% IV. Alignment Settings

AlignWaveFormsAtFirstPeak = 0;                  % 1 = align waveforms at first peak, 0 = keep original
Align_Between_Time_Range  = [0.3e-6 0.4e-6];   % Find first peak between these time ranges.
AlignPeakType = 'positive';                     % Options: 'positive', 'negative', or 'both'

% 3D Peak Alignment Settings (for create3DPeakPlots)
FirstPeakAmplitudeThreshold = 1.0;  % Minimum amplitude threshold for first peak alignment
                                    % Only peaks/valleys with |amplitude| >= this value will be
                                    % considered for alignment reference. Set to 0 to disable.
                                    % Recommended: 0.1 (10% of max amplitude)

AlignDiagnostics = 0;               % Enable diagnostic logging for alignment
AlignDiagnosticSamples = 100;       % Number of random waveforms to sample for diagnostics

%% V. Pre-Computing Settings
savePlot        = 1;         % 1 = savee plots, 0 = display only
SkipOverWriteRequests = 1;   % 1 = skip overwrite prompts AND enable peak data caching
saveData = 1;

%% VI. Statistical Analysis Options
ComputeMaxAmplitude     = 1;
ComputeRMS              = 0;
ComputeVariance         = 0;
ComputeSkewness         = 0;
ComputeKurtosis         = 0;

%% VII. Time-Space Visualizations
EnableMainPlottingApplication   = 1; % Enable Main Plotting Application

%% VIII. Peak Detection and Analysis Settings
EnablePeakExtraction    = 1; % Master control switch for peak extraction and analysis
Enable3DPeakVallyPlotting = 1; % 1 = Enable 3D peak/valley plot generation, 0 = Disable


% Peak Detection Settings
PeakDetectionType = 'peaks'; % Options: 'peaks', 'valleys', 'both'
minPeakHeight = 0.1;        % Minimum peak height as fraction of max amplitude
minPeakDistance = 2;        % Minimum distance between peaks in samples (must be >= 1)
minPeakProminence = 0.2;    % Minimum peak prominence as fraction of max amplitude
useSlopeDetection = true;   % Use slope-based detection for better accuracy
slopeThreshold = 0.05;      % Slope threshold for detecting transitions

% Peak 3D Visualization Options
align3DPlotsWithTime = true;  % Align 3D plots using time values (not shifted)
show3DPeaksAndValleys = true; % Show both peaks and valleys in 3D representation

%% 1. Check for Precomputed Data First
% Add alignment flag to FileNamingArray if it exists
if exist('AlignWaveFormsAtFirstPeak', 'var')
    alignFlag = AlignWaveFormsAtFirstPeak;
else
    alignFlag = 0;
end

% Use DATASET as case number for consistency
caseNumber = DATASET;
FileNamingArray = buildFileNamingArray('DATASET', DATASET, ...
    'CaseNumber', caseNumber, ...
    'TrimWaveData', TrimWaveData, 'TrimTimeRange', TrimTimeRange, ...
    'xRange', xRange, 'yRange', yRange, 'TimeRange', TimeRange, ...
    'AlignAtFirstPeak', alignFlag);
folderPath = fullfile(pwd, 'Saved Wave Forms');
fileExtension = '.mat';
fileType = 'Waveform';

[dataFile, precomputedExists] = buildAndCheckFile(fileType, FileNamingArray, folderPath, fileExtension);

if precomputedExists && SkipOverWriteRequests
    disp('Precomputed waveform file already exists. Skipping all data loading and processing (SkipOverWriteRequests enabled).');
    % Load essential parameters from the precomputed file
    loadedData = load(dataFile);
    t = loadedData.savedData.t;
    X_Coordinates = loadedData.savedData.X_Coordinates;
    Y_Coordinates = loadedData.savedData.Y_Coordinates;

    % Calculate derived parameters
    xlength = max(X_Coordinates) - min(X_Coordinates);
    ylength = max(Y_Coordinates) - min(Y_Coordinates);
    aspectRatio = [xlength ylength 1];
    numY_sub = length(Y_Coordinates);
    numX_sub = length(X_Coordinates);

    % Reconstruct 3D matrix from saved 2D waveform array for statistical computations
    waveformArray = loadedData.savedData.waveformArray;
    numTimePoints = length(t);

    % Vectorized 2D to 3D conversion (replaces slow nested loop)
    % waveformArray is [numWaveforms, numTimePoints] with Y-major ordering
    % Target: waveformFull is [Y, X, T]
    waveformFull = reshape(waveformArray, [numY_sub, numX_sub, numTimePoints]);

    disp(['Loaded precomputed data with ', num2str(length(t)), ' time points.']);

    % Skip to segmentation section
    skipDataProcessing = true;
else
    % Need to process data from scratch
    skipDataProcessing = false;

    disp('Loading dataset...');
    [DataStructure] = loadData(DATASET);

    % Use raw data only
    waveform3DMatrix = DataStructure.waveform3DMatrix;
    waveformFull = waveform3DMatrix;

    % Generate time vector for alignment
    SampRate = DataStructure.sampRate * 1e6;  % Convert to Hz
    dt = 1 / SampRate;
    tSize = DataStructure.tSize;
    t_full = dt:dt:(tSize * dt);

    % Align waveforms at first peak if enabled (before trimming)
    if exist('AlignWaveFormsAtFirstPeak', 'var') && AlignWaveFormsAtFirstPeak == 1
    % Check if time range for alignment is specified
    if exist('Align_Between_Time_Range', 'var') && ~isempty(Align_Between_Time_Range) && length(Align_Between_Time_Range) == 2
        % Check if peak type and diagnostics are specified
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
                [waveformFull, shiftIndices] = AlignWaveformsAtFirstPeak(waveformFull, t_full, 0, Align_Between_Time_Range, 'positive', AlignDiagnostics);
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

        % Update the original matrices with aligned data
        waveform3DMatrix = waveformFull;
        DataStructure.waveform3DMatrix = waveformFull;
    end

    %% 2. Data Extraction and Trimming (only if not using precomputed data)
    disp('Extracting and trimming data...');
    [SampRate, ScanRes, indexRes, xSize, ySize, tSize, waveform3DMatrix, waveformEnvelope3DMatrix, t, X_Coordinates, Y_Coordinates] = ...
        ExtractAndTrimData(DataStructure, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange);

    % Use raw data only (after trimming)
    waveformFull = waveform3DMatrix;

    % Calculate aspect ratio
    xlength = max(X_Coordinates) - min(X_Coordinates);
    ylength = max(Y_Coordinates) - min(Y_Coordinates);
    aspectRatio = [xlength ylength 1];

    % Calculate grid dimensions for peak extraction
    numY_sub = length(Y_Coordinates);
    numX_sub = length(X_Coordinates);

    disp(['Number of time points after trimming: ', num2str(length(t))]);

    % Precompute and Save Waveforms
    disp('Precomputing and saving waveforms...');
    PreComputeWaveformData(waveform3DMatrix, waveformEnvelope3DMatrix, t, X_Coordinates, Y_Coordinates, FileNamingArray);
end



%% 3. Waveform Segmentation
% Check if segmented data already exists
[~, caseNumber, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange, ~] = unpackFileNamingArray(FileNamingArray);

% Generate range strings based on trimming flags (consistent with buildAndCheckFile)
if TrimWaveData == 0
    xStr = 'FullX';
    yStr = 'FullY';
else
    xStr = sprintf('X%.0f-%.0f', xRange(1), xRange(2));
    yStr = sprintf('Y%.0f-%.0f', yRange(1), yRange(2));
end

if TrimTimeRange == 0
    tStr = 'FullT';
else
    tStr = sprintf('T%.1e-%.1e', TimeRange(1), TimeRange(2));
end

segmentedFileName = sprintf('SegmentedWaveformCase%d_Raw_%s_%s_%s_TotalSlices_%d.mat', ...
    caseNumber, xStr, yStr, tStr, TotalSlices);
segmentedFilePath = fullfile(pwd, 'Saved Wave Forms', segmentedFileName);

if exist(segmentedFilePath, 'file') && SkipOverWriteRequests
    disp('Segmented waveform file already exists. Skipping segmentation (SkipOverWriteRequests enabled).');
    % Load existing segmented data for later use
    loadedSegmented = load(segmentedFilePath);
    segmentedDataTotalSlices = loadedSegmented.segmentedData;
else
    disp('Performing total slices segmentation...');
    segmentedDataTotalSlices = SegmentPrecomputedWaveform(dataFile, 'totalSlices', [], TotalSlices);
    MasterSave('SegmentedWaveform', segmentedDataTotalSlices, FileNamingArray, 'totalSlices', TotalSlices);
end



%% 4. Peak Extraction and Analysis
if EnablePeakExtraction == 1
    disp('Performing peak extraction and analysis...');

    % Create peak options structure
    peakOptions = struct();
    peakOptions.minPeakHeight = minPeakHeight;
    peakOptions.minPeakDistance = minPeakDistance;
    peakOptions.minPeakProminence = minPeakProminence;
    peakOptions.useSlopeDetection = useSlopeDetection;
    peakOptions.slopeThreshold = slopeThreshold;

    % Create 3D visualization options structure
    viz3DOptions = struct();
    viz3DOptions.Enable3DPeakVallyPlotting = Enable3DPeakVallyPlotting;
    viz3DOptions.align3DPlotsWithTime = align3DPlotsWithTime;
    viz3DOptions.show3DPeaksAndValleys = show3DPeaksAndValleys;

    % Create comprehensive cache key for peak extraction results
    peakCacheKey = struct();
    peakCacheKey.FileNamingArray = FileNamingArray;
    peakCacheKey.PeakDetectionType = PeakDetectionType;
    peakCacheKey.peakOptions = peakOptions;
    peakCacheKey.viz3DOptions = viz3DOptions;

    % Generate cache filename based on all relevant parameters
    cacheHash = generatePeakCacheHash(peakCacheKey);
    peakCacheFolder = fullfile(pwd, 'Peak Cache');
    if ~exist(peakCacheFolder, 'dir')
        mkdir(peakCacheFolder);
    end
    peakCacheFilename = fullfile(peakCacheFolder, sprintf('PeakCache_Dataset%d_%s_%s.mat', DATASET, PeakDetectionType, cacheHash));

    % Check if cached peak data exists and is valid
    fprintf('Checking for cached peak data...\n');
    fprintf('Cache file: %s\n', peakCacheFilename);
    fprintf('File exists: %s\n', mat2str(exist(peakCacheFilename, 'file')));
    fprintf('SkipOverWriteRequests: %d\n', SkipOverWriteRequests);

    if exist(peakCacheFilename, 'file') && SkipOverWriteRequests == 1
        fprintf('Cached peak extraction results found. Loading from cache...\n');
        try
            cachedData = load(peakCacheFilename);
            peakData = cachedData.peakData;
            fprintf('Successfully loaded cached peak data from: %s\n', peakCacheFilename);
            fprintf('Loaded %d waveforms of peak data from cache.\n', length(peakData));

            % Check if spatial coordinates are cached
            if isfield(cachedData, 't') && isfield(cachedData, 'X_Coordinates')
                % Use cached spatial coordinates
                t = cachedData.t;
                X_Coordinates = cachedData.X_Coordinates;
                Y_Coordinates = cachedData.Y_Coordinates;
                numY_sub = cachedData.numY_sub;
                numX_sub = cachedData.numX_sub;
                fprintf('Using cached spatial coordinates.\n');

                % Only load waveform data if 3D visualization is enabled
                if viz3DOptions.Enable3DPeakVallyPlotting == 1
                    fprintf('Loading waveform data for 3D visualization...\n');
                    [waveformArray, ~, ~, ~, ~, ~] = SmartDataLoader(FileNamingArray);
                end
            else
                % Fallback: load all data including spatial coordinates
                fprintf('Cached spatial coordinates not found. Loading all data...\n');
                [waveformArray, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub] = SmartDataLoader(FileNamingArray);

                % Update the cache with spatial coordinates for next time
                fprintf('Updating cache with spatial coordinates...\n');
                try
                    save(peakCacheFilename, 'peakData', 'peakCacheKey', 'PeakDetectionType', 'DATASET', ...
                         't', 'X_Coordinates', 'Y_Coordinates', 'numY_sub', 'numX_sub', '-append');
                    fprintf('Cache updated with spatial coordinates.\n');
                catch ME
                    fprintf('Warning: Failed to update cache with spatial coordinates (%s)\n', ME.message);
                end
            end

            % Create 3D visualization if enabled
            if viz3DOptions.Enable3DPeakVallyPlotting == 1
                fprintf('Creating 3D peak/valley representations from cached data...\n');
                create3DPeakPlots([], [], peakData, waveformArray, t, ...
                    X_Coordinates, Y_Coordinates, numY_sub, numX_sub, PeakDetectionType, ...
                    viz3DOptions.align3DPlotsWithTime, viz3DOptions.show3DPeaksAndValleys, ...
                    FirstPeakAmplitudeThreshold);
            end
        catch ME
            fprintf('Warning: Failed to load cached peak data (%s). Recomputing...\n', ME.message);
            peakData = []; % Force recomputation
        end
    else
        if ~exist(peakCacheFilename, 'file')
            fprintf('No cached peak data found. Will compute from scratch.\n');
        elseif SkipOverWriteRequests ~= 1
            fprintf('Peak data caching disabled (SkipOverWriteRequests = %d). Will compute from scratch.\n', SkipOverWriteRequests);
        end
        peakData = []; % No cache or cache disabled
    end

    % Perform peak extraction if no valid cached data
    if isempty(peakData)
        fprintf('Computing peak extraction...\n');

        % Use SmartDataLoader for optimized data loading and caching
        fprintf('Loading waveform data using SmartDataLoader...\n');
        [waveformArray, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub] = SmartDataLoader(FileNamingArray);

        % Check if data was loaded successfully
        if isempty(t) || size(waveformArray, 2) == 0
            error('No time data loaded. Check your TimeRange parameter or set TrimTimeRange = 0 to see the full data range.');
        end

        fprintf('Loaded %d waveforms with %d time points each.\n', size(waveformArray, 1), size(waveformArray, 2));

        % Reconstruct 3D matrix for statistical computations (needed later)
        numTimePoints = size(waveformArray, 2);
        % waveformArray rows were built with j-fast ordering (X varying fastest) via permute([2 1 3])
        % To rebuild [Y X T], we need to inverse that ordering: reshape then ipermute
        waveformFull = ipermute(reshape(waveformArray, [numX_sub, numY_sub, numTimePoints]), [2 1 3]);

        % Process peak extraction using consistent variable names
        peakData = PeakExtractionProcessor(...
            waveformArray, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, ...
            PeakDetectionType, peakOptions, viz3DOptions);

        % Save peak extraction results to cache with spatial coordinates
        try
            save(peakCacheFilename, 'peakData', 'peakCacheKey', 'PeakDetectionType', 'DATASET', ...
                 't', 'X_Coordinates', 'Y_Coordinates', 'numY_sub', 'numX_sub');
            fprintf('Peak extraction results cached to: %s\n', peakCacheFilename);
        catch ME
            fprintf('Warning: Failed to save peak cache (%s)\n', ME.message);
        end

        % Also save to a dedicated folder (instead of main directory) for backward compatibility
        legacyFolder = fullfile(pwd, 'Peak Extraction Results');
        if ~exist(legacyFolder, 'dir')
            mkdir(legacyFolder);
        end
        legacyFilename = fullfile(legacyFolder, sprintf('PeakExtractionResults_Dataset%d_%s.mat', DATASET, PeakDetectionType));
        try
            save(legacyFilename, 'peakData', 'PeakDetectionType', 'DATASET');
            fprintf('Peak extraction results also saved to: %s\n', legacyFilename);
        catch ME
            fprintf('Warning: Failed to save legacy peak results (%s)\n', ME.message);
        end
    end
end

%% 5. Statistical Computations and Statistical Plotting
% Pre-allocate StatsToCompute array to avoid dynamic arrays
maxStats = 5; % Maximum possible number of stats
StatsToCompute = cell(maxStats, 1);
statCount = 0;

if ComputeMaxAmplitude == 1
    statCount = statCount + 1;
    StatsToCompute{statCount} = 'MaxAmplitude';
end
if ComputeVariance == 1
    statCount = statCount + 1;
    StatsToCompute{statCount} = 'Variance';
end
if ComputeSkewness == 1
    statCount = statCount + 1;
    StatsToCompute{statCount} = 'Skewness';
end
if ComputeKurtosis == 1
    statCount = statCount + 1;
    StatsToCompute{statCount} = 'Kurtosis';
end
if ComputeRMS == 1
    statCount = statCount + 1;
    StatsToCompute{statCount} = 'RMS';
end

% Trim to actual size
StatsToCompute = StatsToCompute(1:statCount);


if ~isempty(StatsToCompute)
    disp('Computing raw stats for total slices...');
    % Pass precomputed segmentation if available to avoid re-segmentation
    if exist('segmentedDataTotalSlices','var') && ~isempty(segmentedDataTotalSlices)
        ComputeAndTransformStats(waveformFull, t, X_Coordinates, Y_Coordinates, FileNamingArray, StatsToCompute, TotalSlices, 'raw', segmentedDataTotalSlices);
    else
        ComputeAndTransformStats(waveformFull, t, X_Coordinates, Y_Coordinates, FileNamingArray, StatsToCompute, TotalSlices, 'raw');
    end
end

%% 6. Time-Space Visualizations
if EnableMainPlottingApplication == 1
    disp('Generating time-space visualizations...');
    % Call XtVsYPlot.m directly instead of the refactored Main_Plotting_Application
    XtVsYPlot(FileNamingArray, StatsToCompute, 'totalSlices', TotalSlices, aspectRatio, savePlot);
end



%% 8. Final Message
disp('Analysis complete!');