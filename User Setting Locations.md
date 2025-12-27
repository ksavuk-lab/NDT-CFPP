# User Setting Locations Guide
**NDE Analysis Program - Configuration Reference**

---

## Purpose
This document provides a comprehensive guide to all user-configurable settings within the NDE Analysis Program. It directs users to the specific files and line numbers where settings can be modified to customize program behavior, data processing, visualization, and analysis parameters.

---

## Table of Contents
1. [Critical Settings - Start Here](#critical-settings---start-here)
2. [Main.m - Primary User Settings](#mainm---primary-user-settings)
3. [Data File Locations](#data-file-locations)
4. [Processing & Analysis Settings](#processing--analysis-settings)
5. [Visualization Settings](#visualization-settings)
6. [Advanced Settings](#advanced-settings)

---

## Critical Settings - Start Here

### 1. **Data File Locations**
**File:** `Data and Computation Scripts/loadData.m`  
**Lines:** 9-37

This is where the program searches for your data files. The program tries multiple directory paths in order:

1. **Environment Variable (Recommended for new users):**
   - Set environment variable `NDE_DATA_ROOT` to point to your data directory
   - Example (Mac/Linux): `export NDE_DATA_ROOT=/path/to/your/Data`
   - Example (Windows): `set NDE_DATA_ROOT=C:\path\to\your\Data`

2. **Relative Path:** `../../Data/LaminateBeam/L0/` (relative to code directory)

3. **User HOME Path:** `~/Documents/X_ University/Uni Research/NDT/NDE_SDSU-USD/Data/LaminateBeam/L0/`

4. **Hardcoded Paths (Lines 35-37):** Modify these if your data is in a different location

**Dataset Filenames (Lines 41-54):**
- Dataset 1: `L0P5S25_10Mhz_TimingCorrected.mat`
- Dataset 2: `L8P12S5_10Mhz.mat`
- Dataset 3: `L16P1S9_10MHz.mat`
- Dataset 4: `L16P1S10_10MHz.mat`
- Dataset 5: `L16P34S3_10MHz.mat`

### 2. **Scan Step Resolution**
**File:** `Data and Computation Scripts/ExtractAndTrimData.m`  
**Lines:** 7-9, 28-29

The scan resolution is extracted from the loaded data file:
- `ScanRes = DataStructure.scanRes` (Line 8)
- Used to generate spatial coordinates: `X_Coordinates = (0:xSize-1) * ScanRes` (Line 28)

**Note:** This is read from the data file, not directly user-configurable. To change it, you would need to modify the source data or the DataStructure.

### 3. **Frequency Settings (Sampling Rate)**
**File:** `Data and Computation Scripts/ExtractAndTrimData.m`  
**Line:** 7

The sampling rate is extracted from the loaded data:
- `SampRate = DataStructure.sampRate * 1e6` (converts to Hz)

**Note:** This is read from the data file. The original data was collected at a specific sampling rate (typically 100 MHz for 10 MHz transducers).

**Alternative Location - Wave Velocity Calculator:**  
**File:** `Tools/WaveVelocityCalculator.m`  
**Lines:** 5-7
- `OutputFreq_Mhz = 10` - Transducer output frequency
- `SamplingRate_Mhz = 100` - Sampling rate (10x the output frequency)

---

## Main.m - Primary User Settings

**File:** `Main.m`

This is the main control script where most user settings are configured.

### I. Dataset Selection (Lines 21-28)
```matlab
DATASET = 4;  % Line 22
```
Options:
- 1: L0P5S25_10Mhz (No Defect, Sample 25)
- 2: L8P12S5_10Mhz (8 Degree Defect, Sample 5)
- 3: L16P1S9_10MHz (16 Degree Defect, Sample 9)
- 4: L16P1S10_10MHz (16 Degree Defect, Sample 10)
- 5: L16P34S3_10MHz (16 Degree Defect, Sample 3)

### II. Data Processing Options (Lines 30-38)

**Spatial and Time Trimming:**
- `TrimWaveData = 1` (Line 32) - 1 = trim spatial data, 0 = keep all
- `TrimTimeRange = 1` (Line 33) - 1 = trim time data, 0 = keep all

**Trim Ranges:**
- `xRange = [10 70]` (Line 35) - X range in mm
- `yRange = [5 20]` (Line 36) - Y range in mm
- `TimeRange = [0.0e-6 1.60e-6]` (Line 37) - Time range in seconds
- `sampleThickness = 1.67` (Line 38) - Sample Z-height in mm

### III. Waveform Segmentation Settings (Lines 40-41)
- `TotalSlices = 200` (Line 41) - Number of time slices for segmentation

### IV. Alignment Settings (Lines 43-58)

**Basic Alignment:**
- `AlignWaveFormsAtFirstPeak = 1` (Line 45) - 1 = align, 0 = no alignment
- `AlignWaveFormsAtFirstPeak_V2 = 1` (Line 46) - Version 2 alignment flag
- `Align_Between_Time_Range = [0.3e-6 0.35e-6]` (Line 47) - Time range to search for first peak
- `AlignPeakType = 'positive'` (Line 48) - Options: 'positive', 'negative', 'both'

**3D Peak Alignment:**
- `FirstPeakAmplitudeThreshold = 1.0` (Line 51) - Minimum amplitude for alignment reference
  - Set to 0 to disable threshold
  - Recommended: 0.1 (10% of max amplitude)

**Diagnostics:**
- `AlignDiagnostics = 1` (Line 56) - Enable diagnostic logging
- `AlignDiagnosticSamples = 100` (Line 57) - Number of waveforms to sample for diagnostics

### V. Pre-Computing Settings (Lines 64-67)
- `savePlot = 1` (Line 65) - 1 = save plots, 0 = display only
- `SkipOverWriteRequests = 1` (Line 66) - 1 = skip overwrite prompts AND enable caching
- `saveData = 1` (Line 67) - 1 = save processed data

### VI. Statistical Analysis Options (Lines 69-74)
- `ComputeMaxAmplitude = 1` (Line 70) - Compute maximum amplitude statistic
- `ComputeRMS = 0` (Line 71) - Compute RMS statistic
- `ComputeVariance = 0` (Line 72) - Compute variance statistic
- `ComputeSkewness = 0` (Line 73) - Compute skewness statistic
- `ComputeKurtosis = 0` (Line 74) - Compute kurtosis statistic

### VII. Time-Space Visualizations (Lines 76-77)
- `EnableMainPlottingApplication = 1` (Line 77) - Enable main plotting interface

### VIII. Peak Detection and Analysis Settings (Lines 79-94)

**Master Controls:**
- `EnablePeakExtraction = 1` (Line 80) - Master switch for peak extraction
- `Enable3DPeakVallyPlotting = 1` (Line 81) - Enable 3D peak/valley plots

**Peak Detection Parameters:**
- `PeakDetectionType = 'both'` (Line 85) - Options: 'peaks', 'valleys', 'both'
- `minPeakHeight = 0.2` (Line 86) - Minimum peak height (fraction of max amplitude)
- `minPeakDistance = 2` (Line 87) - Minimum distance between peaks (samples, must be >= 1)
- `minPeakProminence = 0.2` (Line 88) - Minimum peak prominence (fraction of max amplitude)
- `useSlopeDetection = true` (Line 89) - Use slope-based detection for better accuracy
- `slopeThreshold = 0.05` (Line 90) - Slope threshold for detecting transitions

**3D Visualization Options:**
- `align3DPlotsWithTime = true` (Line 93) - Align 3D plots using time values
- `show3DPeaksAndValleys = true` (Line 94) - Show both peaks and valleys in 3D

---

## Data File Locations

### Primary Data Loading
**File:** `Data and Computation Scripts/loadData.m`

**Environment Variable Setup (Lines 12-17):**
- Checks for `NDE_DATA_ROOT` environment variable
- Automatically constructs paths: `$NDE_DATA_ROOT/LaminateBeam/L0/`

**Fallback Paths (Lines 19-37):**
1. Relative path from code directory
2. User HOME-based path
3. Hardcoded legacy paths (modify these for your system)

### Output Directories
**File:** `Main.m`

**Saved Waveforms:**
- Directory: `Saved Wave Forms/` (Line 111)
- Created automatically if doesn't exist

**Statistical Analysis:**
- **File:** `Data and Computation Scripts/ComputeAndTransformStats.m`
- Directory: `Statistical Analysis/` (Line 30)

**Saved Plots:**
- Directory: `Saved Plots/`
- Configured in `XtVsYPlot.m`

**Peak Cache:**
- Directory: `Peak Cache/` (Line 297 in Main.m)
- Stores computed peak extraction results for faster reloading

**Diagnostic Logs:**
- **File:** `Data and Computation Scripts/AlignWaveformsAtFirstPeak.m`
- Directory: `Logs/` (Line 84)
- Created when `AlignDiagnostics = 1`

---

## Processing & Analysis Settings

### Wave Speed Calculation
**File:** `Data and Computation Scripts/ComputeAndTransformStats.m`  
**Line:** 27
- `Speed_of_Wave = 2968.67` (m/s)
- Used for depth calculations: `depth = (Speed_of_Wave * tof / 2) * 1000` (mm)

### Alignment Algorithm
**File:** `Data and Computation Scripts/AlignWaveformsAtFirstPeak.m`

**Default Parameters (Lines 36-65):**
- Default time range: entire waveform if not specified
- Default peak type: 'positive' (Line 54)
- Default diagnostics: disabled (Line 59)
- Default diagnostic samples: 100 (Line 64)

**Peak Detection Logic (Lines 68-77):**
- Supports positive peaks, negative peaks (valleys), or both
- Uses magnitude comparison when 'both' is selected

### Peak Extraction Processor
**File:** `Data and Computation Scripts/PeakExtractionProcessor.m`

**Processing Logic (Lines 44-58):**
- Automatically sets `findValleys` and `findPeaks` flags based on `PeakDetectionType`
- Processes each waveform individually
- Counts peaks, valleys, and total transitions

---

## Visualization Settings

### XtVsYPlot - Main Visualization Interface
**File:** `XtVsYPlot.m`

**Initial View Settings (Lines 180-191):**
- `currentView = 'XtVsY'` - Default view (X,t vs Y)
- `currentSliceIndex = 1` - Starting slice
- `useZScore = false` - Z-score normalization toggle
- `showAllWaveformsInView = false` - Waveform overlay toggle
- `showLayerPaths = false` - Layer path overlay toggle

**Layer Detection Parameters (Lines 2975-2986):**
- Peak prominence slider (adjustable in GUI)
- Minimum distance input (adjustable in GUI)
- Expected layers input (adjustable in GUI)

**Plot Settings (Lines 3708-3716):**
- `useGlobalScale` - Global vs. local color scaling
- `plotType` - Type of plot visualization
- `selectedColormap` - Color scheme
- `contrastEnhancement` - Contrast adjustment
- `useTimeRange` - Time range filtering
- `timeRangeMin`, `timeRangeMax` - Time range bounds

### 3D Peak Visualization
**File:** `Data and Computation Scripts/create3DPeakPlots.m`

**Plate Generation Options (Lines 600-606):**
- `enablePlateGeneration = false` (Line 601) - DISABLED by default to prevent hanging
- `amplitudeTolerance = 0.1` (Line 602) - Amplitude grouping tolerance
- `timeTolerance = 1` (Line 603) - Time grouping tolerance (microseconds)
- `minPointsPerPlate = 10` (Line 604) - Minimum points to form a plate
- `plateType = 'both'` (Line 605) - 'peaks', 'valleys', or 'both'

**Visualization Mode (Line 609):**
- `visualizationMode = 'points'` - Options: 'points' or 'plates'

**Plate Generation Algorithms (Lines 611-620):**
- `advancedAlgorithm = 'amplitude_time_spatial'` (Line 613) - Default algorithm
  - Algorithm #1: Amplitude + Time + Spatial
  - Algorithm #2: Adaptive Cube
  - Algorithm #3: Group Averaging
  - Algorithm #4: Topological Plates

**Visual Rendering Methods (Lines 622-628):**
- `method = 'scatter_points'` (Line 624) - Rendering method
  - Options: 'scatter_points', 'triangulated_surfaces', 'convex_hull', 'voxel_blocks', 'smooth_surfaces', 'wire_frames'
- `transparency = 5` (Line 625) - Transparency level (1-10 scale)
- `edgeStyle = 'smooth'` (Line 626) - Options: 'show_edges', 'smooth', 'wireframe_only'
- `colorCoding = 'amplitude'` (Line 627) - Options: 'amplitude', 'depth', 'plate_id', 'defect_severity'
- `thicknessThreshold = 3` (Line 628) - Threshold in microseconds for thick vs thin plates

**Additional Grouping (Lines 630-635):**
- `method = 'none'` (Line 632) - Post-processing grouping method
  - Options: 'none', 'spatial_proximity', 'depth_layers', 'grid_zones', 'defect_severity'
- `spatialRadius = 5.0` (Line 633) - Spatial proximity radius (mm)
- `gridZones = 3` (Line 634) - NxN grid zones
- `depthLayers = 3` (Line 635) - Number of depth layers

**Surface Options (Lines 639-652):**
- `creationMethod = 'triangulation'` (Line 641) - Surface creation method
- `hideThreshold = 8` (Line 648) - Threshold for hiding surrounded points
- `hideTolerance = 0.5` (Line 649) - Multiplier for hiding amplitude tolerance
- `cubeSize = 2.0` (Line 650) - Default cube size in mm

**3D Visual Settings (Lines 4700-4731):**
- Slice planes: alpha = 0.7, shading = 'interp'
- Smooth surfaces: alpha = 0.8, shading = 'interp'
- Pseudocolor: alpha = 0.9, shading = 'flat'
- Filled contour: alpha = 0.6
- Heatmap: alpha = 0.8, shading = 'interp'

---

## Advanced Settings

### Error Handling and Diagnostics
**File:** `Data and Computation Scripts/ErrorHandler.m`

**Verbosity Levels (Lines 6-11):**
- `SILENT = 0` - No output
- `WARNING = 1` - Warning messages only
- `ERROR = 2` - Error messages with stack trace (default, Line 14)
- `VERBOSE = 3` - Detailed debugging information

**Logging (Lines 15-16):**
- `logToFile = false` - Enable file logging
- `logFileName = ''` - Log file name

### Performance Monitoring
**File:** `Tools/run_with_profiling.m`

**Profiling Setup (Lines 10-20):**
- Sets error handler to VERBOSE mode
- Enables `SkipOverWriteRequests` for non-interactive behavior
- Starts performance monitor with MATLAB profiler

### Wave Velocity Calculator
**File:** `Tools/WaveVelocityCalculator.m`

**Transducer Properties (Lines 5-7):**
- `OutputFreq_Mhz = 10` - Output frequency in MHz
- `SamplingRate_Mhz = 100` - Sampling rate in MHz

**Material Properties (Line 10):**
- `Velocity_Range = [2000, 3500]` - Expected velocity range (m/s)

**Front Wall Detection (Line 14):**
- `FrontWall_Range_Samples = [4, 15]` - Sample range for front wall detection

**Data Quality Control (Lines 17-20):**
- `Enable_Data_Filtering = true` - Enable filtering of problematic samples
- `Problematic_Samples = [3, 4, 5]` - Samples with known defects

### Caching System
**File:** `Main.m`

**Peak Cache (Lines 289-301):**
- Cache hash generation based on all relevant parameters
- Cache directory: `Peak Cache/`
- Cache filename format: `PeakCache_Dataset{N}_{type}_{hash}.mat`

**Waveform Cache (Lines 117-149):**
- Checks for precomputed waveform files
- Loads from cache when `SkipOverWriteRequests = 1`
- Skips all data loading and processing if cache exists

---

## Quick Reference Summary

### Most Commonly Modified Settings:

1. **Dataset Selection:** `Main.m`, Line 22
2. **Spatial Range:** `Main.m`, Lines 35-36
3. **Time Range:** `Main.m`, Line 37
4. **Number of Slices:** `Main.m`, Line 41
5. **Alignment Enable:** `Main.m`, Line 45
6. **Peak Detection Type:** `Main.m`, Line 85
7. **Enable Visualizations:** `Main.m`, Lines 77, 80, 81
8. **Data File Path:** `Data and Computation Scripts/loadData.m`, Lines 35-37 (or set `NDE_DATA_ROOT` environment variable)

### Settings That Are Read from Data Files (Not Directly Configurable):
- Sampling Rate (`SampRate`)
- Scan Resolution (`ScanRes`)
- Index Resolution (`indexRes`)

### Settings Hardcoded in Processing Scripts:
- Wave Speed: `Data and Computation Scripts/ComputeAndTransformStats.m`, Line 27

---

## Notes for New Users

1. **Start with Main.m** - This is your primary configuration file
2. **Use Environment Variables** - Set `NDE_DATA_ROOT` instead of modifying hardcoded paths
3. **Enable Caching** - Keep `SkipOverWriteRequests = 1` to avoid recomputing data
4. **Check Logs** - Enable `AlignDiagnostics = 1` to troubleshoot alignment issues
5. **Visualization Settings** - Most visualization settings can be adjusted in the GUI after running the program

---

**Document Version:** 1.0  
**Last Updated:** 2025-10-15  
**Program Version:** Current Revision

