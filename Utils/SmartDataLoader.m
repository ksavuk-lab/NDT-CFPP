function [waveformArray, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub] = SmartDataLoader(FileNamingArray)
    % SMARTDATALOADER - Load and process waveform data with smart caching
    %
    % This function loads waveform data using the standard data loading pipeline
    % and converts it to a 2D array format suitable for peak analysis.
    %
    % Inputs:
    %   FileNamingArray - Array with naming parameters [DATASET, caseNumber, ...]
    %
    % Outputs:
    %   waveformArray   - 2D array [numWaveforms, numTimePoints]
    %   t               - Time vector
    %   X_Coordinates   - X spatial coordinates
    %   Y_Coordinates   - Y spatial coordinates
    %   numY_sub        - Number of Y points (after trimming)
    %   numX_sub        - Number of X points (after trimming)

    % Handle cache clearing request
    if isempty(FileNamingArray)
        % Clear persistent variables
        clear cachedData cachedFileNamingArray
        waveformArray = [];
        t = [];
        X_Coordinates = [];
        Y_Coordinates = [];
        numY_sub = 0;
        numX_sub = 0;
        return;
    end

    % Unpack the FileNamingArray to get parameters
    [DATASET, caseNumber, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange, alignAtFirstPeak] = unpackFileNamingArray(FileNamingArray);

    % Check for cached waveform data to avoid reloading
    persistent cachedData cachedFileNamingArray

    % Check if we can use cached data
    if ~isempty(cachedData) && ~isempty(cachedFileNamingArray) && ...
       isequal(FileNamingArray, cachedFileNamingArray)
        fprintf('SmartDataLoader: Using cached waveform data for Dataset %d...\n', DATASET);
        waveformArray = cachedData.waveformArray;
        t = cachedData.t;
        X_Coordinates = cachedData.X_Coordinates;
        Y_Coordinates = cachedData.Y_Coordinates;
        numY_sub = cachedData.numY_sub;
        numX_sub = cachedData.numX_sub;
        return;
    end

    fprintf('SmartDataLoader: Loading Dataset %d...\n', DATASET);
    
    % Step 1: Load raw data using the standard loadData function
    try
        [DataStructure] = loadData(DATASET);
        fprintf('Raw data loaded successfully.\n');
    catch ME
        error('Failed to load dataset %d: %s', DATASET, ME.message);
    end
    
    % Step 2: Apply alignment if requested (with caching)
    if alignAtFirstPeak == 1
        fprintf('Applying waveform alignment...\n');

        % Get raw data for alignment
        waveform3DMatrix = DataStructure.waveform3DMatrix;

        % Generate time vector for alignment
        SampRate = DataStructure.sampRate * 1e6;  % Convert to Hz
        dt = 1 / SampRate;
        tSize = DataStructure.tSize;
        t_full = dt:dt:(tSize * dt);

        % Check for cached alignment
        alignmentCacheKey = struct('dataset', DATASET, 'timeRange', [], ...
            'peakType', 'positive', 'dataSize', size(waveform3DMatrix));
        alignmentHash = generateAlignmentCacheHash(alignmentCacheKey);

        alignmentCacheFolder = fullfile(pwd, 'Alignment Cache');
        if ~exist(alignmentCacheFolder, 'dir'), mkdir(alignmentCacheFolder); end
        alignmentCacheFile = fullfile(alignmentCacheFolder, ...
            sprintf('AlignmentCache_Dataset%d_%s.mat', DATASET, alignmentHash));

        if exist(alignmentCacheFile, 'file')
            % Load and apply cached alignment
            fprintf('Loading cached alignment: %s\n', alignmentCacheFile);
            cachedAlign = load(alignmentCacheFile);
            waveform3DMatrix = applyCachedAlignment(waveform3DMatrix, cachedAlign.shiftIndices);
            DataStructure.waveform3DMatrix = waveform3DMatrix;
            fprintf('Cached alignment applied.\n');
        elseif exist('AlignWaveformsAtFirstPeak', 'file')
            try
                % Compute alignment fresh
                [waveform3DMatrix, shiftIndices] = AlignWaveformsAtFirstPeak(waveform3DMatrix, t_full, 0);
                DataStructure.waveform3DMatrix = waveform3DMatrix;
                fprintf('Waveform alignment completed.\n');

                % Save to cache
                fprintf('Saving alignment cache: %s\n', alignmentCacheFile);
                save(alignmentCacheFile, 'shiftIndices', 'alignmentCacheKey', '-v7.3');
            catch ME
                fprintf('Warning: Alignment failed (%s). Continuing without alignment.\n', ME.message);
            end
        else
            fprintf('Warning: AlignWaveformsAtFirstPeak function not found. Skipping alignment.\n');
        end
    end
    
    % Step 3: Extract and trim data using the standard pipeline
    try
        fprintf('Extracting and trimming data...\n');

        % First, let's check the original data time range
        SampRate_orig = DataStructure.sampRate * 1e6;  % Convert to Hz
        dt_orig = 1 / SampRate_orig;
        tSize_orig = DataStructure.tSize;
        t_orig = dt_orig:dt_orig:(tSize_orig * dt_orig);

        fprintf('Original data time range: %.2e to %.2e seconds (%d points)\n', ...
                t_orig(1), t_orig(end), length(t_orig));
        fprintf('Requested TimeRange: [%.2e, %.2e] seconds\n', TimeRange(1), TimeRange(2));

        [SampRate, ScanRes, indexRes, xSize, ySize, tSize, waveform3DMatrix, ~, t, X_Coordinates, Y_Coordinates] = ...
            ExtractAndTrimData(DataStructure, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange);
        fprintf('Data extraction completed.\n');
    catch ME
        error('Failed to extract and trim data: %s', ME.message);
    end
    
    % Step 4: Convert 3D matrix to 2D waveform array
    fprintf('Converting 3D matrix to 2D waveform array...\n');
    
    % Get dimensions after trimming
    numY_sub = length(Y_Coordinates);
    numX_sub = length(X_Coordinates);
    numTimePoints = length(t);
    
    % Initialize waveform array
    totalWaveforms = numY_sub * numX_sub;
    waveformArray = zeros(totalWaveforms, numTimePoints);
    
    % Convert 3D to 2D by flattening spatial dimensions (vectorized) preserving j-fast ordering
    waveformArray = reshape(permute(waveform3DMatrix, [2 1 3]), totalWaveforms, numTimePoints);

    fprintf('SmartDataLoader: Successfully loaded %d waveforms with %d time points each.\n', ...
            totalWaveforms, numTimePoints);
    fprintf('Spatial dimensions: %d x %d (Y x X)\n', numY_sub, numX_sub);

    % Check if time vector is valid before printing
    if ~isempty(t) && length(t) > 0
        fprintf('Time range: %.2e to %.2e seconds\n', t(1), t(end));
    else
        fprintf('Warning: Time vector is empty! Check TimeRange parameter.\n');
    end

    if ~isempty(X_Coordinates) && ~isempty(Y_Coordinates)
        fprintf('Spatial ranges: X=[%.1f, %.1f] mm, Y=[%.1f, %.1f] mm\n', ...
                X_Coordinates(1), X_Coordinates(end), Y_Coordinates(1), Y_Coordinates(end));
    end

    % Error check for empty time data
    if isempty(t) || length(t) == 0
        error('Time vector is empty. This usually means the TimeRange parameter [%.2e, %.2e] is outside the actual data time range. Try using a larger time range or set TrimTimeRange = 0.', ...
              TimeRange(1), TimeRange(2));
    end

    % Cache the loaded data for future use
    cachedData = struct();
    cachedData.waveformArray = waveformArray;
    cachedData.t = t;
    cachedData.X_Coordinates = X_Coordinates;
    cachedData.Y_Coordinates = Y_Coordinates;
    cachedData.numY_sub = numY_sub;
    cachedData.numX_sub = numX_sub;
    cachedFileNamingArray = FileNamingArray;

    fprintf('SmartDataLoader: Data cached for future use.\n');
end

function clearSmartDataLoaderCache()
    % CLEARSMARTDATALOADERCACHE - Clear the persistent cache in SmartDataLoader
    % Call this function if you need to force reloading of data

    % Clear persistent variables by calling with empty inputs
    SmartDataLoader([]);
    fprintf('SmartDataLoader cache cleared.\n');
end
