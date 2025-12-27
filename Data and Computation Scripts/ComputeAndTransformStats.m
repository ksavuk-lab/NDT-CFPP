function ComputeAndTransformStats(waveformFull, t, X_Coordinates, Y_Coordinates, FileNamingArray, statsToCompute, numSlices, transformType, varargin)
    % Inputs:
    %   waveformFull    - Full 3D waveform matrix [y, x, t]
    %   t               - Time vector
    %   X_Coordinates   - Full X coordinates
    %   Y_Coordinates   - Full Y coordinates
    %   FileNamingArray - For saving
    %   statsToCompute  - Cell array of statistics to compute
    %   numSlices       - Number of slices for 'totalSlices' segmentation
    %   transformType   - Transformation type ('raw', 'zscore_global', 'scale_layer', 'zscore_layer')

    % Optional precomputed segmentation to avoid re-segmentation
    segmentedData = [];
    if ~isempty(varargin)
        for k = 1:numel(varargin)
            if iscell(varargin{k})
                segmentedData = varargin{k};
                break;
            end
        end
    end

    % Hardcode method for consistency
    method = 'totalSlices';

    % Hardcode the speed of wave
    Speed_of_Wave = 2968.67; % m/s

    % Determine output folder
    folderPath = fullfile(pwd, 'Statistical Analysis');
    fileExtension = '.mat';

    % Segment the full waveform data using 'totalSlices' (only if not precomputed)
    if isempty(segmentedData)
        segmentedData = segmentFullWaveform(waveformFull, t, numSlices);
    else
        % Validate precomputed segmentation matches requested slices
        if numel(segmentedData) ~= numSlices
            fprintf('Warning: Precomputed segmentation has %d slices, expected %d. Proceeding with provided data.\n', numel(segmentedData), numSlices);
        end
    end

    % Compute statistics on full grid
    numY = size(waveformFull, 1);
    numX = size(waveformFull, 2);
    numSegments = length(segmentedData);

    % Compute TOF and Depth on the full waveform once (with caching)
    tofMap = nan(numY, numX);
    depthMap = nan(numY, numX);
    if any(strcmp(statsToCompute, 'TOF')) || any(strcmp(statsToCompute, 'Depth'))
        % Check for cached TOF/Depth maps
        caseNumber = FileNamingArray(2);
        tofCacheFile = fullfile(folderPath, sprintf('TOFDepthCache_Case%d_%dx%d.mat', caseNumber, numY, numX));

        if exist(tofCacheFile, 'file')
            % Load cached TOF/Depth maps
            fprintf('Loading cached TOF/Depth maps: %s\n', tofCacheFile);
            cachedTOF = load(tofCacheFile);
            tofMap = cachedTOF.tofMap;
            depthMap = cachedTOF.depthMap;
        else
            % Compute TOF/Depth (vectorized)
            fprintf('Computing TOF/Depth maps...\n');
            absWaveform = abs(waveformFull); % [Y, X, T]
            [~, maxIdx] = max(absWaveform, [], 3); % Index of max along time dimension

            % Convert indices to TOF values [s]
            tofMap = t(maxIdx);

            % Depth = (Speed * TOF / 2) * 1000 to convert to mm
            % d = v * t / 2 (divide by 2 for round-trip)
            depthMap = (Speed_of_Wave * tofMap / 2) * 1000; % [mm]

            % Save to cache
            fprintf('Saving TOF/Depth cache: %s\n', tofCacheFile);
            save(tofCacheFile, 'tofMap', 'depthMap', 'Speed_of_Wave', '-v7.3');
        end
    end

    % Compute RelativeTOF per segment (vectorized)
    relativeTofMaps = cell(numSegments, 1);
    if any(strcmp(statsToCompute, 'RelativeTOF'))
        for seg = 1:numSegments
            segmentTime = segmentedData{seg}.time;
            segmentWaveform = segmentedData{seg}.waveform; % [Y, X, T_seg]

            % Vectorized: find max absolute amplitude along time dimension
            absSegment = abs(segmentWaveform);
            [~, maxIdx] = max(absSegment, [], 3); % [Y, X] indices into T_seg

            % Convert indices to TOF values
            segmentTofMap = segmentTime(maxIdx); % [Y, X] in seconds

            % Compute the mean TOF for this segment (excluding NaNs)
            meanTof = mean(segmentTofMap(:), 'omitnan');
            if isnan(meanTof)
                meanTof = 0; % If all values are NaN, set mean to 0
            end
            % Compute RelativeTOF as deviation from the mean
            relativeTofMaps{seg} = segmentTofMap - meanTof; % [s]
        end
    end

    % Helper function for statistic computation - computes ONLY the requested stat
    % (Optimized: previous version computed ALL stats then selected one)
    function statResult = computeStat(data, statType, seg)
        switch statType
            case 'RMS'
                statResult = sqrt(mean(data.^2, 3, 'omitnan'));
            case 'MaxAmplitude'
                % Delegate to Utils/computeSignedMaxAmplitudeVec for performance
                statResult = computeSignedMaxAmplitudeVec(data);
            case 'Variance'
                statResult = var(data, 0, 3, 'omitnan');
            case 'Skewness'
                statResult = skewness(data, 1, 3);
            case 'Kurtosis'
                statResult = kurtosis(data, 1, 3);
            case 'TOF'
                statResult = tofMap; % Precomputed static TOF map
            case 'Depth'
                statResult = depthMap; % Precomputed static Depth map
            case 'RelativeTOF'
                statResult = relativeTofMaps{seg}; % Varies per segment
            otherwise
                error('ComputeAndTransformStats:UnknownStat', ...
                    'Unknown statistic type: %s', statType);
        end
    end

    h = waitbar(0, sprintf('Computing %s stats...', transformType));
    for statIdx = 1:length(statsToCompute)
        statType = statsToCompute{statIdx};

        % Define statName based on transformType
        switch transformType
            case 'raw', statName = statType;
            case 'zscore_global', statName = ['ZScore_' statType];
            case 'scale_layer', statName = [statType '_Scaled'];
            case 'zscore_layer', statName = [statType '_ZScorePerLayer'];
            otherwise, error('Invalid transformType: %s', transformType);
        end

        % Extract case number, X, Y, T ranges from FileNamingArray
        caseNumber = FileNamingArray(2);
        xRange = ['X' num2str(FileNamingArray(5)) '-' num2str(FileNamingArray(6))];
        yRange = ['Y' num2str(FileNamingArray(7)) '-' num2str(FileNamingArray(8))];
        tRange = ['T' sprintf('%.1e', FileNamingArray(9)) '-' sprintf('%.1e', FileNamingArray(10))];

        % Use raw data only (envelope processing has been removed)
        dataType = 'Raw';

        % Determine Z-score and scaling parameters based on transformType
        if strcmp(transformType, 'zscore_global') || strcmp(transformType, 'zscore_layer')
            zScore = 'Yes_Z-Score';
        else
            zScore = 'No_Z-Score';
        end

        if strcmp(transformType, 'scale_layer') || strcmp(transformType, 'zscore_layer')
            scaling = 'PerLayerScaling';
        else
            scaling = 'GlobalScaling';
        end

        % Construct file name according to the desired format:
        % (Raw or Envelope)_TotalSlices(##)_(Statistical Analysis Method)_(PerLayerScaling or GlobalScaling)_(Yes/No_Z-Score)_Case(##)_X(range)_Y(range)_T(range)
        fileName = [dataType '_TotalSlices(' num2str(numSlices) ')_' statType '_' scaling '_' zScore '_Case' num2str(caseNumber) '_' xRange '_' yRange '_' tRange];
        [filePath, fileExists] = buildAndCheckFile(fileName, FileNamingArray, folderPath, fileExtension);

        if fileExists
            fprintf('File for %s already exists: %s. Skipping computation.\n', statName, filePath);
            waitbar(statIdx / length(statsToCompute), h);
            continue;
        end

        % Compute the statistic
        statMaps = cell(numSegments, 1);

        % Compute global statistics if needed
        if strcmp(transformType, 'zscore_global')
            globalStat = computeStat(waveformFull, statType, 1); % seg=1 for global stat
            globalMean = mean(globalStat(:), 'omitnan');
            globalStd = std(globalStat(:), 'omitnan');
        end

        % Compute and transform statistics per segment
        for seg = 1:numSegments
            waveformSegment = segmentedData{seg}.waveform; % [y, x, t_seg]
            statVal = computeStat(waveformSegment, statType, seg); % Pass seg explicitly

            % Safeguard: Ensure map shape matches grid size
            if ~isequal(size(statVal), [numY, numX])
                error('ComputeAndTransformStats:ShapeMismatch', ...
                    'Stat map size mismatch for %s (segment %d): got [%s], expected [%d %d]', ...
                    statType, seg, num2str(size(statVal)), numY, numX);
            end

            % Apply transformation
            switch transformType
                case 'raw'
                    transformedStat = statVal;
                case 'zscore_global'
                    transformedStat = (statVal - globalMean) / globalStd;
                case 'scale_layer'
                    if all(isnan(statVal(:)))
                        transformedStat = nan(size(statVal));
                    else
                        % For MaxAmplitude, use max absolute value for scaling to preserve sign
                        if strcmp(statType, 'MaxAmplitude')
                            maxAbsStat = max(abs(statVal(:)), [], 'omitnan');
                            transformedStat = statVal ./ (maxAbsStat + (maxAbsStat == 0));
                        else
                            % For other statistics, use the original scaling method
                            maxStat = max(statVal(:), [], 'omitnan');
                            transformedStat = statVal ./ (maxStat + (maxStat == 0));
                        end
                    end
                case 'zscore_layer'
                    meanStat = mean(statVal(:), 'omitnan');
                    stdStat = std(statVal(:), 'omitnan');
                    transformedStat = (statVal - meanStat) / (stdStat + (stdStat == 0));
            end

            statMaps{seg} = transformedStat;
        end

        % Prepare and save data
        statData.maps = statMaps;
        statData.X_sub = X_Coordinates;
        statData.Y_sub = Y_Coordinates;
        statData.segmentTimes = cellfun(@(x) x.time, segmentedData, 'UniformOutput', false);
        statData.statType = statName;
        statData.method = method;
        statData.param = numSlices;
        statData.transform_type = transformType; % Add transform_type for consistency with Python
        statData.FileNamingArray = FileNamingArray;
        MasterSave('StatData', statData, FileNamingArray, statName);

        waitbar(statIdx / length(statsToCompute), h);
    end
    close(h);
end

% Helper function for conditional string assignment
function str = ifelse(condition, trueStr, falseStr)
    if condition
        str = trueStr;
    else
        str = falseStr;
    end
end

% Helper function to segment full waveform data using 'totalSlices'
function segmentedData = segmentFullWaveform(waveformFull, t, numSlices)
    % Use the updated SegmentPrecomputedWaveform function instead of reimplementing
    % This ensures consistent behavior between the two segmentation methods
    segmentedData = SegmentPrecomputedWaveform(waveformFull, t, 'totalSlices', numSlices);
end