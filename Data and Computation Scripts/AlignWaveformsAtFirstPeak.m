function [alignedWaveform, shiftIndices] = AlignWaveformsAtFirstPeak(waveform3DMatrix, t, isEnvelope, timeRange, peakType, enableDiagnostics, numDiagnosticSamples)
    % AlignWaveformsAtFirstPeak - Aligns waveforms at their first significant peak or valley
    %
    % This function finds the first significant peak or valley in each waveform within a specified time range
    % and shifts the data to align these features at the same time point (effectively setting a new zero time).
    % This function only processes raw waveforms using square wave peak detection algorithms.
    %
    % Syntax:
    %   [alignedWaveform, shiftIndices] = AlignWaveformsAtFirstPeak(waveform3DMatrix, t, isEnvelope, timeRange, peakType, enableDiagnostics, numDiagnosticSamples)
    %
    % Inputs:
    %   waveform3DMatrix - 3D matrix of waveforms [y, x, time]
    %   t                - Time vector
    %   isEnvelope       - Boolean flag (ignored - always treated as raw data)
    %   timeRange        - Optional [min max] time range to search for peaks (in seconds)
    %                      If not provided, the entire waveform is searched
    %   peakType         - Optional parameter to specify which type of peaks to detect:
    %                      'positive' - Only detect positive peaks (default)
    %                      'negative' - Only detect negative peaks (valleys)
    %                      'both'     - Detect both and use the one with larger magnitude
    %   enableDiagnostics - Optional flag to enable diagnostic logging (default: false)
    %                      When true, logs details about peak selection for a random subset of waveforms
    %   numDiagnosticSamples - Optional number of random waveforms to sample for diagnostics (default: 100)
    %
    % Outputs:
    %   alignedWaveform - 3D matrix of aligned waveforms
    %   shiftIndices    - Matrix of indices by which each waveform was shifted

    % Start timing
    tic;

    % Force envelope flag to false (always process as raw data)
    isEnvelope = false;

    % Handle optional parameters
    if nargin < 4 || isempty(timeRange)
        % If timeRange not provided, search the entire waveform
        timeRangeIndices = [1, length(t)];
        searchRangeMsg = 'searching entire waveform';
    else
        % Find indices corresponding to the specified time range
        timeRangeIndices = [find(t >= timeRange(1), 1, 'first'), find(t <= timeRange(2), 1, 'last')];
        if isempty(timeRangeIndices(1))
            timeRangeIndices(1) = 1;
        end
        if isempty(timeRangeIndices(2))
            timeRangeIndices(2) = length(t);
        end
        searchRangeMsg = ['searching between ', num2str(t(timeRangeIndices(1))*1e6), ' and ', num2str(t(timeRangeIndices(2))*1e6), ' μs'];
    end

    % Handle peak type parameter
    if nargin < 5 || isempty(peakType)
        peakType = 'positive'; % Default to positive peaks only
    end

    % Handle diagnostics parameter
    if nargin < 6 || isempty(enableDiagnostics)
        enableDiagnostics = false; % Default to no diagnostics
    end

    % Handle number of diagnostic samples parameter
    if nargin < 7 || isempty(numDiagnosticSamples)
        numDiagnosticSamples = 100; % Default to 100 samples
    end

    % Display start message based on peak type
    switch lower(peakType)
        case 'positive'
            disp(['Starting waveform alignment at first positive peak (', searchRangeMsg, ')...']);
        case 'negative'
            disp(['Starting waveform alignment at first negative peak (', searchRangeMsg, ')...']);
        case 'both'
            disp(['Starting waveform alignment at first significant peak (positive or negative) (', searchRangeMsg, ')...']);
        otherwise
            error('Invalid peakType parameter. Must be ''positive'', ''negative'', or ''both''.');
    end

    % Initialize diagnostics if enabled
    if enableDiagnostics
        disp(['Diagnostic logging enabled - will sample ', num2str(numDiagnosticSamples), ' random waveforms']);

        % Create Logs folder if it doesn't exist
        logsFolder = fullfile(pwd, 'Logs');
        if ~exist(logsFolder, 'dir')
            mkdir(logsFolder);
            disp(['Created Logs folder: ', logsFolder]);
        end

        % Create a timestamp for the log file
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        logFileName = ['alignment_diagnostics_', timestamp, '.txt'];
        logFilePath = fullfile(logsFolder, logFileName);
        logFile = fopen(logFilePath, 'w');
        fprintf(logFile, 'Alignment Diagnostics Log - %s\n\n', datestr(now));
        fprintf(logFile, 'Peak Type: %s\n', peakType);
        fprintf(logFile, 'Time Range: %s\n\n', searchRangeMsg);

        % Create a formatted table header
        fprintf(logFile, '%-10s | %-15s | %-15s | %-15s | %s\n', 'Position', 'Peak Index', 'Time (μs)', 'Peak Value', 'Context Values');
        fprintf(logFile, '%-10s-|-%-15s-|-%-15s-|-%-15s-|-%s\n', '----------', '---------------', '---------------', '---------------', '--------------------');
    end

    % Get dimensions
    [numY, numX, numT] = size(waveform3DMatrix);
    alignedWaveform = zeros(size(waveform3DMatrix));
    shiftIndices = zeros(numY, numX);

    % Define threshold for peak detection (can be adjusted)
    threshold = 1.0; % Amplitude threshold for considering a peak

    % For diagnostics: randomly select waveforms to log
    if enableDiagnostics
        % Calculate total number of waveforms
        totalWaveforms = numY * numX;
        % Generate random indices (or fewer if there are fewer waveforms)
        numSamples = min(numDiagnosticSamples, totalWaveforms);
        % Create a linear index array and shuffle it
        linearIndices = randperm(totalWaveforms, numSamples);
        % Convert to subscript indices
        [diagY, diagX] = ind2sub([numY, numX], linearIndices);
        % Create a logical matrix to quickly check if a waveform should be logged
        diagMask = false(numY, numX);
        for i = 1:numSamples
            diagMask(diagY(i), diagX(i)) = true;
        end
    end

    % Process each waveform
    for y = 1:numY
        for x = 1:numX
            % Extract current waveform
            currentWaveform = squeeze(waveform3DMatrix(y, x, :));

            % Skip if waveform contains NaN values
            if any(isnan(currentWaveform))
                alignedWaveform(y, x, :) = currentWaveform;
                continue;
            end

            % Find the first peak based on data type within the specified time range
            % Extract the portion of the waveform within the time range
            rangeWaveform = currentWaveform(timeRangeIndices(1):timeRangeIndices(2));

            % Always process as raw data - handle square waves with the specified peak type
            peakIdx = findFirstPeakInRawWaveform(currentWaveform, threshold, timeRangeIndices, peakType);

            % If a peak was found, shift the waveform
            if ~isempty(peakIdx) && peakIdx > 0
                shiftIndices(y, x) = peakIdx;

                % Shift the waveform to align the peak at the beginning
                shiftedWaveform = zeros(size(currentWaveform));

                % Copy the data from the peak onwards
                shiftedWaveform(1:numT-peakIdx+1) = currentWaveform(peakIdx:end);

                % Store the shifted waveform
                alignedWaveform(y, x, :) = shiftedWaveform;

                % Log diagnostic information for selected waveforms
                if enableDiagnostics && diagMask(y, x)
                    % Get context values (5 points before and after the peak)
                    contextStart = max(1, peakIdx - 5);
                    contextEnd = min(length(currentWaveform), peakIdx + 5);
                    contextValues = currentWaveform(contextStart:contextEnd);

                    % Format context values as a string (vectorized - avoids O(n²) string growth)
                    contextStr = sprintf('%.4f ', contextValues);

                    % Log the information in table format
                    fprintf(logFile, '[%3d,%3d]  | %-15d | %-15.2f | %-15.4f | %s\n', ...
                        y, x, peakIdx, t(peakIdx)*1e6, currentWaveform(peakIdx), contextStr);
                end
            else
                % If no peak found, keep the original waveform
                alignedWaveform(y, x, :) = currentWaveform;

                % Log diagnostic information for selected waveforms with no peak
                if enableDiagnostics && diagMask(y, x)
                    fprintf(logFile, '[%3d,%3d]  | %-15s | %-15s | %-15s | %s\n', y, x, 'N/A', 'N/A', 'N/A', 'No peak found - Using original waveform');
                end
            end
        end
    end

    % Report completion time
    disp(['Waveform alignment completed in ', num2str(toc), ' seconds']);

    % Close diagnostic log file if it was opened
    if enableDiagnostics
        fprintf(logFile, '\nAlignment completed in %.2f seconds\n', toc);
        fclose(logFile);
        disp(['Diagnostic information saved to ', logFilePath]);
    end
end

function peakIdx = findFirstPeakInRawWaveform(waveform, threshold, timeRangeIndices, peakType)
    % Helper function to find the first peak or valley in a raw waveform within a specified time range
    % Handles square waves by looking for a rise/fall followed by a change in direction
    % Updated to handle square wave peaks by finding the center point when multiple time values have the same max amplitude
    %
    % Inputs:
    %   waveform         - The waveform to analyze
    %   threshold        - Threshold for peak detection as a fraction of max amplitude
    %   timeRangeIndices - [start, end] indices defining the range to search
    %   peakType         - 'positive', 'negative', or 'both'

    % Initialize
    peakIdx = [];

    % Extract the portion of the waveform within the time range
    rangeStart = timeRangeIndices(1);
    rangeEnd = timeRangeIndices(2);
    rangeWaveform = waveform(rangeStart:rangeEnd);
    numPoints = length(rangeWaveform);

    % Process based on peak type
    switch lower(peakType)
        case 'positive'
            % Find positive peaks
            % Adjust threshold based on the range
            adjustedThreshold = threshold * max(rangeWaveform);

            % Find where the signal first crosses the threshold within the range
            thresholdCrossings = find(rangeWaveform >= adjustedThreshold, 1, 'first');

            if isempty(thresholdCrossings)
                % No threshold crossing found in the range
                % Use the maximum value in the range, but find center for square waves
                maxValue = max(rangeWaveform);
                peakIdx = findSquareWavePeakCenter(rangeWaveform, maxValue, rangeStart - 1);
                return;
            end

            % Start from the first threshold crossing
            startIdx = thresholdCrossings;

            % Look for the peak (where signal starts decreasing after rising)
            for i = startIdx:numPoints-1
                if rangeWaveform(i) > rangeWaveform(i+1)
                    % Found a potential peak

                    % Verify it's a true peak by checking a few points ahead
                    % to ensure it's not just noise
                    verificationWindow = min(5, numPoints-i); % Look ahead up to 5 points
                    if all(rangeWaveform(i) >= rangeWaveform(i+1:i+verificationWindow))
                        peakIdx = i + rangeStart - 1; % Convert to global index
                        break;
                    end
                end
            end

            % If no peak found using the above method, use the maximum value
            if isempty(peakIdx)
                subRangeWaveform = rangeWaveform(startIdx:end);
                maxValue = max(subRangeWaveform);
                localPeakIdx = findSquareWavePeakCenter(subRangeWaveform, maxValue, 0);
                peakIdx = localPeakIdx + startIdx + rangeStart - 2; % Convert to global index
            end

        case 'negative'
            % Find negative peaks (valleys)
            % Adjust threshold based on the range
            adjustedThreshold = threshold * abs(min(rangeWaveform));

            % Find where the signal first crosses the negative threshold within the range
            thresholdCrossings = find(rangeWaveform <= -adjustedThreshold, 1, 'first');

            if isempty(thresholdCrossings)
                % No threshold crossing found in the range
                % Use the minimum value in the range, but find center for square waves
                minValue = min(rangeWaveform);
                peakIdx = findSquareWavePeakCenter(rangeWaveform, minValue, rangeStart - 1);
                return;
            end

            % Start from the first threshold crossing
            startIdx = thresholdCrossings;

            % Look for the valley (where signal starts increasing after falling)
            for i = startIdx:numPoints-1
                if rangeWaveform(i) < rangeWaveform(i+1)
                    % Found a potential valley

                    % Verify it's a true valley by checking a few points ahead
                    % to ensure it's not just noise
                    verificationWindow = min(5, numPoints-i); % Look ahead up to 5 points
                    if all(rangeWaveform(i) <= rangeWaveform(i+1:i+verificationWindow))
                        peakIdx = i + rangeStart - 1; % Convert to global index
                        break;
                    end
                end
            end

            % If no valley found using the above method, use the minimum value
            if isempty(peakIdx)
                subRangeWaveform = rangeWaveform(startIdx:end);
                minValue = min(subRangeWaveform);
                localPeakIdx = findSquareWavePeakCenter(subRangeWaveform, minValue, 0);
                peakIdx = localPeakIdx + startIdx + rangeStart - 2; % Convert to global index
            end

        case 'both'
            % Find both positive and negative peaks and use the one with larger magnitude
            % First find positive peak
            posPeakIdx = findFirstPeakInRawWaveform(waveform, threshold, timeRangeIndices, 'positive');

            % Then find negative peak
            negPeakIdx = findFirstPeakInRawWaveform(waveform, threshold, timeRangeIndices, 'negative');

            % Compare magnitudes and choose the one with larger absolute value
            if ~isempty(posPeakIdx) && ~isempty(negPeakIdx)
                posValue = waveform(posPeakIdx);
                negValue = waveform(negPeakIdx);

                if abs(posValue) >= abs(negValue)
                    peakIdx = posPeakIdx;
                else
                    peakIdx = negPeakIdx;
                end
            elseif ~isempty(posPeakIdx)
                peakIdx = posPeakIdx;
            elseif ~isempty(negPeakIdx)
                peakIdx = negPeakIdx;
            else
                % If neither found, use the point with maximum absolute value
                absRangeWaveform = abs(rangeWaveform);
                maxAbsValue = max(absRangeWaveform);
                peakIdx = findSquareWavePeakCenter(absRangeWaveform, maxAbsValue, rangeStart - 1);
            end
    end
end

function centerIdx = findSquareWavePeakCenter(data, targetValue, startOffset)
    % Helper function to find the center point of a square wave peak/valley
    % When multiple time values have the same max/min amplitude, this finds the center
    %
    % Inputs:
    %   data         - The data array to search in
    %   targetValue  - The target value to find (max or min value)
    %   startOffset  - Offset to add to the returned index (default: 0)
    %
    % Outputs:
    %   centerIdx    - Index of the center point of the square wave peak

    if nargin < 3
        startOffset = 0;
    end

    % Find all indices with the target value
    targetIndices = find(data == targetValue);

    if isempty(targetIndices)
        centerIdx = [];
        return;
    end

    % If only one point, return it
    if length(targetIndices) == 1
        centerIdx = targetIndices(1) + startOffset;
        return;
    end

    % For multiple points, find the center
    % Divide the count by 2 to get the center index
    centerPosition = ceil(length(targetIndices) / 2);
    centerIdx = targetIndices(centerPosition) + startOffset;
end
