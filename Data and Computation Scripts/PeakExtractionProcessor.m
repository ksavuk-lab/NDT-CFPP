function peakData = PeakExtractionProcessor(...
    waveformArray, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, ...
    peakDetectionType, peakOptions, viz3DOptions)
    % PEAKEXTRACTIONPROCESSOR - Main function for peak extraction and 3D visualization
    %
    % This function handles the peak extraction workflow:
    % 1. Extract peaks/valleys from waveforms
    % 2. Generate 3D visualization if enabled
    %
    % Inputs:
    %   waveformArray       - 2D array of waveform data [numWaveforms x numTimePoints]
    %   t                   - Time vector
    %   X_Coordinates       - X spatial coordinates
    %   Y_Coordinates       - Y spatial coordinates
    %   numY_sub            - Number of Y points
    %   numX_sub            - Number of X points
    %   peakDetectionType   - 'peaks', 'valleys', or 'both'
    %   peakOptions         - Structure with peak detection parameters
    %   viz3DOptions        - Structure with 3D visualization options
    %
    % Outputs:
    %   peakData            - Cell array of peak data for each waveform

    fprintf('\n=== Peak Extraction and 3D Visualization ===\n');
    fprintf('Peak detection type: %s\n', peakDetectionType);

    %% Step 1: Peak Extraction
    fprintf('\nExtracting peaks from waveforms...\n');

    % SmartDataLoader guarantees 2D format: [numWaveforms, numTimePoints]
    [numWaveforms, numTimePoints] = size(waveformArray);
    fprintf('Processing %d waveforms with %d time points each\n', numWaveforms, numTimePoints);

    % Validate dimensions match
    if numTimePoints ~= length(t)
        error('Time dimension mismatch: waveform has %d time points but time vector has %d points', numTimePoints, length(t));
    end

    peakData = cell(numWaveforms, 1);
    peakCounts = zeros(numWaveforms, 1);
    valleyCounts = zeros(numWaveforms, 1);
    totalTransitionCounts = zeros(numWaveforms, 1);
    
    % Set what to find based on user selection
    switch lower(peakDetectionType)
        case 'peaks'
            peakOptions.findValleys = false;
            fprintf('Detecting peaks only...\n');
        case 'valleys'
            peakOptions.findValleys = true;
            peakOptions.findPeaks = false;
            fprintf('Detecting valleys only...\n');
        case 'both'
            peakOptions.findValleys = true;
            fprintf('Detecting both peaks and valleys...\n');
        otherwise
            error('Invalid PeakDetectionType. Use ''peaks'', ''valleys'', or ''both''.');
    end
    
    % Process each waveform
    fprintf('Processing waveforms: ');
    progressStep = max(1, floor(numWaveforms/10));
    
    for i = 1:numWaveforms
        % Show progress
        if mod(i, progressStep) == 0
            fprintf('%d%% ', round(i/numWaveforms*100));
        end
        
        % Extract peaks and valleys for this waveform
        peakData{i} = detectPeaksAndValleys(waveformArray(i,:), t, peakOptions);
        
        % Count peaks and valleys
        if ~isempty(peakData{i})
            if strcmp(lower(peakDetectionType), 'valleys')
                % For valleys only, count valleys as peaks
                peakCounts(i) = sum(peakData{i}.TransitionType == -1);
                valleyCounts(i) = 0;
            else
                peakCounts(i) = sum(peakData{i}.TransitionType == 1);
                valleyCounts(i) = sum(peakData{i}.TransitionType == -1);
            end
            totalTransitionCounts(i) = height(peakData{i});
        end
    end
    
    fprintf('\nPeak extraction complete.\n');
    
    % Calculate statistics
    fprintf('\n=== Peak Extraction Statistics ===\n');
    if strcmp(lower(peakDetectionType), 'valleys')
        fprintf('Average number of valleys per waveform: %.2f\n', mean(peakCounts));
        fprintf('Min valleys: %d, Max valleys: %d\n', min(peakCounts), max(peakCounts));
    else
        fprintf('Average number of peaks per waveform: %.2f\n', mean(peakCounts));
        if strcmp(lower(peakDetectionType), 'both')
            fprintf('Average number of valleys per waveform: %.2f\n', mean(valleyCounts));
            fprintf('Average total transitions per waveform: %.2f\n', mean(totalTransitionCounts));
        end
        fprintf('Min peaks: %d, Max peaks: %d\n', min(peakCounts), max(peakCounts));
    end

    %% Step 2: 3D Visualization
    if viz3DOptions.Enable3DPeakVallyPlotting == 1
        fprintf('Creating 3D peak/valley representations...\n');
        create3DPeakPlots([], [], peakData, waveformArray, t, ...
            X_Coordinates, Y_Coordinates, numY_sub, numX_sub, peakDetectionType, ...
            viz3DOptions.align3DPlotsWithTime, viz3DOptions.show3DPeaksAndValleys);
    else
        fprintf('3D peak/valley plotting disabled (Enable3DPeakVallyPlotting = 0).\n');
    end
    
    fprintf('\nPeak Extraction and 3D Visualization completed successfully.\n');
end

%% Helper Functions

function transitionData = detectPeaksAndValleys(waveform, t, options)
    % DETECTPEAKSANDVALLEYS - Detect peaks and valleys in a waveform using advanced methods

    % Validate input dimensions
    if length(waveform) ~= length(t)
        error('Waveform and time vector must have the same length. Waveform: %d, Time: %d', length(waveform), length(t));
    end

    % Set default options if fields are missing
    if ~isfield(options, 'minPeakHeight')
        options.minPeakHeight = 0.1;
    end
    if ~isfield(options, 'minPeakDistance')
        options.minPeakDistance = 5;
    end
    if ~isfield(options, 'minPeakProminence')
        options.minPeakProminence = 0.1;
    end
    if ~isfield(options, 'useSlopeDetection')
        options.useSlopeDetection = true;
    end
    if ~isfield(options, 'slopeThreshold')
        options.slopeThreshold = 0.1;
    end
    if ~isfield(options, 'findValleys')
        options.findValleys = true;
    end

    % Ensure minPeakDistance is at least 1 sample
    if options.minPeakDistance < 1
        options.minPeakDistance = 1;
    end

    % Initialize cell arrays to collect data, then convert to arrays
    transitionData_temp = struct();
    transitionData_temp.indices = {};
    transitionData_temp.amplitudes = {};
    transitionData_temp.types = {};
    transitionData_temp.slopesBefore = {};
    transitionData_temp.slopesAfter = {};
    dataCount = 0;

    % Calculate absolute maximum for thresholds
    absMax = max(abs(waveform));

    % Calculate the slope of the waveform
    dt = t(2) - t(1); % Time step
    slope = diff(waveform) / dt;
    slope = [slope(1), slope]; % Pad to match original length

    % Find peaks using MATLAB's findpeaks function
    minHeightThreshold = options.minPeakHeight * absMax;
    minProminence = options.minPeakProminence * absMax;

    % Find positive peaks
    [peakAmplitudes, peakIndices] = findpeaks(waveform, ...
        'MinPeakDistance', options.minPeakDistance, ...
        'MinPeakProminence', minProminence);

    % Filter peaks by relative height and bounds check
    if ~isempty(peakAmplitudes)
        validPeaks = abs(peakAmplitudes) >= minHeightThreshold;
        peakAmplitudes = peakAmplitudes(validPeaks);
        peakIndices = peakIndices(validPeaks);

        % Ensure indices are within waveform bounds
        validBounds = peakIndices >= 1 & peakIndices <= length(waveform);
        peakAmplitudes = peakAmplitudes(validBounds);
        peakIndices = peakIndices(validBounds);
    end

    % Store peaks
    numPeaks = length(peakIndices);
    if numPeaks > 0
        dataCount = dataCount + 1;
        transitionData_temp.indices{dataCount} = peakIndices;
        transitionData_temp.amplitudes{dataCount} = peakAmplitudes;
        transitionData_temp.types{dataCount} = ones(1, numPeaks);

        % Calculate slopes before and after peaks
        peakSlopesBefore = zeros(1, numPeaks);
        peakSlopesAfter = zeros(1, numPeaks);
        for i = 1:numPeaks
            idx = peakIndices(i);
            if idx > 1 && idx < length(waveform)
                peakSlopesBefore(i) = slope(idx-1);
                peakSlopesAfter(i) = slope(idx+1);
            else
                peakSlopesBefore(i) = 0;
                peakSlopesAfter(i) = 0;
            end
        end
        transitionData_temp.slopesBefore{dataCount} = peakSlopesBefore;
        transitionData_temp.slopesAfter{dataCount} = peakSlopesAfter;
    end

    % Find valleys if requested
    if options.findValleys
        [valleyAmplitudes, valleyIndices] = findpeaks(-waveform, ...
            'MinPeakDistance', options.minPeakDistance, ...
            'MinPeakProminence', minProminence);

        % Convert valley amplitudes back to original scale
        valleyAmplitudes = -valleyAmplitudes;

        % Filter valleys by relative height and bounds check
        if ~isempty(valleyAmplitudes)
            validValleys = abs(valleyAmplitudes) >= minHeightThreshold;
            valleyAmplitudes = valleyAmplitudes(validValleys);
            valleyIndices = valleyIndices(validValleys);

            % Ensure indices are within waveform bounds
            validBounds = valleyIndices >= 1 & valleyIndices <= length(waveform);
            valleyAmplitudes = valleyAmplitudes(validBounds);
            valleyIndices = valleyIndices(validBounds);
        end

        % Store valleys
        numValleys = length(valleyIndices);
        if numValleys > 0
            dataCount = dataCount + 1;
            transitionData_temp.indices{dataCount} = valleyIndices;
            transitionData_temp.amplitudes{dataCount} = valleyAmplitudes;
            transitionData_temp.types{dataCount} = -ones(1, numValleys);

            % Calculate slopes before and after valleys
            valleySlopesBefore = zeros(1, numValleys);
            valleySlopesAfter = zeros(1, numValleys);
            for i = 1:numValleys
                idx = valleyIndices(i);
                if idx > 1 && idx < length(waveform)
                    valleySlopesBefore(i) = slope(idx-1);
                    valleySlopesAfter(i) = slope(idx+1);
                else
                    valleySlopesBefore(i) = 0;
                    valleySlopesAfter(i) = 0;
                end
            end
            transitionData_temp.slopesBefore{dataCount} = valleySlopesBefore;
            transitionData_temp.slopesAfter{dataCount} = valleySlopesAfter;
        end
    end

    % Check if any transitions were found
    if dataCount == 0
        % Create an empty table with the correct variable names
        transitionData = table([], [], [], [], [], [], ...
            'VariableNames', {'TransitionNumber', 'TransitionTime', 'TransitionAmplitude', ...
            'TransitionType', 'SlopeBeforeTransition', 'SlopeAfterTransition'});
        return;
    end

    % Combine all collected data - pre-allocate to avoid dynamic arrays
    % First count total elements
    totalElements = 0;
    for i = 1:dataCount
        totalElements = totalElements + length(transitionData_temp.indices{i});
    end

    % Pre-allocate arrays
    allTransitionIndices = zeros(1, totalElements);
    allTransitionAmplitudes = zeros(1, totalElements);
    transitionTypes = zeros(1, totalElements);
    slopesBefore = zeros(1, totalElements);
    slopesAfter = zeros(1, totalElements);

    % Fill arrays using indexing instead of concatenation
    currentIdx = 0;
    for i = 1:dataCount
        numElements = length(transitionData_temp.indices{i});
        if numElements > 0
            indices = currentIdx+1:currentIdx+numElements;
            allTransitionIndices(indices) = transitionData_temp.indices{i};
            allTransitionAmplitudes(indices) = transitionData_temp.amplitudes{i};
            transitionTypes(indices) = transitionData_temp.types{i};
            slopesBefore(indices) = transitionData_temp.slopesBefore{i};
            slopesAfter(indices) = transitionData_temp.slopesAfter{i};
            currentIdx = currentIdx + numElements;
        end
    end

    % Sort transitions by their position in the waveform
    [sortedIndices, sortOrder] = sort(allTransitionIndices);
    sortedAmplitudes = allTransitionAmplitudes(sortOrder);
    sortedTypes = transitionTypes(sortOrder);
    sortedSlopesBefore = slopesBefore(sortOrder);
    sortedSlopesAfter = slopesAfter(sortOrder);

    % Validate indices are within bounds of time vector
    validIndices = sortedIndices >= 1 & sortedIndices <= length(t);
    if ~all(validIndices)
        fprintf('Warning: Found %d peak indices outside time vector bounds. Removing invalid peaks.\n', sum(~validIndices));
        sortedIndices = sortedIndices(validIndices);
        sortedAmplitudes = sortedAmplitudes(validIndices);
        sortedTypes = sortedTypes(validIndices);
        sortedSlopesBefore = sortedSlopesBefore(validIndices);
        sortedSlopesAfter = sortedSlopesAfter(validIndices);
    end

    % Check if any valid transitions remain after filtering
    if isempty(sortedIndices)
        % Create an empty table with the correct variable names
        transitionData = table([], [], [], [], [], [], ...
            'VariableNames', {'TransitionNumber', 'TransitionTime', 'TransitionAmplitude', ...
            'TransitionType', 'SlopeBeforeTransition', 'SlopeAfterTransition'});
        return;
    end

    % Convert indices to time
    transitionTimes = t(sortedIndices);

    % Create a table with the results
    transitionData = table((1:length(sortedIndices))', transitionTimes(:), sortedAmplitudes(:), ...
        sortedTypes(:), sortedSlopesBefore(:), sortedSlopesAfter(:), ...
        'VariableNames', {'TransitionNumber', 'TransitionTime', 'TransitionAmplitude', ...
        'TransitionType', 'SlopeBeforeTransition', 'SlopeAfterTransition'});
end
