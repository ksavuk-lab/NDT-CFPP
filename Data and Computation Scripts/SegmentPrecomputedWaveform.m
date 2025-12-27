% SEGMENTPRECOMPUTEDWAVEFORM - Segments precomputed waveform data
%
% ## WHAT THIS FUNCTION DOES
% This function provides two ways to segment waveform data:
% 1. From a .mat file: Loads waveform data saved by PreComputeWaveformData and segments it
% 2. Directly from a 3D matrix: Segments a provided waveform matrix
%
% The function segments each waveform into parts based on the specified method, and returns
% the results in a cell array for further analysis.
%
% ## USAGE OPTIONS
% Option 1: SegmentPrecomputedWaveform(filePath, method, interval, numSlices)
%   - Loads data from a file and segments it
%
% Option 2: SegmentPrecomputedWaveform(waveformFull, t, method, param)
%   - Directly segments a provided 3D waveform matrix
%   - This replaces the functionality of SegmentPrecomputedWaveformHelper
%
% ## INPUTS FOR OPTION 1
% 1. filePath   - Path to the .mat file containing precomputed waveform data
% 2. method     - Segmentation method: 'equalSpacing' or 'totalSlices'
% 3. interval   - Time interval in seconds (for 'equalSpacing'; ignored otherwise)
% 4. numSlices  - Number of slices (for 'totalSlices'; ignored otherwise)
%
% ## INPUTS FOR OPTION 2
% 1. waveformFull - 3D matrix of waveforms [y, x, time] or 2D matrix [numWaveforms, time]
% 2. t            - Time vector corresponding to the waveform data
% 3. method       - Segmentation method: 'equalSpacing' or 'totalSlices'
% 4. param        - Parameter value (interval for 'equalSpacing', numSlices for 'totalSlices')
%
% ## OUTPUT
% segmentedData - Cell array where each element is a struct with:
%                 .waveform - Matrix of waveform segments
%                 .time     - Corresponding time vector for the segment

function segmentedData = SegmentPrecomputedWaveform(varargin)
    % Determine which function mode to use based on input arguments
    if nargin < 3 || nargin > 4
        error('SegmentPrecomputedWaveform requires 3 or 4 input arguments');
    end

    % Check if first argument is a text path or a numeric array (waveform data)
    arg1 = varargin{1};
    isStrSupported = exist('string','class') == 8; % robust check for older MATLAB
    isText = ischar(arg1) || (isStrSupported && isa(arg1, 'string'));

    if isText
        % Option 1: Load from file and segment
        segmentedData = segmentFromFile(varargin{:});
    elseif isnumeric(arg1)
        % Option 2: Segment directly from provided waveform data
        segmentedData = segmentFromMatrix(varargin{:});
    else
        error('First argument must be either a file path or a waveform matrix');
    end
end

% Function to segment waveform data loaded from a file
function segmentedData = segmentFromFile(filePath, method, interval, numSlices)
    % Load the precomputed data
    loaded = load(filePath);
    if ~isfield(loaded, 'savedData')
        error('Invalid file: "savedData" structure not found in %s', filePath);
    end
    data = loaded.savedData;
    waveformArray = data.waveformArray;  % [numWaveforms, numTimePoints]
    t = data.t;                          % Time vector

    % Check if we have spatial information to reconstruct 3D structure
    is3D = isfield(data, 'numY_sub') && isfield(data, 'numX_sub');
    if is3D
        ySize = data.numY_sub;
        xSize = data.numX_sub;
    end

    % Get dimensions and time step
    M = size(waveformArray, 2);  % Number of time points
    dt = t(2) - t(1);            % Time step (assumes uniform sampling)

    % Non-blocking progress dialog (fallback to waitbar if not available)
    useProgress = exist('ProgressDialog','file') == 2;
    progress = [];

    % Segment based on method
    if strcmp(method, 'equalSpacing')
        if isempty(interval) || interval <= 0
            error('Interval must be a positive number for equalSpacing method');
        end
        numPointsPerSegment = round(interval / dt);
        if numPointsPerSegment < 1
            error('Interval too small for the sampling rate (dt = %.6f s)', dt);
        end
        numSegments = ceil(M / numPointsPerSegment);
        segmentedData = cell(numSegments, 1);

        % Console-based progress updates every 10%
        lastPctReported = -1;
        fprintf('Segmenting waveforms (Equal Spacing): ');
        % Will print 0%,10%,...,100% as progress advances


        for i = 1:numSegments
            start_idx = (i-1) * numPointsPerSegment + 1;
            end_idx = min(i * numPointsPerSegment, M);
            time_indices = start_idx:end_idx;

            % Store the segmented waveform
            segment_waveform = waveformArray(:, time_indices);

            % If we have spatial information, reshape to 3D
            if is3D
                segmentedData{i}.waveform = reshape(segment_waveform, ySize, xSize, []);
                segmentedData{i}.original_shape = [ySize, xSize, length(time_indices)];
            else
                segmentedData{i}.waveform = segment_waveform;
            end

            segmentedData{i}.time = t(time_indices);
            % Update progress to console every 10%
            pct = floor((i / numSegments) * 10) * 10;
            if pct ~= lastPctReported
                fprintf('%d%% ', pct);
                lastPctReported = pct;
            end
        end
    elseif strcmp(method, 'totalSlices')
        if isempty(numSlices) || numSlices <= 0 || floor(numSlices) ~= numSlices
            error('numSlices must be a positive integer for totalSlices method');
        end
        % Check if requested slices exceed available time points
        if numSlices > M
            originalNumSlices = numSlices;
            numSlices = M;
            fprintf('Requested %d slices, but only %d time points available. Adjusting to %d slices.\n', ...
                    originalNumSlices, M, numSlices);
        end
        edges = round(linspace(1, M+1, numSlices+1));
        segmentedData = cell(numSlices, 1);

        % Console-based progress updates every 10%
        lastPctReported = -1;
        fprintf('Segmenting waveforms (Total Slices): ');
        % Will print 0%,10%,...,100% as progress advances

        for i = 1:numSlices
            time_indices = edges(i):edges(i+1)-1;

            % Store the segmented waveform
            segment_waveform = waveformArray(:, time_indices);

            % If we have spatial information, reshape to 3D
            if is3D
                segmentedData{i}.waveform = reshape(segment_waveform, ySize, xSize, []);
                segmentedData{i}.original_shape = [ySize, xSize, length(time_indices)];
            else
                segmentedData{i}.waveform = segment_waveform;
            end

            segmentedData{i}.time = t(time_indices);
            % Update progress to console every 10%
            pct = floor((i / numSlices) * 10) * 10;
            if pct ~= lastPctReported
                fprintf('%d%% ', pct);
                lastPctReported = pct;
            end
        end
    else
        error('Invalid method: use ''equalSpacing'' or ''totalSlices''');
    end

    % Finalize console progress
    fprintf('100%%\n');
end

% Function to segment waveform data directly from a provided matrix
% This replaces the functionality of SegmentPrecomputedWaveformHelper
function segmentedData = segmentFromMatrix(waveformFull, t, method, param)
    % Check dimensions of waveformFull
    dims = size(waveformFull);

    % Store original dimensions for later use
    is3D = (length(dims) == 3);

    % Handle both 3D [y, x, time] and 2D [numWaveforms, time] inputs
    if is3D
        % 3D matrix [y, x, time]
        ySize = dims(1);
        xSize = dims(2);
        tSize = dims(3);

        % Check for empty data
        if tSize == 0 || ySize == 0 || xSize == 0
            error('waveformFull is empty. Cannot segment.');
        end

        % Reshape to 2D [numWaveforms, time] for processing
        % We'll reshape back to 3D later
        waveformArray = reshape(waveformFull, ySize * xSize, tSize);
    elseif length(dims) == 2
        % Already in 2D format [numWaveforms, time]
        waveformArray = waveformFull;
        tSize = dims(2);

        % Check for empty data
        if tSize == 0 || dims(1) == 0
            error('waveformFull is empty. Cannot segment.');
        end
    else
        error('waveformFull must be either a 3D matrix [y, x, time] or a 2D matrix [numWaveforms, time]');
    end

    % Process based on segmentation method
    if strcmp(method, 'equalSpacing')
        interval = param;
        if isempty(interval) || interval <= 0
            error('interval must be a positive number for equalSpacing method');
        end

        % Calculate time step
        dt = t(2) - t(1);  % Time step (assumes uniform sampling)

        % Calculate number of points per segment
        numPointsPerSegment = round(interval / dt);
        if numPointsPerSegment < 1
            error('Interval too small for the sampling rate (dt = %.6f s)', dt);
        end

        % Calculate number of segments
        numSegments = ceil(tSize / numPointsPerSegment);
        segmentedData = cell(numSegments, 1);

        for i = 1:numSegments
            start_idx = (i-1) * numPointsPerSegment + 1;
            end_idx = min(i * numPointsPerSegment, tSize);
            time_indices = start_idx:end_idx;

            % Store the segmented waveform
            segment_waveform = waveformArray(:, time_indices);

            % If original data was 3D, reshape back to preserve spatial structure
            if is3D
                segmentedData{i}.waveform = reshape(segment_waveform, ySize, xSize, []);
                segmentedData{i}.original_shape = [ySize, xSize, length(time_indices)];
            else
                segmentedData{i}.waveform = segment_waveform;
            end

            segmentedData{i}.time = t(time_indices);
        end
    elseif strcmp(method, 'totalSlices')
        numSlices = param;
        if isempty(numSlices) || numSlices <= 0 || floor(numSlices) ~= numSlices
            error('numSlices must be a positive integer for totalSlices method');
        end

        % Check if requested slices exceed available time points
        if numSlices > tSize
            originalNumSlices = numSlices;
            numSlices = tSize;
            fprintf('Requested %d slices, but only %d time points available. Adjusting to %d slices.\n', ...
                    originalNumSlices, tSize, numSlices);
        end

        % Calculate edges for each slice
        edges = round(linspace(1, tSize+1, numSlices+1));
        segmentedData = cell(numSlices, 1);

        for i = 1:numSlices
            idx = edges(i):edges(i+1)-1;
            if isempty(idx)
                % Handle empty segments (can happen with very small time ranges)
                if is3D
                    segmentedData{i}.waveform = nan(ySize, xSize, 1); % Fallback for empty segment
                    segmentedData{i}.original_shape = [ySize, xSize, 1];
                else
                    segmentedData{i}.waveform = nan(size(waveformArray, 1), 1); % Fallback for empty segment
                end
                segmentedData{i}.time = t(1); % Single time point
            else
                % Store the segmented waveform
                segment_waveform = waveformArray(:, idx);

                % If original data was 3D, reshape back to preserve spatial structure
                if is3D
                    segmentedData{i}.waveform = reshape(segment_waveform, ySize, xSize, []);
                    segmentedData{i}.original_shape = [ySize, xSize, length(idx)];
                else
                    segmentedData{i}.waveform = segment_waveform;
                end
                segmentedData{i}.time = t(idx);
            end
        end
    else
        error('Invalid method: use ''equalSpacing'' or ''totalSlices''');
    end
end