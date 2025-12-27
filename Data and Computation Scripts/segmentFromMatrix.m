% SEGMENTFROMMATRIX - Extracted from Data and Computation Scripts/SegmentPrecomputedWaveform.m
%
% This function was automatically extracted from a nested function.

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
