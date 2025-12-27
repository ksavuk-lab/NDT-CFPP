% SEGMENTFROMFILE - Extracted from Data and Computation Scripts/SegmentPrecomputedWaveform.m
%
% This function was automatically extracted from a nested function.

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
    % Initialize waitbar
    h = waitbar(0, 'Segmenting waveforms...');
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
            % Update waitbar
            waitbar(i / numSegments, h, sprintf('Segmenting (Equal Spacing): %d of %d', i, numSegments));
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
            % Update waitbar
            waitbar(i / numSlices, h, sprintf('Segmenting (Total Slices): %d of %d', i, numSlices));
        end
    else
        error('Invalid method: use ''equalSpacing'' or ''totalSlices''');
    end
    % Close waitbar
    close(h);
end
% Function to segment waveform data directly from a provided matrix
% This replaces the functionality of SegmentPrecomputedWaveformHelper
