function [SampRate, ScanRes, indexRes, xSize, ySize, tSize, waveform3DMatrix, waveformEnvelope3DMatrix, t, X_Coordinates, Y_Coordinates] = ExtractAndTrimData(DataStructure, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange)
    % Start timing
    tic;
    disp('Start Extract and Trim Data');
    
    % Extract fields from DataStructure
    SampRate = DataStructure.sampRate * 1e6;  % Convert to Hz
    ScanRes = DataStructure.scanRes;
    indexRes = DataStructure.indexRes;
    xSize = DataStructure.xSize;
    ySize = DataStructure.ySize;
    tSize = DataStructure.tSize;
    waveform3DMatrix = DataStructure.waveform3DMatrix;
    
    % Handle optional envelope data
    if isfield(DataStructure, 'waveformEnvelope3DMatrix')
        waveformEnvelope3DMatrix = DataStructure.waveformEnvelope3DMatrix;
    else
        waveformEnvelope3DMatrix = [];
        warning('Envelope data not found in DataStructure!');
    end
    
    % Generate time vector
    dt = 1 / SampRate;
    t = dt:dt:(tSize * dt);
    
    % Generate spatial coordinates
    X_Coordinates = (0:xSize-1) * ScanRes;
    Y_Coordinates = (0:ySize-1) * ScanRes;
    
    % Trim spatial data if enabled
    if TrimWaveData == 1
        xIndices = find(X_Coordinates >= xRange(1) & X_Coordinates <= xRange(2));
        yIndices = find(Y_Coordinates >= yRange(1) & Y_Coordinates <= yRange(2));
        X_Coordinates = X_Coordinates(xIndices);
        Y_Coordinates = Y_Coordinates(yIndices);
    else
        xIndices = 1:length(X_Coordinates);
        yIndices = 1:length(Y_Coordinates);
    end
    
    % Trim time data if enabled
    if TrimTimeRange == 1
        timeIndices = find(t >= TimeRange(1) & t <= TimeRange(2));
        t = t(timeIndices);
    else
        timeIndices = 1:length(t);
    end
    
    % Trim waveform matrices
    waveform3DMatrix = waveform3DMatrix(yIndices, xIndices, timeIndices);
    if ~isempty(waveformEnvelope3DMatrix)
        waveformEnvelope3DMatrix = waveformEnvelope3DMatrix(yIndices, xIndices, timeIndices);
    end
    
    % Replace entire zero waveforms with NaN
    zeroWaveformMask = all(waveform3DMatrix == 0, 3);  % Check along time dimension
    waveform3DMatrix(repmat(zeroWaveformMask, [1, 1, size(waveform3DMatrix, 3)])) = NaN;
    
    if ~isempty(waveformEnvelope3DMatrix)
        zeroEnvelopeMask = all(waveformEnvelope3DMatrix == 0, 3);
        waveformEnvelope3DMatrix(repmat(zeroEnvelopeMask, [1, 1, size(waveformEnvelope3DMatrix, 3)])) = NaN;
    end
    
    % Report completion time
    disp(['Extract and Trim Data: ', num2str(toc), ' seconds']);
end
