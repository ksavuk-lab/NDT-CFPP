function PreComputeWaveformData(waveform3DMatrix, waveformEnvelope3DMatrix, t, X_Coordinates, Y_Coordinates, FileNamingArray, varargin)
    % Optional parameters for scaling
    p = inputParser;
    addParameter(p, 'ScaleAmp', false, @islogical);  % Flag to enable scaling
    addParameter(p, 'maxAmp', 1.0, @isnumeric);     % Target maximum amplitude
    parse(p, varargin{:});

    tic;  % Start timing

    % File setup
    folderPath = fullfile(pwd, 'Saved Wave Forms');

    % Check if the directory exists, and create it if it doesn't
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
        fprintf('Created directory: %s\n', folderPath);
    end

    fileType = 'Waveform';
    fileExtension = '.mat';
    [dataFile, fileExists] = buildAndCheckFile(fileType, FileNamingArray, folderPath, fileExtension);

    % Check file existence and prompt user
    if fileExists
        % Check if SkipOverWriteRequests is defined in the base workspace
        skipOverwrite = 0;
        if evalin('base', 'exist(''SkipOverWriteRequests'', ''var'')')
            skipOverwrite = evalin('base', 'SkipOverWriteRequests');
        end

        if skipOverwrite
            % Skip the overwrite without prompting
            fprintf('File already exists. Skipping extraction and saving (SkipOverWriteRequests enabled).\n');
            fprintf('PreComputeWaveformData completed in %.2f seconds.\n', toc);
            return;
        else
            % Prompt the user
            choice = input('File already exists. Overwrite? (y/n): ', 's');
            if ~strcmpi(choice, 'y')
                fprintf('Skipping extraction and saving.\n');
                fprintf('PreComputeWaveformData completed in %.2f seconds.\n', toc);
                return;
            end
        end
    end

    % Unpack metadata and select data source
    [~, caseNumber, ~, ~, xRange, yRange, TimeRange, alignAtFirstPeak] = unpackFileNamingArray(FileNamingArray);

    % Use raw waveform data (Main.m only uses raw data)
    dataToPlot = waveform3DMatrix;

    % Always use raw data
    dataType = 'Raw';

    % Extract all waveforms (no subsampling) - vectorized reshape preserving original row ordering (j fast)
    numY = length(Y_Coordinates);
    numX = length(X_Coordinates);
    totalWaveforms = numY * numX;
    % dataToPlot is [y,x,t]; original code iterated i=Y (outer), j=X (inner, fastest).
    % This creates Y-major order: [all X at Y=1, all X at Y=2, ...]
    % Reshape directly to maintain Y-major order
    waveformArray = reshape(dataToPlot, totalWaveforms, size(dataToPlot, 3));

    % Optional amplitude scaling
    if p.Results.ScaleAmp
        currentMax = max(abs(waveformArray(:)));
        if currentMax > 0
            waveformArray = (waveformArray / currentMax) * p.Results.maxAmp;
        else
            warning('All amplitudes are zero. No scaling applied.');
        end
    end

    % Save data
    savedData.waveformArray = waveformArray;
    savedData.t = t;
    savedData.X_Coordinates = X_Coordinates;
    savedData.Y_Coordinates = Y_Coordinates;

    savedData.dataType = dataType;
    savedData.caseNumber = caseNumber;
    savedData.xRange = xRange;
    savedData.yRange = yRange;
    savedData.TimeRange = TimeRange;
    savedData.numY_sub = numY;
    savedData.numX_sub = numX;
    savedData.alignAtFirstPeak = alignAtFirstPeak;
    if p.Results.ScaleAmp
        savedData.maxAmp = p.Results.maxAmp;
        savedData.originalMax = currentMax;
    end
    save(dataFile, 'savedData', '-v7.3');
    fprintf('Waveforms saved in %s\n', dataFile);

    fprintf('PreComputeWaveformData completed in %.2f seconds.\n', toc);
end