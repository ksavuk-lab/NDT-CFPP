function MasterSave(saveType, data, FileNamingArray, varargin)
    % MASTERSAVE - Centralized saving functionality for plots and data
    %
    % Inputs:
    %   saveType - Type of save operation ('StatData', 'SegmentedWaveform', 'StatPlot')
    %   data - Data to save
    %   FileNamingArray - Array with naming parameters
    %   varargin - Additional parameters specific to save type
    
    try
        switch saveType
            case 'StatData'
                % Save statistical data
                statName = varargin{1};
                saveStatisticalData(data, FileNamingArray, statName);
                
            case 'SegmentedWaveform'
                % Save segmented waveform data
                method = varargin{1};
                param = varargin{2};
                saveSegmentedWaveform(data, FileNamingArray, method, param);
                
            case 'StatPlot'
                % Save statistical plot
                plotType = varargin{1};
                savePlotFigure(data, FileNamingArray, plotType);
                
            otherwise
                fprintf('Warning: Unknown save type: %s\n', saveType);
        end
        
    catch ME
        fprintf('Error in MasterSave: %s\n', ME.message);
    end
end

function saveStatisticalData(statData, FileNamingArray, statName)
    % Save statistical data to Statistical Analysis folder
    
    % Unpack FileNamingArray
    [~, caseNumber, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange, ~] = unpackFileNamingArray(FileNamingArray);
    
    % Data type is always Raw (envelope handling is legacy and no longer used)
    dataType = 'Raw';
    
    % Create folder if it doesn't exist
    folderPath = fullfile(pwd, 'Statistical Analysis');
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    
    % Build filename (match XtVsYPlot expected format)
    xStr = sprintf('X%.0f-%.0f', xRange(1), xRange(2));
    yStr = sprintf('Y%.0f-%.0f', yRange(1), yRange(2));
    tStr = sprintf('T%.1e-%.1e', TimeRange(1), TimeRange(2));
    
    if isfield(statData, 'method') && isfield(statData, 'param')
        method = statData.method;
        param = statData.param;
        if strcmp(method, 'totalSlices')
            methodStr = sprintf('TotalSlices(%d)', param);
        else
            methodStr = sprintf('%s_%.1f', method, param);
        end
    else
        methodStr = 'TotalSlices(400)'; % Default
    end
    
    fileName = sprintf('%s_%s_%s_GlobalScaling_No_Z-Score_Case%d_%s_%s_%s.mat', ...
        dataType, methodStr, statName, caseNumber, xStr, yStr, tStr);
    
    filePath = fullfile(folderPath, fileName);
    
    % Check if file already exists
    if exist(filePath, 'file')
        fprintf('Statistical data file already exists: %s\n', filePath);
        fprintf('Skipping save for this file.\n');
        return;
    end
    
    % Save the data
    save(filePath, 'statData', '-v7.3');
    fprintf('Statistical data saved: %s\n', filePath);
end

function saveSegmentedWaveform(segmentedData, FileNamingArray, method, param)
    % Save segmented waveform data to Saved Wave Forms folder
    
    % Unpack FileNamingArray
    [~, caseNumber, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange, ~] = unpackFileNamingArray(FileNamingArray);
    
    % Data type is always Raw (envelope handling is legacy and no longer used)
    dataType = 'Raw';
    
    % Create folder if it doesn't exist
    folderPath = fullfile(pwd, 'Saved Wave Forms');
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    
    % Build filename (match XtVsYPlot expected format)
    xStr = sprintf('X%.0f-%.0f', xRange(1), xRange(2));
    yStr = sprintf('Y%.0f-%.0f', yRange(1), yRange(2));
    tStr = sprintf('T%.1e-%.1e', TimeRange(1), TimeRange(2));
    
    if strcmp(method, 'totalSlices')
        methodStr = sprintf('TotalSlices_%d', param);
    else
        methodStr = sprintf('%s_%.1f', method, param);
    end
    
    fileName = sprintf('SegmentedWaveformCase%d_%s_%s_%s_%s_%s.mat', ...
        caseNumber, dataType, xStr, yStr, tStr, methodStr);
    
    filePath = fullfile(folderPath, fileName);
    
    % Check if file already exists
    if exist(filePath, 'file')
        fprintf('Segmented data file already exists: %s\n', filePath);
        fprintf('Skipping save for this file.\n');
        return;
    end
    
    % Save the data
    save(filePath, 'segmentedData', '-v7.3');
    fprintf('Segmented waveform data saved: %s\n', filePath);
end

function savePlotFigure(figHandle, FileNamingArray, plotType)
    % Save plot figure
    
    % Unpack FileNamingArray
    [~, caseNumber, ~, ~, xRange, yRange, TimeRange, ~] = unpackFileNamingArray(FileNamingArray);
    
    % Create folder if it doesn't exist
    folderPath = fullfile(pwd, 'Saved Plots');
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    
    % Build filename (match XtVsYPlot expected format)
    xStr = sprintf('X%.0f-%.0f', xRange(1), xRange(2));
    yStr = sprintf('Y%.0f-%.0f', yRange(1), yRange(2));
    tStr = sprintf('T%.1e-%.1e', TimeRange(1), TimeRange(2));
    
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    fileName = sprintf('%s_Case%d_%s_%s_%s_%s', ...
        plotType, caseNumber, xStr, yStr, tStr, timestamp);
    
    % Save as both .fig and .png (with headless/UI-safe fallbacks)
    figPath = fullfile(folderPath, [fileName '.fig']);
    pngPath = fullfile(folderPath, [fileName '.png']);

    try
        % In headless or when UI components cause save issues, prefer exportgraphics/exportapp
        if ~usejava('desktop')
            ax = findall(figHandle, 'Type', 'axes');
            if ~isempty(ax)
                exportgraphics(ax(1), pngPath);
                fprintf('Plot saved (axes only, headless): %s\n', pngPath);
            else
                fprintf('No axes found to export in headless mode. Skipping PNG export.\n');
            end
            % Skip saving .fig in headless to avoid UI component errors
            return;
        end

        % Interactive mode: try standard save first
        savefig(figHandle, figPath);
        saveas(figHandle, pngPath);
        fprintf('Plot saved: %s\n', figPath);
        fprintf('Plot saved: %s\n', pngPath);
    catch ME
        % Fallbacks: exportapp for uifigure or exportgraphics for axes
        try
            if exist('exportapp','file') == 2
                exportapp(figHandle, pngPath);
                fprintf('Plot saved via exportapp: %s\n', pngPath);
            else
                ax = findall(figHandle, 'Type', 'axes');
                if ~isempty(ax)
                    exportgraphics(ax(1), pngPath);
                    fprintf('Plot saved via exportgraphics: %s\n', pngPath);
                else
                    fprintf('Error saving plot (no axes available): %s\n', ME.message);
                end
            end
        catch ME2
            fprintf('Error saving plot: %s\n', ME2.message);
        end
    end
end
