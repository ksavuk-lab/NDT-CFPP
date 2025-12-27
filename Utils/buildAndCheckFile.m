    % buildAndCheckFile - Builds a file name and checks if it exists
    %
    % Inputs:
    %   fileType      - String indicating the file type (e.g., 'Waveform', 'Plot2D', 'FFT')
    %   FileNamingArray - Array with naming parameters [~, caseNumber, dataType, ...]
    %   folderPath    - Directory where the file is stored (e.g., 'Saved Wave Forms')
    %   fileExtension - File extension (e.g., '.mat', '.fig')
    %
    % Outputs:
    %   filePath     - Full path to the file
    %   fileExists   - Boolean (true if file exists, false otherwise)

    function [filePath, fileExists] = buildAndCheckFile(fileType, FileNamingArray, folderPath, fileExtension)

    % Unpack the FileNamingArray
    [~, caseNumber, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange, alignAtFirstPeak] = unpackFileNamingArray(FileNamingArray);

    % Use Raw data type (current Main.m only uses raw data)
    dataType = 'Raw';

    % Generate range strings based on trimming flags
    if TrimWaveData == 0
        xStr = 'FullX';
        yStr = 'FullY';
    else
        xStr = sprintf('X%.1f-%.1f', xRange(1), xRange(2));
        yStr = sprintf('Y%.1f-%.1f', yRange(1), yRange(2));
    end
    if TrimTimeRange == 0
        tStr = 'FullT';
    else
        tStr = sprintf('T%.1e-%.1e', TimeRange(1), TimeRange(2));
    end

    % Add alignment and zero-shifting info to the file name if enabled
    alignmentStr = '';
    if alignAtFirstPeak == 1
        alignmentStr = 'AlignedAt1stPeak_';
    end



    % Build the file name using a consistent format (simplified for current usage)
    filename = sprintf('%sCase%d_%s_%s%s_%s_%s%s', ...
                       fileType, caseNumber, dataType, alignmentStr, xStr, yStr, tStr, fileExtension);
    filePath = fullfile(folderPath, filename);

    % Check if the file exists
    fileExists = exist(filePath, 'file') == 2;
end