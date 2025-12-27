%% Helper Function %%
% Helps other functions unpack naming variables and keeping code clean.

function [DATASET, caseNumber, TrimWaveData, TrimTimeRange, xRange, yRange, TimeRange, alignAtFirstPeak] = unpackFileNamingArray(FileNamingArray)
    DATASET = FileNamingArray(1);
    caseNumber = FileNamingArray(2);
    TrimWaveData = FileNamingArray(3);
    TrimTimeRange = FileNamingArray(4);
    xRange = FileNamingArray(5:6);       % 2-element array
    yRange = FileNamingArray(7:8);       % 2-element array
    TimeRange = FileNamingArray(9:10);   % 2-element array

    % Check if alignment flag is included
    if length(FileNamingArray) >= 11
        alignAtFirstPeak = FileNamingArray(11);
    else
        alignAtFirstPeak = 0; % Default to no alignment
    end
end