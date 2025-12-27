% SEGMENTFULLWAVEFORM - Extracted from Data and Computation Scripts/ComputeAndTransformStats.m
%
% This function was automatically extracted from a nested function.

function segmentedData = segmentFullWaveform(waveformFull, t, numSlices)
    % Use the updated SegmentPrecomputedWaveform function instead of reimplementing
    % This ensures consistent behavior between the two segmentation methods
    segmentedData = SegmentPrecomputedWaveform(waveformFull, t, 'totalSlices', numSlices);
end
