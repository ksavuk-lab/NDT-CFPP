function result = computeSignedMaxAmplitudeVec(data)
% computeSignedMaxAmplitudeVec - Vectorized signed max amplitude across 3rd dim
% data: [Y X T]
% result: [Y X] selecting the sample with maximum absolute amplitude, preserving sign

[numY, numX, numT] = size(data);
if numT == 0
    result = zeros(numY, numX);
    return;
end

[~, idxMax] = max(abs(data), [], 3);
rowIdx = repmat((1:numY)', 1, numX);
colIdx = repmat(1:numX, numY, 1);
linIdx = sub2ind([numY, numX, numT], rowIdx, colIdx, idxMax);
result = data(linIdx);
end

