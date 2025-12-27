function alignedWaveform = applyCachedAlignment(waveform3DMatrix, shiftIndices)
    % APPLYCACHEDALIGNMENT - Apply pre-computed shift indices to align waveforms
    %
    % This function applies cached alignment shifts to waveform data, avoiding
    % the need to recompute peak detection when alignment parameters haven't changed.
    %
    % Inputs:
    %   waveform3DMatrix - 3D matrix of waveforms [Y, X, T]
    %   shiftIndices     - 2D matrix of shift amounts [Y, X] (in samples)
    %
    % Output:
    %   alignedWaveform  - 3D matrix of aligned waveforms [Y, X, T]
    
    [numY, numX, numT] = size(waveform3DMatrix);
    alignedWaveform = zeros(numY, numX, numT);
    
    for y = 1:numY
        for x = 1:numX
            shift = shiftIndices(y, x);
            currentWaveform = squeeze(waveform3DMatrix(y, x, :));
            
            if shift > 0
                % Shift waveform left (remove beginning, pad end with zeros)
                alignedWaveform(y, x, :) = [currentWaveform(shift+1:end); zeros(shift, 1)];
            elseif shift < 0
                % Shift waveform right (pad beginning with zeros, remove end)
                alignedWaveform(y, x, :) = [zeros(-shift, 1); currentWaveform(1:end+shift)];
            else
                % No shift needed
                alignedWaveform(y, x, :) = currentWaveform;
            end
        end
    end
    
    fprintf('Applied cached alignment shifts to %d x %d waveforms.\n', numY, numX);
end

