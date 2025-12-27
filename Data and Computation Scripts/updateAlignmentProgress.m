function updateAlignmentProgress(mainFigure, alignmentPhase, sliceIndex, numSlices, iteration, columnProgress, overallProgress)
    % UPDATEALIGNMENTPROGRESS - Real-time visual updates during alignment computation
    % Provides continuous feedback showing which slice and column are being processed
    %
    % Inputs:
    %   mainFigure - Handle to the 3D plot figure
    %   alignmentPhase - 'Y-slices' or 'X-slices'
    %   sliceIndex - Current slice being processed
    %   numSlices - Total number of slices
    %   iteration - Current alignment iteration
    %   columnProgress - Progress within current column (0-1)
    %   overallProgress - Overall progress of alignment (0-1)
    
    persistent lastUpdateTime;
    
    try
        % Throttle updates reasonably (max 10 FPS - enough to see changes, won't slow computation)
        currentTime = now;
        if ~isempty(lastUpdateTime) && (currentTime - lastUpdateTime) < (1/10)/86400
            return; % Skip this update
        end
        lastUpdateTime = currentTime;
        
        % Validate figure handle
        if ~ErrorHandler.validateFigureHandle(mainFigure, 'updateAlignmentProgress')
            return;
        end
        
        % Get figure data
        userData = get(mainFigure, 'UserData');
        if isstruct(userData) && isfield(userData, 'originalFigData')
            figData = userData.originalFigData;
        else
            figData = userData;
        end
        
        % Get main axes
        mainAx = DataAccessUtils.getMainAxes(figData);
        if isempty(mainAx)
            return;
        end
        
        % Update the annotation with detailed progress
        try
            prevAnnotation = findobj(mainAx, 'Tag', 'ActiveSliceAnnotation');
            if ~isempty(prevAnnotation)
                % Create detailed progress string
                progressStr = sprintf('ALIGNING: %s %d/%d (Iter %d, %.0f%%)', ...
                    alignmentPhase, sliceIndex, numSlices, iteration, columnProgress*100);
                
                set(prevAnnotation, 'String', progressStr);
                
                % Change color based on progress for visual feedback
                if columnProgress < 0.3
                    bgColor = [1, 1, 0]; % Yellow - starting
                    textColor = [0, 0, 0]; % Black text
                elseif columnProgress < 0.7
                    bgColor = [1, 0.5, 0]; % Orange - in progress
                    textColor = [1, 1, 1]; % White text
                else
                    bgColor = [0, 1, 0]; % Green - nearly done
                    textColor = [0, 0, 0]; % Black text
                end
                
                set(prevAnnotation, 'BackgroundColor', bgColor, 'Color', textColor);
            end
        catch
            % Ignore annotation errors
        end
        
        % Update figure title with detailed progress
        try
            titleStr = sprintf('3D Peak Visualization - %s %d/%d (Iter %d, %.1f%% complete)', ...
                alignmentPhase, sliceIndex, numSlices, iteration, overallProgress*100);
            set(mainFigure, 'Name', titleStr);
        catch
            % Ignore title errors
        end
        
        % Update status text if available
        try
            statusText = findobj(mainFigure, 'Tag', 'AlignmentStatusText');
            if ~isempty(statusText)
                statusStr = sprintf('Status: %s %d/%d - Iteration %d (%.1f%%)', ...
                    alignmentPhase, sliceIndex, numSlices, iteration, columnProgress*100);
                set(statusText, 'String', statusStr);
            end
        catch
            % Ignore status errors
        end
        
        % EFFICIENT UPDATE: Non-blocking refresh that doesn't slow computation
        drawnow limitrate;  % Rate-limited refresh - won't block computation
        
    catch ME
        % Use error handler but don't disrupt alignment
        ErrorHandler.handleError(ME, 'updateAlignmentProgress', ErrorHandler.SILENT);
    end
end
