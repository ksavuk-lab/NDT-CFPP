function restoreNormalVisualization(mainFig)
    % RESTORENORMALVISUALIZATION - Restore normal visibility to all points
    % Called when alignment completes to show all data points normally
    
    try
        % Standardized figure validation
        if ~ErrorHandler.validateFigureHandle(mainFig, 'restoreNormalVisualization')
            return;
        end
        
        % Get the current figure data
        userData = get(mainFig, 'UserData');

        % CRITICAL FIX: During alignment, the UserData structure changes
        fprintf('DEBUG: Checking data structure for restoration\n');

        if isstruct(userData) && isfield(userData, 'originalFigData')
            % We're in alignment mode - use the original figData
            figData = userData.originalFigData;
            fprintf('DEBUG: Using originalFigData for restoration\n');
        else
            % Regular mode - userData IS the figData
            figData = userData;
            fprintf('DEBUG: Using regular figData for restoration\n');
        end

        % Get main axes using unified access
        mainAx = DataAccessUtils.getMainAxes(figData);
        if isempty(mainAx)
            fprintf('DEBUG: No valid main axes found for restoration\n');
            return;
        end
        
        % Find all scatter plot objects in the main axes
        scatterObjects = findobj(mainAx, 'Type', 'scatter');
        if isempty(scatterObjects)
            fprintf('DEBUG: No scatter objects found for restore\n');
            return;
        end
        
        fprintf('DEBUG: Restoring normal visualization for %d scatter objects\n', length(scatterObjects));
        
        % Restore full opacity and normal size to all points
        for i = 1:length(scatterObjects)
            scatter_obj = scatterObjects(i);
            
            % Reset transparency to full opacity
            set(scatter_obj, 'MarkerFaceAlpha', 1.0);
            set(scatter_obj, 'MarkerEdgeAlpha', 1.0);
            set(scatter_obj, 'AlphaDataMode', 'auto');
            
            % Reset marker size to default
            xData = get(scatter_obj, 'XData');
            if ~isempty(xData)
                defaultSize = 20; % Default marker size
                set(scatter_obj, 'SizeData', defaultSize);
            end
        end
        
        % Update status
        statusText = findobj(mainFig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', 'Status: Alignment Complete - All data visible');
        end

        % CLEANUP: Remove active slice annotation
        try
            prevAnnotation = findobj(mainAx, 'Tag', 'ActiveSliceAnnotation');
            if ~isempty(prevAnnotation)
                delete(prevAnnotation);
            end
        catch
            % Ignore cleanup errors
        end

        % CLEANUP: Reset figure title
        set(mainFig, 'Name', '3D Peak Visualization - Alignment Complete');

        % Force visual update with multiple refreshes
        drawnow;
        pause(0.01);
        drawnow;

        fprintf('DEBUG: Normal visualization restored with cleanup\n');
        
    catch ME
        ErrorHandler.handleError(ME, 'restoreNormalVisualization', ErrorHandler.WARNING);
    end
end
