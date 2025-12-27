function highlightActiveSlice(mainFig, viewType, sliceIndex, totalSlices, alignmentPhase)
    % HIGHLIGHTACTIVESLICE - Highlight the currently active slice being aligned
    % Make non-active points almost invisible so users can see alignment changes in real-time
    %
    % Inputs:
    %   mainFig - Main 3D figure handle
    %   viewType - 'XtVsY' or 'YtVsX' 
    %   sliceIndex - Current slice being processed (1-based)
    %   totalSlices - Total number of slices
    %   alignmentPhase - 'Y-slices' or 'X-slices' for user feedback
    
    try
        % Standardized figure validation
        if ~ErrorHandler.validateFigureHandle(mainFig, 'highlightActiveSlice')
            return;
        end
        
        % Check if slice highlighting is enabled
        showActiveSliceCheckbox = findobj(mainFig, 'Tag', 'ShowActiveSliceCheckbox');
        if ~isempty(showActiveSliceCheckbox) && ~get(showActiveSliceCheckbox, 'Value')
            fprintf('DEBUG: Slice highlighting disabled by checkbox\n');
            return; % Feature is disabled
        end
        
        fprintf('DEBUG: Highlighting %s slice %d/%d in %s\n', alignmentPhase, sliceIndex, totalSlices, viewType);

        % Get the current figure data
        userData = get(mainFig, 'UserData');

        % CRITICAL FIX: During alignment, the UserData structure changes
        % The original figData is stored in userData.originalFigData
        fprintf('DEBUG: Checking data structure during alignment\n');

        if isstruct(userData) && isfield(userData, 'originalFigData')
            % We're in alignment mode - use the original figData
            figData = userData.originalFigData;
            fprintf('DEBUG: Using originalFigData from alignment mode\n');
        else
            % Regular mode - userData IS the figData
            figData = userData;
            fprintf('DEBUG: Using regular figData mode\n');
        end

        % DEBUG: Show what fields are available in the actual figData
        if isstruct(figData)
            fields = fieldnames(figData);
            fprintf('DEBUG: Available figData fields: %s\n', strjoin(fields, ', '));

            % Check for specific fields we need
            hasAllPeakX = isfield(figData, 'allPeakX');
            hasAllPeakY = isfield(figData, 'allPeakY');
            hasMainAx = isfield(figData, 'mainAx');
            hasXCoords = isfield(figData, 'X_Coordinates');
            hasYCoords = isfield(figData, 'Y_Coordinates');

            fprintf('DEBUG: Key fields - allPeakX:%d, allPeakY:%d, mainAx:%d, X_Coordinates:%d, Y_Coordinates:%d\n', ...
                hasAllPeakX, hasAllPeakY, hasMainAx, hasXCoords, hasYCoords);
        else
            fprintf('DEBUG: figData is not a structure! Type: %s\n', class(figData));
            return;
        end

        % Validate data structure
        if ~DataAccessUtils.validateDataStructure(figData, 'highlightActiveSlice')
            fprintf('DEBUG: Data structure validation failed\n');
            return;
        end

        % Get main axes using unified access
        mainAx = DataAccessUtils.getMainAxes(figData);
        if isempty(mainAx)
            fprintf('DEBUG: No valid main axes found\n');
            return;
        end

        % Get peak coordinates using unified access
        coords = DataAccessUtils.getPeakCoordinates(figData);
        if isempty(coords.X) || isempty(coords.Y)
            fprintf('DEBUG: Peak coordinate data not available\n');
            return;
        end

        allPeakX = coords.X;
        allPeakY = coords.Y;

        % Get grid coordinates using unified access
        gridCoords = DataAccessUtils.getGridCoordinates(figData);
        if isempty(gridCoords.X_Coordinates) || isempty(gridCoords.Y_Coordinates)
            fprintf('DEBUG: Grid coordinate data not available\n');
            return;
        end

        X_Coordinates = gridCoords.X_Coordinates;
        Y_Coordinates = gridCoords.Y_Coordinates;

        % Find all scatter plot objects in the main axes
        scatterObjects = findobj(mainAx, 'Type', 'scatter');
        if isempty(scatterObjects)
            fprintf('DEBUG: No scatter objects found in main axes\n');
            return;
        end

        fprintf('DEBUG: Found %d scatter objects\n', length(scatterObjects));
        
        % Determine which points belong to the active slice
        switch viewType
            case 'XtVsY'
                % Aligning Y-slices: highlight points at specific Y coordinate
                if sliceIndex <= length(Y_Coordinates)
                    activeY = Y_Coordinates(sliceIndex);
                    tolerance = (max(Y_Coordinates) - min(Y_Coordinates)) / (2 * length(Y_Coordinates));
                    activePointMask = abs(allPeakY - activeY) <= tolerance;
                    fprintf('DEBUG: Y-slice %d, activeY=%.3f, tolerance=%.3f, active points=%d\n', ...
                            sliceIndex, activeY, tolerance, sum(activePointMask));
                else
                    activePointMask = false(size(allPeakY));
                end
                
            case 'YtVsX'
                % Aligning X-slices: highlight points at specific X coordinate
                if sliceIndex <= length(X_Coordinates)
                    activeX = X_Coordinates(sliceIndex);
                    tolerance = (max(X_Coordinates) - min(X_Coordinates)) / (2 * length(X_Coordinates));
                    fprintf('DEBUG: X-slice %d, activeX=%.3f, tolerance=%.3f\n', ...
                            sliceIndex, activeX, tolerance);
                else
                    fprintf('DEBUG: X-slice index %d out of range\n', sliceIndex);
                    return;
                end
                
            otherwise
                fprintf('DEBUG: Unknown viewType: %s\n', viewType);
                return;
        end
        
        % CRITICAL FIX: Handle separate peak/valley scatter objects
        % The 3D plot creates separate scatter objects for peaks and valleys
        % We need to create separate masks for each scatter object

        fprintf('DEBUG: Creating separate masks for peak/valley scatter objects\n');

        % Get peak types to separate peaks from valleys
        if isfield(figData, 'allPeakTypes')
            allPeakTypes = figData.allPeakTypes;
            fprintf('DEBUG: Found peak types data with %d entries\n', length(allPeakTypes));
        else
            fprintf('DEBUG: No peak types found - using combined approach\n');
            allPeakTypes = [];
        end

        % Update scatter plot visualization
        for i = 1:length(scatterObjects)
            scatter_obj = scatterObjects(i);

            % Get current data from scatter object
            xData = get(scatter_obj, 'XData');
            yData = get(scatter_obj, 'YData');
            zData = get(scatter_obj, 'ZData');

            fprintf('DEBUG: Scatter object %d has %d points\n', i, length(xData));

            if isempty(xData)
                fprintf('DEBUG: Scatter object %d is empty, skipping\n', i);
                continue;
            end

            % Create mask for this specific scatter object's data
            % Match the scatter object's coordinates with our active slice coordinates
            scatterMask = false(size(xData));

            % Find which points in this scatter object belong to the active slice
            for j = 1:length(xData)
                switch viewType
                    case 'XtVsY'
                        % Check if this point's Y coordinate matches the active slice
                        if abs(yData(j) - activeY) <= tolerance
                            scatterMask(j) = true;
                        end
                    case 'YtVsX'
                        % Check if this point's X coordinate matches the active slice
                        if abs(xData(j) - activeX) <= tolerance
                            scatterMask(j) = true;
                        end
                end
            end

            activePointsInScatter = sum(scatterMask);
            fprintf('DEBUG: Scatter object %d: %d active points out of %d total\n', ...
                i, activePointsInScatter, length(xData));
            
            % EXTREME VISUAL CHANGES: Make the difference very obvious
            alphaValues = ones(size(xData)) * 0.01; % Nearly invisible (1% opacity)
            alphaValues(scatterMask) = 1.0;         % Full opacity for active slice

            % Dramatic marker size changes
            markerSizes = ones(size(xData)) * 2;    % Tiny markers for inactive
            markerSizes(scatterMask) = 50;          % Very large markers for active slice

            % Apply extreme visual changes
            set(scatter_obj, 'MarkerFaceAlpha', 'flat', 'AlphaData', alphaValues);
            set(scatter_obj, 'MarkerEdgeAlpha', 'flat', 'AlphaDataMode', 'manual');
            set(scatter_obj, 'SizeData', markerSizes);

            % ADDITIONAL: Change color for active slice
            if any(scatterMask)
                % Get current color data
                currentColors = get(scatter_obj, 'CData');
                if ~isempty(currentColors) && size(currentColors, 1) == length(xData)
                    % Make active points bright red for maximum visibility
                    newColors = currentColors;
                    if size(newColors, 2) == 3 % RGB colors
                        newColors(scatterMask, :) = repmat([1, 0, 0], sum(scatterMask), 1); % Bright red
                    else % Single color values
                        newColors(scatterMask) = max(newColors(:)); % Maximum color value
                    end
                    set(scatter_obj, 'CData', newColors);
                end
            end

            fprintf('DEBUG: Applied highlighting to scatter object %d (%d active points)\n', i, activePointsInScatter);
        end
        
        % Update status to show which slice is active
        statusText = findobj(mainFig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', sprintf('Status: %s - Processing slice %d/%d', ...
                alignmentPhase, sliceIndex, totalSlices));
        end

        % ADDITIONAL VISUAL INDICATOR: Update figure title
        titleStr = sprintf('3D Peak Visualization - ACTIVE: %s Slice %d/%d', ...
            alignmentPhase, sliceIndex, totalSlices);
        set(mainFig, 'Name', titleStr);

        % ADDITIONAL: Add text annotation on the plot
        try
            % Remove previous annotation if it exists
            prevAnnotation = findobj(mainAx, 'Tag', 'ActiveSliceAnnotation');
            if ~isempty(prevAnnotation)
                delete(prevAnnotation);
            end

            % Add new annotation
            annotationText = sprintf('ACTIVE: %s %d/%d', alignmentPhase, sliceIndex, totalSlices);
            text(mainAx, 0.02, 0.98, annotationText, ...
                'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', ...
                'Color', 'red', 'BackgroundColor', 'yellow', ...
                'Tag', 'ActiveSliceAnnotation', 'Interpreter', 'none');
        catch
            % Ignore annotation errors
        end
        
        % EFFICIENT VISUAL UPDATE: Single refresh without slowing computation
        drawnow limitrate;  % Rate-limited refresh - won't slow down computation

        fprintf('DEBUG: Highlighting complete for slice %d - VISUAL UPDATE FORCED\n', sliceIndex);
        
    catch ME
        ErrorHandler.handleError(ME, 'highlightActiveSlice', ErrorHandler.WARNING);
    end
end
