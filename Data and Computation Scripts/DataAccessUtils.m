classdef DataAccessUtils < handle
    % DATAACCESSUTILS - Unified data access utilities for 3D Peak Visualization
    % Provides consistent access to figData structures across different modes
    
    methods (Static)
        function data = getAlignmentData(figData)
            % GETALIGNMENTDATA - Get alignment data with consistent structure handling
            % Handles both cell array and direct structure access
            
            data = [];
            
            try
                if ~isstruct(figData)
                    ErrorHandler.handleError(MException('DataAccess:InvalidInput', ...
                        'figData must be a structure'), 'getAlignmentData', ErrorHandler.ERROR);
                    return;
                end
                
                % Check for alignment data
                if isfield(figData, 'alignmentData')
                    alignmentData = figData.alignmentData;
                    
                    % Handle different alignment data formats
                    if isfield(alignmentData, 'alignedStatMaps')
                        data = alignmentData.alignedStatMaps;
                    elseif isfield(alignmentData, 'originalStatMaps')
                        data = alignmentData.originalStatMaps;
                    end
                    
                    % Handle cell array wrapping
                    if iscell(data) && length(data) == 1
                        data = data{1};
                    end
                end
                
            catch ME
                ErrorHandler.handleError(ME, 'getAlignmentData', ErrorHandler.ERROR);
            end
        end
        
        function success = setAlignmentData(figData, data, isAligned)
            % SETALIGNMENTDATA - Set alignment data with consistent structure
            success = false;
            
            try
                if ~isstruct(figData)
                    ErrorHandler.handleError(MException('DataAccess:InvalidInput', ...
                        'figData must be a structure'), 'setAlignmentData', ErrorHandler.ERROR);
                    return;
                end
                
                % Ensure alignmentData structure exists
                if ~isfield(figData, 'alignmentData')
                    figData.alignmentData = struct();
                end
                
                % Store data consistently
                if nargin >= 3 && isAligned
                    figData.alignmentData.alignedStatMaps = data;
                    figData.alignmentData.isAligned = true;
                else
                    figData.alignmentData.originalStatMaps = data;
                    figData.alignmentData.isAligned = false;
                end
                
                success = true;
                
            catch ME
                ErrorHandler.handleError(ME, 'setAlignmentData', ErrorHandler.ERROR);
            end
        end
        
        function axes_handle = getMainAxes(figData)
            % GETMAINAXES - Get main axes handle with fallback logic
            axes_handle = [];
            
            try
                % Check for alignment mode first
                if isfield(figData, 'alignmentData') && isfield(figData, 'originalFigData')
                    % In alignment mode, use original figData
                    originalFigData = figData.originalFigData;
                    if isfield(originalFigData, 'mainAx') && ishghandle(originalFigData.mainAx)
                        axes_handle = originalFigData.mainAx;
                        return;
                    end
                end
                
                % Regular mode
                if isfield(figData, 'mainAx') && ishghandle(figData.mainAx)
                    axes_handle = figData.mainAx;
                    return;
                end
                
                % Fallback: look for any valid axes
                fields = fieldnames(figData);
                for i = 1:length(fields)
                    field = fields{i};
                    if contains(field, 'Ax') && isfield(figData, field)
                        candidate = figData.(field);
                        if ishghandle(candidate) && strcmp(get(candidate, 'Type'), 'axes')
                            axes_handle = candidate;
                            return;
                        end
                    end
                end
                
            catch ME
                ErrorHandler.handleError(ME, 'getMainAxes', ErrorHandler.WARNING);
            end
        end
        
        function coords = getPeakCoordinates(figData)
            % GETPEAKCOORDINATES - Get peak coordinate data with fallback
            coords = struct('X', [], 'Y', [], 'Z', [], 'Amps', []);

            try
                % CRITICAL FIX: The figData passed here is from the 3D plot figure itself
                % During alignment, this should be the original 3D plot data
                coords = DataAccessUtils.extractCoordinatesFromStruct(figData);

                if ~isempty(coords.X)
                    return;
                end

                % If that fails, try to find the data in alignment structure
                if isfield(figData, 'alignmentData') && isfield(figData, 'originalFigData')
                    originalFigData = figData.originalFigData;
                    coords = DataAccessUtils.extractCoordinatesFromStruct(originalFigData);
                end

            catch ME
                ErrorHandler.handleError(ME, 'getPeakCoordinates', ErrorHandler.WARNING);
            end
        end
        
        function coords = extractCoordinatesFromStruct(figData)
            % EXTRACTCOORDINATESFROMSTRUCT - Extract coordinates from structure
            coords = struct('X', [], 'Y', [], 'Z', [], 'Amps', []);

            % Try different field name patterns
            coordFields = {
                {'allPeakX', 'allPeakY', 'allPeakZ', 'allPeakAmps'};
                {'peakX', 'peakY', 'peakZ', 'peakAmps'};
                {'X', 'Y', 'Z', 'Amps'};
                {'x', 'y', 'z', 'amps'}
            };

            fprintf('DEBUG: Trying to extract coordinates from structure\n');

            for i = 1:length(coordFields)
                fields = coordFields{i};
                hasAllFields = all(isfield(figData, fields));
                fprintf('DEBUG: Pattern %d (%s): %d\n', i, strjoin(fields, ', '), hasAllFields);

                if hasAllFields
                    coords.X = figData.(fields{1});
                    coords.Y = figData.(fields{2});
                    coords.Z = figData.(fields{3});
                    coords.Amps = figData.(fields{4});

                    fprintf('DEBUG: Successfully extracted coordinates - X:%d, Y:%d, Z:%d, Amps:%d points\n', ...
                        length(coords.X), length(coords.Y), length(coords.Z), length(coords.Amps));
                    break;
                end
            end

            if isempty(coords.X)
                fprintf('DEBUG: No coordinate patterns matched!\n');
            end
        end
        
        function gridCoords = getGridCoordinates(figData)
            % GETGRIDCOORDINATES - Get spatial grid coordinates
            gridCoords = struct('X_Coordinates', [], 'Y_Coordinates', []);
            
            try
                % Check alignment mode first
                if isfield(figData, 'alignmentData') && isfield(figData, 'originalFigData')
                    originalFigData = figData.originalFigData;
                    gridCoords = DataAccessUtils.extractGridFromStruct(originalFigData);
                    if ~isempty(gridCoords.X_Coordinates)
                        return;
                    end
                end
                
                % Regular mode
                gridCoords = DataAccessUtils.extractGridFromStruct(figData);
                
            catch ME
                ErrorHandler.handleError(ME, 'getGridCoordinates', ErrorHandler.WARNING);
            end
        end
        
        function gridCoords = extractGridFromStruct(figData)
            % EXTRACTGRIDFROMSTRUCT - Extract grid coordinates from structure
            gridCoords = struct('X_Coordinates', [], 'Y_Coordinates', []);
            
            % Try different field name patterns
            if isfield(figData, 'X_Coordinates') && isfield(figData, 'Y_Coordinates')
                gridCoords.X_Coordinates = figData.X_Coordinates;
                gridCoords.Y_Coordinates = figData.Y_Coordinates;
            elseif isfield(figData, 'xCoords') && isfield(figData, 'yCoords')
                gridCoords.X_Coordinates = figData.xCoords;
                gridCoords.Y_Coordinates = figData.yCoords;
            end
        end
        
        function success = validateDataStructure(figData, context)
            % VALIDATEDATASTRUCTURE - Comprehensive data structure validation
            success = false;
            
            try
                % Basic structure validation
                if ~ErrorHandler.validateFigData(figData, {}, context)
                    return;
                end
                
                % Check for required visualization data
                coords = DataAccessUtils.getPeakCoordinates(figData);
                if isempty(coords.X) || isempty(coords.Y) || isempty(coords.Z)
                    ErrorHandler.handleError(MException('DataAccess:MissingCoords', ...
                        'Missing peak coordinate data'), context, ErrorHandler.ERROR);
                    return;
                end
                
                % Check for main axes
                mainAx = DataAccessUtils.getMainAxes(figData);
                if isempty(mainAx)
                    ErrorHandler.handleError(MException('DataAccess:MissingAxes', ...
                        'No valid axes found'), context, ErrorHandler.WARNING);
                    % Don't return false - axes might be created later
                end
                
                success = true;
                
            catch ME
                ErrorHandler.handleError(ME, context, ErrorHandler.ERROR);
            end
        end
    end
end
