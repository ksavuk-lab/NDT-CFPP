classdef ErrorHandler < handle
    % ERRORHANDLER - Unified error handling system for 3D Peak Visualization
    % Provides consistent error reporting, logging, and recovery mechanisms
    
    properties (Constant)
        % Error severity levels
        SILENT = 0;     % No output
        WARNING = 1;    % Warning messages only
        ERROR = 2;      % Error messages with stack trace
        VERBOSE = 3;    % Detailed debugging information
    end
    
    properties
        verbosityLevel = 2;  % Default to ERROR level
        logToFile = false;
        logFileName = '';
    end
    
    methods (Static)
        function instance = getInstance()
            % Singleton pattern for global error handler
            persistent errorHandlerInstance;
            if isempty(errorHandlerInstance)
                errorHandlerInstance = ErrorHandler();
            end
            instance = errorHandlerInstance;
        end
        
        function handleError(ME, context, severity)
            % HANDLEERROR - Unified error handling with context
            % 
            % Inputs:
            %   ME - MException object
            %   context - String describing where error occurred
            %   severity - Error severity level (optional, default: ERROR)
            
            if nargin < 3
                severity = ErrorHandler.ERROR;
            end
            
            handler = ErrorHandler.getInstance();
            
            if severity <= handler.verbosityLevel
                switch severity
                    case ErrorHandler.SILENT
                        % No output
                        return;
                        
                    case ErrorHandler.WARNING
                        fprintf('Warning in %s: %s\n', context, ME.message);
                        
                    case ErrorHandler.ERROR
                        fprintf('Error in %s: %s\n', context, ME.message);
                        if handler.verbosityLevel >= ErrorHandler.VERBOSE
                            fprintf('Stack trace:\n');
                            for i = 1:length(ME.stack)
                                fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
                            end
                        end
                        
                    case ErrorHandler.VERBOSE
                        fprintf('DEBUG Error in %s: %s\n', context, ME.message);
                        fprintf('Full stack trace:\n');
                        for i = 1:length(ME.stack)
                            fprintf('  %s (line %d) in %s\n', ME.stack(i).name, ME.stack(i).line, ME.stack(i).file);
                        end
                end
            end
            
            % Log to file if enabled
            if handler.logToFile && ~isempty(handler.logFileName)
                handler.logErrorToFile(ME, context, severity);
            end
        end
        
        function success = validateFigureHandle(fig, context)
            % VALIDATEFIGUREHANDLE - Standardized figure handle validation
            success = false;
            
            try
                if isempty(fig) || ~ishghandle(fig)
                    ErrorHandler.handleError(MException('ErrorHandler:InvalidHandle', ...
                        'Invalid figure handle'), context, ErrorHandler.WARNING);
                    return;
                end
                success = true;
            catch ME
                ErrorHandler.handleError(ME, context, ErrorHandler.ERROR);
            end
        end
        
        function success = validateFigData(figData, requiredFields, context)
            % VALIDATEFIGDATA - Standardized figData validation
            success = false;
            
            try
                if ~isstruct(figData)
                    ErrorHandler.handleError(MException('ErrorHandler:InvalidData', ...
                        'figData is not a structure'), context, ErrorHandler.ERROR);
                    return;
                end
                
                for i = 1:length(requiredFields)
                    if ~isfield(figData, requiredFields{i})
                        ErrorHandler.handleError(MException('ErrorHandler:MissingField', ...
                            sprintf('Missing required field: %s', requiredFields{i})), ...
                            context, ErrorHandler.ERROR);
                        return;
                    end
                end
                
                success = true;
            catch ME
                ErrorHandler.handleError(ME, context, ErrorHandler.ERROR);
            end
        end
    end
    
    methods (Access = private)
        function logErrorToFile(obj, ME, context, severity)
            % Log error to file
            try
                fid = fopen(obj.logFileName, 'a');
                if fid ~= -1
                    fprintf(fid, '[%s] %s in %s: %s\n', ...
                        datestr(now), obj.getSeverityString(severity), context, ME.message);
                    fclose(fid);
                end
            catch
                % Ignore file logging errors
            end
        end
        
        function str = getSeverityString(~, severity)
            switch severity
                case ErrorHandler.SILENT
                    str = 'SILENT';
                case ErrorHandler.WARNING
                    str = 'WARNING';
                case ErrorHandler.ERROR
                    str = 'ERROR';
                case ErrorHandler.VERBOSE
                    str = 'VERBOSE';
                otherwise
                    str = 'UNKNOWN';
            end
        end
    end
end
