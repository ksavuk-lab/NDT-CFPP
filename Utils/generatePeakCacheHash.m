function hashStr = generatePeakCacheHash(peakCacheKey)
    % GENERATEPEAKCACHEHASH - Generate a unique hash string for peak cache identification
    %
    % This function creates a unique hash string based on all parameters that affect
    % peak extraction results, ensuring that cached data is only used when all
    % relevant parameters match exactly.
    %
    % Input:
    %   peakCacheKey - Structure containing all parameters affecting peak extraction
    %
    % Output:
    %   hashStr - Short hash string for filename use
    
    try
        % Convert the entire structure to a string representation
        keyStr = struct2str(peakCacheKey);
        
        % Create a simple hash using MATLAB's built-in functions
        % Use Java's hashCode for consistency across sessions
        javaStr = java.lang.String(keyStr);
        hashCode = javaStr.hashCode();
        
        % Convert to positive hex string (8 characters)
        hashStr = sprintf('%08X', typecast(int32(hashCode), 'uint32'));
        
    catch ME
        % Fallback: use current time if hashing fails
        fprintf('Warning: Hash generation failed (%s). Using timestamp fallback.\n', ME.message);
        hashStr = datestr(now, 'yyyymmdd_HHMMSS');
    end
end

function str = struct2str(s)
    % Convert structure to string representation for hashing
    % This function recursively converts all fields to strings
    
    if isstruct(s)
        fields = fieldnames(s);
        str = '';
        for i = 1:length(fields)
            fieldName = fields{i};
            fieldValue = s.(fieldName);
            str = [str, fieldName, ':', struct2str(fieldValue), ';'];
        end
    elseif isnumeric(s)
        if isscalar(s)
            str = sprintf('%.6g', s);
        else
            str = sprintf('%.6g,', s(:));
            str = str(1:end-1); % Remove trailing comma
        end
    elseif islogical(s)
        if s
            str = 'true';
        else
            str = 'false';
        end
    elseif ischar(s) || isstring(s)
        str = char(s);
    elseif iscell(s)
        str = '';
        for i = 1:length(s)
            str = [str, struct2str(s{i}), ','];
        end
        if ~isempty(str)
            str = str(1:end-1); % Remove trailing comma
        end
    else
        % For other types, convert to string
        str = char(string(s));
    end
end
