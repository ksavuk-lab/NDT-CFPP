function hashStr = generateAlignmentCacheHash(alignmentCacheKey)
    % GENERATEALIGNMENTCACHEHASH - Generate a unique hash for alignment cache
    %
    % Creates a unique hash based on all parameters affecting waveform alignment,
    % ensuring cached alignment data is only used when parameters match exactly.
    %
    % Input:
    %   alignmentCacheKey - Structure containing:
    %       .dataset        - Dataset number
    %       .timeRange      - [min max] time range for peak search
    %       .peakType       - 'positive', 'negative', or 'both'
    %       .dataSize       - [Y, X, T] size of waveform data
    %
    % Output:
    %   hashStr - 8-character hex hash string for filename use
    
    try
        % Convert the structure to a string representation
        keyStr = struct2str(alignmentCacheKey);
        
        % Create hash using Java's hashCode for consistency
        javaStr = java.lang.String(keyStr);
        hashCode = javaStr.hashCode();
        
        % Convert to positive hex string (8 characters)
        hashStr = sprintf('%08X', typecast(int32(hashCode), 'uint32'));
        
    catch ME
        % Fallback: use timestamp if hashing fails
        fprintf('Warning: Alignment hash generation failed (%s). Using timestamp.\n', ME.message);
        hashStr = datestr(now, 'yyyymmdd_HHMMSS');
    end
end

function str = struct2str(s)
    % Convert structure to string representation for hashing
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
        str = char(string(s));
    elseif ischar(s) || isstring(s)
        str = char(s);
    elseif iscell(s)
        str = '';
        for i = 1:length(s)
            str = [str, struct2str(s{i}), ','];
        end
        if ~isempty(str)
            str = str(1:end-1);
        end
    else
        str = char(string(s));
    end
end

