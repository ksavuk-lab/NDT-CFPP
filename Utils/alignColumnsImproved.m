function [alignedData, convergenceInfo] = alignColumnsImproved(inputData, varargin)
    % ALIGNCOLUMNSIMPROVED - Enhanced column alignment with distance-weighted row alignment
    %
    % This function aligns columns in a 2D matrix using a distance-weighted approach
    % where each column is aligned to the entire row, with closer columns having
    % more influence. Designed for carbon fiber plate scanning with natural warping.
    %
    % Note: This utility has been moved to the Utils directory for better organization
    %
    % Syntax:
    %   [alignedData, convergenceInfo] = alignColumnsImproved(inputData)
    %   [alignedData, convergenceInfo] = alignColumnsImproved(inputData, 'Parameter', Value, ...)
    %
    % Inputs:
    %   inputData - 2D matrix (M x N) where each column will be aligned
    %
    % Optional Parameters:
    %   'MaxShift'         - Maximum shift range in pixels (default: 10)
    %   'CostFunction'     - Cost function: 'mse', 'correlation', 'ncc' (default: 'mse')
    %   'AlignmentMethod'  - Alignment method: 'full', 'local', 'average', 'optimizer' (default: 'optimizer')
    %   'LocalScope'       - Number of columns to compare on each side for 'local' method (default: 5)
    %   'OptimizerType'    - Optimizer: 'fminbnd', 'fminsearch', 'patternsearch' (default: 'fminbnd')
    %   'PadMethod'        - Padding method: 'circular', 'zeros', 'replicate' (default: 'circular')
    %   'Verbose'          - Display progress information (default: false)
    %   'ConvergenceThreshold' - Relative change threshold for convergence (default: 1e-4)
    %   'MaxIterations'    - Maximum iterations for convergence (default: 50)
    %   'WeightingFunction' - Distance weighting: 'exponential', 'gaussian', 'linear' (default: 'exponential')
    %   'WeightingScale'   - Scale parameter for weighting function (default: 3.0)
    %   'ProgressCallback' - Function handle for progress updates (default: [])
    %
    % Outputs:
    %   alignedData     - 2D matrix with aligned columns
    %   convergenceInfo - Structure with convergence information
    %
    % Example:
    %   [aligned, info] = alignColumnsImproved(data, 'MaxShift', 15, 'WeightingScale', 2.5);

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'inputData', @(x) isnumeric(x) && ismatrix(x));
    addParameter(p, 'MaxShift', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'CostFunction', 'mse', @(x) ischar(x) && ismember(lower(x), {'mse', 'correlation', 'ncc'}));
    addParameter(p, 'AlignmentMethod', 'average', @(x) ischar(x) && ismember(lower(x), {'full', 'average'}));
    addParameter(p, 'LocalScope', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'OptimizerType', 'fminbnd', @(x) ischar(x) && ismember(lower(x), {'fminbnd', 'fminsearch', 'patternsearch'}));
    addParameter(p, 'PadMethod', 'circular', @(x) ischar(x) && ismember(lower(x), {'circular', 'zeros', 'replicate'}));
    addParameter(p, 'Verbose', false, @(x) islogical(x) || (isnumeric(x) && ismember(x, [0, 1])));
    addParameter(p, 'ConvergenceThreshold', 0.01, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MaxIterations', 50, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'WeightingFunction', 'exponential', @(x) ischar(x) && ismember(lower(x), {'exponential', 'gaussian', 'linear'}));
    addParameter(p, 'WeightingScale', 3.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'ProgressCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
    
    parse(p, inputData, varargin{:});
    
    % Extract parameters
    maxShift = p.Results.MaxShift;
    costFunction = lower(p.Results.CostFunction);
    alignmentMethod = lower(p.Results.AlignmentMethod);
    localScope = p.Results.LocalScope;
    optimizerType = lower(p.Results.OptimizerType);
    padMethod = lower(p.Results.PadMethod);
    verbose = logical(p.Results.Verbose);
    convergenceThreshold = p.Results.ConvergenceThreshold;
    maxIterations = p.Results.MaxIterations;
    weightingFunction = lower(p.Results.WeightingFunction);
    weightingScale = p.Results.WeightingScale;
    progressCallback = p.Results.ProgressCallback;
    
    % Get matrix dimensions
    [numRows, numCols] = size(inputData);
    
    if verbose
        fprintf('Enhanced column alignment: %d columns, max shift ±%d pixels\n', numCols, maxShift);
        fprintf('Method: %s, convergence threshold: %.2e\n', alignmentMethod, convergenceThreshold);
        if strcmp(alignmentMethod, 'full')
            fprintf('Weighting: %s (scale=%.2f)\n', weightingFunction, weightingScale);
        end
    end
    
    % Initialize output matrix and tracking
    alignedData = inputData;
    previousShifts = zeros(1, numCols);
    convergenceHistory = zeros(1, maxIterations); % Pre-allocate for performance
    convergenceHistoryIndex = 0; % Track actual length

    % Calculate total columns to process for progress tracking
    totalColumnsToProcess = numCols; % Process all columns including edges


    % Initialize cost cache for performance (only for methods that can benefit)
    costCache = containers.Map('KeyType', 'char', 'ValueType', 'double');

    % Pre-compute row average for 'average' method
    if strcmp(alignmentMethod, 'average')
        rowAverage = mean(alignedData, 2); % Average across all columns for each row
        initialRowAverage = rowAverage; % Store initial average for convergence stability
    else
        rowAverage = [];
        initialRowAverage = [];
    end

    % Convergence loop
    for iteration = 1:maxIterations
        iterationStartTime = tic; % Time this iteration

        if verbose
            fprintf('Iteration %d/%d: ', iteration, maxIterations);
        end

        % Store shifts from this iteration
        currentShifts = zeros(1, numCols);
        totalCostImprovement = 0;

        % Process each column (including first and last columns)
        % Process middle columns first, then edge columns last for better convergence
        columnOrder = [2:(numCols-1), 1, numCols]; % Middle columns first, then edges

        for colIdx = 1:length(columnOrder)
            col = columnOrder(colIdx);
            columnStartTime = tic; % Time this column

            % Calculate column progress within iteration
            columnProgress = colIdx / totalColumnsToProcess;
            overallProgress = (iteration - 1 + columnProgress) / maxIterations;

            % Column-level progress callback
            if ~isempty(progressCallback)
                progressCallback(iteration, columnProgress, alignedData, col, numCols, overallProgress);
            end

            % Get current column
            currentCol = alignedData(:, col);

            % Find optimal shift using selected alignment method
            % For edge columns, use neighbor-only comparison for better local alignment
            if col == 1 || col == numCols
                % Edge columns: use single neighbor method for local alignment
                [bestShift, costImprovement] = findOptimalShiftEdge(currentCol, alignedData, col, ...
                    maxShift, costFunction, padMethod, costCache);
            else
                % Middle columns: use selected method
                switch alignmentMethod
                    case 'full'
                        % Original method: compare against all columns with weighting
                        [bestShift, costImprovement] = findOptimalShiftWeighted(currentCol, alignedData, col, ...
                            maxShift, costFunction, padMethod, weightingFunction, weightingScale, costCache);
                    case 'average'
                        % Row average method: compare against row average (fastest)
                        [bestShift, costImprovement] = findOptimalShiftAverage(currentCol, rowAverage, ...
                            maxShift, costFunction, padMethod, costCache);
                    otherwise
                        error('Unknown alignment method: %s', alignmentMethod);
                end
            end

            % Apply the shift (only integer shifts now - no subpixel)
            if abs(bestShift) > 0.5  % Use 0.5 threshold for integer shifts
                alignedData(:, col) = applyIntegerShift(currentCol, round(bestShift), padMethod);
                currentShifts(col) = round(bestShift);
                totalCostImprovement = totalCostImprovement + costImprovement;

                % Update row average if using average method (but less frequently for stability)
                if strcmp(alignmentMethod, 'average')
                    % Only update row average every few columns to improve convergence stability
                    if mod(col, 5) == 0 || col == numCols-1
                        rowAverage = mean(alignedData, 2); % Recalculate average after several shifts
                    end
                end
            end

            if verbose && mod(col, max(1, round(totalColumnsToProcess/10))) == 0
                columnTime = toc(columnStartTime);
                fprintf('  Column %d/%d (%.1f%%) - %.3fs\n', col, numCols, columnProgress*100, columnTime);
            end
        end

        % Calculate convergence metrics
        shiftChange = norm(currentShifts - previousShifts) / max(norm(currentShifts), 1e-10);
        convergenceHistoryIndex = convergenceHistoryIndex + 1;
        convergenceHistory(convergenceHistoryIndex) = shiftChange; % Use indexing instead of end+1

        iterationTime = toc(iterationStartTime);

        if verbose
            fprintf('shift change: %.6f, cost improvement: %.6f, time: %.3fs\n', shiftChange, totalCostImprovement, iterationTime);
        end

        % Check for convergence
        if shiftChange < convergenceThreshold
            if verbose
                fprintf('Converged after %d iterations (threshold: %.2e)\n', iteration, convergenceThreshold);
            end
            break;
        end

        % Update previous shifts
        previousShifts = currentShifts;

        % Progress callback for end of iteration
        if ~isempty(progressCallback)
            progressCallback(iteration, 1.0, alignedData, numCols+1, numCols, iteration/maxIterations);
        end
    end
    
    % Prepare convergence information
    convergenceInfo = struct();
    convergenceInfo.iterations = convergenceHistoryIndex; % Use actual length
    convergenceInfo.converged = (convergenceHistory(convergenceHistoryIndex) < convergenceThreshold);
    convergenceInfo.finalChange = convergenceHistory(convergenceHistoryIndex);
    convergenceInfo.convergenceHistory = convergenceHistory(1:convergenceHistoryIndex); % Trim to actual size
    convergenceInfo.threshold = convergenceThreshold;
    
    if verbose
        if convergenceInfo.converged
            fprintf('Alignment converged successfully\n');
        else
            fprintf('Alignment reached maximum iterations without full convergence\n');
        end
    end
end

function [bestShift, costImprovement] = findOptimalShiftWeighted(currentCol, fullData, colIndex, ...
    maxShift, costFunction, padMethod, weightingFunction, weightingScale, costCache)
    % Find optimal shift using distance-weighted comparison to entire row (original method)

    [numRows, numCols] = size(fullData);

    % Only integer shifts now (removed subpixel for performance)
    shiftRange = -maxShift:maxShift;

    bestCost = inf;
    bestShift = 0;
    originalCost = calculateWeightedRowCost(currentCol, fullData, colIndex, costFunction, weightingFunction, weightingScale, costCache);

    % Test each possible shift
    for shift = shiftRange
        % Apply shift to current column (integer only)
        shiftedCol = applyIntegerShift(currentCol, shift, padMethod);

        % Calculate weighted cost against entire row
        cost = calculateWeightedRowCost(shiftedCol, fullData, colIndex, costFunction, weightingFunction, weightingScale, costCache);

        % Update best shift if this is better
        if cost < bestCost
            bestCost = cost;
            bestShift = shift;
        end
    end

    costImprovement = originalCost - bestCost;
end

function [bestShift, costImprovement] = findOptimalShiftEdge(currentCol, fullData, colIndex, ...
    maxShift, costFunction, padMethod, costCache)
    % Find optimal shift for edge columns using only immediate neighbor

    [numRows, numCols] = size(fullData);

    % Only integer shifts (removed subpixel for performance)
    shiftRange = -maxShift:maxShift;

    bestCost = inf;
    bestShift = 0;

    % Determine neighbor column
    if colIndex == 1
        % First column: compare only with column 2
        neighborCol = fullData(:, 2);
    else
        % Last column: compare only with second-to-last column
        neighborCol = fullData(:, numCols - 1);
    end

    % Calculate original cost
    originalCost = calculateCost(currentCol, neighborCol, costFunction);

    % Test each possible shift
    for shift = shiftRange
        % Apply shift to current column (integer only)
        shiftedCol = applyIntegerShift(currentCol, shift, padMethod);

        % Calculate cost against single neighbor
        cost = calculateCost(shiftedCol, neighborCol, costFunction);

        % Update best shift if this is better
        if cost < bestCost
            bestCost = cost;
            bestShift = shift;
        end
    end

    costImprovement = originalCost - bestCost;
end

function [bestShift, costImprovement] = findOptimalShiftAverage(currentCol, rowAverage, ...
    maxShift, costFunction, padMethod, costCache)
    % Find optimal shift using comparison to row average - FASTEST METHOD

    shiftRange = -maxShift:maxShift;

    bestCost = inf;
    bestShift = 0;
    originalCost = calculateCost(currentCol, rowAverage, costFunction);

    % Test each possible shift
    for shift = shiftRange
        % Apply shift to current column (integer only)
        shiftedCol = applyIntegerShift(currentCol, shift, padMethod);

        % Calculate cost against row average
        cost = calculateCost(shiftedCol, rowAverage, costFunction);

        % Update best shift if this is better
        if cost < bestCost
            bestCost = cost;
            bestShift = shift;
        end
    end

    costImprovement = originalCost - bestCost;
end

function cost = calculateWeightedRowCost(column, fullData, colIndex, costFunction, weightingFunction, weightingScale, costCache)
    % Calculate distance-weighted cost against entire row (with caching)

    [numRows, numCols] = size(fullData);
    totalCost = 0;
    totalWeight = 0;

    % Compare against all other columns with distance weighting
    for otherCol = 1:numCols
        if otherCol == colIndex
            continue; % Skip self-comparison
        end

        % Calculate distance weight
        distance = abs(otherCol - colIndex);
        weight = calculateDistanceWeight(distance, weightingFunction, weightingScale);

        % Calculate cost between columns (with caching)
        columnCost = calculateCostCached(column, fullData(:, otherCol), costFunction, costCache, colIndex, otherCol);

        % Accumulate weighted cost
        totalCost = totalCost + weight * columnCost;
        totalWeight = totalWeight + weight;
    end

    % Return normalized weighted cost
    if totalWeight > 0
        cost = totalCost / totalWeight;
    else
        cost = inf;
    end
end

function weight = calculateDistanceWeight(distance, weightingFunction, scale)
    % Calculate distance-based weight (higher weight for closer columns)

    switch weightingFunction
        case 'exponential'
            weight = exp(-distance / scale);
        case 'gaussian'
            weight = exp(-(distance^2) / (2 * scale^2));
        case 'linear'
            weight = max(0, 1 - distance / scale);
        otherwise
            weight = 1; % Uniform weighting as fallback
    end
end

function cost = calculateCostCached(col1, col2, costFunction, costCache, colIndex1, colIndex2)
    % Calculate cost between two columns with caching for performance

    % Create cache key (order-independent)
    if colIndex1 < colIndex2
        cacheKey = sprintf('%d_%d_%s', colIndex1, colIndex2, costFunction);
    else
        cacheKey = sprintf('%d_%d_%s', colIndex2, colIndex1, costFunction);
    end

    % Check cache first
    if isKey(costCache, cacheKey)
        cost = costCache(cacheKey);
        return;
    end

    % Calculate cost and cache it
    cost = calculateCost(col1, col2, costFunction);
    costCache(cacheKey) = cost;
end

% Optimized helper functions (removed subpixel interpolation for performance)

function shiftedCol = applyIntegerShift(column, shift, padMethod)
    % Apply integer pixel shift with zero padding only (no wrapping)

    n = length(column);
    shiftedCol = zeros(size(column));

    if shift == 0
        shiftedCol = column;
        return;
    end

    % Only use zero padding - no circular wrapping allowed
    if shift > 0  % Shift down
        shiftedCol(shift+1:end) = column(1:end-shift);
        % shiftedCol(1:shift) remains zero (zero padding at top)
    else  % Shift up
        shift = abs(shift);
        shiftedCol(1:end-shift) = column(shift+1:end);
        % shiftedCol(end-shift+1:end) remains zero (zero padding at bottom)
    end
end

% Subpixel interpolation removed for performance optimization

function cost = calculateCost(col1, col2, costFunction)
    % Calculate cost between two columns, excluding edge rows to avoid padding effects

    % Remove NaN values for cost calculation
    validIdx = ~isnan(col1) & ~isnan(col2);
    if sum(validIdx) < 2
        cost = inf;  % Not enough valid data points
        return;
    end

    % Exclude edge rows (5 from top and bottom) to avoid padding effects
    edgeRows = 5;
    dataLength = length(col1);
    if dataLength > 2 * edgeRows
        % Create mask to exclude edge rows
        excludeIdx = false(size(col1));
        excludeIdx(1:edgeRows) = true;  % Top edge
        excludeIdx(end-edgeRows+1:end) = true;  % Bottom edge

        % Combine with NaN exclusion
        validIdx = validIdx & ~excludeIdx;
    end

    if sum(validIdx) < 2
        cost = inf;  % Not enough valid data points after excluding edges
        return;
    end

    col1_valid = col1(validIdx);
    col2_valid = col2(validIdx);

    switch costFunction
        case 'mse'
            % Mean Squared Error (lower is better)
            cost = mean((col1_valid - col2_valid).^2);

        case 'correlation'
            % Negative correlation (so lower is better for minimization)
            corrCoeff = corrcoef(col1_valid, col2_valid);
            if size(corrCoeff, 1) > 1
                cost = -corrCoeff(1, 2);  % Negative correlation
            else
                cost = -1;  % Perfect correlation if only one unique value
            end

        case 'ncc'
            % Normalized Cross-Correlation (negative for minimization)
            % Normalize the vectors
            col1_norm = (col1_valid - mean(col1_valid)) / (std(col1_valid) + eps);
            col2_norm = (col2_valid - mean(col2_valid)) / (std(col2_valid) + eps);

            % Calculate NCC
            ncc = mean(col1_norm .* col2_norm);
            cost = -ncc;  % Negative for minimization

        otherwise
            error('Unknown cost function: %s', costFunction);
    end
end
