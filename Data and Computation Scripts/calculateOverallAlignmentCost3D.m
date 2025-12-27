function totalCost = calculateOverallAlignmentCost3D(inputData)
    % CALCULATEOVERALLALIGNMENTCOST3D - Calculate overall alignment cost for convergence checking
    % Imported from XtVsYPlot.m and adapted for 3D peak plot compatibility
    %
    % This function calculates the overall alignment cost by computing the mean squared error
    % between all pairs of time segments across all statistics. For sparse peak data, it only
    % considers non-zero elements to avoid bias from empty regions.
    %
    % Input:
    %   inputData - Cell array of statistical data structures, each containing:
    %               - maps: Cell array of 2D matrices (one per time segment)
    %
    % Output:
    %   totalCost - Scalar value representing overall alignment cost
    %               Lower values indicate better alignment
    %               Returns inf if calculation fails
    %
    % Algorithm:
    %   For each statistic:
    %     For each pair of time segments:
    %       Calculate MSE between segments (only non-zero elements for peak data)
    %       Add to total cost
    %
    % Optimization for Peak Data:
    %   - Only considers non-zero elements to avoid bias from sparse data
    %   - Handles empty segments gracefully
    %   - Provides robust error handling

    totalCost = 0;
    nStats = length(inputData);

    try
        for statIdx = 1:nStats
            numSegments = length(inputData{statIdx}.maps);
            [numY, numX] = size(inputData{statIdx}.maps{1});

            % Calculate cost for each segment pair
            for seg1 = 1:numSegments-1
                for seg2 = seg1+1:numSegments
                    map1 = inputData{statIdx}.maps{seg1};
                    map2 = inputData{statIdx}.maps{seg2};

                    % Calculate MSE between segments (measure of alignment quality)
                    % CRITICAL: Only consider non-zero elements - zeros are not real data points
                    % Zero values represent "no peak detected" and should not influence cost calculation

                    % Only compare positions where BOTH maps have non-zero values
                    % This ensures we're only measuring alignment quality of actual peak data
                    nonZeroMask = (map1 ~= 0) & (map2 ~= 0);

                    if any(nonZeroMask(:))
                        % Calculate MSE only for positions with actual peak data in both segments
                        segmentCost = mean((map1(nonZeroMask) - map2(nonZeroMask)).^2);
                    else
                        % No overlapping peak data between segments - no cost contribution
                        segmentCost = 0;
                    end
                    totalCost = totalCost + segmentCost;
                end
            end
        end

    catch ME
        fprintf('Error calculating alignment cost: %s\n', ME.message);
        totalCost = inf;
    end
end
