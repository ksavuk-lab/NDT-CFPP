function alignedData = alignAllSlicesInView3D(inputData, viewType, progressDialog, mainFigure)
    % ALIGNALLSLICESINVIEW3D - Align all slices in a specific view (XtVsY or YtVsX)
    % Imported from XtVsYPlot.m and adapted for 3D peak plot compatibility
    %
    % Inputs:
    %   inputData - Cell array of statistical data structures
    %   viewType - 'XtVsY' or 'YtVsX'
    %   progressDialog - Optional progress dialog handle for updates
    %   mainFigure - Optional main figure handle for live visualization updates
    %
    % This function aligns all slices in either XtVsY or YtVsX view:
    % - XtVsY: Aligns Y-slices (rows) across time segments
    % - YtVsX: Aligns X-slices (columns) across time segments
    %
    % Inputs:
    %   inputData - Cell array of statistical data structures
    %   viewType  - 'XtVsY' or 'YtVsX'
    %
    % Output:
    %   alignedData - Cell array with aligned statistical data
    %
    % The function uses alignColumnsImproved with parameters optimized for peak data:
    % - Reduced MaxShift (10 vs 15) for less noisy peak data
    % - Reduced LocalScope (3 vs 5) for faster processing
    % - Relaxed convergence (0.01 vs 0.001) appropriate for sparse data
    % - Reduced iterations (10 vs 20) for faster processing

    alignedData = inputData;
    [numY, numX] = size(inputData{1}.maps{1});
    numSegments = length(inputData{1}.maps);
    nStats = length(inputData);

    % Determine slice range based on view
    switch viewType
        case 'XtVsY'
            numSlices = numY;
            sliceType = 'Y';
        case 'YtVsX'
            numSlices = numX;
            sliceType = 'X';
        otherwise
            error('Unknown view type: %s', viewType);
    end

    fprintf('    Aligning %d %s-slices in %s view...\n', numSlices, sliceType, viewType);

    % Start timing for this alignment phase
    phaseStartTime = tic;

    % Alignment parameters (optimized for peak data)
    alignmentMethod = 'full'; % Use 'full' method for better alignment quality
    convergenceThreshold = 0.001;

    % Progress tracking
    processedSlices = 0;
    skippedSlices = 0;
    alignedSlices = 0;
    lastProgressTime = tic;
    progressInterval = 2; % Report every 2 slices for reasonable visual updates
    timeInterval = 1.0;   % Report every 1 second - enough to see progress without slowing computation

    % Process each slice
    for sliceIndex = 1:numSlices

        % Progress reporting every N slices or every N seconds
        if mod(sliceIndex, progressInterval) == 0 || toc(lastProgressTime) > timeInterval
            fprintf('      Processing %s slice %d/%d (%.1f%%) - Aligned: %d, Skipped: %d\n', ...
                    sliceType, sliceIndex, numSlices, (sliceIndex/numSlices)*100, alignedSlices, skippedSlices);

            % Update progress dialog with slice-level progress
            if nargin >= 3 && ~isempty(progressDialog) && isstruct(progressDialog) && ...
               isfield(progressDialog, 'Figure') && isvalid(progressDialog.Figure)

                progressMsg = sprintf(['    Aligning %d %s-slices in %s view...\n' ...
                    '      Processing %s slice %d/%d (%.1f%%)\n' ...
                    '      Aligned: %d, Skipped: %d\n' ...
                    '      Zero values excluded from alignment'], ...
                    numSlices, sliceType, viewType, sliceType, sliceIndex, numSlices, ...
                    (sliceIndex/numSlices)*100, alignedSlices, skippedSlices);

                set(progressDialog.Text, 'String', progressMsg);

                % Update progress bar for slice progress
                sliceProgress = sliceIndex / numSlices;
                progressBarMsg = sprintf('Processing %s slice %d/%d (%.1f%%)', ...
                    sliceType, sliceIndex, numSlices, sliceProgress*100);
                set(progressDialog.Bar, 'String', progressBarMsg, ...
                    'BackgroundColor', [0.2 + 0.6*sliceProgress, 0.8, 0.2]);
            end

            % HIGHLIGHT ACTIVE SLICE: Show user which slice is being processed
            if nargin >= 4 && ~isempty(mainFigure) && ishghandle(mainFigure)
                alignmentPhase = sprintf('%s-slices', sliceType);
                highlightActiveSlice(mainFigure, viewType, sliceIndex, numSlices, alignmentPhase);
            end

            lastProgressTime = tic;
        end
        % Process each statistic
        for statIdx = 1:nStats
            % Extract slice data across all time segments
            sliceData = zeros(numSegments, size(inputData{statIdx}.maps{1}, 2));
            
            switch viewType
                case 'XtVsY'
                    for seg = 1:numSegments
                        sliceData(seg, :) = inputData{statIdx}.maps{seg}(sliceIndex, :);
                    end
                case 'YtVsX'
                    for seg = 1:numSegments
                        sliceData(seg, :) = inputData{statIdx}.maps{seg}(:, sliceIndex)';
                    end
            end

            % CRITICAL: Skip if slice has no non-zero data or insufficient variation
            % Zero values represent "no peak detected" and should not participate in alignment
            nonZeroMask = sliceData ~= 0;
            nonZeroCount = sum(nonZeroMask(:));

            if nonZeroCount == 0
                % No actual peak data in this slice - skip alignment
                skippedSlices = skippedSlices + 1;
                continue;
            end

            if nonZeroCount < 3
                % Too few non-zero points for meaningful alignment - skip
                skippedSlices = skippedSlices + 1;
                continue;
            end

            % Check variation only among non-zero values
            nonZeroData = sliceData(nonZeroMask);
            if std(nonZeroData) < 1e-6
                % Insufficient variation in actual peak data - skip
                skippedSlices = skippedSlices + 1;
                continue;
            end

            try
                % Note: Removed real-time computation updates to avoid slowing down math

                % REAL-TIME CALLBACK: Create progress callback for continuous updates
                progressCallback = [];
                if nargin >= 4 && ~isempty(mainFigure) && ishghandle(mainFigure)
                    progressCallback = @(iteration, columnProgress, alignedData, col, numCols, overallProgress) ...
                        updateAlignmentProgress(mainFigure, alignmentPhase, sliceIndex, numSlices, ...
                                              iteration, columnProgress, overallProgress);
                end

                % Apply alignment using alignColumnsImproved with parameters optimized for peak data
                % IMPORTANT: Zero values will be preserved as zeros (no interpolation)
                [alignedSliceData, ~] = alignColumnsImproved(sliceData, ...
                    'MaxShift', 10, ...        % Reduced for peak data
                    'CostFunction', 'mse', ...
                    'AlignmentMethod', alignmentMethod, ...
                    'LocalScope', 3, ...       % Reduced for speed
                    'PadMethod', 'zeros', ...  % Maintains zeros as zeros
                    'Verbose', false, ...
                    'ConvergenceThreshold', 0.01, ... % Relaxed for peak data
                    'MaxIterations', 10, ...   % Reduced for speed
                    'WeightingFunction', 'exponential', ...
                    'WeightingScale', 2.0, ... % Reduced for peak data
                    'ProgressCallback', progressCallback); % REAL-TIME UPDATES!

                % Put aligned data back
                switch viewType
                    case 'XtVsY'
                        for seg = 1:numSegments
                            alignedData{statIdx}.maps{seg}(sliceIndex, :) = alignedSliceData(seg, :);
                        end
                    case 'YtVsX'
                        for seg = 1:numSegments
                            alignedData{statIdx}.maps{seg}(:, sliceIndex) = alignedSliceData(seg, :)';
                        end
                end

                alignedSlices = alignedSlices + 1;

            catch ME
                fprintf('Warning: Alignment failed for %s slice %d: %s\n', sliceType, sliceIndex, ME.message);
                % Continue with original data for this slice
                skippedSlices = skippedSlices + 1;
            end
        end

        processedSlices = processedSlices + 1;
    end

    % Calculate and display alignment improvement for this phase
    alignmentEfficiency = (alignedSlices / numSlices) * 100;
    avgTimePerSlice = (toc(phaseStartTime)) / numSlices;

    fprintf('    ✅ Completed %s alignment: %d total slices (%d aligned, %d skipped)\n', ...
            viewType, numSlices, alignedSlices, skippedSlices);
    fprintf('       • Alignment efficiency: %.1f%% (%d/%d slices successfully aligned)\n', ...
            alignmentEfficiency, alignedSlices, numSlices);
    fprintf('       • Average time per slice: %.3f seconds\n', avgTimePerSlice);

    if alignmentEfficiency >= 95
        fprintf('       • Quality assessment: 🟢 EXCELLENT (≥95%% success rate)\n');
    elseif alignmentEfficiency >= 85
        fprintf('       • Quality assessment: 🟡 GOOD (85-95%% success rate)\n');
    elseif alignmentEfficiency >= 70
        fprintf('       • Quality assessment: 🟠 MODERATE (70-85%% success rate)\n');
    else
        fprintf('       • Quality assessment: 🔴 POOR (<70%% success rate)\n');
    end

    % Show completion of this alignment phase
    if nargin >= 4 && ~isempty(mainFigure) && isvalid(mainFigure)
        statusText = findobj(mainFigure, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', sprintf('Status: %s alignment complete (%d/%d aligned)', ...
                viewType, alignedSlices, numSlices));
        end
        drawnow;
    end
end
