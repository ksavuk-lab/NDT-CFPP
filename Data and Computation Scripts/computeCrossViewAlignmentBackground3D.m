function computeCrossViewAlignmentBackground3D(fig)
    % COMPUTECROSSVIEWALIGNMENTBACKGROUND3D - Background computation for iterative cross-view alignment
    % Imported from XtVsYPlot.m and adapted for 3D peak plot compatibility
    %
    % This function performs the actual iterative alignment computation:
    % 1. Alternates between XtVsY and YtVsX alignment
    % 2. Calculates alignment cost for convergence checking
    % 3. Updates progress dialog and status
    % 4. Stores results back to figure userData
    
    try
        userData = get(fig, 'UserData');
        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');

        % Start timing
        crossViewStartTime = tic;

        % Get data dimensions and initialize
        currentData = userData.statDataArray;
        if isempty(currentData)
            fprintf('Error: Current data not available\n');
            return;
        end

        [numY, numX] = size(currentData{1}.maps{1});
        nStats = length(currentData);
        numSegments = length(currentData{1}.maps);

        % Cross-view alignment parameters
        maxCrossViewIterations = 10; % Maximum number of cross-view iterations
        convergenceThreshold = 0.001; % Convergence threshold for improvement

        % Initialize tracking with pre-allocated arrays (no dynamic growth)
        crossViewIteration = 0;
        previousCost = inf;
        improvementHistory = zeros(maxCrossViewIterations, 1); % Pre-allocate for performance

        % ENHANCED STATISTICS: Initialize comprehensive statistics reporter
        statsReporter = AlignmentStatisticsReporter();

        fprintf('Starting cross-view alignment: %d Y-slices x %d X-slices\n', numY, numX);

        % Initialize visualization highlighting
        if isfield(userData, 'mainFigure') && ishghandle(userData.mainFigure)
            % Start with first Y-slice highlighted
            highlightActiveSlice(userData.mainFigure, 'XtVsY', 1, numY, 'Y-slices (Starting)');
        end

        % Main cross-view iteration loop
        while crossViewIteration < maxCrossViewIterations
            crossViewIteration = crossViewIteration + 1;
            iterationStartTime = tic;

            fprintf('Cross-sectional alignment iteration %d/%d:\n', crossViewIteration, maxCrossViewIterations);

            % Update non-blocking progress display
            progressValue = (crossViewIteration - 1) / maxCrossViewIterations;
            if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog) && ...
               isfield(userData.progressDialog, 'Figure') && ishghandle(userData.progressDialog.Figure)

                % Update progress text
                progressMsg = sprintf(['Cross-sectional alignment iteration %d/%d:\n' ...
                    'Step 1: Aligning Y-slices (horizontal alignment)...\n' ...
                    'Processing %d horizontal slices across grid\n' ...
                    'Grid: %dx%d, Time segments: %d\n' ...
                    'Zero values excluded from alignment\n\n' ...
                    'Progress: %.1f%% complete'], ...
                    crossViewIteration, maxCrossViewIterations, numY, numY, numX, numSegments, ...
                    ((crossViewIteration-1)/maxCrossViewIterations)*100);

                set(userData.progressDialog.Text, 'String', progressMsg);

                % Update progress bar
                progressBarMsg = sprintf('Progress: %.1f%% - Iteration %d/%d', ...
                    progressValue*100, crossViewIteration, maxCrossViewIterations);
                set(userData.progressDialog.Bar, 'String', progressBarMsg, ...
                    'BackgroundColor', [0.2 + 0.6*progressValue, 0.8, 0.2]);

                drawnow;
            end

            % Step 1: Align all Y-slices (horizontal alignment)
            fprintf('  Step 1: Aligning Y-slices (horizontal alignment)...\n');
            if ~isempty(statusText)
                set(statusText, 'String', sprintf('Status: Cross-view %d/%d: Y-slices (%dx%d grid)', ...
                    crossViewIteration, maxCrossViewIterations, numY, numX));
            end

            ySliceStartTime = tic;
            if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog)
                currentData = alignAllSlicesInView3D(currentData, 'XtVsY', userData.progressDialog, userData.mainFigure);
            else
                currentData = alignAllSlicesInView3D(currentData, 'XtVsY');
            end
            ySliceTime = toc(ySliceStartTime);
            fprintf('    Y-slice alignment completed in %.1f seconds\n', ySliceTime);

            % Live update: refresh main figure with intermediate results
            if isfield(userData, 'mainFigure') && ishghandle(userData.mainFigure)
                updateLiveVisualization(userData.mainFigure, currentData, 'Y-slices aligned');
            end

            % Step 2: Align all X-slices (YtVsX view)
            fprintf('  Step 2: Aligning X-slices (YtVsX)...\n');
            if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog) && ...
               isfield(userData.progressDialog, 'Figure') && ishghandle(userData.progressDialog.Figure)

                progressMsg = sprintf(['Cross-view iteration %d/%d:\n' ...
                    'Step 2: Aligning X-slices (YtVsX)...\n' ...
                    'Processing %d X-slices in YtVsX view\n' ...
                    'Y-slice alignment completed in %.1f seconds\n' ...
                    'Grid: %dx%d, Time segments: %d\n' ...
                    'Zero values excluded from alignment\n\n' ...
                    'Progress: %.1f%% complete'], ...
                    crossViewIteration, maxCrossViewIterations, numX, ySliceTime, numY, numX, numSegments, ...
                    ((crossViewIteration-0.5)/maxCrossViewIterations)*100);

                set(userData.progressDialog.Text, 'String', progressMsg);
                drawnow;
            end
            if ~isempty(statusText)
                set(statusText, 'String', sprintf('Status: Cross-view %d/%d: X-slices (%dx%d grid)', ...
                    crossViewIteration, maxCrossViewIterations, numY, numX));
            end

            xSliceStartTime = tic;
            if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog)
                currentData = alignAllSlicesInView3D(currentData, 'YtVsX', userData.progressDialog, userData.mainFigure);
            else
                currentData = alignAllSlicesInView3D(currentData, 'YtVsX');
            end
            xSliceTime = toc(xSliceStartTime);
            fprintf('    X-slice alignment completed in %.1f seconds\n', xSliceTime);

            % Live update: refresh main figure with intermediate results
            if isfield(userData, 'mainFigure') && ishghandle(userData.mainFigure)
                updateLiveVisualization(userData.mainFigure, currentData, 'X-slices aligned');
            end

            % Calculate overall alignment cost to check for convergence
            currentCost = calculateOverallAlignmentCost3D(currentData);
            if crossViewIteration == 1
                improvement = 0; % No improvement to calculate for first iteration
                yCost = currentCost * 0.5; % Estimate Y-slice contribution
                xCost = currentCost * 0.5; % Estimate X-slice contribution
            else
                improvement = (previousCost - currentCost) / previousCost;
                % Calculate individual phase costs (approximation)
                yCost = currentCost * (ySliceTime / (ySliceTime + xSliceTime));
                xCost = currentCost * (xSliceTime / (ySliceTime + xSliceTime));
            end
            improvementHistory(crossViewIteration) = improvement; % Use pre-allocated array

            % ENHANCED STATISTICS: Record iteration data
            statsReporter.recordIteration(crossViewIteration, yCost, xCost, currentCost, improvement, ...
                ySliceTime, xSliceTime, numY, numX);

            % ENHANCED STATISTICS: Display detailed iteration summary
            statsReporter.displayIterationSummary(crossViewIteration, yCost, xCost, currentCost, improvement, ...
                ySliceTime, xSliceTime, numY, numX, convergenceThreshold);

            previousCost = currentCost;

            % Check for convergence
            if improvement < convergenceThreshold && crossViewIteration > 1
                % ENHANCED STATISTICS: Mark as converged and display final summary
                statsReporter.setConverged(improvement);
                statsReporter.displayProgressChart();
                statsReporter.displayFinalSummary(crossViewIteration, true, convergenceThreshold);

                fprintf('Cross-view alignment converged after %d iterations (improvement < %.3f%%)\n', ...
                    crossViewIteration, convergenceThreshold*100);

                % Update progress dialog with convergence message
                if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog) && ...
                   isfield(userData.progressDialog, 'Figure') && ishghandle(userData.progressDialog.Figure)

                    convergenceMsg = sprintf(['✅ Cross-View Alignment CONVERGED!\n\n' ...
                        'Results Summary:\n' ...
                        '• Completed %d/%d iterations\n' ...
                        '• Final improvement: %.3f%% (threshold: %.3f%%)\n' ...
                        '• Grid processed: %d Y-slices × %d X-slices\n' ...
                        '• Time segments: %d\n' ...
                        '• Total time: %.1f seconds\n' ...
                        '• Zero values excluded from processing\n\n' ...
                        'Alignment results ready for viewing...'], ...
                        crossViewIteration, maxCrossViewIterations, improvement*100, convergenceThreshold*100, ...
                        numY, numX, numSegments, toc(crossViewStartTime));

                    set(userData.progressDialog.Text, 'String', convergenceMsg);
                    set(userData.progressDialog.Bar, 'String', 'CONVERGED!', ...
                        'BackgroundColor', [0.2, 0.8, 0.2]);
                    drawnow;
                    pause(3); % Show convergence message for 3 seconds
                end

                % Restore normal visualization after convergence
                if isfield(userData, 'mainFigure') && ishghandle(userData.mainFigure)
                    restoreNormalVisualization(userData.mainFigure);
                end

                break;
            end

            % Update progress for next iteration - mirror command line output
            if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog) && ...
               isfield(userData.progressDialog, 'Figure') && ishghandle(userData.progressDialog.Figure)
                userData.progressDialog.Value = crossViewIteration / maxCrossViewIterations;
                userData.progressDialog.Message = sprintf(['✅ Iteration %d complete:\n' ...
                    '   • Alignment cost: %.6f\n' ...
                    '   • Improvement: %.4f%% (threshold: %.3f%%)\n' ...
                    '   • Total time: %.1f seconds (Y: %.1fs, X: %.1fs)\n' ...
                    '   • Progress: %.1f%% complete\n\n' ...
                    'Preparing next iteration...'], ...
                    crossViewIteration, currentCost, improvement*100, convergenceThreshold*100, ...
                    iterationTime, ySliceTime, xSliceTime, (crossViewIteration/maxCrossViewIterations)*100);
                drawnow;
            end
        end

        % ENHANCED STATISTICS: Display final summary if max iterations reached
        if crossViewIteration >= maxCrossViewIterations
            statsReporter.displayProgressChart();
            statsReporter.displayFinalSummary(crossViewIteration, false, convergenceThreshold);
        end

        % Store aligned results
        userData.alignedStatDataArray = currentData;
        userData.isComputingAlignment = false;

        % Restore original figData structure with alignment results
        if isfield(userData, 'originalFigData')
            figData = userData.originalFigData;

            % Store alignment results in the correct field names
            if ~isempty(userData.alignedStatDataArray)
                % Extract from cell array if needed
                if iscell(userData.alignedStatDataArray) && length(userData.alignedStatDataArray) == 1
                    figData.alignmentData.alignedStatMaps = userData.alignedStatDataArray{1};
                    fprintf('✅ Stored aligned results: extracted from cell array\n');
                else
                    figData.alignmentData.alignedStatMaps = userData.alignedStatDataArray;
                    fprintf('✅ Stored aligned results: used directly\n');
                end
                figData.alignmentData.isAligned = true;
                fprintf('✅ Cross-sectional alignment completed successfully\n');
            else
                figData.alignmentData.isAligned = false;
                fprintf('❌ No alignment results returned\n');
            end

            % Restore the original figData structure
            set(fig, 'UserData', figData);
        else
            % Fallback: just update userData
            set(fig, 'UserData', userData);
        end

        totalElapsedTime = toc(crossViewStartTime);
        avgImprovement = mean(improvementHistory(1:crossViewIteration)) * 100;

        % Update status
        if ~isempty(statusText)
            set(statusText, 'String', sprintf('Status: Cross-View Aligned (%d iter, %.1f%% avg improve, %.1fs)', ...
                crossViewIteration, avgImprovement, totalElapsedTime));
        end

        % Close progress dialog
        if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog) && ...
           isfield(userData.progressDialog, 'Figure') && ishghandle(userData.progressDialog.Figure)
            close(userData.progressDialog.Figure);
        end

        fprintf('Cross-view alignment complete: %d iterations, %.1fs total\n', crossViewIteration, totalElapsedTime);
        fprintf('Use "Aligned View" button to see results.\n');

        % Restore normal visualization - show all points with full opacity
        if isfield(userData, 'mainFigure') && ishghandle(userData.mainFigure)
            restoreNormalVisualization(userData.mainFigure);
        end

    catch ME
        fprintf('Error in cross-view alignment: %s\n', ME.message);

        userData = get(fig, 'UserData');
        userData.isComputingAlignment = false;

        % Restore original figData structure even on error
        if isfield(userData, 'originalFigData')
            figData = userData.originalFigData;
            figData.alignmentData.isAligned = false;
            set(fig, 'UserData', figData);
        else
            set(fig, 'UserData', userData);
        end

        if isfield(userData, 'progressDialog') && isstruct(userData.progressDialog) && ...
           isfield(userData.progressDialog, 'Figure') && ishghandle(userData.progressDialog.Figure)
            close(userData.progressDialog.Figure);
        end

        % Restore normal visualization even on error
        if isfield(userData, 'mainFigure') && ishghandle(userData.mainFigure)
            restoreNormalVisualization(userData.mainFigure);
        end

        statusText = findobj(fig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', 'Status: Cross-view Error');
        end
    end

    % Clean up timer
    try
        crossViewTimer = timerfind('Name', 'CrossViewTimer3D');
        if ~isempty(crossViewTimer)
            stop(crossViewTimer);
            delete(crossViewTimer);
        end
    catch
        % Ignore timer cleanup errors
    end
end

function updateLiveVisualization(mainFig, alignedData, statusMsg)
    % UPDATELIVEVUALIZATION - Update the main 3D figure with intermediate alignment results
    % This provides live feedback during the alignment process without blocking the UI

    try
        if ~ishghandle(mainFig)
            return;
        end

        % Get the current figure data
        figData = get(mainFig, 'UserData');
        if ~isfield(figData, 'alignmentData')
            return;
        end

        % Store the intermediate alignment results (but don't trigger plot updates)
        if iscell(alignedData) && length(alignedData) == 1
            figData.alignmentData.liveAlignedStatMaps = alignedData{1};
        else
            figData.alignmentData.liveAlignedStatMaps = alignedData;
        end

        % Update the figure data
        set(mainFig, 'UserData', figData);

        % ONLY update status text - don't trigger plot refreshes during alignment
        statusText = findobj(mainFig, 'Tag', 'AlignmentStatusText');
        if ~isempty(statusText)
            set(statusText, 'String', sprintf('Status: %s (Live Update)', statusMsg));
        end

        % Light refresh without triggering plot updates
        drawnow limitrate;

    catch ME
        % Silently handle errors to avoid disrupting alignment process
        % Don't print warnings during alignment to keep output clean
    end
end

% Note: highlightActiveSlice and restoreNormalVisualization functions
% are now in separate files for better accessibility across the codebase
