function create3DPeakPlots(~, ~, peakData, ~, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, peakDetectionType, ~, showPeaksAndValleys, firstPeakAmplitudeThreshold)
    % CREATE3DPEAKPLOTS - Create 3D scatter plot of all peaks/valleys with alignment
    %
    % This function creates a single 3D visualization showing all detected peaks and valleys.
    % All waveforms are aligned by their first peak/valley (based on detection type), then
    % all discovered peaks/valleys are plotted with:
    % - X, Y: Spatial coordinates
    % - Z: Time values (depth)
    % - Color: Amplitude values
    % Grouping structure is ignored for this visualization.
    %
    % CRITICAL ZERO VALUE POLICY:
    % - Zero amplitude values represent "no peak detected" and are excluded from:
    %   * Cross-sectional alignment algorithms
    %   * Plate generation and grouping processes
    %   * Neighbor finding and binding operations
    %   * Cost calculations and averaging operations
    % - Zero values are preserved in visualizations but do not participate in processing
    %
    % Inputs:
    %   ~                   - Unused: groupedPeakData (ignored)
    %   ~                   - Unused: groupNames (ignored)
    %   peakData            - Cell array of peak data for each waveform
    %   ~                   - Unused: waveformArray (ignored)
    %   t                   - Time vector
    %   X_Coordinates       - X spatial coordinates
    %   Y_Coordinates       - Y spatial coordinates
    %   numY_sub            - Number of Y points
    %   numX_sub            - Number of X points
    %   peakDetectionType   - Type of peaks detected ('peaks', 'valleys', 'both')
    %   ~                   - Unused: alignWithTime (ignored - always uses alignment)
    %   showPeaksAndValleys - Boolean, show both peaks and valleys (true) or separate (false)

    fprintf('Generating 3D peak/valley visualization (alignment-based)...\n');

    % Check if we have valid data
    if isempty(peakData)
        fprintf('No peak data available for 3D plotting.\n');
        return;
    end

    % Validate peak data structure
    fprintf('Validating peak data structure...\n');
    fprintf('Peak data contains %d waveforms\n', length(peakData));

    % Check a few sample waveforms to ensure data structure is correct
    validCount = 0;
    requiredColumns = {'TransitionTime', 'TransitionAmplitude', 'TransitionType'};

    for i = 1:min(5, length(peakData))
        if ~isempty(peakData{i})
            if istable(peakData{i})
                % Check for required columns
                if all(ismember(requiredColumns, peakData{i}.Properties.VariableNames))
                    % Additional validation: check data types and ranges
                    try
                        times = peakData{i}.TransitionTime;
                        amps = peakData{i}.TransitionAmplitude;
                        types = peakData{i}.TransitionType;

                        % Validate data types
                        if isnumeric(times) && isnumeric(amps) && isnumeric(types)
                            % Validate ranges
                            if all(isfinite(times)) && all(isfinite(amps)) && all(ismember(types, [-1, 1]))
                                validCount = validCount + 1;
                            else
                                fprintf('Warning: Waveform %d has invalid data values (non-finite or invalid types).\n', i);
                            end
                        else
                            fprintf('Warning: Waveform %d has non-numeric data columns.\n', i);
                        end
                    catch ME
                        fprintf('Warning: Error validating waveform %d data: %s\n', i, ME.message);
                    end
                else
                    fprintf('Warning: Waveform %d missing required columns. Has: %s\n', i, strjoin(peakData{i}.Properties.VariableNames, ', '));
                end
            else
                fprintf('Warning: Waveform %d is not a table (type: %s).\n', i, class(peakData{i}));
            end
        end
    end

    if validCount == 0
        fprintf('Error: Peak data structure is invalid. Expected table with columns: %s\n', strjoin(requiredColumns, ', '));
        fprintf('TransitionType should contain values: 1 (peaks) or -1 (valleys)\n');
        fprintf('TransitionTime and TransitionAmplitude should contain finite numeric values\n');
        return;
    end

    fprintf('Peak data structure validation passed.\n');

    % Set default amplitude threshold if not provided
    if nargin < 13 || isempty(firstPeakAmplitudeThreshold)
        firstPeakAmplitudeThreshold = 0.1; % Default 10% of max amplitude
    end

    if firstPeakAmplitudeThreshold > 0
        fprintf('Using first peak amplitude threshold: %.3f (only peaks/valleys with |amplitude| >= %.3f will be considered for alignment)\n', ...
                firstPeakAmplitudeThreshold, firstPeakAmplitudeThreshold);
    else
        fprintf('No amplitude threshold applied (all peaks/valleys considered for alignment)\n');
    end

    % Create aligned 3D plot of all peaks/valleys
    createAligned3DPlot(peakData, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, peakDetectionType, showPeaksAndValleys, firstPeakAmplitudeThreshold);

    fprintf('3D peak/valley visualization complete.\n');
end

function createAligned3DPlot(peakData, t, X_Coordinates, Y_Coordinates, numY_sub, numX_sub, peakDetectionType, showPeaksAndValleys, firstPeakAmplitudeThreshold)
    % Create aligned 3D plot of all peaks/valleys
    % t parameter: time vector for cross-sectional alignment

    fprintf('Aligning waveforms and collecting peak/valley data...\n');

    numWaveforms = length(peakData);

    % Step 1: ROBUST PEAK ALIGNMENT ANALYSIS (with caching)
    fprintf('=== ROBUST PEAK ALIGNMENT ANALYSIS ===\n');

    try
        % Check for cached alignment analysis results
        alignCacheFolder = fullfile(pwd, 'Peak Cache');
        if ~exist(alignCacheFolder, 'dir'), mkdir(alignCacheFolder); end

        % Create cache key based on peak data characteristics
        alignCacheKey = struct('numWaveforms', numWaveforms, 'peakType', peakDetectionType, ...
            'gridSize', [numY_sub, numX_sub], 'timeRange', [t(1), t(end)]);
        alignCacheHash = generatePeakCacheHash(alignCacheKey);
        alignCacheFile = fullfile(alignCacheFolder, sprintf('RobustAlignCache_%s.mat', alignCacheHash));

        if exist(alignCacheFile, 'file')
            % Load cached alignment analysis
            fprintf('Loading cached robust alignment analysis: %s\n', alignCacheFile);
            cachedAlign = load(alignCacheFile);
            alignmentReference = cachedAlign.alignmentReference;
            analysisResults = cachedAlign.analysisResults;
        else
            % DEBUG: Check peak data structure before analysis
            fprintf('DEBUG: Checking peak data structure...\n');
            validPeakCount = 0;
            for i = 1:min(10, length(peakData))
                if ~isempty(peakData{i}) && isstruct(peakData{i})
                    transitions = peakData{i};
                    if isfield(transitions, 'TransitionType')
                        peakIndices = find(transitions.TransitionType == 1);
                        fprintf('  Waveform %d: %d peaks detected\n', i, length(peakIndices));
                        if length(peakIndices) >= 2
                            validPeakCount = validPeakCount + 1;
                        end
                    end
                end
            end
            fprintf('DEBUG: %d waveforms have sufficient peaks for analysis\n', validPeakCount);

            % Use statistical analysis to determine optimal alignment reference
            [alignmentReference, analysisResults] = RobustPeakAlignmentAnalyzer.determineOptimalAlignment(peakData, t, true);

            % Save to cache
            fprintf('Saving robust alignment cache: %s\n', alignCacheFile);
            save(alignCacheFile, 'alignmentReference', 'analysisResults', 'alignCacheKey', '-v7.3');
        end

        % Extract robust alignment parameters
        referenceTime = alignmentReference.referenceTime;
        optimalPeakPosition = alignmentReference.peakPosition;

        % Determine alignment type based on analysis results
        if strcmp(lower(peakDetectionType), 'peaks')
            alignmentType = sprintf('robust_peak_%d', optimalPeakPosition);
            fprintf('✅ ROBUST ALIGNMENT: Using peak position %d (statistically optimal)\n', optimalPeakPosition);
        elseif strcmp(lower(peakDetectionType), 'valleys')
            alignmentType = sprintf('robust_valley_%d', optimalPeakPosition);
            fprintf('✅ ROBUST ALIGNMENT: Using valley position %d (statistically optimal)\n', optimalPeakPosition);
        else % 'both'
            alignmentType = sprintf('robust_peak_%d', optimalPeakPosition);
            fprintf('✅ ROBUST ALIGNMENT: Using peak position %d from combined peak/valley data\n', optimalPeakPosition);
        end

        fprintf('  → Reference time: %.3f μs (consensus: %.1f%%)\n', ...
            referenceTime * 1e6, alignmentReference.consensusStrength * 100);
        fprintf('  → Quality score: %.3f/1.0\n', analysisResults.qualityMetrics.qualityScore);

        robustAlignmentUsed = true;

    catch ME
        % Fallback to original method if robust analysis fails
        fprintf('⚠️  Robust alignment analysis failed: %s\n', ME.message);
        fprintf('⚠️  Error occurred in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        fprintf('⚠️  This explains why you''re still seeing misaligned peaks!\n');
        fprintf('Falling back to original first-feature alignment method...\n');

        if strcmp(lower(peakDetectionType), 'peaks')
            alignmentType = 'first_peak';
            fprintf('Aligning all waveforms by their first peak...\n');
        elseif strcmp(lower(peakDetectionType), 'valleys')
            alignmentType = 'first_valley';
            fprintf('Aligning all waveforms by their first valley...\n');
        else % 'both'
            alignmentType = 'first_peak'; % Default to first peak when both are detected
            fprintf('Aligning all waveforms by their first peak (both peaks and valleys detected)...\n');
        end

        optimalPeakPosition = 1; % Use first feature as fallback
        robustAlignmentUsed = false;
    end

    % Step 2: Find alignment reference time for each waveform
    alignmentTimes = zeros(numWaveforms, 1);
    validWaveforms = false(numWaveforms, 1);

    for i = 1:numWaveforms
        if ~isempty(peakData{i})
            transitions = peakData{i};

            if robustAlignmentUsed
                % ROBUST ALIGNMENT: Use statistically-determined optimal peak position with amplitude threshold
                if contains(alignmentType, 'robust_peak')
                    % Find peaks (positive transitions) above amplitude threshold
                    peakIndices = find(transitions.TransitionType == 1);
                    if ~isempty(peakIndices)
                        % Apply amplitude threshold filter with debugging (NO abs() for peaks - they should be positive)
                        peakAmplitudes = transitions.TransitionAmplitude(peakIndices);
                        validPeakIndices = peakIndices(peakAmplitudes >= firstPeakAmplitudeThreshold);

                        % Debug output for first few waveforms
                        if i <= 3
                            fprintf('    Robust Peak Waveform %d: threshold=%.3f, %d/%d peaks pass\n', ...
                                    i, firstPeakAmplitudeThreshold, length(validPeakIndices), length(peakIndices));
                        end

                        if length(validPeakIndices) >= optimalPeakPosition
                            % Use the optimal peak position (not always first) from filtered peaks
                            transitionIdx = validPeakIndices(optimalPeakPosition);
                            alignmentTimes(i) = transitions.TransitionTime(transitionIdx);
                            validWaveforms(i) = true;

                            % Debug: Show which amplitude was used
                            if i <= 3
                                usedAmplitude = transitions.TransitionAmplitude(transitionIdx);
                                fprintf('      → Using robust peak #%d with amplitude %.3f\n', optimalPeakPosition, usedAmplitude);
                            end
                        end
                    end
                elseif contains(alignmentType, 'robust_valley')
                    % Find valleys (negative transitions) above amplitude threshold
                    valleyIndices = find(transitions.TransitionType == -1);
                    if ~isempty(valleyIndices)
                        % Apply amplitude threshold filter with debugging (abs() for valleys - they should be negative)
                        valleyAmplitudes = transitions.TransitionAmplitude(valleyIndices);
                        validValleyIndices = valleyIndices(abs(valleyAmplitudes) >= firstPeakAmplitudeThreshold);

                        % Debug output for first few waveforms
                        if i <= 3
                            fprintf('    Robust Valley Waveform %d: threshold=%.3f, %d/%d valleys pass\n', ...
                                    i, firstPeakAmplitudeThreshold, length(validValleyIndices), length(valleyIndices));
                        end

                        if length(validValleyIndices) >= optimalPeakPosition
                            % Use the optimal valley position from filtered valleys
                            transitionIdx = validValleyIndices(optimalPeakPosition);
                            alignmentTimes(i) = transitions.TransitionTime(transitionIdx);
                            validWaveforms(i) = true;

                            % Debug: Show which amplitude was used
                            if i <= 3
                                usedAmplitude = transitions.TransitionAmplitude(transitionIdx);
                                fprintf('      → Using robust valley #%d with amplitude %.3f\n', optimalPeakPosition, usedAmplitude);
                            end
                        end
                    end
                end
            else
                % FALLBACK: Original first-feature alignment with amplitude threshold
                if strcmp(alignmentType, 'first_peak')
                    % Find first peak (positive transition) above amplitude threshold
                    peakIndices = find(transitions.TransitionType == 1);
                    if ~isempty(peakIndices)
                        % Apply amplitude threshold filter with detailed debugging (NO abs() for peaks)
                        peakAmplitudes = transitions.TransitionAmplitude(peakIndices);
                        thresholdMask = peakAmplitudes >= firstPeakAmplitudeThreshold;
                        validPeakIndices = peakIndices(thresholdMask);

                        % Debug output for first few waveforms
                        if i <= 3
                            fprintf('    Waveform %d PEAK DEBUG:\n', i);
                            fprintf('      Threshold: %.3f (for positive peaks)\n', firstPeakAmplitudeThreshold);
                            fprintf('      Peak amplitudes: [%s]\n', sprintf('%.3f ', peakAmplitudes));
                            fprintf('      Threshold mask: [%s]\n', sprintf('%d ', thresholdMask));
                            fprintf('      Valid indices: %d peaks pass threshold\n', length(validPeakIndices));
                        end

                        if ~isempty(validPeakIndices)
                            alignmentTimes(i) = transitions.TransitionTime(validPeakIndices(1));
                            validWaveforms(i) = true;

                            % Debug: Show which amplitude was actually used for alignment
                            if i <= 3
                                usedAmplitude = transitions.TransitionAmplitude(validPeakIndices(1));
                                fprintf('      → Using peak with amplitude %.3f for alignment\n', usedAmplitude);
                            end
                        else
                            % Debug: Show why this waveform was rejected
                            if i <= 3
                                fprintf('      → No peaks meet threshold %.3f, waveform rejected\n', firstPeakAmplitudeThreshold);
                            end
                        end
                    end
                else % 'first_valley'
                    % Find first valley (negative transition) above amplitude threshold
                    valleyIndices = find(transitions.TransitionType == -1);
                    if ~isempty(valleyIndices)
                        % Apply amplitude threshold filter with detailed debugging
                        valleyAmplitudes = transitions.TransitionAmplitude(valleyIndices);
                        valleyAbsAmplitudes = abs(valleyAmplitudes);
                        thresholdMask = valleyAbsAmplitudes >= firstPeakAmplitudeThreshold;
                        validValleyIndices = valleyIndices(thresholdMask);

                        % Debug output for first few waveforms
                        if i <= 3
                            fprintf('    Waveform %d VALLEY DEBUG:\n', i);
                            fprintf('      Threshold: %.3f (applied to |valley amplitude|)\n', firstPeakAmplitudeThreshold);
                            fprintf('      Valley amplitudes: [%s]\n', sprintf('%.3f ', valleyAmplitudes));
                            fprintf('      Abs amplitudes: [%s]\n', sprintf('%.3f ', abs(valleyAmplitudes)));
                            fprintf('      Threshold mask: [%s]\n', sprintf('%d ', thresholdMask));
                            fprintf('      Valid indices: %d valleys pass threshold\n', length(validValleyIndices));
                        end

                        if ~isempty(validValleyIndices)
                            alignmentTimes(i) = transitions.TransitionTime(validValleyIndices(1));
                            validWaveforms(i) = true;

                            % Debug: Show which amplitude was actually used for alignment
                            if i <= 3
                                usedAmplitude = transitions.TransitionAmplitude(validValleyIndices(1));
                                fprintf('      → Using valley with amplitude %.3f for alignment\n', usedAmplitude);
                            end
                        else
                            % Debug: Show why this waveform was rejected
                            if i <= 3
                                fprintf('      → No valleys meet threshold %.3f, waveform rejected\n', firstPeakAmplitudeThreshold);
                            end
                        end
                    end
                end
            end
        end
    end

    numValidWaveforms = sum(validWaveforms);
    fprintf('Found %d waveforms with valid alignment features out of %d total.\n', numValidWaveforms, numWaveforms);

    if firstPeakAmplitudeThreshold > 0
        numFilteredOut = numWaveforms - numValidWaveforms;
        if numFilteredOut > 0
            fprintf('  → %d waveforms filtered out due to amplitude threshold (%.3f)\n', numFilteredOut, firstPeakAmplitudeThreshold);
        else
            fprintf('  → All waveforms met the amplitude threshold requirement\n');
        end

        % Debug: Check some sample amplitudes to verify threshold is working
        if numValidWaveforms > 0
            sampleCount = min(5, numValidWaveforms);
            fprintf('  → Debug: Sample first peak amplitudes from %d waveforms:\n', sampleCount);
            debugCount = 0;
            for i = 1:numWaveforms
                if validWaveforms(i) && debugCount < sampleCount
                    debugCount = debugCount + 1;
                    transitions = peakData{i};
                    if strcmp(alignmentType, 'first_peak')
                        peakIndices = find(transitions.TransitionType == 1);
                        if ~isempty(peakIndices)
                            validPeakIndices = peakIndices(abs(transitions.TransitionAmplitude(peakIndices)) >= firstPeakAmplitudeThreshold);
                            if ~isempty(validPeakIndices)
                                firstPeakAmp = transitions.TransitionAmplitude(validPeakIndices(1));
                                fprintf('    Waveform %d: First valid peak amplitude = %.3f\n', i, firstPeakAmp);
                            end
                        end
                    else
                        valleyIndices = find(transitions.TransitionType == -1);
                        if ~isempty(valleyIndices)
                            validValleyIndices = valleyIndices(abs(transitions.TransitionAmplitude(valleyIndices)) >= firstPeakAmplitudeThreshold);
                            if ~isempty(validValleyIndices)
                                firstValleyAmp = transitions.TransitionAmplitude(validValleyIndices(1));
                                fprintf('    Waveform %d: First valid valley amplitude = %.3f\n', i, firstValleyAmp);
                            end
                        end
                    end
                end
            end
        end
    end
    fprintf('  → Alignment synchronizes waveforms to a common time reference for better comparison.\n');

    if numValidWaveforms == 0
        fprintf('No waveforms found with alignment features. Cannot create 3D plot.\n');
        return;
    end

    % Validate spatial coordinate arrays
    if length(X_Coordinates) ~= numX_sub || length(Y_Coordinates) ~= numY_sub
        fprintf('Error: Spatial coordinate array sizes do not match expected dimensions.\n');
        fprintf('X_Coordinates length: %d, expected: %d\n', length(X_Coordinates), numX_sub);
        fprintf('Y_Coordinates length: %d, expected: %d\n', length(Y_Coordinates), numY_sub);
        return;
    end


    % Step 3: Calculate alignment offset
    if robustAlignmentUsed
        % Use the statistically-determined reference time from robust analysis
        % (referenceTime was already set during robust analysis)
        fprintf('Using robust alignment reference time: %.3f μs\n', referenceTime * 1e6);
    else
        % Fallback: align to earliest first peak/valley
        referenceTime = min(alignmentTimes(validWaveforms));
        fprintf('Using earliest feature time as reference: %.3f μs\n', referenceTime * 1e6);
    end
    if robustAlignmentUsed
        fprintf('  → Robust alignment avoids spurious early peaks by using statistical analysis.\n');
        fprintf('  → Peak position %d was determined to be the most consistent alignment reference.\n', optimalPeakPosition);
        fprintf('  → All waveforms will be time-shifted so their optimal peak aligns to this reference.\n');
        fprintf('  → This creates a physically-meaningful baseline for comparing peak patterns.\n');
    else
        fprintf('  → This is the earliest first peak/valley time across all waveforms.\n');
        fprintf('  → All waveforms will be time-shifted so their first peak/valley aligns to this reference.\n');
        fprintf('  → This creates a common time baseline for comparing peak/valley patterns.\n');
    end

    % Step 4: Collect all aligned peak/valley data
    % First, count total number of transitions to pre-allocate accurately
    totalTransitions = 0;
    validWaveformIndices = find(validWaveforms);

    for i = 1:length(validWaveformIndices)
        waveformIdx = validWaveformIndices(i);
        if ~isempty(peakData{waveformIdx}) && istable(peakData{waveformIdx})
            totalTransitions = totalTransitions + height(peakData{waveformIdx});
        end
    end

    if totalTransitions == 0
        fprintf('No valid transitions found for 3D plotting.\n');
        return;
    end

    fprintf('Pre-allocating arrays for %d total transitions...\n', totalTransitions);

    % Pre-allocate arrays with exact size needed
    allPeakX = zeros(totalTransitions, 1);
    allPeakY = zeros(totalTransitions, 1);
    allPeakZ = zeros(totalTransitions, 1); % Time values (aligned)
    allPeakAmps = zeros(totalTransitions, 1);
    allPeakTypes = zeros(totalTransitions, 1);

    currentIndex = 0;

    for i = 1:length(validWaveformIndices)
        waveformIdx = validWaveformIndices(i);
        transitions = peakData{waveformIdx};

        % Calculate spatial coordinates for this waveform with bounds checking
        if waveformIdx < 1 || waveformIdx > (numY_sub * numX_sub)
            fprintf('Warning: Waveform index %d is out of bounds (max: %d). Skipping.\n', waveformIdx, numY_sub * numX_sub);
            continue;
        end

        [yIdx, xIdx] = ind2sub([numY_sub, numX_sub], waveformIdx);

        % Additional validation of calculated indices
        if xIdx < 1 || xIdx > numX_sub || yIdx < 1 || yIdx > numY_sub
            fprintf('Warning: Waveform %d has invalid spatial indices (x=%d, y=%d). Skipping.\n', waveformIdx, xIdx, yIdx);
            continue;
        end

        % Validate coordinate arrays have sufficient elements
        if xIdx > length(X_Coordinates) || yIdx > length(Y_Coordinates)
            fprintf('Warning: Coordinate arrays too small for indices (x=%d, y=%d). Skipping.\n', xIdx, yIdx);
            continue;
        end

        spatialX = X_Coordinates(xIdx);
        spatialY = Y_Coordinates(yIdx);

        % Calculate time offset for alignment
        timeOffset = alignmentTimes(waveformIdx) - referenceTime;

        % Add all peaks/valleys from this waveform
        numTransitions = height(transitions);
        for j = 1:numTransitions
            currentIndex = currentIndex + 1;
            allPeakX(currentIndex) = spatialX;
            allPeakY(currentIndex) = spatialY;

            % Apply alignment offset to time
            alignedTime = transitions.TransitionTime(j) - timeOffset;
            allPeakZ(currentIndex) = alignedTime * 1e6; % Convert to μs for plotting

            allPeakAmps(currentIndex) = transitions.TransitionAmplitude(j);
            allPeakTypes(currentIndex) = transitions.TransitionType(j);
        end
    end

    if currentIndex == 0
        fprintf('No aligned peak data found for 3D plotting.\n');
        return;
    end

    fprintf('Collected %d aligned peaks/valleys for 3D visualization.\n', currentIndex);

    % Step 5: Create enhanced 3D visualization with controls
    create3DVisualizationWithControls(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, ...
                                     peakDetectionType, showPeaksAndValleys, alignmentType, ...
                                     numValidWaveforms, numWaveforms, peakData, ...
                                     X_Coordinates, Y_Coordinates, numY_sub, numX_sub, t);
end

function create3DVisualizationWithControls(allPeakX, allPeakY, allPeakZ, allPeakAmps, allPeakTypes, ...
                                          peakDetectionType, showPeaksAndValleys, alignmentType, ...
                                          numValidWaveforms, numWaveforms, peakData, ...
                                          X_Coordinates, Y_Coordinates, numY_sub, numX_sub, t)
    % Enhanced 3D visualization with interactive controls and comprehensive data validation

    fprintf('Creating enhanced 3D visualization...\n');
    fprintf('Data points: %d peaks/valleys\n', length(allPeakX));

    % Validate input data consistency
    dataLengths = [length(allPeakX), length(allPeakY), length(allPeakZ), ...
                   length(allPeakAmps), length(allPeakTypes)];
    if ~all(dataLengths == dataLengths(1))
        error('Input data arrays have inconsistent lengths: X=%d, Y=%d, Z=%d, Amps=%d, Types=%d', ...
              dataLengths(1), dataLengths(2), dataLengths(3), dataLengths(4), dataLengths(5));
    end

    % Validate spatial coordinate arrays
    if length(X_Coordinates) ~= numX_sub || length(Y_Coordinates) ~= numY_sub
        error('Spatial coordinate array sizes do not match grid dimensions. X: %d vs %d, Y: %d vs %d', ...
              length(X_Coordinates), numX_sub, length(Y_Coordinates), numY_sub);
    end

    % Create main figure with enhanced layout
    fig = figure('Name', sprintf('Enhanced 3D Peak/Valley Plot - %s', peakDetectionType), ...
                 'Position', [100, 100, 1400, 800], ...
                 'MenuBar', 'figure', ...
                 'ToolBar', 'figure', ...
                 'Color', [0.94, 0.94, 0.94]);

    % Add resize callback to maintain UI positions (responsive UI like XtVsYPlot)
    set(fig, 'ResizeFcn', @(src, event) repositionUIControlsOnResize(src));

    % Set minimum figure size to prevent UI overlap
    set(fig, 'Units', 'pixels');
    currentPos = get(fig, 'Position');
    minWidth = 800;  % Minimum width to accommodate UI panel
    minHeight = 600; % Minimum height for proper layout
    if currentPos(3) < minWidth || currentPos(4) < minHeight
        set(fig, 'Position', [currentPos(1), currentPos(2), max(currentPos(3), minWidth), max(currentPos(4), minHeight)]);
    end

    % Create main plot area (left side) - more space for plot, fixed positioning for controls
    mainAx = axes('Parent', fig, 'Position', [0.05, 0.05, 0.72, 0.9], 'Units', 'normalized');

    % Store initial data in figure UserData with comprehensive validation
    figData = struct();

    % Core data arrays (validated above)
    figData.allPeakX = allPeakX;
    figData.allPeakY = allPeakY;
    figData.allPeakZ = allPeakZ;
    figData.allPeakAmps = allPeakAmps;
    figData.allPeakTypes = allPeakTypes;

    % Analysis parameters
    figData.peakDetectionType = peakDetectionType;
    figData.showPeaksAndValleys = showPeaksAndValleys;
    figData.alignmentType = alignmentType;
    figData.numValidWaveforms = numValidWaveforms;
    figData.numWaveforms = numWaveforms;

    % Axes handles
    figData.mainAx = mainAx;
    figData.leftAx = [];  % Initialize for separate mode
    figData.rightAx = []; % Initialize for separate mode

    % Store data arrays
    figData.X_Coordinates = X_Coordinates;
    figData.Y_Coordinates = Y_Coordinates;
    figData.numY_sub = numY_sub;
    figData.numX_sub = numX_sub;
    figData.t = t; % Store time vector for cross-sectional alignment

    % Default visualization settings (ensure all required fields exist)
    figData.currentMarkerSize = 5; % Start very small
    figData.currentPeakMarker = 'o'; % Circles instead of triangles
    figData.currentValleyMarker = 'o'; % Circles instead of triangles
    figData.currentColormap = 'jet'; % Start with jet colormap

    % Color range settings with proper validation
    if ~isempty(allPeakAmps)
        dataRange = [min(allPeakAmps), max(allPeakAmps)];
        % Ensure valid range
        if dataRange(1) == dataRange(2)
            dataRange = dataRange(1) + [-0.1, 0.1]; % Add small range if all values are identical
        end
        figData.currentColorRange = dataRange;
        fprintf('Data amplitude range: [%.3f, %.3f]\n', dataRange(1), dataRange(2));
    else
        figData.currentColorRange = [-1, 1]; % Default range
        fprintf('No amplitude data available, using default range [-1, 1]\n');
    end

    figData.autoColorRange = true;
    figData.showColorbar = false; % Start with colorbar off
    figData.currentView = [45, 30];
    figData.sliceMode = 'none';
    figData.sliceValue = 0;
    figData.plotMode = 'combined'; % 'combined' or 'separate'

    set(fig, 'UserData', figData);

    % Create right-anchored control panel following MainPlottingApplication pattern
    createRightAnchoredControlPanel(fig);

    % Add diagnostic information BEFORE plotting (in case plotting fails)
    fprintf('\n=== 3D UI Diagnostic Information ===\n');
    fprintf('Figure handle: %s\n', class(fig));
    figPos = get(fig, 'Position');
    fprintf('Figure position: [%d, %d, %d, %d]\n', figPos(1), figPos(2), figPos(3), figPos(4));
    axPos = get(mainAx, 'Position');
    fprintf('Main axes position: [%.2f, %.2f, %.2f, %.2f]\n', axPos(1), axPos(2), axPos(3), axPos(4));
    fprintf('Data points to plot: %d\n', length(allPeakX));

    % Display data size information
    fprintf('Dataset contains %d data points for visualization.\n', length(allPeakX));

    fprintf('If you cannot see the UI controls, try:\n');
    fprintf('1. Resize the figure window\n');
    fprintf('2. Check if controls are outside visible area\n');
    fprintf('3. Run test_3d_ui.m to test UI layout\n');
    fprintf('=====================================\n\n');

    % Initial plot with error handling
    try
        fprintf('Starting 3D plot rendering...\n');
        updatePlot(fig);
        fprintf('3D plot rendering completed successfully.\n');
    catch ME
        fprintf('Error during 3D plot rendering: %s\n', ME.message);
        fprintf('Creating fallback simple plot...\n');
        createFallbackPlot(fig);
    end
end

function isValid = validateFigDataIntegrity(figData)
    % VALIDATEFIGDATAINTEGRITY - Comprehensive validation of figData structure
    % This function ensures all required fields exist and data is consistent

    isValid = false;

    try
        % Check if figData is a structure
        if ~isstruct(figData)
            fprintf('Error: figData is not a structure\n');
            return;
        end

        % Required data fields
        requiredDataFields = {'allPeakX', 'allPeakY', 'allPeakZ', 'allPeakAmps', 'allPeakTypes'};
        for i = 1:length(requiredDataFields)
            if ~isfield(figData, requiredDataFields{i})
                fprintf('Error: Missing required data field: %s\n', requiredDataFields{i});
                return;
            end
        end

        % Check data array consistency
        dataLengths = [length(figData.allPeakX), length(figData.allPeakY), ...
                      length(figData.allPeakZ), length(figData.allPeakAmps), ...
                      length(figData.allPeakTypes)];
        if ~all(dataLengths == dataLengths(1))
            fprintf('Error: Data array length mismatch: X=%d, Y=%d, Z=%d, Amps=%d, Types=%d\n', ...
                    dataLengths(1), dataLengths(2), dataLengths(3), dataLengths(4), dataLengths(5));
            return;
        end

        % Required settings fields
        requiredSettingsFields = {'currentMarkerSize', 'currentPeakMarker', 'currentValleyMarker', ...
                                 'currentColormap', 'currentColorRange', 'autoColorRange', ...
                                 'showColorbar', 'currentView', 'plotMode'};
        for i = 1:length(requiredSettingsFields)
            if ~isfield(figData, requiredSettingsFields{i})
                fprintf('Error: Missing required settings field: %s\n', requiredSettingsFields{i});
                return;
            end
        end

        % Validate spatial data if present
        if isfield(figData, 'X_Coordinates') && isfield(figData, 'Y_Coordinates') && ...
           isfield(figData, 'numX_sub') && isfield(figData, 'numY_sub')
            if length(figData.X_Coordinates) ~= figData.numX_sub || ...
               length(figData.Y_Coordinates) ~= figData.numY_sub
                fprintf('Error: Spatial coordinate array size mismatch\n');
                return;
            end
        end

        % Validate color range
        if length(figData.currentColorRange) ~= 2 || ...
           figData.currentColorRange(1) >= figData.currentColorRange(2)
            fprintf('Error: Invalid color range: [%.3f, %.3f]\n', ...
                    figData.currentColorRange(1), figData.currentColorRange(2));
            return;
        end

        % Validate plot mode
        if ~ismember(figData.plotMode, {'combined', 'separate'})
            fprintf('Error: Invalid plot mode: %s\n', figData.plotMode);
            return;
        end

        % All validations passed
        isValid = true;

    catch ME
        fprintf('Error during figData validation: %s\n', ME.message);
        isValid = false;
    end
end

function createFallbackPlot(fig)
    % Create a simple fallback plot when main plotting fails
    figData = get(fig, 'UserData');

    % Clear and setup axes
    cla(figData.mainAx);
    hold(figData.mainAx, 'on');

    % Use all data for fallback plot
    plotX = figData.allPeakX;
    plotY = figData.allPeakY;
    plotZ = figData.allPeakZ;
    plotAmps = figData.allPeakAmps;
    fprintf('Fallback plot using all %d data points.\n', length(plotX));

    % Simple scatter plot
    scatter3(figData.mainAx, plotX, plotY, plotZ, 20, plotAmps, 'filled');

    % Basic formatting
    xlabel(figData.mainAx, 'X Position (mm)');
    ylabel(figData.mainAx, 'Y Position (mm)');
    zlabel(figData.mainAx, 'Aligned Time (μs)');
    title(figData.mainAx, sprintf('3D Peak/Valley Plot (Fallback - %d points)', length(plotX)));
    colormap(figData.mainAx, 'parula');
    colorbar(figData.mainAx);
    grid(figData.mainAx, 'on');
    view(figData.mainAx, [45, 30]);

    hold(figData.mainAx, 'off');
end

function createRightAnchoredControlPanel(fig)
    % Create right-anchored control panel following MainPlottingApplication pattern

    fprintf('Creating right-anchored control panel...\n');

    % Get figure dimensions for right-anchored positioning
    figPos = get(fig, 'Position');
    figWidth = figPos(3);
    figHeight = figPos(4);

    % Define UI panel constants - RIGHT-TOP ANCHORED SYSTEM
    UI_PANEL_WIDTH = 180;           % Fixed width in pixels
    UI_MARGIN_RIGHT = 10;           % Distance from right edge to right edge of UI controls
    UI_MARGIN_TOP = 10;             % Distance from top edge to first UI element
    UI_ELEMENT_HEIGHT = 25;         % Standard height for dropdowns/buttons
    UI_TEXT_HEIGHT = 18;            % Height for text elements
    UI_CHECKBOX_HEIGHT = 20;        % Height for checkbox elements
    UI_SPACING = 5;                 % Spacing between elements

    % Calculate UI_LEFT so that UI_LEFT + UI_PANEL_WIDTH = figWidth - UI_MARGIN_RIGHT
    % This ensures the RIGHT EDGE of UI controls is exactly UI_MARGIN_RIGHT pixels from figure right edge
    UI_RIGHT_EDGE = figWidth - UI_MARGIN_RIGHT;  % Where the right edge of UI should be
    UI_LEFT = UI_RIGHT_EDGE - UI_PANEL_WIDTH;    % Left position to achieve this

    % Ensure UI_LEFT is always positive (minimum 10 pixels from left edge)
    if UI_LEFT < 10
        UI_LEFT = 10;
        UI_PANEL_WIDTH = UI_RIGHT_EDGE - UI_LEFT; % Adjust width if figure is too narrow
    end

    % Calculate starting Y position from TOP edge (MATLAB uses bottom-left origin)
    UI_TOP_Y = figHeight - UI_MARGIN_TOP - UI_ELEMENT_HEIGHT; % Start from top, accounting for element height

    % Debug output for UI positioning
    fprintf('UI Positioning: figWidth=%.0f, figHeight=%.0f\n', figWidth, figHeight);
    fprintf('UI Panel: X=%.0f, Y_start=%.0f, Width=%.0f\n', UI_LEFT, UI_TOP_Y, UI_PANEL_WIDTH);
    fprintf('✓ UI controls anchored %.0f pixels from right edge, %.0f pixels from top edge\n', UI_MARGIN_RIGHT, UI_MARGIN_TOP);

    % Start positioning from calculated top-anchored position and work down
    currentY = UI_TOP_Y; % Start from calculated top-anchored position

    % Plot Mode Dropdown - RIGHT-ANCHORED
    plotModeDropdown = uicontrol('Parent', fig, 'Style', 'popupmenu', ...
        'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'Combined Plot', 'Separate Plots'}, 'Value', 1, ...
        'Tag', 'PlotModeDropdown', 'FontSize', 9, ...
        'Callback', @(src, event) changePlotMode(src, fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Peak Marker Dropdown - RIGHT-ANCHORED
    uicontrol('Parent', fig, 'Style', 'text', ...
              'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT], ...
              'String', 'Peak Markers:', ...
              'HorizontalAlignment', 'center', ...
              'BackgroundColor', [0.94, 0.94, 0.94], ...
              'FontSize', 10, 'FontWeight', 'bold', ...
              'Tag', 'PeakMarkerLabel');
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;

    peakMarkerDropdown = uicontrol('Parent', fig, 'Style', 'popupmenu', ...
        'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'Triangle ▲', 'Circle ●', 'Square ■', 'Diamond ◆'}, 'Value', 2, ...
        'Tag', 'PeakMarkerDropdown', 'FontSize', 9, ...
        'Callback', @(src, event) changePeakMarkerDropdown(src, fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Valley Marker Dropdown - RIGHT-ANCHORED
    uicontrol('Parent', fig, 'Style', 'text', ...
              'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT], ...
              'String', 'Valley Markers:', ...
              'HorizontalAlignment', 'center', ...
              'BackgroundColor', [0.94, 0.94, 0.94], ...
              'FontSize', 10, 'FontWeight', 'bold', ...
              'Tag', 'ValleyMarkerLabel');
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;

    valleyMarkerDropdown = uicontrol('Parent', fig, 'Style', 'popupmenu', ...
        'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'Triangle ▼', 'Circle ●', 'Square ■', 'Diamond ◆'}, 'Value', 2, ...
        'Tag', 'ValleyMarkerDropdown', 'FontSize', 9, ...
        'Callback', @(src, event) changeValleyMarkerDropdown(src, fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Marker size controls - RIGHT-ANCHORED
    uicontrol('Parent', fig, 'Style', 'text', ...
              'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT], ...
              'String', 'Marker Size:', ...
              'HorizontalAlignment', 'center', ...
              'BackgroundColor', [0.94, 0.94, 0.94], ...
              'FontSize', 10, 'FontWeight', 'bold', ...
              'Tag', 'MarkerSizeLabel');
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;

    % Side-by-side size buttons
    buttonWidth = (UI_PANEL_WIDTH - UI_SPACING) / 2;
    uicontrol('Parent', fig, 'Style', 'pushbutton', ...
              'String', 'Small', ...
              'Position', [UI_LEFT, currentY, buttonWidth, UI_ELEMENT_HEIGHT], ...
              'FontSize', 9, ...
              'Tag', 'SmallSizeButton', ...
              'Callback', @(src, event) changeMarkerSize(fig, 5));

    uicontrol('Parent', fig, 'Style', 'pushbutton', ...
              'String', 'Large', ...
              'Position', [UI_LEFT + buttonWidth + UI_SPACING, currentY, buttonWidth, UI_ELEMENT_HEIGHT], ...
              'FontSize', 9, ...
              'Tag', 'LargeSizeButton', ...
              'Callback', @(src, event) changeMarkerSize(fig, 20));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Colormap Dropdown - RIGHT-ANCHORED (from XtVsYPlot.m)
    uicontrol('Parent', fig, 'Style', 'text', ...
              'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT], ...
              'String', 'Colormaps:', ...
              'HorizontalAlignment', 'center', ...
              'BackgroundColor', [0.94, 0.94, 0.94], ...
              'FontSize', 10, 'FontWeight', 'bold', ...
              'Tag', 'ColormapLabel');
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;

    colormapDropdown = uicontrol('Parent', fig, 'Style', 'popupmenu', ...
        'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'jet', 'parula', 'hsv', 'hot', 'gray', 'bone', 'coolwarm'}, 'Value', 1, ...
        'Tag', 'ColormapDropdown', 'FontSize', 9, ...
        'Callback', @(src, event) changeColormapDropdown(src, fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % View Angle Dropdown - RIGHT-ANCHORED
    uicontrol('Parent', fig, 'Style', 'text', ...
              'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT], ...
              'String', 'View Angles:', ...
              'HorizontalAlignment', 'center', ...
              'BackgroundColor', [0.94, 0.94, 0.94], ...
              'FontSize', 10, 'FontWeight', 'bold', ...
              'Tag', 'ViewLabel');
    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;

    viewDropdown = uicontrol('Parent', fig, 'Style', 'popupmenu', ...
        'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
        'String', {'Default (45°, 30°)', 'Top View (0°, 90°)', 'Side View (0°, 0°)', 'Front View (90°, 0°)', 'Isometric (45°, 45°)', 'Back View (-90°, 0°)'}, 'Value', 1, ...
        'Tag', 'ViewDropdown', 'FontSize', 9, ...
        'Callback', @(src, event) changeViewDropdown(src, fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;

    % Toggle colorbar - RIGHT-ANCHORED
    uicontrol('Parent', fig, 'Style', 'pushbutton', ...
              'String', 'Toggle Colorbar', ...
              'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT], ...
              'FontSize', 9, ...
              'Tag', 'ColorbarToggleButton', ...
              'Callback', @(src, event) toggleColorbarSimple(fig));
    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;


    fprintf('Right-anchored control panel created successfully.\n');
    fprintf('UI elements positioned from right edge at %d pixels.\n', UI_MARGIN_RIGHT);

    % Force figure to refresh and bring to front
    drawnow;
    figure(fig);

    % Add a small delay and then reposition to ensure proper layout
    pause(0.1);
    repositionUIControlsOnResize(fig);
end



% Callback functions for dropdown controls
function changePeakMarkerDropdown(src, fig)
    figData = get(fig, 'UserData');
    markerOptions = {'^', 'o', 's', 'd'}; % Triangle, Circle, Square, Diamond
    selectedMarker = markerOptions{src.Value};
    figData.currentPeakMarker = selectedMarker;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    fprintf('Peak marker changed to: %s\n', selectedMarker);
end

function changeValleyMarkerDropdown(src, fig)
    figData = get(fig, 'UserData');
    markerOptions = {'v', 'o', 's', 'd'}; % Triangle down, Circle, Square, Diamond
    selectedMarker = markerOptions{src.Value};
    figData.currentValleyMarker = selectedMarker;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    fprintf('Valley marker changed to: %s\n', selectedMarker);
end

% Legacy callback functions for button controls (kept for compatibility)
function changePeakMarker(fig, markerShape)
    figData = get(fig, 'UserData');
    figData.currentPeakMarker = markerShape;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    fprintf('Peak marker changed to: %s\n', markerShape);
end

function changeValleyMarker(fig, markerShape)
    figData = get(fig, 'UserData');
    figData.currentValleyMarker = markerShape;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    fprintf('Valley marker changed to: %s\n', markerShape);
end

function changeMarkerSize(fig, newSize)
    figData = get(fig, 'UserData');
    figData.currentMarkerSize = newSize;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    fprintf('Marker size changed to: %d\n', newSize);
end

function changeColormapDropdown(src, fig)
    figData = get(fig, 'UserData');
    colormapOptions = {'jet', 'parula', 'hsv', 'hot', 'gray', 'bone', 'coolwarm'};
    selectedColormap = colormapOptions{src.Value};
    figData.currentColormap = selectedColormap;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    fprintf('Colormap changed to: %s\n', selectedColormap);
end

function changeViewDropdown(src, fig)
    figData = get(fig, 'UserData');
    viewOptions = {[45, 30], [0, 90], [0, 0], [90, 0], [45, 45], [-90, 0]};
    selectedView = viewOptions{src.Value};
    figData.currentView = selectedView;
    set(fig, 'UserData', figData);

    % Apply view to appropriate axes based on plot mode
    if isfield(figData, 'plotMode') && strcmp(figData.plotMode, 'separate')
        if isfield(figData, 'leftAx') && all(ishandle(figData.leftAx))
            view(figData.leftAx, selectedView);
        end
        if isfield(figData, 'rightAx') && all(ishandle(figData.rightAx))
            view(figData.rightAx, selectedView);
        end
    else
        if isfield(figData, 'mainAx') && all(ishandle(figData.mainAx))
            view(figData.mainAx, selectedView);
        end
    end

    drawnow;
    viewNames = {'Default (45°, 30°)', 'Top View (0°, 90°)', 'Side View (0°, 0°)', 'Front View (90°, 0°)', 'Isometric (45°, 45°)', 'Back View (-90°, 0°)'};
    fprintf('View changed to: %s\n', viewNames{src.Value});
end

% Legacy callback functions (kept for compatibility)
function changeColormap(fig, colormapName)
    figData = get(fig, 'UserData');
    figData.currentColormap = colormapName;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    fprintf('Colormap changed to: %s\n', colormapName);
end

function changeView(fig, viewAngles)
    figData = get(fig, 'UserData');
    figData.currentView = viewAngles;
    set(fig, 'UserData', figData);

    % Apply view to appropriate axes based on plot mode
    if isfield(figData, 'plotMode') && strcmp(figData.plotMode, 'separate')
        if isfield(figData, 'leftAx') && all(ishandle(figData.leftAx))
            view(figData.leftAx, viewAngles);
        end
        if isfield(figData, 'rightAx') && all(ishandle(figData.rightAx))
            view(figData.rightAx, viewAngles);
        end
    else
        if isfield(figData, 'mainAx') && all(ishandle(figData.mainAx))
            view(figData.mainAx, viewAngles);
        end
    end

    drawnow;
    fprintf('View changed to: [%d, %d]\n', viewAngles(1), viewAngles(2));
end

function toggleColorbarSimple(fig)
    figData = get(fig, 'UserData');
    figData.showColorbar = ~figData.showColorbar;
    set(fig, 'UserData', figData);
    updatePlotLayout(fig);
    drawnow;
    if figData.showColorbar
        fprintf('Colorbar enabled\n');
    else
        fprintf('Colorbar disabled\n');
    end
end











function updatePlot(fig)
    % Main plotting function that updates the 3D visualization with comprehensive validation

    % Validate figure handle
    if ~isvalid(fig) || ~ishghandle(fig)
        fprintf('Error: Invalid figure handle in updatePlot\n');
        return;
    end

    try
        figData = get(fig, 'UserData');
    catch ME
        fprintf('Error: Cannot access figure UserData (%s)\n', ME.message);
        return;
    end

    % Comprehensive figData validation
    if ~validateFigDataIntegrity(figData)
        fprintf('Error: figData validation failed in updatePlot\n');
        return;
    end

    % Validate main axes exists and is valid
    if ~isfield(figData, 'mainAx') || ~isvalid(figData.mainAx) || ~ishghandle(figData.mainAx)
        fprintf('Warning: Main axes invalid. Recreating axes...\n');
        % Recreate the main axes
        figData.mainAx = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.72, 0.90]);
        set(fig, 'UserData', figData);
    end

    % Clear the main axes
    cla(figData.mainAx);
    hold(figData.mainAx, 'on');


    % Point view - show ONLY peaks/valleys (NO plates, even if they exist)
    % Apply slicing if enabled
    if strcmp(figData.sliceMode, 'none')
        % No slicing - use all data
        plotX = figData.allPeakX;
        plotY = figData.allPeakY;
        plotZ = figData.allPeakZ;
        plotAmps = figData.allPeakAmps;
        plotTypes = figData.allPeakTypes;
    else
        % Apply slicing
        tolerance = 0.1; % Tolerance for slice selection

        switch figData.sliceMode
            case 'x'
                validIndices = abs(figData.allPeakX - figData.sliceValue) <= tolerance;
            case 'y'
                validIndices = abs(figData.allPeakY - figData.sliceValue) <= tolerance;
            case 'z'
                validIndices = abs(figData.allPeakZ - figData.sliceValue) <= tolerance;
        end

        plotX = figData.allPeakX(validIndices);
        plotY = figData.allPeakY(validIndices);
        plotZ = figData.allPeakZ(validIndices);
        plotAmps = figData.allPeakAmps(validIndices);
        plotTypes = figData.allPeakTypes(validIndices);
    end

    % Check if we have data to plot
    if isempty(plotX)
        text(figData.mainAx, 0.5, 0.5, 0.5, 'No data in current slice', ...
             'HorizontalAlignment', 'center', 'FontSize', 14);
        return;
    end

    % Keep all data points - no reduction applied
    fprintf('Plotting all %d data points without reduction.\n', length(plotX));

    % Determine color range
    if figData.autoColorRange
        colorRange = [min(plotAmps), max(plotAmps)];
        % Update the manual controls to show current range if they exist
        if isfield(figData, 'controls') && isfield(figData.controls, 'colorMinEdit')
            set(figData.controls.colorMinEdit, 'String', sprintf('%.3f', colorRange(1)));
            set(figData.controls.colorMaxEdit, 'String', sprintf('%.3f', colorRange(2)));
        end
    else
        colorRange = figData.currentColorRange;
    end

    % Plot based on detection type and settings
    if figData.showPeaksAndValleys && strcmp(lower(figData.peakDetectionType), 'both')
        % Plot peaks and valleys with different markers
        peakIndices = plotTypes == 1;
        valleyIndices = plotTypes == -1;

        if any(peakIndices)
            scatter3(figData.mainAx, plotX(peakIndices), plotY(peakIndices), plotZ(peakIndices), ...
                    figData.currentMarkerSize, plotAmps(peakIndices), 'filled', figData.currentPeakMarker, ...
                    'DisplayName', 'Peaks');
        end

        if any(valleyIndices)
            scatter3(figData.mainAx, plotX(valleyIndices), plotY(valleyIndices), plotZ(valleyIndices), ...
                    figData.currentMarkerSize, plotAmps(valleyIndices), 'filled', figData.currentValleyMarker, ...
                    'DisplayName', 'Valleys');
        end

        if any(peakIndices) && any(valleyIndices)
            legend(figData.mainAx, 'show', 'Location', 'best');
        end

    else
        % Plot all with same marker type
        if strcmp(lower(figData.peakDetectionType), 'peaks') || all(plotTypes == 1)
            markerShape = figData.currentPeakMarker;
        else
            markerShape = figData.currentValleyMarker;
        end

        scatter3(figData.mainAx, plotX, plotY, plotZ, ...
                figData.currentMarkerSize, plotAmps, 'filled', markerShape);
    end

    % Set colormap and color range
    if strcmp(figData.currentColormap, 'coolwarm')
        colormap(figData.mainAx, coolwarm());
    else
        colormap(figData.mainAx, figData.currentColormap);
    end
    caxis(figData.mainAx, colorRange);

    % Show/hide colorbar
    if figData.showColorbar
        colorbar(figData.mainAx);
    else
        colorbar(figData.mainAx, 'off');
    end

    % NOTE: Plates are NOT visualized in point view mode
    % Plates are only shown when explicitly in plate view mode via plotPlatesOnly()

    % Set labels and title
    xlabel(figData.mainAx, 'X Position (mm)');
    ylabel(figData.mainAx, 'Y Position (mm)');
    zlabel(figData.mainAx, 'Aligned Time (μs)');

    % Create title based on slicing mode and plates
    if strcmp(figData.sliceMode, 'none')
        titleStr = sprintf('3D Aligned %s in Space-Time', figData.peakDetectionType);
    else
        titleStr = sprintf('3D %s - %s Slice at %.2f', figData.peakDetectionType, ...
                          upper(figData.sliceMode), figData.sliceValue);
    end

    title(figData.mainAx, titleStr);

    % Set view and grid
    view(figData.mainAx, figData.currentView);
    grid(figData.mainAx, 'on');

    hold(figData.mainAx, 'off');
end

function changePlotMode(src, fig)
    % Change between combined and separate plot modes with data validation

    try
        figData = get(fig, 'UserData');
    catch ME
        fprintf('Error accessing figData in changePlotMode: %s\n', ME.message);
        return;
    end

    selectedMode = src.String{src.Value};

    % Validate data availability for separate mode
    if strcmp(selectedMode, 'Separate Plots')
        % Check if we have both peaks and valleys for meaningful separation
        if isfield(figData, 'allPeakTypes') && ~isempty(figData.allPeakTypes)
            numPeaks = sum(figData.allPeakTypes == 1);
            numValleys = sum(figData.allPeakTypes == -1);

            if numPeaks == 0 && numValleys == 0
                fprintf('Cannot switch to Separate Plots: No peak or valley data available.\n');
                % Reset dropdown to combined mode
                set(src, 'Value', 1);
                return;
            elseif numPeaks == 0
                fprintf('Warning: Switching to Separate Plots with no peaks (only %d valleys)\n', numValleys);
            elseif numValleys == 0
                fprintf('Warning: Switching to Separate Plots with no valleys (only %d peaks)\n', numPeaks);
            else
                fprintf('Switching to Separate Plots: %d peaks, %d valleys\n', numPeaks, numValleys);
            end
        else
            fprintf('Warning: No peak/valley type data available for separate plotting\n');
        end

        figData.plotMode = 'separate';
    else
        figData.plotMode = 'combined';
    end

    set(fig, 'UserData', figData);

    % Clear existing axes and recreate layout
    clearExistingAxes(fig);

    % Update the plot layout
    updatePlotLayout(fig);
    drawnow;
    fprintf('Plot mode changed to: %s\n', selectedMode);
end

function clearExistingAxes(fig)
    % Clear all existing axes before creating new layout with comprehensive validation

    try
        figData = get(fig, 'UserData');
    catch ME
        fprintf('Error accessing figData in clearExistingAxes: %s\n', ME.message);
        return;
    end

    % Safely delete existing axes with proper validation
    if isfield(figData, 'mainAx') && ~isempty(figData.mainAx)
        try
            if isvalid(figData.mainAx) && ishghandle(figData.mainAx)
                delete(figData.mainAx);
                fprintf('Deleted main axes\n');
            end
        catch ME
            fprintf('Warning: Could not delete main axes (%s)\n', ME.message);
        end
    end

    if isfield(figData, 'leftAx') && ~isempty(figData.leftAx)
        try
            if isvalid(figData.leftAx) && ishghandle(figData.leftAx)
                delete(figData.leftAx);
                fprintf('Deleted left axes\n');
            end
        catch ME
            fprintf('Warning: Could not delete left axes (%s)\n', ME.message);
        end
    end

    if isfield(figData, 'rightAx') && ~isempty(figData.rightAx)
        try
            if isvalid(figData.rightAx) && ishghandle(figData.rightAx)
                delete(figData.rightAx);
                fprintf('Deleted right axes\n');
            end
        catch ME
            fprintf('Warning: Could not delete right axes (%s)\n', ME.message);
        end
    end

    % Clear axes references
    figData.mainAx = [];
    figData.leftAx = [];
    figData.rightAx = [];

    % Update figData
    try
        set(fig, 'UserData', figData);
    catch ME
        fprintf('Error updating figData in clearExistingAxes: %s\n', ME.message);
    end
end

function updatePlotLayout(fig)
    % Update plot layout based on plot mode with comprehensive validation

    % Validate figure handle
    if ~isvalid(fig) || ~ishghandle(fig)
        fprintf('Error: Invalid figure handle in updatePlotLayout\n');
        return;
    end

    try
        figData = get(fig, 'UserData');
    catch ME
        fprintf('Error: Cannot access figure UserData (%s)\n', ME.message);
        return;
    end

    % Validate figData structure
    if ~isstruct(figData)
        fprintf('Error: figData is not a valid structure\n');
        return;
    end

    % Ensure plotMode field exists
    if ~isfield(figData, 'plotMode')
        fprintf('Warning: plotMode not found in figData. Defaulting to combined.\n');
        figData.plotMode = 'combined';
        set(fig, 'UserData', figData);
    end

    fprintf('Updating plot layout: %s mode\n', figData.plotMode);


    % Execute the appropriate plotting mode with comprehensive validation
    if strcmp(figData.plotMode, 'separate')
        % Create separate plots for peaks and valleys
        needToCreateAxes = true;

        % Check if axes exist and are valid
        if isfield(figData, 'leftAx') && isfield(figData, 'rightAx') && ...
           ~isempty(figData.leftAx) && ~isempty(figData.rightAx)
            try
                if isvalid(figData.leftAx) && ishghandle(figData.leftAx) && ...
                   isvalid(figData.rightAx) && ishghandle(figData.rightAx)
                    needToCreateAxes = false;
                end
            catch
                % Axes are invalid, need to recreate
                needToCreateAxes = true;
            end
        end

        if needToCreateAxes
            fprintf('Creating separate plot axes...\n');
            createSeparatePlots(fig);
            % Refresh figData after axes creation
            try
                figData = get(fig, 'UserData');
            catch ME
                fprintf('Error refreshing figData after axes creation: %s\n', ME.message);
                return;
            end
        end

        % Now plot data on the separate axes with validation
        if isfield(figData, 'leftAx') && isfield(figData, 'rightAx') && ...
           ~isempty(figData.leftAx) && ~isempty(figData.rightAx)
            try
                if isvalid(figData.leftAx) && ishghandle(figData.leftAx) && ...
                   isvalid(figData.rightAx) && ishghandle(figData.rightAx)
                    fprintf('Updating separate plots...\n');
                    plotSeparateData(fig, figData.leftAx, figData.rightAx);
                else
                    fprintf('Error: Separate axes are invalid after creation\n');
                end
            catch ME
                fprintf('Error validating separate axes: %s\n', ME.message);
            end
        else
            fprintf('Error: Separate axes still not available after creation\n');
        end
    else
        % Use combined plot (existing behavior) with validation
        needToCreateMainAx = true;

        if isfield(figData, 'mainAx') && ~isempty(figData.mainAx)
            try
                if isvalid(figData.mainAx) && ishghandle(figData.mainAx)
                    needToCreateMainAx = false;
                end
            catch
                needToCreateMainAx = true;
            end
        end

        if needToCreateMainAx
            fprintf('Creating combined plot axes...\n');
            createCombinedPlot(fig);
        else
            fprintf('Updating existing combined plot...\n');
            updatePlot(fig);
        end
    end

    fprintf('Plot layout update complete.\n');
end

function createCombinedPlot(fig)
    % Create single combined plot axes
    figData = get(fig, 'UserData');

    % Create main plot area (left side) - more space for plot, fixed positioning for controls
    figData.mainAx = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.72, 0.9]);

    set(fig, 'UserData', figData);
    updatePlot(fig);
end

function createSeparatePlots(fig)
    % Create separate side-by-side plots for peaks and valleys with full data validation
    figData = get(fig, 'UserData');

    % Validate figData structure
    if ~isstruct(figData)
        fprintf('Error: Invalid figData structure in createSeparatePlots\n');
        return;
    end

    % Safely clear existing main axes if it exists
    if isfield(figData, 'mainAx') && ~isempty(figData.mainAx) && all(ishandle(figData.mainAx))
        try
            delete(figData.mainAx);
            fprintf('Cleared existing main axes for separate plot mode\n');
        catch ME
            fprintf('Warning: Failed to delete main axes (%s)\n', ME.message);
        end
    end

    % Clear any existing separate axes
    if isfield(figData, 'leftAx') && ~isempty(figData.leftAx) && all(ishandle(figData.leftAx))
        try
            delete(figData.leftAx);
        catch
            % Ignore deletion errors
        end
    end

    if isfield(figData, 'rightAx') && ~isempty(figData.rightAx) && all(ishandle(figData.rightAx))
        try
            delete(figData.rightAx);
        catch
            % Ignore deletion errors
        end
    end

    fprintf('Creating separate side-by-side axes for peaks and valleys...\n');

    % Create two side-by-side axes with proper spacing for UI controls
    % Left plot: peaks
    leftAx = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.05, 0.05, 0.33, 0.9]);
    set(leftAx, 'Tag', 'LeftPeakAxes');

    % Right plot: valleys
    rightAx = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.39, 0.05, 0.33, 0.9]);
    set(rightAx, 'Tag', 'RightValleyAxes');

    % Validate axes were created successfully
    if ~all(ishghandle(leftAx)) || ~all(ishghandle(rightAx))
        fprintf('Error: Failed to create separate plot axes\n');
        return;
    end

    % Update figData with new axes and ensure all required fields exist
    figData.leftAx = leftAx;
    figData.rightAx = rightAx;
    figData.mainAx = leftAx; % Keep mainAx for compatibility with existing functions

    % Ensure autoColorRange field exists for proper color scaling
    if ~isfield(figData, 'autoColorRange')
        figData.autoColorRange = true;
    end

    % Ensure currentColorRange exists
    if ~isfield(figData, 'currentColorRange')
        figData.currentColorRange = [-1, 1];
    end

    % Store updated figData
    set(fig, 'UserData', figData);

    fprintf('Separate axes created successfully. Plotting data...\n');

    % Plot peaks on left, valleys on right with full data validation
    plotSeparateData(fig, leftAx, rightAx);

    fprintf('Separate plots creation complete.\n');
end

function plotSeparateData(fig, leftAx, rightAx)
    % Plot peaks and valleys on separate axes with full data validation and synchronization

    % Validate figure handle
    if ~isvalid(fig) || ~ishghandle(fig)
        fprintf('Error: Invalid figure handle in plotSeparateData\n');
        return;
    end

    try
        figData = get(fig, 'UserData');
    catch ME
        fprintf('Error: Cannot access figure UserData (%s)\n', ME.message);
        return;
    end

    % Use comprehensive figData validation
    if ~validateFigDataIntegrity(figData)
        fprintf('Error: figData validation failed in plotSeparateData\n');
        return;
    end

    % Comprehensive axes validation
    if isempty(leftAx) || isempty(rightAx) || ~all(ishghandle(leftAx)) || ~all(ishghandle(rightAx))
        fprintf('Error: Separate plot axes invalid. Cannot plot separate data.\n');
        fprintf('Left axis valid: %s, Right axis valid: %s\n', ...
                mat2str(~isempty(leftAx) && all(ishghandle(leftAx))), ...
                mat2str(~isempty(rightAx) && all(ishghandle(rightAx))));
        return;
    end

    % Get and validate data arrays
    allPeakX = figData.allPeakX;
    allPeakY = figData.allPeakY;
    allPeakZ = figData.allPeakZ;
    allPeakAmps = figData.allPeakAmps;
    allPeakTypes = figData.allPeakTypes;

    % Validate data array consistency
    dataLengths = [length(allPeakX), length(allPeakY), length(allPeakZ), ...
                   length(allPeakAmps), length(allPeakTypes)];
    if ~all(dataLengths == dataLengths(1))
        fprintf('Error: Data array length mismatch in separate plot data.\n');
        fprintf('Lengths: X=%d, Y=%d, Z=%d, Amps=%d, Types=%d\n', dataLengths);
        return;
    end

    if isempty(allPeakX)
        fprintf('Warning: No data available for separate plotting.\n');
        % Clear both axes and show message
        cla(leftAx);
        cla(rightAx);
        text(leftAx, 0.5, 0.5, 0.5, 'No Peak Data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        text(rightAx, 0.5, 0.5, 0.5, 'No Valley Data', 'HorizontalAlignment', 'center', 'FontSize', 14);
        return;
    end

    fprintf('Plotting separate data: %d total points\n', length(allPeakX));


    % Point visualization mode: separate peaks and valleys with validation
    peakIndices = allPeakTypes == 1;
    valleyIndices = allPeakTypes == -1;

    numPeaks = sum(peakIndices);
    numValleys = sum(valleyIndices);
    fprintf('  Peaks: %d, Valleys: %d\n', numPeaks, numValleys);

    % Determine color range (same logic as combined plot)
    if figData.autoColorRange
        colorRange = [min(allPeakAmps), max(allPeakAmps)];
    else
        colorRange = figData.currentColorRange;
    end

    % Plot peaks on left axes with comprehensive setup
    cla(leftAx);
    hold(leftAx, 'on');

    if numPeaks > 0
        scatter3(leftAx, allPeakX(peakIndices), allPeakY(peakIndices), allPeakZ(peakIndices), ...
                figData.currentMarkerSize, allPeakAmps(peakIndices), 'filled', figData.currentPeakMarker);
        fprintf('  Left axis: Plotted %d peaks\n', numPeaks);
    else
        text(leftAx, 0.5, 0.5, 0.5, 'No Peaks Found', 'HorizontalAlignment', 'center', 'FontSize', 12);
        fprintf('  Left axis: No peaks to plot\n');
    end

    % Apply identical formatting to left axis
    xlabel(leftAx, 'X Position (mm)');
    ylabel(leftAx, 'Y Position (mm)');
    zlabel(leftAx, 'Aligned Time (μs)');
    title(leftAx, sprintf('Peaks (%d points)', numPeaks));

    % Apply colormap and color range
    if strcmp(figData.currentColormap, 'coolwarm')
        colormap(leftAx, coolwarm());
    else
        colormap(leftAx, figData.currentColormap);
    end
    caxis(leftAx, colorRange);

    % Show/hide colorbar
    if figData.showColorbar
        colorbar(leftAx);
    else
        colorbar(leftAx, 'off');
    end

    grid(leftAx, 'on');
    view(leftAx, figData.currentView);
    hold(leftAx, 'off');

    % Plot valleys on right axes with identical setup
    cla(rightAx);
    hold(rightAx, 'on');

    if numValleys > 0
        scatter3(rightAx, allPeakX(valleyIndices), allPeakY(valleyIndices), allPeakZ(valleyIndices), ...
                figData.currentMarkerSize, allPeakAmps(valleyIndices), 'filled', figData.currentValleyMarker);
        fprintf('  Right axis: Plotted %d valleys\n', numValleys);
    else
        text(rightAx, 0.5, 0.5, 0.5, 'No Valleys Found', 'HorizontalAlignment', 'center', 'FontSize', 12);
        fprintf('  Right axis: No valleys to plot\n');
    end

    % Apply identical formatting to right axis
    xlabel(rightAx, 'X Position (mm)');
    ylabel(rightAx, 'Y Position (mm)');
    zlabel(rightAx, 'Aligned Time (μs)');
    title(rightAx, sprintf('Valleys (%d points)', numValleys));

    % Apply identical colormap and color range
    if strcmp(figData.currentColormap, 'coolwarm')
        colormap(rightAx, coolwarm());
    else
        colormap(rightAx, figData.currentColormap);
    end
    caxis(rightAx, colorRange);

    % Show/hide colorbar (identical to left)
    if figData.showColorbar
        colorbar(rightAx);
    else
        colorbar(rightAx, 'off');
    end

    grid(rightAx, 'on');
    view(rightAx, figData.currentView);
    hold(rightAx, 'off');

    fprintf('Separate plot data rendering complete.\n');
end



function maintain3DUIPositions(fig)
    % Maintain UI positions when figure is resized (following MainPlottingApplication pattern)

    % Get current figure dimensions
    figPos = get(fig, 'Position');
    figWidth = figPos(3);
    figHeight = figPos(4);

    % Define UI panel constants - same as in creation (RIGHT-TOP ANCHORED)
    UI_PANEL_WIDTH = 180;           % Fixed width in pixels
    UI_MARGIN_RIGHT = 10;           % Distance from right edge to right edge of UI controls
    UI_MARGIN_TOP = 10;             % Distance from top edge to first UI element
    UI_ELEMENT_HEIGHT = 25;         % Standard height for dropdowns/buttons
    UI_TEXT_HEIGHT = 18;            % Height for text elements
    UI_CHECKBOX_HEIGHT = 20;        % Height for checkbox elements
    UI_SPACING = 5;                 % Spacing between elements

    % Calculate UI_LEFT so that UI_LEFT + UI_PANEL_WIDTH = figWidth - UI_MARGIN_RIGHT
    % This ensures the RIGHT EDGE of UI controls is exactly UI_MARGIN_RIGHT pixels from figure right edge
    UI_RIGHT_EDGE = figWidth - UI_MARGIN_RIGHT;  % Where the right edge of UI should be
    UI_LEFT = UI_RIGHT_EDGE - UI_PANEL_WIDTH;    % Left position to achieve this

    % Ensure UI_LEFT is always positive (minimum 10 pixels from left edge)
    if UI_LEFT < 10
        UI_LEFT = 10;
        UI_PANEL_WIDTH = UI_RIGHT_EDGE - UI_LEFT; % Adjust width if figure is too narrow
    end

    % Calculate starting Y position from TOP edge (MATLAB uses bottom-left origin)
    UI_TOP_Y = figHeight - UI_MARGIN_TOP - UI_ELEMENT_HEIGHT; % Start from top, accounting for element height

    % Start positioning from calculated top-anchored position and work down
    currentY = UI_TOP_Y; % Start from calculated top-anchored position

    % Find and reposition all UI controls in order
    plotModeDropdown = findobj(fig, 'Tag', 'PlotModeDropdown');
    if ~isempty(plotModeDropdown)
        set(plotModeDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % Peak marker label and dropdown
    peakMarkerLabel = findobj(fig, 'Tag', 'PeakMarkerLabel');
    if ~isempty(peakMarkerLabel)
        set(peakMarkerLabel, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    peakMarkerDropdown = findobj(fig, 'Tag', 'PeakMarkerDropdown');
    if ~isempty(peakMarkerDropdown)
        set(peakMarkerDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % Valley marker label and dropdown
    valleyMarkerLabel = findobj(fig, 'Tag', 'ValleyMarkerLabel');
    if ~isempty(valleyMarkerLabel)
        set(valleyMarkerLabel, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    valleyMarkerDropdown = findobj(fig, 'Tag', 'ValleyMarkerDropdown');
    if ~isempty(valleyMarkerDropdown)
        set(valleyMarkerDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % Marker size controls
    markerSizeLabel = findobj(fig, 'Tag', 'MarkerSizeLabel');
    if ~isempty(markerSizeLabel)
        set(markerSizeLabel, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    buttonWidth = (UI_PANEL_WIDTH - UI_SPACING) / 2;
    smallSizeButton = findobj(fig, 'Tag', 'SmallSizeButton');
    largeSizeButton = findobj(fig, 'Tag', 'LargeSizeButton');
    if ~isempty(smallSizeButton)
        set(smallSizeButton, 'Position', [UI_LEFT, currentY, buttonWidth, UI_ELEMENT_HEIGHT]);
    end
    if ~isempty(largeSizeButton)
        set(largeSizeButton, 'Position', [UI_LEFT + buttonWidth + UI_SPACING, currentY, buttonWidth, UI_ELEMENT_HEIGHT]);
    end
    if ~isempty(smallSizeButton) || ~isempty(largeSizeButton)
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % Colormap controls
    colormapLabel = findobj(fig, 'Tag', 'ColormapLabel');
    if ~isempty(colormapLabel)
        set(colormapLabel, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    colormapDropdown = findobj(fig, 'Tag', 'ColormapDropdown');
    if ~isempty(colormapDropdown)
        set(colormapDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % View controls
    viewLabel = findobj(fig, 'Tag', 'ViewLabel');
    if ~isempty(viewLabel)
        set(viewLabel, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
        currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
    end

    viewDropdown = findobj(fig, 'Tag', 'ViewDropdown');
    if ~isempty(viewDropdown)
        set(viewDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end

    % Colorbar toggle
    colorbarToggleButton = findobj(fig, 'Tag', 'ColorbarToggleButton');
    if ~isempty(colorbarToggleButton)
        set(colorbarToggleButton, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
        currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
    end


    % Update plot area based on mode
    figData = get(fig, 'UserData');
    if isfield(figData, 'plotMode') && strcmp(figData.plotMode, 'separate')
        % Adjust separate plot positions
        if isfield(figData, 'leftAx') && all(ishandle(figData.leftAx))
            set(figData.leftAx, 'Position', [0.05, 0.05, 0.33, 0.9]);
        end
        if isfield(figData, 'rightAx') && all(ishandle(figData.rightAx))
            set(figData.rightAx, 'Position', [0.39, 0.05, 0.33, 0.9]);
        end
    else
        % Adjust single plot position
        if isfield(figData, 'mainAx') && all(ishandle(figData.mainAx))
            set(figData.mainAx, 'Position', [0.05, 0.05, 0.72, 0.9]);
        end
    end
end


% Helper function to create a coolwarm colormap (from XtVsYPlot.m)



















function t = reconstructTimeVector(peakData)
% Attempt to reconstruct time vector from peak data
t = [];

try
    % Find all unique times from peak data
    allTimes = [];
    for i = 1:length(peakData)
        if ~isempty(peakData{i}) && istable(peakData{i})
            allTimes = [allTimes; peakData{i}.TransitionTime];
        end
    end

    if isempty(allTimes)
        return;
    end

    % Create time vector spanning the range with reasonable resolution
    minTime = min(allTimes);
    maxTime = max(allTimes);

    % Estimate time step from peak data density
    numSamples = 1000; % Default reasonable number of samples
    t = linspace(minTime, maxTime, numSamples);

    fprintf('Reconstructed time vector: %.2f to %.2f μs (%d samples)\n', ...
            minTime*1e6, maxTime*1e6, length(t));

catch ME
    fprintf('Failed to reconstruct time vector: %s\n', ME.message);
    t = [];
end
end



function repositionUIControlsOnResize(fig)
    % REPOSITIONUICONTROLSONRESIZE - Reposition UI controls when figure is resized
    % This ensures the right-anchored UI panel stays fixed to the right edge

    try
        % Validate figure handle before proceeding
        if ~ishandle(fig) || ~ishghandle(fig, 'figure')
            return; % Figure has been deleted, exit silently
        end

        % Get current figure dimensions
        figPos = get(fig, 'Position');
        if length(figPos) < 4
            fprintf('Warning: Invalid figure position data\n');
            return;
        end
        figWidth = figPos(3);
        figHeight = figPos(4);

        % Validate dimensions are reasonable
        if figWidth < 100 || figHeight < 100
            return; % Figure too small, skip repositioning
        end

        % Define UI panel constants - same as in creation (RIGHT-TOP ANCHORED)
        UI_PANEL_WIDTH = 180;           % Fixed width in pixels
        UI_MARGIN_RIGHT = 10;           % Distance from right edge to right edge of UI controls
        UI_MARGIN_TOP = 10;             % Distance from top edge to first UI element
        UI_ELEMENT_HEIGHT = 25;         % Standard height for dropdowns/buttons
        UI_TEXT_HEIGHT = 18;            % Height for text elements
        UI_CHECKBOX_HEIGHT = 20;        % Height for checkbox elements
        UI_SPACING = 5;                 % Spacing between elements

        % Calculate UI_LEFT so that UI_LEFT + UI_PANEL_WIDTH = figWidth - UI_MARGIN_RIGHT
        % This ensures the RIGHT EDGE of UI controls is exactly UI_MARGIN_RIGHT pixels from figure right edge
        UI_RIGHT_EDGE = figWidth - UI_MARGIN_RIGHT;  % Where the right edge of UI should be
        UI_LEFT = UI_RIGHT_EDGE - UI_PANEL_WIDTH;    % Left position to achieve this

        % Ensure UI_LEFT is always positive (minimum 10 pixels from left edge)
        if UI_LEFT < 10
            UI_LEFT = 10;
            UI_PANEL_WIDTH = UI_RIGHT_EDGE - UI_LEFT; % Adjust width if figure is too narrow
        end

        % Calculate starting Y position from TOP edge (MATLAB uses bottom-left origin)
        UI_TOP_Y = figHeight - UI_MARGIN_TOP - UI_ELEMENT_HEIGHT; % Start from top, accounting for element height

        % Start positioning from calculated top-anchored position and work down
        currentY = UI_TOP_Y; % Start from calculated top-anchored position

        % Find and reposition all UI controls in order
        plotModeDropdown = findobj(fig, 'Tag', 'PlotModeDropdown');
        if ~isempty(plotModeDropdown)
            set(plotModeDropdown, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
            currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
        end

        % CRITICAL FIX: Complete controls list following exact creation order
        % This matches the order in createRightAnchoredControlPanel function
        controlsInOrder = {
            'PeakMarkerLabel', 'PeakMarkerDropdown',
            'ValleyMarkerLabel', 'ValleyMarkerDropdown',
            'MarkerSizeLabel', % Note: Size buttons handled separately
            'ColormapLabel', 'ColormapDropdown',
            'ViewLabel', 'ViewDropdown',
            'ColorbarToggleButton'
            % Note: Tolerance inputs and buttons handled separately below
        };

        % Process controls in creation order
        for i = 1:length(controlsInOrder)
            control = findobj(fig, 'Tag', controlsInOrder{i});
            if ~isempty(control)
                % Ensure pixel units for precise right-anchored placement
                set(control, 'Units', 'pixels');
                if contains(controlsInOrder{i}, 'Label')
                    set(control, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_TEXT_HEIGHT]);
                    currentY = currentY - UI_TEXT_HEIGHT - UI_SPACING;
                else
                    set(control, 'Position', [UI_LEFT, currentY, UI_PANEL_WIDTH, UI_ELEMENT_HEIGHT]);
                    currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
                end
            end
        end

        % Handle special cases for side-by-side buttons and input fields
        buttonWidth = (UI_PANEL_WIDTH - UI_SPACING) / 2;
        smallSizeButton = findobj(fig, 'Tag', 'SmallSizeButton');
        largeSizeButton = findobj(fig, 'Tag', 'LargeSizeButton');
        if ~isempty(smallSizeButton)
            set(smallSizeButton, 'Units','pixels','Position', [UI_LEFT, currentY, buttonWidth, UI_ELEMENT_HEIGHT]);
        end
        if ~isempty(largeSizeButton)
            set(largeSizeButton, 'Units','pixels','Position', [UI_LEFT + buttonWidth + UI_SPACING, currentY, buttonWidth, UI_ELEMENT_HEIGHT]);
        end
        if ~isempty(smallSizeButton) || ~isempty(largeSizeButton)
            currentY = currentY - UI_ELEMENT_HEIGHT - UI_SPACING;
        end


        % Update plot area to maintain proper spacing from UI controls
        figData = get(fig, 'UserData');
        % Compute normalized right boundary based on UI panel width in pixels
        contentLeftNorm = 0.05; contentBottomNorm = 0.05; contentHeightNorm = 0.9;
        contentRightNorm = (figWidth - UI_PANEL_WIDTH - UI_MARGIN_RIGHT) / figWidth;
        contentWidthNorm = max(0.1, contentRightNorm - contentLeftNorm);
        gapNorm = 0.02; % gap between separate axes

        if isfield(figData, 'plotMode') && strcmp(figData.plotMode, 'separate')
            % Adjust separate plot positions responsively
            leftWidth = max(0.1, (contentWidthNorm - gapNorm) / 2);
            rightWidth = leftWidth;
            leftPos = [contentLeftNorm, contentBottomNorm, leftWidth, contentHeightNorm];
            rightPos = [contentLeftNorm + leftWidth + gapNorm, contentBottomNorm, rightWidth, contentHeightNorm];
            if isfield(figData, 'leftAx') && ~isempty(figData.leftAx) && ishandle(figData.leftAx)
                set(figData.leftAx, 'Units', 'normalized', 'Position', leftPos);
            end
            if isfield(figData, 'rightAx') && ~isempty(figData.rightAx) && ishandle(figData.rightAx)
                set(figData.rightAx, 'Units', 'normalized', 'Position', rightPos);
            end
            % Keep UI panel at fixed pixel width on the right by ensuring axes' right edge
            % respects UI_PANEL_WIDTH. This is already handled by contentRightNorm above.
        else
            % Adjust single plot position responsively
            if isfield(figData, 'mainAx') && ~isempty(figData.mainAx) && ishandle(figData.mainAx)
                set(figData.mainAx, 'Units', 'normalized', 'Position', [contentLeftNorm, contentBottomNorm, contentWidthNorm, contentHeightNorm]);
            end
        end

        % Force UI refresh
        drawnow;

    catch ME
        % Silently handle resize errors to avoid disrupting user experience
        fprintf('Warning: UI resize failed: %s\n', ME.message);
    end
end

function fixUIPositioning(fig)
    % FIXUIPOSITIONING - Manual function to fix UI positioning issues
    % Call this function if UI controls are not properly positioned

    if nargin < 1
        fig = gcf; % Use current figure if none specified
    end

    fprintf('Manually fixing UI positioning...\n');

    % Force repositioning
    repositionUIControlsOnResize(fig);

    % Force refresh
    drawnow;

    fprintf('UI positioning fixed. Controls should now be properly anchored to right edge.\n');
end

function verifyUIPositioning(fig)
    % VERIFYUIPOSITIONING - Verify that UI controls are properly positioned
    % Call this function to check if UI positioning is working correctly

    if nargin < 1
        fig = gcf; % Use current figure if none specified
    end

    fprintf('\n=== UI Positioning Verification ===\n');

    % Get figure dimensions
    figPos = get(fig, 'Position');
    figWidth = figPos(3);

    % Expected positioning
    UI_MARGIN_RIGHT = 10;
    expectedRightEdge = figWidth - UI_MARGIN_RIGHT;

    fprintf('Figure width: %.0f pixels\n', figWidth);
    fprintf('Expected UI right edge: %.0f pixels from left (%.0f pixels from right)\n', ...
            expectedRightEdge, UI_MARGIN_RIGHT);

    % Check a few key UI controls
    controls = {'PlotModeDropdown'};

    for i = 1:length(controls)
        control = findobj(fig, 'Tag', controls{i});
        if ~isempty(control)
            pos = get(control, 'Position');
            actualRightEdge = pos(1) + pos(3);
            distanceFromRight = figWidth - actualRightEdge;

            fprintf('Control "%s": Right edge at %.0f pixels (%.0f pixels from right edge)\n', ...
                    controls{i}, actualRightEdge, distanceFromRight);

            if abs(distanceFromRight - UI_MARGIN_RIGHT) < 2 % Allow 2 pixel tolerance
                fprintf('  ✓ Correctly positioned\n');
            else
                fprintf('  ❌ Incorrectly positioned (should be %.0f pixels from right)\n', UI_MARGIN_RIGHT);
            end
        else
            fprintf('Control "%s": Not found\n', controls{i});
        end
    end

    fprintf('=================================\n\n');
end

% REMOVED: createPointVisualization - consolidated with createScatterPointsVisualization
% All references updated to use createScatterPointsVisualization for consistency

function surfaceHandle = createTriangulatedSurface(ax, X, Y, Z, color, alpha, edgeStyle)
% CREATETRIANGULATEDSURFACE - Create surface using Delaunay triangulation
% FIXED: Handles duplicate points and provides robust triangulation

surfaceHandle = [];

try
    if length(X) < 3
        return;
    end

    % CRITICAL FIX: Remove duplicate points before triangulation
    [uniquePoints, uniqueIdx] = removeDuplicatePoints(X(:), Y(:), Z(:));

    if size(uniquePoints, 1) < 3
        fprintf('  Warning: Less than 3 unique points after duplicate removal. Using point visualization.\n');
        surfaceHandle = createPointVisualization(ax, X, Y, Z, color, alpha);
        return;
    end

    X_unique = uniquePoints(:, 1);
    Y_unique = uniquePoints(:, 2);
    Z_unique = uniquePoints(:, 3);

    fprintf('  Triangulation: %d points → %d unique points\n', length(X), length(X_unique));

    % Create 2D Delaunay triangulation in XY plane with unique points
    try
        DT = delaunayTriangulation(X_unique, Y_unique);
    catch ME
        fprintf('  Warning: Delaunay triangulation failed (%s). Using convex hull.\n', ME.message);
        surfaceHandle = createConvexHullSurface(ax, X, Y, Z, color, alpha, edgeStyle);
        return;
    end

    if isempty(DT.ConnectivityList)
        fprintf('  Warning: No triangulation connectivity. Using scatter points.\n');
        surfaceHandle = createScatterPointsVisualization(ax, X, Y, Z, color, alpha);
        return;
    end

    % Create triangulated surface using unique points
    if strcmp(edgeStyle, 'smooth')
        surfaceHandle = trisurf(DT.ConnectivityList, X_unique, Y_unique, Z_unique, ...
                               'FaceColor', color, 'FaceAlpha', alpha, ...
                               'EdgeColor', 'none', 'Parent', ax);
    else
        surfaceHandle = trisurf(DT.ConnectivityList, X_unique, Y_unique, Z_unique, ...
                               'FaceColor', color, 'FaceAlpha', alpha, ...
                               'EdgeColor', 'k', 'LineWidth', 0.5, 'Parent', ax);
    end

    fprintf('  ✓ Triangulated surface created with %d triangles\n', size(DT.ConnectivityList, 1));

catch ME
    fprintf('Triangulation failed: %s\n', ME.message);
end
end

function surfaceHandle = createGridBasedSurface(ax, X, Y, Z, color, alpha, edgeStyle)
% CREATEGRIDBASEDSURFACE - Create surface by interpolating onto regular grid

surfaceHandle = [];

try
    if length(X) < 4
        return;
    end

    % Create regular grid
    gridSize = min(20, ceil(sqrt(length(X))));
    xi = linspace(min(X), max(X), gridSize);
    yi = linspace(min(Y), max(Y), gridSize);
    [Xi, Yi] = meshgrid(xi, yi);

    % Interpolate Z values onto grid
    F = scatteredInterpolant(X(:), Y(:), Z(:), 'linear', 'none');
    Zi = F(Xi, Yi);

    % Remove NaN values
    validMask = ~isnan(Zi);
    if sum(validMask(:)) < 3
        return;
    end

    % Create surface
    if strcmp(edgeStyle, 'smooth')
        surfaceHandle = surf(ax, Xi, Yi, Zi, 'FaceColor', color, 'FaceAlpha', alpha, ...
                            'EdgeColor', 'none');
    else
        surfaceHandle = surf(ax, Xi, Yi, Zi, 'FaceColor', color, 'FaceAlpha', alpha, ...
                            'EdgeColor', 'k', 'LineWidth', 0.5);
    end

catch ME
    fprintf('Grid-based surface failed: %s\n', ME.message);
end
end

function [uniquePoints, uniqueIdx] = removeDuplicatePoints(X, Y, Z)
% REMOVEDUPLICATEPOINTS - Remove duplicate XY coordinates, keeping average Z
%
% For points with identical X,Y coordinates, this function:
% 1. Groups them together
% 2. Computes the average Z value
% 3. Returns unique XY points with averaged Z values
%
% This prevents Delaunay triangulation warnings about duplicate points

% Input validation
if length(X) ~= length(Y) || length(Y) ~= length(Z)
    error('Input arrays X, Y, Z must have the same length. Got lengths: X=%d, Y=%d, Z=%d', ...
          length(X), length(Y), length(Z));
end

if isempty(X)
    uniquePoints = [];
    uniqueIdx = [];
    return;
end

% Convert to column vectors for consistency
X = X(:);
Y = Y(:);
Z = Z(:);

% Check for non-finite values
finiteIdx = isfinite(X) & isfinite(Y) & isfinite(Z);
if ~all(finiteIdx)
    fprintf('Warning: Removing %d non-finite points before duplicate removal\n', sum(~finiteIdx));
    X = X(finiteIdx);
    Y = Y(finiteIdx);
    Z = Z(finiteIdx);
end

if isempty(X)
    uniquePoints = [];
    uniqueIdx = [];
    return;
end

try
    % Find unique XY coordinates with tolerance for floating point precision
    tolerance = 1e-10;
    [~, uniqueIdx, groupIdx] = uniquetol([X, Y], tolerance, 'ByRows', true);

    % For each unique XY location, average the Z values
    uniquePoints = zeros(length(uniqueIdx), 3);
    for i = 1:length(uniqueIdx)
        groupMembers = (groupIdx == i);
        uniquePoints(i, 1) = X(uniqueIdx(i)); % X coordinate
        uniquePoints(i, 2) = Y(uniqueIdx(i)); % Y coordinate
        uniquePoints(i, 3) = mean(Z(groupMembers)); % Average Z value
    end

    % Report duplicate removal statistics
    numOriginal = length(X);
    numUnique = size(uniquePoints, 1);
    if numOriginal > numUnique
        fprintf('  Removed %d duplicate points (%d → %d unique)\n', ...
                numOriginal - numUnique, numOriginal, numUnique);
    end

catch ME
    fprintf('Error in duplicate point removal: %s\n', ME.message);
    % Fallback: return original points
    uniquePoints = [X, Y, Z];
    uniqueIdx = (1:length(X))';
end
end











function computeCrossViewAlignmentDirect(fig)
% COMPUTECROSSVIEWALIGNMENTDIRECT - Direct computation without separate initialization
%
% This function directly calls the cross-sectional alignment computation
% without going through the separate initialization step.

try
    figData = get(fig, 'UserData');

    % Check if alignment data exists
    if ~isfield(figData, 'alignmentData') || isempty(figData.alignmentData)
        fprintf('Error: Alignment data not initialized. Cannot compute alignment.\n');
        return;
    end

    % Convert our data structure to match XtVsYPlot.m format
    fprintf('Converting data structure for XtVsYPlot.m compatibility...\n');

    % CRITICAL FIX: Access the correct field name from crossSectionalAlignment.m
    % The data is stored as 'originalStatMaps', not 'statDataArray'
    if isfield(figData.alignmentData, 'originalStatMaps')
        originalStatMaps = figData.alignmentData.originalStatMaps;
        fprintf('Successfully accessed originalStatMaps from alignmentData\n');
    else
        fprintf('Error: originalStatMaps not found in alignmentData\n');
        fprintf('Available fields in alignmentData: %s\n', strjoin(fieldnames(figData.alignmentData), ', '));
        return;
    end

    % CRITICAL FIX: Ensure statDataArray is a cell array (XtVsYPlot.m expects this)
    % The originalStatMaps is a single struct, but XtVsYPlot expects cell array of structs
    if ~iscell(originalStatMaps)
        statDataArray = {originalStatMaps}; % Wrap in cell array
        fprintf('Wrapped originalStatMaps in cell array for XtVsYPlot.m compatibility\n');
    else
        statDataArray = originalStatMaps;
        fprintf('originalStatMaps is already a cell array\n');
    end

    % Validate the data structure thoroughly
    if isempty(statDataArray)
        fprintf('Error: statDataArray is empty after conversion\n');
        return;
    end

    if length(statDataArray) < 1 || isempty(statDataArray{1})
        fprintf('Error: statDataArray{1} is empty\n');
        return;
    end

    if ~isstruct(statDataArray{1})
        fprintf('Error: statDataArray{1} is not a struct\n');
        return;
    end

    if ~isfield(statDataArray{1}, 'maps') || isempty(statDataArray{1}.maps)
        fprintf('Error: statDataArray{1} does not contain maps field\n');
        fprintf('Available fields in statDataArray{1}: %s\n', strjoin(fieldnames(statDataArray{1}), ', '));
        return;
    end

    if ~iscell(statDataArray{1}.maps) || isempty(statDataArray{1}.maps{1})
        fprintf('Error: statDataArray{1}.maps is not a valid cell array\n');
        return;
    end

    fprintf('✅ Data structure validation passed:\n');
    fprintf('   - Number of statistics: %d\n', length(statDataArray));
    fprintf('   - Number of time segments: %d\n', length(statDataArray{1}.maps));
    fprintf('   - Map dimensions: %dx%d\n', size(statDataArray{1}.maps{1}));

    % Create userData structure compatible with XtVsYPlot.m
    userData = struct();
    userData.statDataArray = statDataArray;
    userData.originalStatDataArray = statDataArray; % Keep original
    userData.alignedStatDataArray = []; % Will be populated by alignment
    userData.isComputingAlignment = false;
    userData.progressDialog = [];

    % Store the converted userData in the figure
    set(fig, 'UserData', userData);

    % Store original figData in userData for later restoration
    userData.originalFigData = figData;
    set(fig, 'UserData', userData);

    % Use the imported function from XtVsYPlot.m (defined at end of this file)
    fprintf('Calling computeCrossViewAlignment3D function...\n');
    computeCrossViewAlignment3D([], fig);
    fprintf('computeCrossViewAlignment3D function completed.\n');

    % NOTE: Do NOT restore figData here! The timer callback needs userData.statDataArray
    % The restoration will happen in the timer callback completion in computeCrossViewAlignmentBackground3D.m


catch ME
    fprintf('Error in direct cross-view alignment computation: %s\n', ME.message);
end
end

function isValid = validateSystemIntegrity()
% VALIDATESYSTEMINTEGRITY - Comprehensive system validation
%
% This function checks for all the critical dependencies and validates
% the system is ready for operation.
%
% Output:
%   isValid - Boolean indicating if system passed all validation checks

isValid = true;
fprintf('\n=== System Integrity Validation ===\n');

% Check for required external functions
requiredFunctions = {
    'PlateGenerationProcessor', 'Plate generation algorithm #1';
    'AdaptiveCubeProcessor', 'Plate generation algorithm #2';
    'GroupAveragingProcessor', 'Plate generation algorithm #3';
    'TopologicalPlateGenerator', 'Plate generation algorithm #4'
};

fprintf('Checking external function dependencies:\n');
for i = 1:size(requiredFunctions, 1)
    funcName = requiredFunctions{i, 1};
    description = requiredFunctions{i, 2};

    if exist(funcName, 'file')
        fprintf('  ✓ %s - %s\n', funcName, description);
    else
        fprintf('  ⚠ %s - %s (fallback will be used)\n', funcName, description);
    end
end

% Check MATLAB version compatibility
matlabVersion = version('-release');
fprintf('MATLAB version: %s\n', matlabVersion);

% Validate critical functions exist in this file
internalFunctions = {
    'generatePlateCacheHash', 'Cache management';
    'clearPlateCache', 'Cache clearing';
    'createFallbackPlateData', 'Fallback plate generation';
    'removeDuplicatePoints', 'Triangulation preprocessing';
    'validateFigDataIntegrity', 'Data validation'
};

fprintf('Checking internal function integrity:\n');
for i = 1:size(internalFunctions, 1)
    funcName = internalFunctions{i, 1};
    purpose = internalFunctions{i, 2};

    if exist(funcName, 'file') == 2 % Function exists in this file
        fprintf('  ✓ %s - %s\n', funcName, purpose);
    else
        fprintf('  ✗ %s - %s (MISSING - CRITICAL ERROR)\n', funcName, purpose);
        isValid = false;
    end
end

if isValid
    fprintf('✓ System integrity validation PASSED\n');
else
    fprintf('✗ System integrity validation FAILED - some critical functions missing\n');
end

fprintf('=====================================\n\n');
end

% Removed: generatePlatesWithLayerDetection - algorithm removed

% Removed: generatePlatesWithFiberDirection - algorithm removed

function cacheHash = generatePlateCacheHash(plateCacheKey)
% GENERATEPLATECACHEHASH - Generate unique hash for plate cache key
%
% This function creates a unique hash string from the plate cache key structure
% to identify cached plate data files.
%
% Input:
%   plateCacheKey - Structure containing all parameters that affect plate generation
%
% Output:
%   cacheHash - Unique hash string for cache identification

try
    % Convert structure to string representation
    keyString = struct2string(plateCacheKey);

    % Create hash using built-in MATLAB hash function (if available)
    if exist('hash', 'builtin') || exist('hash', 'file')
        cacheHash = hash(keyString, 'MD5');
    else
        % Fallback: create simple hash from string
        cacheHash = createSimpleHash(keyString);
    end

    % Ensure hash is valid filename (remove special characters)
    cacheHash = regexprep(cacheHash, '[^a-zA-Z0-9]', '');

    % Limit length to reasonable filename size
    if length(cacheHash) > 32
        cacheHash = cacheHash(1:32);
    end

catch ME
    fprintf('Warning: Failed to generate cache hash (%s). Using timestamp.\n', ME.message);
    % Fallback: use timestamp-based hash
    cacheHash = sprintf('Cache_%s', datestr(now, 'yyyymmdd_HHMMSS'));
end
end

function hashStr = createSimpleHash(inputStr)
% CREATESIMPLEHASH - Create simple hash from input string
% Fallback hash function when built-in hash is not available

hashValue = 0;
for i = 1:length(inputStr)
    hashValue = mod(hashValue * 31 + double(inputStr(i)), 2^32);
end
hashStr = sprintf('%08X', hashValue);
end

function str = struct2string(s)
% STRUCT2STRING - Convert structure to string representation
% Recursively converts structure to deterministic string

if isstruct(s)
    fields = sort(fieldnames(s)); % Sort for consistency
    str = '{';
    for i = 1:length(fields)
        field = fields{i};
        value = s.(field);
        str = [str, field, ':', struct2string(value), ';'];
    end
    str = [str, '}'];
elseif iscell(s)
    str = '[';
    for i = 1:length(s)
        str = [str, struct2string(s{i}), ','];
    end
    str = [str, ']'];
elseif isnumeric(s)
    if length(s) == 1
        str = sprintf('%.6f', s);
    else
        str = sprintf('%.6f,', s);
        str = ['[', str(1:end-1), ']']; % Remove trailing comma
    end
elseif ischar(s) || isstring(s)
    str = char(s);
elseif islogical(s)
    str = sprintf('%d', s);
else
    str = class(s);
end
end



% ========================================================================
% CROSS-VIEW ALIGNMENT FUNCTIONS (Now in separate files)
% ========================================================================
% The alignment functions have been moved to separate files:
% - computeCrossViewAlignment3D.m
% - computeCrossViewAlignmentBackground3D.m
% - alignAllSlicesInView3D.m
% - calculateOverallAlignmentCost3D.m


