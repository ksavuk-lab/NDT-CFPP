function testRobustAlignment()
    % TESTROBUSTALIGNMENT - Test function for the robust peak alignment system
    % Creates synthetic peak data to validate the robust alignment analyzer
    
    fprintf('=== TESTING ROBUST PEAK ALIGNMENT SYSTEM ===\n');
    
    % Create synthetic test data
    fprintf('Creating synthetic test data...\n');
    
    % Time vector
    t = linspace(0, 2e-6, 200); % 0 to 2 microseconds, 200 points
    
    % Create test peak data with spurious early peaks in some waveforms
    numWaveforms = 50;
    peakData = cell(numWaveforms, 1);
    
    for i = 1:numWaveforms
        transitions = struct();
        
        if i <= 10
            % First 10 waveforms: Have spurious early peak + main peaks
            transitionTimes = [0.3e-6, 0.8e-6, 1.2e-6, 1.6e-6]; % Early spurious + 3 main peaks
            transitionAmps = [0.2, 0.8, 0.6, 0.4]; % Spurious is smaller
            transitionTypes = [1, 1, 1, 1]; % All peaks
            
        elseif i <= 40
            % Most waveforms: Normal pattern starting at 0.8μs (main peak)
            transitionTimes = [0.8e-6, 1.2e-6, 1.6e-6]; % 3 main peaks
            transitionAmps = [0.8, 0.6, 0.4];
            transitionTypes = [1, 1, 1]; % All peaks
            
        else
            % Last 10 waveforms: Different pattern (should be minority)
            transitionTimes = [0.9e-6, 1.3e-6, 1.7e-6]; % Slightly shifted
            transitionAmps = [0.7, 0.5, 0.3];
            transitionTypes = [1, 1, 1]; % All peaks
        end
        
        % Add some noise to make it realistic
        transitionTimes = transitionTimes + (rand(size(transitionTimes)) - 0.5) * 0.05e-6;
        transitionAmps = transitionAmps + (rand(size(transitionAmps)) - 0.5) * 0.1;
        
        % Create transitions structure
        transitions.TransitionTime = transitionTimes';
        transitions.TransitionAmplitude = transitionAmps';
        transitions.TransitionType = transitionTypes';
        
        peakData{i} = transitions;
    end
    
    fprintf('Created %d synthetic waveforms:\n', numWaveforms);
    fprintf('  - 10 waveforms with spurious early peaks at 0.3μs\n');
    fprintf('  - 30 waveforms with main peak at 0.8μs (majority)\n');
    fprintf('  - 10 waveforms with shifted peak at 0.9μs (minority)\n');
    
    % Test the robust alignment analyzer
    fprintf('\nTesting RobustPeakAlignmentAnalyzer...\n');
    
    try
        [alignmentReference, analysisResults] = RobustPeakAlignmentAnalyzer.determineOptimalAlignment(peakData, t, true);
        
        fprintf('\n=== TEST RESULTS ===\n');
        fprintf('✅ Robust alignment analysis completed successfully!\n');
        fprintf('Selected peak position: %d\n', alignmentReference.peakPosition);
        fprintf('Reference time: %.3f μs\n', alignmentReference.referenceTime * 1e6);
        fprintf('Consensus strength: %.1f%%\n', alignmentReference.consensusStrength * 100);
        fprintf('Quality score: %.3f/1.0\n', analysisResults.qualityMetrics.qualityScore);
        
        % Validate results
        expectedPeakPosition = 2; % Should choose 2nd peak (0.8μs) over 1st peak (0.3μs spurious)
        expectedReferenceTime = 0.8e-6; % Should be around 0.8 microseconds
        
        if alignmentReference.peakPosition == expectedPeakPosition
            fprintf('✅ PASS: Correctly identified peak position %d (avoiding spurious early peaks)\n', expectedPeakPosition);
        else
            fprintf('❌ FAIL: Expected peak position %d, got %d\n', expectedPeakPosition, alignmentReference.peakPosition);
        end
        
        if abs(alignmentReference.referenceTime - expectedReferenceTime) < 0.1e-6
            fprintf('✅ PASS: Reference time %.3f μs is close to expected %.3f μs\n', ...
                alignmentReference.referenceTime * 1e6, expectedReferenceTime * 1e6);
        else
            fprintf('❌ FAIL: Reference time %.3f μs is far from expected %.3f μs\n', ...
                alignmentReference.referenceTime * 1e6, expectedReferenceTime * 1e6);
        end
        
        if alignmentReference.consensusStrength >= 0.6
            fprintf('✅ PASS: Strong consensus achieved (%.1f%% >= 60%%)\n', alignmentReference.consensusStrength * 100);
        else
            fprintf('⚠️  WARN: Weak consensus (%.1f%% < 60%%)\n', alignmentReference.consensusStrength * 100);
        end
        
        fprintf('\n=== COMPARISON WITH NAIVE ALIGNMENT ===\n');
        
        % Compare with naive first-peak alignment
        naiveReferenceTimes = [];
        for i = 1:numWaveforms
            if ~isempty(peakData{i})
                transitions = peakData{i};
                peakIndices = find(transitions.TransitionType == 1);
                if ~isempty(peakIndices)
                    naiveReferenceTimes = [naiveReferenceTimes; transitions.TransitionTime(peakIndices(1))];
                end
            end
        end
        
        naiveReferenceTime = min(naiveReferenceTimes);
        
        fprintf('Naive first-peak alignment would use: %.3f μs\n', naiveReferenceTime * 1e6);
        fprintf('Robust alignment uses: %.3f μs\n', alignmentReference.referenceTime * 1e6);
        fprintf('Improvement: %.3f μs shift to avoid spurious peaks\n', ...
            (alignmentReference.referenceTime - naiveReferenceTime) * 1e6);
        
        if alignmentReference.referenceTime > naiveReferenceTime
            fprintf('✅ SUCCESS: Robust alignment correctly avoided spurious early peaks!\n');
        else
            fprintf('❌ ISSUE: Robust alignment did not improve over naive method\n');
        end
        
    catch ME
        fprintf('❌ ERROR: Robust alignment analysis failed: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
    
    fprintf('\n=== TEST COMPLETE ===\n');
end
