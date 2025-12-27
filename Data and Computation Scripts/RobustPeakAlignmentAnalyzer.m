classdef RobustPeakAlignmentAnalyzer < handle
    % ROBUSTPEAKALIGNMENTANALYZER - Statistical analysis for robust peak alignment reference
    % Analyzes first 3 peaks in each waveform to identify the true alignment reference
    % that represents the same physical feature across the majority of waveforms
    
    properties (Constant)
        MAX_ANALYSIS_PEAKS = 3;     % Analyze first 3 peaks for pattern recognition
        MIN_PEAK_COUNT = 2;         % Minimum peaks required for analysis
        SIMILARITY_THRESHOLD = 0.15; % Time similarity threshold (in microseconds)
        MAJORITY_THRESHOLD = 0.6;   % 60% of waveforms must agree for consensus
    end
    
    methods (Static)
        function [alignmentReference, analysisResults] = determineOptimalAlignment(peakData, t, verbose)
            % DETERMINEOPTIMALALIGNMENT - Main function to determine robust alignment reference
            %
            % Inputs:
            %   peakData - Cell array of peak data for each waveform
            %   t - Time vector
            %   verbose - Display detailed analysis (optional, default: false)
            %
            % Outputs:
            %   alignmentReference - Structure with optimal alignment parameters
            %   analysisResults - Detailed analysis results for debugging

            if nargin < 3
                verbose = false;
            end
            
            fprintf('\n=== ROBUST PEAK ALIGNMENT ANALYSIS ===\n');
            
            % Step 1: Extract first 3 peaks from each waveform
            [peakPatterns, validWaveforms] = RobustPeakAlignmentAnalyzer.extractPeakPatterns(peakData, t, verbose);
            
            if isempty(peakPatterns)
                error('No valid peak patterns found for alignment analysis');
            end
            
            % Step 2: Analyze statistical properties of peak positions
            [peakStats, groupings] = RobustPeakAlignmentAnalyzer.analyzePeakStatistics(peakPatterns, verbose);
            
            % Step 3: Identify consensus alignment reference
            [alignmentReference, consensus] = RobustPeakAlignmentAnalyzer.findConsensusReference(peakStats, groupings, verbose);
            
            % Step 4: Validate alignment quality
            [qualityMetrics] = RobustPeakAlignmentAnalyzer.validateAlignmentQuality(peakPatterns, alignmentReference, verbose);
            
            % Compile results
            analysisResults = struct();
            analysisResults.peakPatterns = peakPatterns;
            analysisResults.validWaveforms = validWaveforms;
            analysisResults.peakStats = peakStats;
            analysisResults.groupings = groupings;
            analysisResults.consensus = consensus;
            analysisResults.qualityMetrics = qualityMetrics;
            analysisResults.totalWaveforms = length(peakData);
            analysisResults.validWaveforms = length(validWaveforms);
            
            if verbose
                RobustPeakAlignmentAnalyzer.displayAnalysisResults(alignmentReference, analysisResults);
            end
            
            fprintf('=== ROBUST ALIGNMENT ANALYSIS COMPLETE ===\n\n');
        end
        
        function [peakPatterns, validWaveforms] = extractPeakPatterns(peakData, t, verbose)
            % EXTRACTPEAKPATTERNS - Extract first 3 peaks from each waveform
            % Updated to work with transitions data structure (TransitionType, TransitionTime)

            peakPatterns = [];
            validWaveforms = [];

            if verbose
                fprintf('Extracting peak patterns from %d waveforms...\n', length(peakData));
            end

            for i = 1:length(peakData)
                if isempty(peakData{i})
                    continue;
                end

                transitions = peakData{i};

                % DEBUG: Check what fields are actually available
                if verbose && i <= 3  % Only show for first 3 waveforms
                    fprintf('DEBUG: Waveform %d fields: %s\n', i, strjoin(fieldnames(transitions), ', '));
                end

                % Check for required fields in transitions structure
                if ~isfield(transitions, 'TransitionType') || ~isfield(transitions, 'TransitionTime') || ~isfield(transitions, 'TransitionAmplitude')
                    if verbose && i <= 3
                        fprintf('DEBUG: Waveform %d missing required fields\n', i);
                    end
                    continue;
                end

                % Extract peaks (positive transitions, TransitionType == 1)
                peakIndices = find(transitions.TransitionType == 1);

                if verbose && i <= 3
                    fprintf('DEBUG: Waveform %d has %d peaks\n', i, length(peakIndices));
                end

                % Need at least MIN_PEAK_COUNT peaks
                if length(peakIndices) < RobustPeakAlignmentAnalyzer.MIN_PEAK_COUNT
                    continue;
                end

                % Extract first 3 peaks (or fewer if not available)
                numPeaksToAnalyze = min(length(peakIndices), RobustPeakAlignmentAnalyzer.MAX_ANALYSIS_PEAKS);

                pattern = struct();
                pattern.waveformIndex = i;
                pattern.numPeaks = numPeaksToAnalyze;

                for p = 1:numPeaksToAnalyze
                    transitionIdx = peakIndices(p);
                    pattern.peakTimes(p) = transitions.TransitionTime(transitionIdx);
                    pattern.peakAmplitudes(p) = transitions.TransitionAmplitude(transitionIdx);
                    pattern.peakIndices(p) = transitionIdx; % Index in transitions array
                end

                % Calculate inter-peak intervals
                if numPeaksToAnalyze > 1
                    for p = 1:numPeaksToAnalyze-1
                        pattern.intervals(p) = pattern.peakTimes(p+1) - pattern.peakTimes(p);
                    end
                else
                    pattern.intervals = [];
                end

                peakPatterns = [peakPatterns; pattern];
                validWaveforms = [validWaveforms; i];
            end

            if verbose
                fprintf('Extracted patterns from %d/%d waveforms (%.1f%%)\n', ...
                    length(peakPatterns), length(peakData), ...
                    100*length(peakPatterns)/length(peakData));
            end
        end
        
        function [peakStats, groupings] = analyzePeakStatistics(peakPatterns, verbose)
            % ANALYZEPEAKSTATISTICS - Analyze statistical properties of peak positions
            
            if verbose
                fprintf('Analyzing statistical properties of peak patterns...\n');
            end
            
            numPatterns = length(peakPatterns);
            maxPeaks = RobustPeakAlignmentAnalyzer.MAX_ANALYSIS_PEAKS;
            
            % Initialize statistics for each peak position
            peakStats = struct();
            for p = 1:maxPeaks
                peakStats(p).position = p;
                peakStats(p).times = [];
                peakStats(p).amplitudes = [];
                peakStats(p).count = 0;
            end
            
            % Collect data for each peak position
            for i = 1:numPatterns
                pattern = peakPatterns(i);
                for p = 1:min(pattern.numPeaks, maxPeaks)
                    peakStats(p).times = [peakStats(p).times; pattern.peakTimes(p)];
                    peakStats(p).amplitudes = [peakStats(p).amplitudes; pattern.peakAmplitudes(p)];
                    peakStats(p).count = peakStats(p).count + 1;
                end
            end
            
            % Calculate statistics for each peak position
            for p = 1:maxPeaks
                if peakStats(p).count > 0
                    peakStats(p).meanTime = mean(peakStats(p).times);
                    peakStats(p).stdTime = std(peakStats(p).times);
                    peakStats(p).medianTime = median(peakStats(p).times);
                    peakStats(p).meanAmplitude = mean(peakStats(p).amplitudes);
                    peakStats(p).stdAmplitude = std(peakStats(p).amplitudes);
                    peakStats(p).coverage = peakStats(p).count / numPatterns;
                    
                    % Calculate consistency metric (lower std = more consistent)
                    peakStats(p).consistency = 1 / (1 + peakStats(p).stdTime);
                else
                    peakStats(p).meanTime = NaN;
                    peakStats(p).stdTime = NaN;
                    peakStats(p).medianTime = NaN;
                    peakStats(p).meanAmplitude = NaN;
                    peakStats(p).stdAmplitude = NaN;
                    peakStats(p).coverage = 0;
                    peakStats(p).consistency = 0;
                end
            end
            
            % Group waveforms by peak pattern similarity
            groupings = RobustPeakAlignmentAnalyzer.groupBySimilarity(peakPatterns, peakStats, verbose);
            
            if verbose
                fprintf('Peak position statistics:\n');
                for p = 1:maxPeaks
                    if peakStats(p).count > 0
                        fprintf('  Peak %d: Mean=%.3f±%.3f μs, Coverage=%.1f%%, Consistency=%.3f\n', ...
                            p, peakStats(p).meanTime*1e6, peakStats(p).stdTime*1e6, ...
                            peakStats(p).coverage*100, peakStats(p).consistency);
                    end
                end
            end
        end
        
        function groupings = groupBySimilarity(peakPatterns, peakStats, verbose)
            % GROUPBYSIMILARITY - Group waveforms by peak pattern similarity
            
            groupings = struct();
            groupings.groups = {};
            groupings.groupSizes = [];
            groupings.groupRepresentatives = [];
            
            % For now, implement a simple grouping based on first peak timing
            % More sophisticated clustering can be added later
            
            if isempty(peakPatterns) || peakStats(1).count == 0
                return;
            end
            
            firstPeakTimes = [peakPatterns.peakTimes];
            firstPeakTimes = firstPeakTimes(1:3:end); % Extract first peak from each pattern
            
            % Simple clustering based on time similarity
            threshold = RobustPeakAlignmentAnalyzer.SIMILARITY_THRESHOLD * 1e-6; % Convert to seconds
            
            groups = {};
            groupSizes = [];
            
            for i = 1:length(firstPeakTimes)
                assigned = false;
                
                % Check if this pattern fits in an existing group
                for g = 1:length(groups)
                    groupMean = mean([peakPatterns(groups{g}).peakTimes]);
                    groupMean = groupMean(1:3:end); % First peaks only
                    
                    if abs(firstPeakTimes(i) - mean(groupMean)) <= threshold
                        groups{g} = [groups{g}, i];
                        groupSizes(g) = groupSizes(g) + 1;
                        assigned = true;
                        break;
                    end
                end
                
                % Create new group if not assigned
                if ~assigned
                    groups{end+1} = i;
                    groupSizes(end+1) = 1;
                end
            end
            
            groupings.groups = groups;
            groupings.groupSizes = groupSizes;
            
            % Find representative for each group (pattern closest to group mean)
            for g = 1:length(groups)
                groupIndices = groups{g};
                groupTimes = firstPeakTimes(groupIndices);
                groupMean = mean(groupTimes);
                
                [~, closestIdx] = min(abs(groupTimes - groupMean));
                groupings.groupRepresentatives(g) = groupIndices(closestIdx);
            end
            
            if verbose
                fprintf('Found %d groups based on peak similarity:\n', length(groups));
                for g = 1:length(groups)
                    fprintf('  Group %d: %d waveforms (%.1f%%)\n', ...
                        g, groupSizes(g), 100*groupSizes(g)/length(peakPatterns));
                end
            end
        end

        function [alignmentReference, consensus] = findConsensusReference(peakStats, groupings, verbose)
            % FINDCONSENSUSREFERENCE - Identify the optimal alignment reference

            if verbose
                fprintf('Determining consensus alignment reference...\n');
            end

            % Find the largest group (majority consensus)
            [maxGroupSize, maxGroupIdx] = max(groupings.groupSizes);
            majorityFraction = maxGroupSize / sum(groupings.groupSizes);

            consensus = struct();
            consensus.majorityGroupIndex = maxGroupIdx;
            consensus.majorityGroupSize = maxGroupSize;
            consensus.majorityFraction = majorityFraction;
            consensus.meetsThreshold = majorityFraction >= RobustPeakAlignmentAnalyzer.MAJORITY_THRESHOLD;

            if ~consensus.meetsThreshold
                warning('Majority consensus not reached (%.1f%% < %.1f%%). Using largest group.', ...
                    majorityFraction*100, RobustPeakAlignmentAnalyzer.MAJORITY_THRESHOLD*100);
            end

            % Determine which peak position to use for alignment
            % Use the peak position with highest consistency in the majority group
            bestPeakPosition = 1; % Default to first peak
            bestConsistency = 0;

            for p = 1:RobustPeakAlignmentAnalyzer.MAX_ANALYSIS_PEAKS
                if peakStats(p).count > 0 && peakStats(p).consistency > bestConsistency
                    bestPeakPosition = p;
                    bestConsistency = peakStats(p).consistency;
                end
            end

            alignmentReference = struct();
            alignmentReference.peakPosition = bestPeakPosition;
            alignmentReference.referenceTime = peakStats(bestPeakPosition).meanTime;
            alignmentReference.timeStd = peakStats(bestPeakPosition).stdTime;
            alignmentReference.consistency = peakStats(bestPeakPosition).consistency;
            alignmentReference.coverage = peakStats(bestPeakPosition).coverage;
            alignmentReference.majorityGroup = maxGroupIdx;
            alignmentReference.consensusStrength = majorityFraction;

            if verbose
                fprintf('Selected alignment reference:\n');
                fprintf('  Peak position: %d (of first %d peaks)\n', bestPeakPosition, RobustPeakAlignmentAnalyzer.MAX_ANALYSIS_PEAKS);
                fprintf('  Reference time: %.3f μs\n', alignmentReference.referenceTime*1e6);
                fprintf('  Time consistency: ±%.3f μs (std)\n', alignmentReference.timeStd*1e6);
                fprintf('  Coverage: %.1f%% of waveforms\n', alignmentReference.coverage*100);
                fprintf('  Consensus strength: %.1f%%\n', alignmentReference.consensusStrength*100);
            end
        end

        function qualityMetrics = validateAlignmentQuality(peakPatterns, alignmentReference, verbose)
            % VALIDATEALIGNMENTQUALITY - Validate the quality of the alignment reference

            if verbose
                fprintf('Validating alignment quality...\n');
            end

            peakPos = alignmentReference.peakPosition;
            refTime = alignmentReference.referenceTime;

            % Calculate alignment errors
            alignmentErrors = [];
            validAlignments = 0;

            for i = 1:length(peakPatterns)
                pattern = peakPatterns(i);
                if pattern.numPeaks >= peakPos
                    error = abs(pattern.peakTimes(peakPos) - refTime);
                    alignmentErrors = [alignmentErrors; error];
                    validAlignments = validAlignments + 1;
                end
            end

            qualityMetrics = struct();
            if ~isempty(alignmentErrors)
                qualityMetrics.meanError = mean(alignmentErrors);
                qualityMetrics.stdError = std(alignmentErrors);
                qualityMetrics.maxError = max(alignmentErrors);
                qualityMetrics.medianError = median(alignmentErrors);
                qualityMetrics.validAlignments = validAlignments;
                qualityMetrics.alignmentRate = validAlignments / length(peakPatterns);

                % Quality score (0-1, higher is better)
                qualityMetrics.qualityScore = alignmentReference.consistency * alignmentReference.coverage * alignmentReference.consensusStrength;
            else
                qualityMetrics.meanError = Inf;
                qualityMetrics.stdError = Inf;
                qualityMetrics.maxError = Inf;
                qualityMetrics.medianError = Inf;
                qualityMetrics.validAlignments = 0;
                qualityMetrics.alignmentRate = 0;
                qualityMetrics.qualityScore = 0;
            end

            if verbose
                fprintf('Alignment quality metrics:\n');
                fprintf('  Mean error: %.3f μs\n', qualityMetrics.meanError*1e6);
                fprintf('  Std error: ±%.3f μs\n', qualityMetrics.stdError*1e6);
                fprintf('  Max error: %.3f μs\n', qualityMetrics.maxError*1e6);
                fprintf('  Valid alignments: %d/%d (%.1f%%)\n', ...
                    qualityMetrics.validAlignments, length(peakPatterns), qualityMetrics.alignmentRate*100);
                fprintf('  Overall quality score: %.3f\n', qualityMetrics.qualityScore);
            end
        end

        function displayAnalysisResults(alignmentReference, analysisResults)
            % DISPLAYANALYSISRESULTS - Display comprehensive analysis results

            fprintf('\n=== ROBUST PEAK ALIGNMENT ANALYSIS RESULTS ===\n');
            fprintf('Dataset Overview:\n');
            fprintf('  Total waveforms: %d\n', analysisResults.totalWaveforms);
            fprintf('  Valid for analysis: %d (%.1f%%)\n', ...
                analysisResults.validWaveforms, 100*analysisResults.validWaveforms/analysisResults.totalWaveforms);

            fprintf('\nAlignment Reference Selected:\n');
            fprintf('  Peak position: %d\n', alignmentReference.peakPosition);
            fprintf('  Reference time: %.3f μs\n', alignmentReference.referenceTime*1e6);
            fprintf('  Consensus strength: %.1f%%\n', alignmentReference.consensusStrength*100);

            fprintf('\nQuality Assessment:\n');
            fprintf('  Alignment quality score: %.3f/1.0\n', analysisResults.qualityMetrics.qualityScore);
            fprintf('  Mean alignment error: %.3f μs\n', analysisResults.qualityMetrics.meanError*1e6);

            if alignmentReference.consensusStrength >= RobustPeakAlignmentAnalyzer.MAJORITY_THRESHOLD
                fprintf('  Status: ✅ ROBUST CONSENSUS ACHIEVED\n');
            else
                fprintf('  Status: ⚠️  WEAK CONSENSUS - Consider manual review\n');
            end

            fprintf('===============================================\n\n');
        end
    end
end
