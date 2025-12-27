classdef AlignmentStatisticsReporter < handle
    % ALIGNMENTSTATISTICSREPORTER - Comprehensive statistics reporting for cross-sectional alignment
    % Provides detailed cost improvement tracking, convergence analysis, and performance metrics
    
    properties
        iterationHistory = [];
        startTime = [];
        totalSlicesProcessed = 0;
        totalShiftsApplied = 0;
        convergenceAchieved = false;
        finalImprovement = 0;
    end
    
    methods
        function obj = AlignmentStatisticsReporter()
            % Constructor - initialize timing
            obj.startTime = tic;
            obj.iterationHistory = struct('iteration', {}, 'yCost', {}, 'xCost', {}, ...
                'totalCost', {}, 'improvement', {}, 'yTime', {}, 'xTime', {}, ...
                'totalTime', {}, 'ySlicesAligned', {}, 'xSlicesAligned', {});
        end
        
        function recordIteration(obj, iteration, yCost, xCost, totalCost, improvement, ...
                yTime, xTime, ySlicesAligned, xSlicesAligned)
            % Record statistics for a single iteration
            
            totalTime = toc(obj.startTime);
            
            % Store iteration data
            obj.iterationHistory(end+1) = struct(...
                'iteration', iteration, ...
                'yCost', yCost, ...
                'xCost', xCost, ...
                'totalCost', totalCost, ...
                'improvement', improvement, ...
                'yTime', yTime, ...
                'xTime', xTime, ...
                'totalTime', totalTime, ...
                'ySlicesAligned', ySlicesAligned, ...
                'xSlicesAligned', xSlicesAligned);
            
            obj.totalSlicesProcessed = obj.totalSlicesProcessed + ySlicesAligned + xSlicesAligned;
        end
        
        function displayIterationSummary(obj, iteration, yCost, xCost, totalCost, improvement, ...
                yTime, xTime, ySlicesAligned, xSlicesAligned, convergenceThreshold)
            % Display detailed iteration summary
            
            fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
            fprintf('║                    ITERATION %2d SUMMARY                      ║\n', iteration);
            fprintf('╠══════════════════════════════════════════════════════════════╣\n');
            
            % Cost breakdown
            fprintf('║ ALIGNMENT COSTS:                                             ║\n');
            fprintf('║   Y-slice alignment cost: %12.6f                    ║\n', yCost);
            fprintf('║   X-slice alignment cost: %12.6f                    ║\n', xCost);
            fprintf('║   Total alignment cost:   %12.6f                    ║\n', totalCost);
            
            % Improvement metrics
            if iteration > 1
                fprintf('║                                                              ║\n');
                fprintf('║ IMPROVEMENT METRICS:                                         ║\n');
                fprintf('║   Cost improvement:       %8.4f%% ', improvement*100);
                if improvement >= convergenceThreshold
                    fprintf('(GOOD)              ║\n');
                else
                    fprintf('(CONVERGED)         ║\n');
                end
                fprintf('║   Convergence threshold:  %8.4f%%                     ║\n', convergenceThreshold*100);
                
                % Improvement trend
                if length(obj.iterationHistory) >= 2
                    prevImprovement = obj.iterationHistory(end-1).improvement;
                    trend = improvement - prevImprovement;
                    if trend > 0
                        fprintf('║   Improvement trend:      ↗ INCREASING (+%.4f%%)          ║\n', trend*100);
                    elseif trend < 0
                        fprintf('║   Improvement trend:      ↘ DECREASING (%.4f%%)          ║\n', trend*100);
                    else
                        fprintf('║   Improvement trend:      → STABLE                        ║\n');
                    end
                end
            end
            
            % Processing statistics
            fprintf('║                                                              ║\n');
            fprintf('║ PROCESSING STATISTICS:                                       ║\n');
            fprintf('║   Y-slices processed:     %3d/101 (%5.1f%%)                ║\n', ...
                ySlicesAligned, 101, (ySlicesAligned/101)*100);
            fprintf('║   X-slices processed:     %3d/101 (%5.1f%%)                ║\n', ...
                xSlicesAligned, 101, (xSlicesAligned/101)*100);
            fprintf('║   Total slices this iter: %3d                               ║\n', ...
                ySlicesAligned + xSlicesAligned);
            
            % Timing information
            fprintf('║                                                              ║\n');
            fprintf('║ TIMING BREAKDOWN:                                            ║\n');
            fprintf('║   Y-slice alignment:      %8.1f seconds                   ║\n', yTime);
            fprintf('║   X-slice alignment:      %8.1f seconds                   ║\n', xTime);
            fprintf('║   Iteration total:        %8.1f seconds                   ║\n', yTime + xTime);
            fprintf('║   Cumulative time:        %8.1f seconds                   ║\n', toc(obj.startTime));
            
            % Performance metrics
            totalSlicesThisIter = ySlicesAligned + xSlicesAligned;
            if totalSlicesThisIter > 0
                avgTimePerSlice = (yTime + xTime) / totalSlicesThisIter;
                fprintf('║   Avg time per slice:     %8.3f seconds                   ║\n', avgTimePerSlice);
                
                % Estimate remaining time
                if iteration < 10 % Assuming max 10 iterations
                    remainingIters = 10 - iteration;
                    estimatedRemaining = remainingIters * (yTime + xTime);
                    fprintf('║   Est. remaining time:    %8.1f seconds                   ║\n', estimatedRemaining);
                end
            end
            
            fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');
        end
        
        function displayFinalSummary(obj, finalIteration, converged, convergenceThreshold)
            % Display comprehensive final summary
            
            totalTime = toc(obj.startTime);
            
            fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
            fprintf('║                    FINAL ALIGNMENT SUMMARY                   ║\n');
            fprintf('╠══════════════════════════════════════════════════════════════╣\n');
            
            % Convergence status
            fprintf('║ CONVERGENCE STATUS:                                          ║\n');
            if converged
                fprintf('║   Status: ✅ CONVERGED after %d iterations                 ║\n', finalIteration);
                fprintf('║   Final improvement: %.4f%% (< %.4f%% threshold)        ║\n', ...
                    obj.finalImprovement*100, convergenceThreshold*100);
            else
                fprintf('║   Status: ⚠️  MAXIMUM ITERATIONS REACHED (%d)              ║\n', finalIteration);
                fprintf('║   Final improvement: %.4f%% (≥ %.4f%% threshold)        ║\n', ...
                    obj.finalImprovement*100, convergenceThreshold*100);
            end
            
            % Overall improvement
            if length(obj.iterationHistory) >= 2
                initialCost = obj.iterationHistory(1).totalCost;
                finalCost = obj.iterationHistory(end).totalCost;
                totalImprovement = (initialCost - finalCost) / initialCost;
                
                fprintf('║                                                              ║\n');
                fprintf('║ OVERALL IMPROVEMENT:                                         ║\n');
                fprintf('║   Initial cost:           %12.6f                    ║\n', initialCost);
                fprintf('║   Final cost:             %12.6f                    ║\n', finalCost);
                fprintf('║   Total improvement:      %8.4f%% reduction             ║\n', totalImprovement*100);
                
                if totalImprovement > 0.1
                    fprintf('║   Quality assessment:     🟢 EXCELLENT (>10%% improvement)  ║\n');
                elseif totalImprovement > 0.05
                    fprintf('║   Quality assessment:     🟡 GOOD (5-10%% improvement)     ║\n');
                elseif totalImprovement > 0.01
                    fprintf('║   Quality assessment:     🟠 MODERATE (1-5%% improvement)  ║\n');
                else
                    fprintf('║   Quality assessment:     🔴 MINIMAL (<1%% improvement)    ║\n');
                end
            end
            
            % Processing statistics
            fprintf('║                                                              ║\n');
            fprintf('║ PROCESSING STATISTICS:                                       ║\n');
            fprintf('║   Total iterations:       %3d                               ║\n', finalIteration);
            fprintf('║   Total slices processed: %5d                             ║\n', obj.totalSlicesProcessed);
            fprintf('║   Avg slices per iter:    %5.1f                             ║\n', obj.totalSlicesProcessed/finalIteration);
            
            % Timing summary
            fprintf('║                                                              ║\n');
            fprintf('║ TIMING SUMMARY:                                              ║\n');
            fprintf('║   Total alignment time:   %8.1f seconds                   ║\n', totalTime);
            fprintf('║   Average per iteration:  %8.1f seconds                   ║\n', totalTime/finalIteration);
            if obj.totalSlicesProcessed > 0
                fprintf('║   Average per slice:      %8.3f seconds                   ║\n', totalTime/obj.totalSlicesProcessed);
            end
            
            % Performance rating
            avgTimePerSlice = totalTime / obj.totalSlicesProcessed;
            if avgTimePerSlice < 0.5
                fprintf('║   Performance rating:     🚀 FAST (<0.5s per slice)        ║\n');
            elseif avgTimePerSlice < 2.0
                fprintf('║   Performance rating:     ⚡ GOOD (0.5-2s per slice)       ║\n');
            elseif avgTimePerSlice < 5.0
                fprintf('║   Performance rating:     🐌 SLOW (2-5s per slice)         ║\n');
            else
                fprintf('║   Performance rating:     🐢 VERY SLOW (>5s per slice)     ║\n');
            end
            
            fprintf('╚══════════════════════════════════════════════════════════════╝\n\n');
        end
        
        function displayProgressChart(obj)
            % Display ASCII progress chart of cost improvements
            
            if length(obj.iterationHistory) < 2
                return;
            end
            
            fprintf('COST IMPROVEMENT PROGRESS CHART:\n');
            fprintf('═══════════════════════════════════════════════════════════════\n');
            
            % Extract costs and improvements
            costs = [obj.iterationHistory.totalCost];
            improvements = [obj.iterationHistory.improvement];
            
            % Normalize for display (0-50 characters wide)
            maxCost = max(costs);
            minCost = min(costs);
            costRange = maxCost - minCost;
            
            for i = 1:length(costs)
                % Cost bar
                if costRange > 0
                    barLength = round(40 * (maxCost - costs(i)) / costRange);
                else
                    barLength = 20;
                end
                
                costBar = repmat('█', 1, barLength);
                costBar = [costBar, repmat('░', 1, 40-barLength)];
                
                % Improvement indicator
                if i > 1
                    if improvements(i) > 0.01
                        indicator = '📈';
                    elseif improvements(i) > 0.001
                        indicator = '📊';
                    else
                        indicator = '📉';
                    end
                else
                    indicator = '🎯';
                end
                
                fprintf('Iter %2d: %s %s Cost: %.6f', i, indicator, costBar, costs(i));
                if i > 1
                    fprintf(' (%.3f%% improvement)', improvements(i)*100);
                end
                fprintf('\n');
            end
            
            fprintf('═══════════════════════════════════════════════════════════════\n\n');
        end
        
        function setConverged(obj, finalImprovement)
            % Mark alignment as converged
            obj.convergenceAchieved = true;
            obj.finalImprovement = finalImprovement;
        end
    end
end
