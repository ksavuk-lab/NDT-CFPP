function updateProgress(progressHandle, currentValue, message, currentFunction)
% UPDATEPROGRESS - Enhanced progress dialog with detailed information
%
% Standalone version so other files (e.g., SegmentPrecomputedWaveform) can
% call it directly.
%
% INPUTS:
%   progressHandle  - Handle returned by ProgressDialog
%   currentValue    - Current progress value
%   message         - Status message to display (optional)
%   currentFunction - Current function/task being executed (optional)

if nargin < 1 || ~isstruct(progressHandle)
    return; % Nothing to update
end

if ~isfield(progressHandle, 'fig') || ~ishghandle(progressHandle.fig)
    return; % Progress dialog was closed
end

try
    % Update current value
    if isfield(progressHandle,'maxValue')
        progressHandle.currentValue = min(currentValue, progressHandle.maxValue);
        maxValue = progressHandle.maxValue;
    else
        progressHandle.currentValue = currentValue;
        maxValue = max(1, currentValue);
    end

    % Calculate percentages
    percentComplete = (progressHandle.currentValue / maxValue) * 100;
    percentRemaining = 100 - percentComplete;

    % Calculate elapsed time and estimate remaining time
    if ~isfield(progressHandle,'startTime') || isempty(progressHandle.startTime)
        progressHandle.startTime = tic;
    end
    elapsedTime = toc(progressHandle.startTime);
    if progressHandle.currentValue > 0
        estimatedTotalTime = elapsedTime * (maxValue / progressHandle.currentValue);
        estimatedRemainingTime = estimatedTotalTime - elapsedTime;
    else
        estimatedRemainingTime = 0;
    end

    % Update progress bar width
    barWidth = (progressHandle.currentValue / maxValue) * 360;
    if isfield(progressHandle,'progressBar') && ishghandle(progressHandle.progressBar)
        set(progressHandle.progressBar, 'Position', [20, 70, barWidth, 20]);
    end

    % Update percentage text with both completed and remaining
    if isfield(progressHandle,'percentText') && ishghandle(progressHandle.percentText)
        percentString = sprintf('%.1f%% Complete | %.1f%% Remaining', percentComplete, percentRemaining);
        set(progressHandle.percentText, 'String', percentString);
    end

    % Update current function if provided
    if nargin >= 4 && ~isempty(currentFunction)
        progressHandle.currentFunction = currentFunction;
    end

    % Update task text
    taskString = sprintf('Function: %s', getfield(progressHandle,'currentFunction','Working...')); %#ok<GFLD>
    taskHandle = findobj(progressHandle.fig, 'Tag', 'TaskText');
    if ~isempty(taskHandle)
        set(taskHandle, 'String', taskString);
    end

    % Update time estimation
    timeHandle = findobj(progressHandle.fig, 'Tag', 'TimeText');
    if ~isempty(timeHandle)
        if estimatedRemainingTime > 0
            if estimatedRemainingTime < 60
                timeString = sprintf('Time: ~%.0f seconds remaining', estimatedRemainingTime);
            elseif estimatedRemainingTime < 3600
                timeString = sprintf('Time: ~%.1f minutes remaining', estimatedRemainingTime / 60);
            else
                timeString = sprintf('Time: ~%.1f hours remaining', estimatedRemainingTime / 3600);
            end
        else
            timeString = 'Time: Calculating...';
        end
        set(timeHandle, 'String', timeString);
    end

    % Update main message if provided
    if nargin >= 3 && ~isempty(message)
        messageHandle = findobj(progressHandle.fig, 'Tag', 'MessageText');
        if ~isempty(messageHandle)
            set(messageHandle, 'String', message);
        end
    end

    % Force display update
    drawnow;

catch ME %#ok<NASGU>
    % Silently ignore progress update errors
end
end

