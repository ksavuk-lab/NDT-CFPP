function progressHandle = ProgressDialog(title, message, maxValue)
% PROGRESSDIALOG - Create a non-blocking progress dialog window
%
% This function creates a progress dialog that doesn't gray out the main window
% and provides real-time updates during plate generation processing.
%
% INPUTS:
%   title    - Title for the progress window
%   message  - Initial message to display
%   maxValue - Maximum value for progress bar (optional, default 100)
%
% OUTPUTS:
%   progressHandle - Handle to progress dialog structure
%
% USAGE:
%   % Create progress dialog
%   progress = ProgressDialog('Processing', 'Starting plate generation...', 1000);
%   
%   % Update progress
%   updateProgress(progress, 250, 'Processing layer 1...');
%   
%   % Close when done
%   closeProgress(progress);

if nargin < 3
    maxValue = 100;
end

% Create figure for progress dialog
progressHandle = struct();
progressHandle.maxValue = maxValue;
progressHandle.currentValue = 0;
progressHandle.startTime = tic; % Track start time for ETA calculation
progressHandle.lastUpdateTime = tic;
progressHandle.currentFunction = 'Initializing';

% Calculate position (center of screen)
screenSize = get(0, 'ScreenSize');
figWidth = 400;
figHeight = 150;
figX = (screenSize(3) - figWidth) / 2;
figY = (screenSize(4) - figHeight) / 2;

% Create progress figure
progressHandle.fig = figure('Name', title, ...
                           'NumberTitle', 'off', ...
                           'MenuBar', 'none', ...
                           'ToolBar', 'none', ...
                           'Resize', 'off', ...
                           'Position', [figX, figY, figWidth, figHeight], ...
                           'WindowStyle', 'normal', ... % Non-modal
                           'CloseRequestFcn', @(src,evt) []); % Prevent closing

% Create UI elements
uicontrol('Parent', progressHandle.fig, ...
          'Style', 'text', ...
          'String', message, ...
          'Position', [20, 100, 360, 30], ...
          'FontSize', 10, ...
          'HorizontalAlignment', 'left', ...
          'Tag', 'MessageText');

% Progress bar background
uicontrol('Parent', progressHandle.fig, ...
          'Style', 'text', ...
          'String', '', ...
          'Position', [20, 70, 360, 20], ...
          'BackgroundColor', [0.9, 0.9, 0.9], ...
          'Tag', 'ProgressBackground');

% Progress bar foreground
progressHandle.progressBar = uicontrol('Parent', progressHandle.fig, ...
                                      'Style', 'text', ...
                                      'String', '', ...
                                      'Position', [20, 70, 1, 20], ...
                                      'BackgroundColor', [0.2, 0.6, 0.2], ...
                                      'Tag', 'ProgressBar');

% Progress percentage text (completed and remaining)
progressHandle.percentText = uicontrol('Parent', progressHandle.fig, ...
                                      'Style', 'text', ...
                                      'String', '0% Complete | 100% Remaining', ...
                                      'Position', [20, 45, 360, 15], ...
                                      'FontSize', 9, ...
                                      'HorizontalAlignment', 'center', ...
                                      'Tag', 'PercentText');

% Current function/task text
progressHandle.taskText = uicontrol('Parent', progressHandle.fig, ...
                                   'Style', 'text', ...
                                   'String', 'Function: Initializing...', ...
                                   'Position', [20, 30, 360, 15], ...
                                   'FontSize', 8, ...
                                   'HorizontalAlignment', 'left', ...
                                   'ForegroundColor', [0.2, 0.2, 0.6], ...
                                   'Tag', 'TaskText');

% Estimated time remaining text
progressHandle.timeText = uicontrol('Parent', progressHandle.fig, ...
                                   'Style', 'text', ...
                                   'String', 'Time: Calculating...', ...
                                   'Position', [20, 15, 360, 15], ...
                                   'FontSize', 8, ...
                                   'HorizontalAlignment', 'left', ...
                                   'ForegroundColor', [0.6, 0.2, 0.2], ...
                                   'Tag', 'TimeText');

% Force display
drawnow;

end

function updateProgress(progressHandle, currentValue, message, currentFunction)
% UPDATEPROGRESS - Enhanced progress dialog with detailed information
%
% INPUTS:
%   progressHandle  - Handle returned by ProgressDialog
%   currentValue    - Current progress value
%   message         - Status message to display (optional)
%   currentFunction - Current function/task being executed (optional)

if ~isfield(progressHandle, 'fig') || ~isvalid(progressHandle.fig)
    return; % Progress dialog was closed
end

try
    % Update current value
    progressHandle.currentValue = min(currentValue, progressHandle.maxValue);

    % Calculate percentages
    percentComplete = (progressHandle.currentValue / progressHandle.maxValue) * 100;
    percentRemaining = 100 - percentComplete;

    % Calculate elapsed time and estimate remaining time
    elapsedTime = toc(progressHandle.startTime);
    if progressHandle.currentValue > 0
        estimatedTotalTime = elapsedTime * (progressHandle.maxValue / progressHandle.currentValue);
        estimatedRemainingTime = estimatedTotalTime - elapsedTime;
    else
        estimatedRemainingTime = 0;
    end

    % Update progress bar width
    barWidth = (progressHandle.currentValue / progressHandle.maxValue) * 360;
    set(progressHandle.progressBar, 'Position', [20, 70, barWidth, 20]);

    % Update percentage text with both completed and remaining
    percentString = sprintf('%.1f%% Complete | %.1f%% Remaining', percentComplete, percentRemaining);
    set(progressHandle.percentText, 'String', percentString);

    % Update current function if provided
    if nargin >= 4 && ~isempty(currentFunction)
        progressHandle.currentFunction = currentFunction;
    end

    % Update task text
    taskString = sprintf('Function: %s', progressHandle.currentFunction);
    taskHandle = findobj(progressHandle.fig, 'Tag', 'TaskText');
    if ~isempty(taskHandle)
        set(taskHandle, 'String', taskString);
    end

    % Update time estimation
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

    timeHandle = findobj(progressHandle.fig, 'Tag', 'TimeText');
    if ~isempty(timeHandle)
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

catch ME
    % Handle errors gracefully
    fprintf('Progress update error: %s\n', ME.message);
end

end

function closeProgress(progressHandle)
% CLOSEPROGRESS - Close progress dialog
%
% INPUTS:
%   progressHandle - Handle returned by ProgressDialog

try
    if isfield(progressHandle, 'fig') && isvalid(progressHandle.fig)
        close(progressHandle.fig);
    end
catch ME
    fprintf('Progress close error: %s\n', ME.message);
end

end
