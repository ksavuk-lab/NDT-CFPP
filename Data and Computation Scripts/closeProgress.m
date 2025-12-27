function closeProgress(progressHandle)
% CLOSEPROGRESS - Safely close a ProgressDialog-like handle
% Accepts the struct returned by ProgressDialog and closes its figure if it exists.

if nargin < 1 || ~isstruct(progressHandle)
    return;
end

try
    if isfield(progressHandle,'fig') && ishghandle(progressHandle.fig)
        close(progressHandle.fig);
    end
catch %#ok<CTCH>
    % Swallow any UI errors silently
end
end

