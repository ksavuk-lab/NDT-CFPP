function clearPlateCache(varargin)
% CLEARPLATECACHE - Clear cached plate generation data
%
% This function removes cached plate data files to force recomputation.
% Useful when you want to ensure fresh calculations or when debugging.
%
% Usage:
%   clearPlateCache()           - Clear all plate cache files
%   clearPlateCache('all')      - Clear all plate cache files
%   clearPlateCache('pattern')  - Clear files matching pattern (e.g., 'Amp0.100*')
%
% Examples:
%   clearPlateCache()                    % Clear all
%   clearPlateCache('Amp0.100*')         % Clear specific amplitude tolerance
%   clearPlateCache('*Time1.0*')         % Clear specific time tolerance
%   clearPlateCache('*both*')            % Clear all 'both' type plates

% Default to clearing all if no arguments
if nargin == 0 || (nargin == 1 && strcmp(varargin{1}, 'all'))
    pattern = 'PlateCache_*.mat';
    clearAll = true;
else
    pattern = ['PlateCache_*' varargin{1} '*.mat'];
    clearAll = false;
end

% Define cache folder
plateCacheFolder = fullfile(pwd, 'Plate Cache');

if ~exist(plateCacheFolder, 'dir')
    fprintf('No plate cache folder found. Nothing to clear.\n');
    return;
end

% Find matching files
cacheFiles = dir(fullfile(plateCacheFolder, pattern));

if isempty(cacheFiles)
    if clearAll
        fprintf('No plate cache files found to clear.\n');
    else
        fprintf('No plate cache files matching pattern "%s" found.\n', varargin{1});
    end
    return;
end

% Display what will be cleared
fprintf('\n=== Clearing Plate Cache ===\n');
if clearAll
    fprintf('Found %d plate cache files to clear:\n', length(cacheFiles));
else
    fprintf('Found %d plate cache files matching pattern "%s":\n', length(cacheFiles), varargin{1});
end

totalSize = 0;
for i = 1:length(cacheFiles)
    filePath = fullfile(plateCacheFolder, cacheFiles(i).name);
    fileSize = cacheFiles(i).bytes;
    totalSize = totalSize + fileSize;
    
    % Display file info
    if fileSize > 1024*1024
        sizeStr = sprintf('%.1f MB', fileSize / (1024*1024));
    elseif fileSize > 1024
        sizeStr = sprintf('%.1f KB', fileSize / 1024);
    else
        sizeStr = sprintf('%d bytes', fileSize);
    end
    
    fprintf('  %s (%s)\n', cacheFiles(i).name, sizeStr);
end

% Display total size
if totalSize > 1024*1024
    totalSizeStr = sprintf('%.1f MB', totalSize / (1024*1024));
elseif totalSize > 1024
    totalSizeStr = sprintf('%.1f KB', totalSize / 1024);
else
    totalSizeStr = sprintf('%d bytes', totalSize);
end
fprintf('Total size: %s\n', totalSizeStr);

% Confirm deletion
if length(cacheFiles) > 5 || totalSize > 10*1024*1024  % More than 5 files or 10MB
    response = input('Are you sure you want to delete these files? (y/N): ', 's');
    if ~strcmpi(response, 'y') && ~strcmpi(response, 'yes')
        fprintf('Cache clearing cancelled.\n');
        return;
    end
end

% Delete files
deletedCount = 0;
failedCount = 0;

for i = 1:length(cacheFiles)
    filePath = fullfile(plateCacheFolder, cacheFiles(i).name);
    try
        delete(filePath);
        deletedCount = deletedCount + 1;
    catch ME
        fprintf('Warning: Failed to delete %s (%s)\n', cacheFiles(i).name, ME.message);
        failedCount = failedCount + 1;
    end
end

% Report results
fprintf('\n=== Cache Clearing Complete ===\n');
fprintf('Successfully deleted: %d files\n', deletedCount);
if failedCount > 0
    fprintf('Failed to delete: %d files\n', failedCount);
end
fprintf('Freed disk space: %s\n', totalSizeStr);
fprintf('=====================================\n\n');

% Also clear any persistent variables in SmartDataLoader if clearing all
if clearAll
    try
        clearSmartDataLoaderCache();
        fprintf('Also cleared SmartDataLoader cache.\n');
    catch
        % Function may not exist, ignore
    end
end

end
