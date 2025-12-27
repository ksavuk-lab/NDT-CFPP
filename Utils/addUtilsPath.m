function addUtilsPath()
    % ADDUTILSPATH - Idempotently add project utility folders to MATLAB path
    %
    % Resolves project root from this file's location to avoid dependence on pwd.

    % This file lives in <projectRoot>/Utils/addUtilsPath.m
    thisFile = mfilename('fullpath');
    utilsDir = fileparts(thisFile);
    projectRoot = fileparts(utilsDir);

    dirsToAdd = {
        fullfile(projectRoot, 'Utils')
        fullfile(projectRoot, 'Tools')
        fullfile(projectRoot, 'Data and Computation Scripts')
    };

    for i = 1:numel(dirsToAdd)
        dirPath = dirsToAdd{i};
        if exist(dirPath, 'dir')
            if ~isOnPath(dirPath)
                addpath(dirPath);
                fprintf('Added directory to MATLAB path: %s\n', dirPath);
            end
        end
    end
end

function tf = isOnPath(dirPath)
    % Check path membership using pathsep boundaries to avoid partial matches
    p = path;
    if isempty(p)
        tf = false; return;
    end
    if endsWith(p, pathsep)
        p = p(1:end-1);
    end
    tf = contains([pathsep p pathsep], [pathsep dirPath pathsep]);
end
