%% LOADDATA Load Data: Puts Raw Data into Matrix
% Initialize DataStructure as an empty structure
% This ensures that DataStructure is defined even if the loading process fails
% Supports multiple directory paths and tries each one until successful

function [DataStructure] = loadData(DATASET)
    DataStructure = struct();

    % Define multiple base directory paths to try (dynamic first, then fallbacks)
    baseDirs = {};

    % 1) Environment override (export NDE_DATA_ROOT=/path/to/Data)
    dataRootEnv = getenv('NDE_DATA_ROOT');
    if ~isempty(dataRootEnv)
        baseDirs{end+1} = fullfile(dataRootEnv, 'LaminateBeam', 'L0', filesep); %#ok<AGROW>
        baseDirs{end+1} = [dataRootEnv filesep]; 
    end

    % 2) Relative to repo root: ../../Data/LaminateBeam/L0 from code root
    try
        repoRelative = fullfile(pwd, '..', '..', 'Data', 'LaminateBeam', 'L0', filesep);
        baseDirs{end+1} = repoRelative; 
    catch
    end

    % 3) User HOME-based typical path
    try
        homeDir = getenv('HOME');
        xUni = fullfile(homeDir, 'Documents', 'X_ University', 'Uni Research', 'NDT', 'NDE_SDSU-USD', 'Data', 'LaminateBeam', 'L0', filesep);
        baseDirs{end+1} = xUni; 
    catch
    end

    % 4) Legacy hardcoded paths (fallbacks)
    baseDirs{end+1} = '/Users/ksavuk/Documents/X_ University/Uni Research/NDT/NDE_SDSU-USD/Data/LaminateBeam/L0/';
    baseDirs{end+1} = '/Users/ksavuk/Documents/University or Personal Education/Uni Research/NDT/NDE_SDSU-USD/Data/LaminateBeam/L0/';
    baseDirs{end+1} = '/Users/kostia/Documents/University or Personal Education/Uni Research/NDT/NDE_SDSU-USD/Data/LaminateBeam/L0/';

    % Define file names for each dataset
    fileNames = {};
    switch DATASET
        case 1
            fileNames{1} = 'L0P5S25_10Mhz_TimingCorrected.mat';
        case 2
            fileNames{1} = 'L8P12S5_10Mhz.mat';
        case 3
            fileNames{1} = 'L16P1S9_10MHz.mat';
        case 4
            fileNames{1} = 'L16P1S10_10MHz.mat';
        case 5
            fileNames{1} = 'L16P34S3_10MHz.mat';
        otherwise
            error('Invalid dataset number'); % Throw an error for any other case
    end

    tic;
    disp('Start Load Data');

    % Try each directory path until successful
    loadSuccess = false;
    errorMessages = {};

    for i = 1:length(baseDirs)
        fullPath = [baseDirs{i}, fileNames{1}];
        try
            disp(['Trying to load from: ', fullPath]);
            temp = load(fullPath);
            % If we get here, the load was successful
            loadSuccess = true;
            disp(['Successfully loaded from: ', fullPath]);
            break;
        catch ME
            % Store the error message and continue to the next directory
            errorMessages{end+1} = ['Failed to load from ', fullPath, ': ', ME.message];
            disp(errorMessages{end});
        end
    end

    % Check if any load attempt was successful
    if ~loadSuccess
        % Combine all error messages
        allErrors = strjoin(errorMessages, '\n');
        error(['Failed to load data from any directory:\n', allErrors]);
    end

    % Extract 'datastructure' from temp var or preprocess if needed
    if isfield(temp, 'datastructure')
        % Standard format - extract datastructure directly
        DataStructure = temp.datastructure;
    else
        % Special case: Check if this is the L0P5S25 dataset that needs preprocessing
        fieldNames = fieldnames(temp);
        if length(fieldNames) == 1 && contains(fieldNames{1}, 'plyWaviness') && istable(temp.(fieldNames{1}))
            disp('Detected L0P5S25 dataset format - applying preprocessing...');
            DataStructure = preprocessL0P5S25Dataset(fullPath);
        else
            % Unknown format
            disp('Available fields in loaded data:');
            for i = 1:length(fieldNames)
                disp(['  ', fieldNames{i}, ' (', class(temp.(fieldNames{i})), ')']);
            end
            error(['Unrecognized data format. Expected ''datastructure'' field but found: ', strjoin(fieldNames, ', ')]);
        end
    end

    disp(['Load Data: ', num2str(toc), ' seconds']); % Time to load data into var.
end
