afunction results = run_all_tests()
% RUN_ALL_TESTS - Convenience entrypoint to execute all unit/smoke tests
% Usage from MATLAB:
%   results = run_all_tests();
% Usage from CLI:
%   matlab -batch "results = run_all_tests(); disp(results)"

% Ensure Utils (where addUtilsPath lives) is on path before calling it
thisDir = fileparts(mfilename('fullpath'));
utilsDir = fullfile(thisDir, 'Utils');
if exist(utilsDir, 'dir') && ~contains(path, utilsDir)
    addpath(utilsDir);
end

% Now centralize all project paths
addUtilsPath();

try
    results = runtests('Tests');
    disp(results);
    % Not all MATLAB releases have 'Skipped'
    hasSkipped = isprop(results, 'Skipped');
    passed = nnz([results.Passed]);
    failed = nnz([results.Failed]);
    incomplete = nnz([results.Incomplete]);
    if hasSkipped
        skipped = nnz([results.Skipped]);
        fprintf('\nSummary: %d Passed, %d Failed, %d Incomplete, %d Skipped\n', passed, failed, incomplete, skipped);
    else
        fprintf('\nSummary: %d Passed, %d Failed, %d Incomplete\n', passed, failed, incomplete);
    end
catch ME
    fprintf('Test run failed: %s\n', ME.message);
    rethrow(ME);
end
end

