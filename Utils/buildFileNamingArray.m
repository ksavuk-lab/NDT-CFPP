function FileNamingArray = buildFileNamingArray(varargin)
% buildFileNamingArray - Constructor helper to create a consistent FileNamingArray
%
% Usage patterns (name-value pairs):
%   FileNamingArray = buildFileNamingArray( ...
%       'DATASET', 3, ...
%       'CaseNumber', 3, ...
%       'TrimWaveData', 1, ...
%       'TrimTimeRange', 1, ...
%       'xRange', [10 70], ...
%       'yRange', [5 20], ...
%       'TimeRange', [0 2e-6], ...
%       'AlignAtFirstPeak', 0);
%
% Returns numeric vector:
%   [DATASET, caseNumber, TrimWaveData, TrimTimeRange, x1, x2, y1, y2, t1, t2, alignAtFirstPeak]
%
% Notes:
% - If CaseNumber is omitted, defaults to DATASET (current convention in Main.m)
% - AlignAtFirstPeak defaults to 0 if omitted

p = inputParser;
addParameter(p, 'DATASET', [], @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'CaseNumber', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'TrimWaveData', 0, @(x) ismember(x,[0 1]));
addParameter(p, 'TrimTimeRange', 0, @(x) ismember(x,[0 1]));
addParameter(p, 'xRange', [NaN NaN], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'yRange', [NaN NaN], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'TimeRange', [NaN NaN], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'AlignAtFirstPeak', 0, @(x) ismember(x,[0 1]));
parse(p, varargin{:});

DATASET = p.Results.DATASET;
caseNumber = p.Results.CaseNumber;
if isempty(caseNumber), caseNumber = DATASET; end
TrimWaveData = p.Results.TrimWaveData;
TrimTimeRange = p.Results.TrimTimeRange;
xRange = p.Results.xRange;
yRange = p.Results.yRange;
TimeRange = p.Results.TimeRange;
alignAtFirstPeak = p.Results.AlignAtFirstPeak;

% Basic validations
assert(~isempty(DATASET), 'DATASET is required');
assert(all(~isnan(xRange)), 'xRange must be [x1 x2]');
assert(all(~isnan(yRange)), 'yRange must be [y1 y2]');
assert(all(~isnan(TimeRange)), 'TimeRange must be [t1 t2]');

FileNamingArray = [DATASET, caseNumber, TrimWaveData, TrimTimeRange, xRange(:).', yRange(:).', TimeRange(:).', alignAtFirstPeak];
end

