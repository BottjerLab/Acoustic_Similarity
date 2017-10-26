function regions = adjustTimeStamps(regions, timeStart, fs)
%ADJUSTTIMESTAMPS moves regions forwards/backwards in time
% 
%   regions = adjustTimeStamps(regions, sampleStart) shifts region 
%   starts/stops by sampleStart (in seconds).
%   If sampleStart is negative, then the regions are moved
%   backwards in time.

if isempty(regions), return; end;
if ~isEvent(regions)
    error('badArgs','First argument should be an event structure');
end;

% warning: this estimation may not work well if regions are very close to 0 time start
if nargin < 3 && numel(regions) >= 1
    fs = regions(1).idxStart / regions(1).start;
end
types = {regions.type};
regions = eventFromStamps([regions.idxStart] + timeStart * fs, ...
    [regions.idxStop] + timeStart * fs, fs);
[regions.type] = types{:};
end
