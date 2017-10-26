function [bouts, nSyllables] = arrangeBouts(events, params, varargin)
%ARRANGEBOUTS group events into longer groups 
%
% arrangeBouts(events) groups events into 'bouts', with a minimum silent period
% between each bout.  Pass with params as 'silentPeriod'

% housekeeping
if nargin < 2 || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});
[bouts,nSyllables] = mergeGaps(events, params.silentPeriod);
% sort events
%{
events = sortBy(events,'start');

boutGap = params.silentPeriod; % in seconds
evStarts = [events.start];
evEnds   = [events.stop ];
interSong = [0 evStarts(2:end) - evEnds(1:end-1)]; % distance between syllables
gaps = find(interSong > boutGap);
boutStarts = [1 gaps];
boutEnds =  [gaps - 1 numel(interSong)];

nBouts = numel(boutStarts);
bouts = initEvents(nBouts);
for ii = 1:nBouts
    bouts(ii) = bookendedClip(events(boutStarts(ii):boutEnds(ii)));
end
nSyllables = boutEnds - boutStarts + 1;
%}
end

