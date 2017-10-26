function [regions,nPer] = mergeGaps(events, gapLength)
% gap length is in the units that events times are in

if isempty(events), regions = []; nPer = []; return; end;

% equivalent to arrangeBouts
events = sortBy(events,'start');

evStarts = [events.start];
evEnds   = [events.stop ];
interSong = [0 (evStarts(2:end) - evEnds(1:end-1))]; % distance between syllables
gaps = find(interSong > gapLength);
mergeStarts = [1 gaps];
mergeStops =  [(gaps - 1) numel(interSong)];

nRegions = numel(mergeStarts);
regions = initEvents(nRegions);
for ii = 1:nRegions
    regions(ii) = bookendedClip(events(mergeStarts(ii):mergeStops(ii)));
end
nPer = mergeStops - mergeStarts + 1;
end