function baseEvents = constructBaseline(events, keyName, startOffsets, stopOffsets, params)
% does automatic exclusion 
% startOffsets is 2 x 1 vector in seconds, and is SIGNED (to go backwards
% in time, go negative)
% stopOffsets is similar
% if we only want one of st artOfsets/stopOffsets, blank the other ([])
fs = params.fs;

% % allow overlap?
% allowSelfOverlap = false; Jenny commented out because not used
% allowOtherOverlap = true;

if ~isempty(keyName)
    isKeyEvent = strcmp(keyName, {events.type}); %events we are going to build baseline around
end

keyStarts = [events(isKeyEvent).start];
keyStops = [events(isKeyEvent).stop];

% construct base events off of key starts
if ~isempty(startOffsets)
    basePreEvents = eventFromTimes(...
        keyStarts + startOffsets(1), ...
        keyStarts + startOffsets(2), fs);
else
    basePreEvents = initEvents(0);
end

% construct base events off of key stops
if ~isempty(stopOffsets)
    basePostEvents = eventFromTimes(...
        keyStops + stopOffsets(1), ... 
        keyStops + stopOffsets(2), fs);
else
    basePostEvents = initEvents(0);
end

% resolve overlaps (method 1): remove all baseline events that intersect
% with 
% warning: may not work for start-derived and stop-derived bases
baseEvents = sortBy([basePreEvents; basePostEvents], 'start'); %Jenny re-named from baseCandidates
%ok so now have a column of pre-motif and post-motif events
%(many of which probably are during song at this point)
overlapBaseWithAny = findOverlapsBi(baseEvents, events);%Jenny re-named from baseCandidates
%baseEvents = baseCandidates; Jenny commented out
baseEvents(overlapBaseWithAny(:,1)) = []; %delete baseline events that overlap with song

% resolve overlaps (method 2): remove only the portions of baseline that
% overlapped with other events
%{
if ~allowSelfOverlap
    baseEvents = findUniqueAreas(baseEvents, events(isKeyEvent));
end
if ~allowOtherOverlap
    baseEvents = findUniqueAreas(baseEvents, events(~isKeyEvent));
end
%}
% merge base events that overlap
baseEvents = findUnion(baseEvents);
