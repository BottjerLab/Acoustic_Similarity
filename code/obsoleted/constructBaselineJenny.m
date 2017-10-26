%removing baseline periods that overlap with events---Jenny
function baselineEvents = constructBaselineJenny(events, keyName, startOffsets, stopOffsets, params)
% does automatic exclusion
% startOffsets is 2 x 1 vector in seconds, and is SIGNED (to go backwards
% in time, go negative)
% stopOffsets is similar
% if we only want one of st artOfsets/stopOffsets, blank the other ([])
fs = params.fs;

% allow overlap?
%allowSelfOverlap = false;
%allowOtherOverlap = true;

% dodgy: default behavior for keyName
isKeyEvent = true(numel(events),1);
if ~isempty(keyName)
    isKeyEvent = strcmp(keyName, {events.type});
end

keyStarts = [events(isKeyEvent).start];
keyStops = [events(isKeyEvent).stop];
nKey = sum(isKeyEvent);

otherEvents = events(~isKeyEvent);


isUsablePreEvent  = false(nKey,1);
isUsablePostEvent = false(nKey,1);

for j = 1:length(events)
    isUsablePreEvent  = isUsablePreEvent | (keyStarts + startOffsets(2) < 
end

for i = 1:length(keyStarts); %go through each event of this type
    for j = 1:length(events); %look for overlap with any event
        preEventsToUse = [];
        if (keyStarts(i) + startOffsets(2)) < events(j).start &...
                (keyStarts(i)+ startOffsets(2)) > events(j).stop;
            %keep track of the baseline periods that do not overlap
            preEventsToUse(end+1) = i;            
        end
        postEventsToUse = [];
        if (keyStops(i)+ stopOffsets(2)) < events(j).start &...
                (keyStops(i)+ stopOffsets(2)) > events(j).stop;
            postEventsToUse(end+1) = i;            
        end
    end
    %save events for which baseline is usable
    correctedOwnStarts = preEventsToUse(keyStarts);
    correctedOwnStops = postEventsToUse(keyStops);
    
    baselinePreEvents = initEvents(0);
    if ~isempty(startOffsets)
        baselinePreEvents = eventFromTimes(...
            correctedOwnStarts + startOffsets(1), ...
            correctedOwnStarts + startOffsets(2), fs);
    end
    baselinePostEvents = initEvents(0);
    if ~isempty(stopOffsets)
        baselinePostEvents = eventFromTimes(...
            correctedOwnStops + stopOffsets(1), ...
            correctedOwnStops + stopOffsets(2), fs);
    end
end
%% Option 1:"timing", removing early
%{
for i = 1:length(ownStarts); %go through each event of this type
    for j = 1:length(events.type); %look for overlap with any event
        k = 1;
        if (ownStarts(i) + startOffsets(2)) < events(j).start &...
                (ownStarts(i)+ startOffsets(2)) > events(j).stop;
            %keep track of the baseline periods that do not overlap
            preEventsToUse(k) = i;
            k = k + 1;
        end
        k = 1;
        if (ownStops(i)+ stopOffsets(2)) < events(j).start &...
                (ownStops(i)+ stopOffsets(2)) > events(j).stop;
            postEventsToUse(k) = i;
            k = k + 1;
        end
    end
    %save events for which baseline is usable
    correctedOwnStarts = preEventsToUse(ownStarts);
    correctedOwnStops = postEventsToUse(ownStops);
    
    baselinePreEvents = initEvents(0);
    if ~isempty(startOffsets)
        baselinePreEvents = eventFromTimes(...
            correctedOwnStarts + startOffsets(1), ...
            correctedOwnStarts + startOffsets(2), fs);
    end
    baselinePostEvents = initEvents(0);
    if ~isempty(stopOffsets)
        baselinePostEvents = eventFromTimes(...
            correctedOwnStops + stopOffsets(1), ...
            correctedOwnStops + stopOffsets(2), fs);
    end
end
%}
%would need to run findUnion for only those baselines that have both
%pre and post?

%%Option 2: "length", removing later
%let everything go as usual through findUniqueAreas, then remove...
%   pre-periods that are shorter in time than 
%   abs(startOffsets(1)) - abs(startOffsets(2))...
%   or post-periods that are shorter than stopOffsets(2) - stopOffsets(1);
%Are the lengths of baselines preserved in findUniqueAreas?
    
    
    