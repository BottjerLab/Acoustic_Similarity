function [evs, idxs] = getSubEvents(mainEvent, subEvents)
%GETSUBEVENTS get all events completely contained within a parent event
% 
%   innerevents = getSubEvents(PARENTEVENT, SUBEVENTS) returns the events
%   completely contained within the parent event.
%   
%   [innerevents, idxs] = getSubEvents(PARENTEVENT, SUBEVENTS) returns the events
%   completely contained within the parent event, and their indexes within 
%   the subEvent structure.

if isempty(subEvents)
    evs = initEvents(0); idxs = [];
    return;
end

assert(numel(mainEvent) == 1)
idxs = find([subEvents.start] >= mainEvent.start & ...
                [subEvents.stop ] <= mainEvent.stop); 
evs = subEvents(idxs); % may look for tolerances