function [parentIdxs, eventsInContext] = findParent(parentEvents, subEvents, params, varargin)
% FUNCTION [PARENTIDXS, EVENTSINCONTEXT] = FINDPARENT(PARENTEVENTS, SUBEVENTS)
% 
% Find the parent of subEvent out of the set of parentEvents,
% returning the index in that event
% also may returns the events adjusted to their parent event 

parentIdxs = NaN(1,numel(subEvents));
if nargout == 2
%    eventsInContext = initEvents(numel(subEvents));
    eventsInContext = initEmptyStructArray(fieldnames(subEvents));
end
if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});

% allow for some soft boundaries in terms of allowing events to contain their children

softMargin = 0.02; % in seconds
starts = [parentEvents.start];
stops  = [parentEvents.stop];
for ii = 1:numel(subEvents)
    currParent = find(subEvents(ii).stop - softMargin < stops & ...
        subEvents(ii).start + softMargin > starts, 1);
    if isempty(currParent) % this event is not contained within any parent
        endStopper = min(stops(stops > subEvents(ii).stop) - subEvents(ii).stop);
        topStopper = min(subEvents(ii).start - starts(starts < subEvents(ii).start));
        
        if params.verbose
            fprintf(['Found event without parent containing event, closest start is %02.f s away, '...
                    'closest end is %0.2f s away\n'], topStopper, endStopper);
        end
        continue;
    end
    parentIdxs(ii) = currParent;
    if nargout == 2
        eventsInContext(ii) = adjustTimeStamps(subEvents(ii),-parentEvents(currParent).start);
    end
end
