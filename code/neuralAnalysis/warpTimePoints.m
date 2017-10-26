function eventSpikeTimes = warpTimePoints(spikeTimes, refEvents, events, prePostParams)

% function eventSpikeTimes(SPIKETIMES, REFEVENTS, EVENTS)
%
% SPIKETIMES is a session recording of spike times
% REFEVENTS and EVENTS is a collection of events (parent and child)
% parent events must have type 'parent'
% or a parent field that is true
% 
% REFEVENTS is the collection of events (1 parent, 0+ child) 
% that serve as the reference to which other events will be conformed
% child events should be positioned w.r.t. the parent event, and 
% completely within the parent events
% note that REFEVENTS can refer to a different file - spiking data is not taken into account for the reference. 

if isfield(refEvents, 'parent')
    isParent = [refEvents.parent];
else
    isParent = strcmpi('parent',{refEvents.type});
end
refParent   = refEvents( isParent);
refChildren = refEvents(~isParent);

if isfield(refEvents, 'parent')
    isParent = [events.parent];
else
    isParent = strcmpi('parent',{events.type});
end
parents  = events( isParent);
children = events(~isParent);

% ensure reference children have unique labels
assert(numel(refChildren) == numel(unique({refChildren.type})))

% sort the children in order
refChildren = sortBy(refChildren, 'start');
refTypes = {refChildren.type};

% ensure all ref children are within reference parent
assert(all(isWithinEvent(refChildren,refParent)))

% get the parent-child assignment
parentIdx = findParent(parents, children);
nParents = numel(parents);

% get the unwarped spikes in the parent events
[~,contextSpikes] = countSpikes(parentEvents, spikeTimes);
warpedSpikes = cell(size(contextSpikes));
for ii = 1:nParents
     itsChildren = refChildren(parentIdx == ii);
     % find the matches 
     refCIdx = NaN(1,numel(refChildren));
     for jj = 1:numel(refChildren)
         foo = find(strcmp({itsChildren.type},refTypes{jj}),1);
         if ~isempty(foo), refCIdx = foo; end;
     end
     
     refMatch = refChildren(refCIdx(~isnan(refCIdx)));
     match    = children(~isnan(refCIdx));
     refFids = [  refParent.start; refMatch.start;   refParent.stop; refMatch.stop];
     fids =    [parents(ii).start;    match.start; parents(ii).stop;    match.stop];
     
     % interpolate and subtract the onset
     warpedSpikes{ii} = interpLinearPoints(refFids, contextSpikes{ii},fids) - refParent.start;
end
    
% NB to john: start writing here

% set up parent-child relationships
% find which subevents belong to which events
% note: if the events and subevents are in one to one correspondence, (i.e.
% one is the pre/post rolled version of another, then the parentage is
% probably just one-to-one
if numel(events) == numel(subEvents)
    parentage = 1:numel(subEvents);
    for ii = 1:numel(subEvents)
        eventsInContext(ii) = adjustTimeStamps(subEvents(ii), ...
            -events(ii).start);
    end
else
    if nargin < 5
        [parentage, eventsInContext] = findParent(events,subEvents);
    else %nargin == 5
        parentage = findParent(events,subEvents);
        events = addPrePost(events, prePostParams);
        % take out orphan events :[
        subEvents(isnan(parentage)) = [];
        parentage(isnan(parentage)) = [];
        for ii = 1:numel(subEvents)
            eventsInContext(ii) = adjustTimeStamps(subEvents(ii), ...
                    -events(parentage(ii)).start);
        end
    end

end


function yy = interpLinearPoints(xx,y1,x1)
% interpolate linearly between each pair of points & interpolate, without
% warping before and after the control points
assert(numel(x1) == numel(y1))
if numel(x1) == 0, yy = xx; return; end;

[x1, idxs] = sort(x1);
y1 = y1(idxs);
intervalPtr = 1;

yy = zeros(size(xx));
for kk = 1:numel(xx)
    intervalPtr = find(xx(kk) <= x1,1);
    if isempty(intervalPtr)
        yy(kk) = xx(kk) - x1(end) + y1(end); % no warping after last control point
    elseif intervalPtr == 1
        yy(kk) = xx(kk) - x1(1) + y1(1); % no warping before first control point
    else
        yy(kk) = interp1(x1(intervalPtr-1:intervalPtr),...
            y1(intervalPtr-1:intervalPtr),xx(kk));
    end
end
end
