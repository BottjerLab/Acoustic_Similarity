function [overlapIndices, overlapAreas] = findOverlaps(events)
% overlapIndices a N x 2 matrix of pairs of indices to the event array
% these pairs of indices indicate events with overlap

% precondition: events are sorted by start?

fs = round(events(end).idxStart/events(end).start);
overlapIndices = zeros(0,2); % sparse lists of overlap
overlapAreas = initEvents;
if ~isEvent(events)
    return;
end

% 0.5 offset to avoid start and stop collision
starts = [events.idxStart] + 0.5;

stops = [events.idxStop];
N = numel(events);
parens = [starts stops];
inds = [1:N 1:N]; 
poles = [ones(N,1) -ones(N,1)];

% this gets a little confusing
[sParens, sIdx] = sort(parens);
sInds = inds(sIdx);
sPoles = poles(sIdx); %if it's a start or stop
[~,ssIdx] = sort(sIdx);
% in which order is the Nth start or stop among the series of all
% starts/stops?
rIdxStart = ssIdx(1:N); rIdxStop = ssIdx((N+1):end);

if any(rIdxStop - rIdxStart ~= 1)
    % we have some closures
    % not the fastest way but simple enough
    
    % step 1: identify the overlaps
    for ii = 1:N
        if (rIdxStop(ii) - rIdxStart(ii) ~= 1)
            intersects = sInds(rIdxStart(ii)+1:rIdxStop(ii)-1)';
            intersects(ii > intersects) = [];
            intersects = unique(intersects); % can this be quicker?
            nI = numel(intersects);
            overlapIndices = [overlapIndices; [ii*ones(nI,1) intersects]];
        end
    end
    
    if nargout == 2
        overlapStarts = max([events(overlapIndices(:,1)).start], ...
            [events(overlapIndices(:,2)).start]);
        overlapStops  = min([events(overlapIndices(:,1)).stop], ...
            [events(overlapIndices(:,2)).stop]);
        overlapAreas = eventFromTimes(overlapStarts, overlapStops, fs);
    end
else
    return;
end