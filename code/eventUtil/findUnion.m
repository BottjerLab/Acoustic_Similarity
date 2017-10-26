function events = findUnion(events)

% precondition: events are sorted by start?

events = sortBy(events, 'start');
if ~isEvent(events)
    return;
end

starts = [events.idxStart] + 0.5;
stops = [events.idxStop];
N = numel(events); 
parens = [starts stops];
inds = [1:N 1:N]; 
poles = [ones(N,1) -ones(N,1)];

% this gets a little confusing
[sParens, sIdx] = sort(parens);
sInds = inds(sIdx);
sPoles = poles(sIdx);

csPoles = cumsum(sPoles);
trueStarts = find(sPoles == 1 & csPoles == 1);
trueStops = find(sPoles == -1 & csPoles == 0);

% transfer the trueStop ends to the trueStart ends
num2cell([events(sInds(trueStops)).stop]); [events(sInds(trueStarts)).stop] = ans{:};
num2cell([events(sInds(trueStops)).idxStop]); [events(sInds(trueStarts)).idxStop] = ans{:};

% make sure that a type is still assigned if one of the events has no type
if ~isempty(trueStarts) 
    for ii = 1:numel(trueStarts)
        if isempty(events(sInds(trueStarts(ii))).type)
            events(sInds(trueStarts(ii))).type = ...
                events(sInds(trueStops(ii))).type;
        end
    end
end

%% take out the truestops
events = events(sInds(trueStarts));