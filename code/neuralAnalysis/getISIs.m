function [ISIs_in, bfEvent] = getISIs(spikeTrain, events) %[ISIs_in, ISIs_out], JMA added bfEvent
if nargin < 2 || numel(events) < 1
    ISIs_in = diff(spikeTrain);
    %     ISIs_out = [];
    bfEvent = NaN; %JMA added
    return;
else
    events = sortBy(events, 'start');
    starts = [events.start]';
    stops = [events.stop]';
    
    isBorder = [ones(numel(events),1); -ones(numel(events),1); false(numel(spikeTrain),1)];
    augSpikeTrain = [starts; stops; spikeTrain];
    [sortedAST, sidx] = sort(augSpikeTrain);
    sortedIsBorder = isBorder(sidx);
    
    inEvent = cumsum(sortedIsBorder); inEvent(sortedIsBorder == 1) = 0;
    
    inSpikes  = sortedAST;  inSpikes(sortedIsBorder ~= 0) = NaN;  inSpikes(inEvent~=1) = NaN;
    %outSpikes = sortedAST; outSpikes(sortedIsBorder ~= 0) = NaN; outSpikes(inEvent==1) = NaN;
    
    ISIs_in  = diff( inSpikes);  
    ISIs_in(isnan(ISIs_in )) = [];
    %ISIs_out = diff(outSpikes); ISIs_out(isnan(ISIs_out)) = [];
    
    %JMA added this part to do get burst fraction per event
    bfEvent = NaN(numel(events),1);
    for i = 1: numel(events)
        eventSpikInd = starts(i) < inSpikes > stops(i);
        if sum(eventSpikInd) > 1 %if there are at least two spikes during event
            spikeTimeInEvent = inSpikes(eventSpikInd);
            eventISIs = diff(spikeTimeInEvent);
            nBins = 400;
            maxLog = 2;
            minLog = -3;
            binsUsed = logspace(minLog,maxLog,nBins);
            binsUsed = binsUsed';
            logISI = log10(eventISIs);
            H = ndhist(logISI', nBins, minLog, maxLog);
            toBurst = find(binsUsed <= 0.01, 1,'last');
            bfEvent(i)= nansum(H(1:toBurst))/nansum(H);
        end
    end
end