function evCommon = findIntersection(evs1, evs2)

if isempty(evs1), evCommon = evs1; return; end;
if isempty(evs2), evCommon = evs2; return; end;
% precondition: evs1 and evs2 have no internal overlap of their own
% note: this is NOT "find all evs1 within evs2".
% this is "find all times t s.t. t exists in both evs1 and evs2".  


    % todo: get the fs reliably
    fs = evs2(1).idxStart/evs2(1).start;
    
    % sort the events independently
    evs1 = sortBy(basicEvents(evs1),'start');  
    evs2 = sortBy(basicEvents(evs2),'start');
    
    N1 = numel(evs1);
    N2 = numel(evs2);

    % flags on a track method: 
    flagTimes = [[evs1.start] [evs1.stop] [evs2.start] [evs2.stop]];
    toggles = [ones(1,N1) -ones(1,N1) 2*ones(1,N2) -2*ones(1,N2)];
    [flagTimes, sortPerm] = sort(flagTimes);
    toggles = toggles(sortPerm);

    % now run along the flags, calculating the sum
    % if the events are proper, sum should only be from 0-3.  
    % 0-> is nothing,
    % 1-> is unique 1,
    % 2-> is unique 2
    % 3-> is both
    
    % the last point should always be zero
    runToggle = cumsum(toggles);
    assert(runToggle(end) == 0);
    
    comStarts = flagTimes(runToggle == 3); 
    comStops  = flagTimes(find(runToggle == 3) + 1); 
    evCommon = eventFromTimes(comStarts, comStops, fs); 
     
    % now remove all events that are zero-measure, i.e. start and stop at
    % the same place
    
    evCommon(comStarts == comStops) = [];
end