function [evUnique1, evUnique2] = findUniqueAreas(evs1, evs2)
% return the areas unique to evs1 and evs2
% precondition: evs1 and evs2 have no internal overlap of their own;

    % todo: get the fs reliably (dodgy)
    evUnique1 = initEvents(0); evUnique2 = initEvents(0);
    if isempty(evs2), evUnique1 = evs1; return; end;
    if isempty(evs1), evUnique2 = evs2; return; end;
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
    
    un1Starts = flagTimes(runToggle == 1); 
    un1Stops  = flagTimes(find(runToggle == 1) + 1); 
    un2Starts = flagTimes(runToggle == 2);
    un2Stops  = flagTimes(find(runToggle == 2) + 1); 
    evUnique1 = eventFromTimes(un1Starts, un1Stops, fs); 
    evUnique2 = eventFromTimes(un2Starts, un2Stops, fs); 
    
    % now remove all events that are zero-measure, i.e. start and stop at
    % the same place
    
    evUnique1(un1Starts == un1Stops) = [];
    evUnique2(un2Starts == un2Stops) = [];
    
end
