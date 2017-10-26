function ret = splitIntoOverlap(event, N, overlapTime)
% SPLITINTOOVERLAP splits an event structure into N smaller events
% 
%   function ret = splitIntoOverlap(event, N) segments an
%   event into N smaller events.
%
%   function ret = splitIntoOverlap(event, N, overlapTime) segments an
%   event into N smaller events, where successive events overlap by the 
%   specified amount of time (in s).


if nargin < 3
    overlapTime = 0;
end

if numel(event) > 1
    %split the events proportionally
    evLengths = [event.stop] - [event.start];
    miniN = max(round(evLengths * N / sum(evLengths)),1);
    
    % correct rounding errors
    [~,indMax] = max(miniN);
    miniN(indMax) = miniN(indMax) + (N - sum(miniN));
    
    ret = initEvents(N);
    intervals = cumsum([1 miniN]);
    for jj = 1:numel(event)
        ret(intervals(jj):intervals(jj+1)-1) = ...
            splitIntoOverlap(event(jj),miniN(jj),overlapTime);
    end
    return
end

fs = event.idxStop / event.stop;
totalLen = event.stop - event.start; % in seconds

% window size
timeLen = (totalLen - overlapTime) / N + overlapTime;
% time between window onsets
timeStep = timeLen - overlapTime;

timeStarts = (0:N-1) * timeStep + event.start;
timeStops = timeStarts + timeLen;
ret = eventFromTimes(timeStarts, timeStops, fs);

end