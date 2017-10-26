function regions = codeEvents(intervals, regions, code, spec, params) 
%CODEEVENTS create events from time points found on a spectrum
%  regions = codeEvents(intervals, regions, code, spec, params) 
%  code song regions
% currently assume everything is correctly sized
% intervals is Nx2, regions is Nx1, code is scalar
% spec must contain the times field
% params must contain fs field

% should be OBSOLETED
    nNewEvents = size(intervals,1);
    for istat = 1:nNewEvents;
        
        %TODO: preallocate additional regions, but also check for overlap
        newEvent.type = code; %song
        newEvent.start = spec.times(intervals(istat,1)); %seconds 
        newEvent.stop  = spec.times(intervals(istat,2)); %seconds
        newEvent.idxStart = floor(newEvent.start * params.fs); % sample indices
        newEvent.idxStop  = floor(newEvent.stop  * params.fs); % sample indices
        regions(end+1) = newEvent; 
    end
end