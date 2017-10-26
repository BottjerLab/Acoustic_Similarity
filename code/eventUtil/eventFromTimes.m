function ev = eventFromTimes(starts, stops, fs)
nR = numel(starts);
ev = initEvents(nR);

assert(numel(starts) == numel(stops));
%#ok<*NOANS> 
%#ok<*USENS>
if isempty(starts), return; end;

num2cell(starts); [ev.start] = ans{:};  
num2cell(stops); [ev.stop] = ans{:}; 
num2cell(floor(starts * fs));    [ev.idxStart] = ans{:}; 
num2cell(floor(stops  * fs)); [ev.idxStop] = ans{:}; 
end