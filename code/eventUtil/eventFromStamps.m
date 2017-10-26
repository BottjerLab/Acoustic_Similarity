function ev = eventFromStamps(idxStarts, idxStops, fs)
nR = numel(idxStarts);
ev = initEvents(nR);

assert(numel(idxStarts) == numel(idxStops));
%#ok<*NOANS> 
%#ok<*USENS>
if isempty(idxStarts), return; end;

num2cell(idxStarts); [ev.idxStart] = ans{:}; 
num2cell(idxStops); [ev.idxStop] = ans{:}; 
num2cell(idxStarts / fs); [ev.start] = ans{:};  
num2cell(idxStops / fs); [ev.stop] = ans{:}; 
end