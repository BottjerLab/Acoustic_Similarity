function boo = isEvent(ev)
if isempty(ev), boo = true; return; end;
if ~isstruct(ev), boo = false; return; end;

reqFlds = {'start','stop','type','idxStart','idxStart'};
if ~all(isfield(ev, reqFlds)), 
    boo = false; return; 
end;

% consistent sampling rate
%if max(fabs([ev.idxStart]./[ev.start] - [ev.idxStop]./[ev.stop])) < 0.0001,
%    boo = false; warning('unexpectedEventFS', ...
%        'sampling rate is not consistent');
%    return;
%end;
boo = true;
end
    