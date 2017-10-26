function events = initEvents(N, exampleEvent)
% function EVENTS = INITEVENTS(N)
%
% initialize sparse coding of status, using events, which will be a struct
% array
%
% all events have the following fields:
% type: a label
% start: the start of the clip (in seconds)
% stop: the end of the clip (in seconds)
% idxStart: the start of the clip (in terms of the index in the waveform)
% idxStop: the end of the clip (in terms of the index in the waveform)
% 'start' and 'stop' should be proportional to their index counterparts, 
% and the ratio should be 1/fs, where fs is the sampling rates.

if nargin == 0
    N = 0;
end
if nargin < 2
    fields = {'type','start','stop','idxStart','idxStop'};
elseif isstruct(exampleEvent)
    fields = fieldnames(exampleEvent);
else
    fields = exampleEvent;
end

events = initEmptyStructArray(fields);
if nargin == 1 && N>0
    events(N) = initEvent;
    for ii = 1:numel(fields)
        events(N).(fields{ii}) = [];
    end
end

function events = initEvent
% flag: isn't that helpful, should be replaced/obsoleted
% initialize sparse coding of status, using events, which will be a struct
% array
innerFields = {'type','start','stop','idxStart','idxStop'};

% fun with cell arrays
foo = cell(2,numel(innerFields));
[foo{1,:}] = innerFields{:};
[foo{2,:}] = deal(NaN);
events = struct(foo{:});