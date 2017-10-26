function [sortedEvents, sortIdx] = sortBy(events, field, mode)
%SORTBY Sort event structure by a given field
%
% function [sortedEvents, sortIdx] = sortBy(events, field) sorts a list of 
% events by the given named (numeric) field.
% The mode parameter can be 'ascend' or 'descend'
if nargin < 3
    mode = 'ascend';
end
data = [events.(field)];

if isnumeric(data), % numeric
    % the 2 might get us in trouble
    [foo, sortIdx] = sort(data, 2, mode);
    sortedEvents = events(sortIdx);    
else    
    data = {events.(field)};
    if ~all(cellfun('isclass', data, 'char')) % must be character
        error('eventUtil:sortBy', ...
                'field %s not pure numeric or pure string in events', field);
    else
        [foo, sortIdx] = sort(data);
        sortedEvents = events(sortIdx);    
    end
end
