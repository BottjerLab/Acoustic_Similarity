function events = basicEvents(events)
basicFields = {'type','start','stop','idxStart','idxStop'};
flds = fieldnames(events);
for ii = 1:numel(flds)
    if ~any(strcmp(flds{ii},basicFields))
        events = rmfield(events, flds{ii});
    end
end