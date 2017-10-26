function childEvents = inheritLabels(childEvents, parentEvents)
parentage = findParent(parentEvents, childEvents);
for ii = 1:numel(childEvents)
    if ~isnan(parentage(ii))
        childEvents(ii).type = parentEvents(parentage(ii)).type;
    end
end