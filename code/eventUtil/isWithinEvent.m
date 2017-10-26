function ans = isWithinEvent(subEvs, parentEv)
    [subEvs.start] >= parentEv.start & [subEvs.stop] <= parentEv.stop;
end