function ret = getWholeSongAsRegion(songStruct)

% instead of constructing with eventFromTimes, just use the common sense
% first/last indices
ret = initEvents(1);
ret.start = songStruct.interval;
ret.stop = ret.start * numel(songStruct.values);
ret.idxStart = 1;
ret.idxStop = numel(songStruct.values);
end