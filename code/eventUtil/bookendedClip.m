function newClip = bookendedClip(clipSet)

% instead of inferring the fs, just inherit from the event
newClip = initEvents(1);
newClip.start = nanmin([clipSet.start]);
newClip.idxStart = nanmin([clipSet.idxStart]);
newClip.stop = nanmax([clipSet.stop]);
newClip.idxStop = nanmax([clipSet.idxStop]);
end