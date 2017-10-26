function plotPSTHWithBaseline(events, spikeTimes, params)
% works with only single events (motifs/individual syllables), 
% NOT with multi-syllable song, etc.
superEvents = addPrePost(events, params);
plotPSTH(superEvents, spikeTimes, events);
