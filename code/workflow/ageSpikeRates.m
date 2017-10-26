[,spikeSummary, IDs] = compileStats;

recordMeta = readDataFile('RECORDINGDATA.txt');
IDs

nNeurons = cellfun(@(x) numel(x.stats), spikeSummary);