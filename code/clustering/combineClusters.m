function combineClusters(songStruct, bouts, vocalizations, boutLinked)
    
% get the longest bout
[foo, imax] = max([bouts.stop] - [bouts.start]);

% get clusters of that bout
if isempty(boutLinked{imax})
    error('combineClusters::monosyllabic', 'no bouts contain songs');
end

% separate syllables into separate clusters,
% either clusters or number of possible clusters
distCutoff = 0.7;
clusterIdxs = cluster(boutLinked{imax},'Cutoff', distCutoff,'criterion','Distance');

% use the length of the events as a major hint
inBout = getSubEvents(bouts(imax), vocalizations);
uClusters = unique(clusterIdxs);

%% get exemplars (representative syllables) from each 
avgSylLength = zeros(1,max(uClusters));
exemplar = initEvents;
for ii = 1:numel(uClusters)
    avgSylLength(uClusters(ii)) = ...
        mean([inBout(clusterIdxs == uClusters(ii)).stop ] - ...
             [inBout(clusterIdxs == uClusters(ii)).start]);     
    % pull an exemplar from each cluster
    exemplar(uClusters(ii)) = inBout(find(clusterIdxs == uClusters(ii),1));
    
    fs = 1/songStruct.interval;
    eve = addPrePost(exemplar(uClusters(ii)),...
        processArgs(defaultParams,'preRoll',0.02,'postRoll',0.02));
    clip = getClip(eve, songStruct);
    %playSound(clip, fs, true); pause(1);
end

keyboard
% now go through each of the other bouts, and match each syllable to other
% syllables in the bout

for ii = 1:numel(bouts)
    inBout = getSubEvents(bouts(imax), vocalizations);
    for jj = 1:numel(exemplar)
        simbatch = similarityRegions(songStruct, exemplar(jj), inBout, [])
        keyboard
    end
end
