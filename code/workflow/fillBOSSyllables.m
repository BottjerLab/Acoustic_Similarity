function fillBOSSyllables(birdID)
% This script takes all the BOS syllables from a bird and parses the
% syllables
if nargin < 1
    birdID = 'Db113'; % bad noise problems, maybe fix in audacity
end
sessionRecords = reportOnData(birdID);

nSessions = numel(sessionRecords);
for ii = nSessions:nSessions
    iRecord = sessionRecords(ii);
    motifs = loadFromManifest(iRecord.manifest,'manualMotifs');
    if isempty(motifs), continue; end;
    
    bosMotifs = motifs(strcmp('BOS',{motifs.type}));
    if isempty(bosMotifs), continue; end;
    
    fprintf('Processing %s...\n',iRecord.sessionID);
    
    noiseMask = loadFromManifest(iRecord.manifest, 'noiseMask');
    [songStruct, matFile] = loadFromManifest(iRecord.manifest, 'songStruct');
    params = processArgs(defaultParams, [birdID '-BOS'], 'fs',1/songStruct.interval, 'preroll',15,'postroll',15);
    bosSyllables = ...
        parseRegionsIntoSyllables(songStruct, bosMotifs, params, ...
        'noiseFilter', noiseMask,'nps.reduction',-20, ...
        'minCenterFreq', 800,...
        'plot', false, 'pause', false, 'verbose', true);   
    
    [bosSyllables.type] = deal('BOSsyll');
    [matpath, matFile] = fileparts(matFile); 
    save([matpath filesep prependForSave('BOSsyllables-',matFile)], 'bosSyllables');
end