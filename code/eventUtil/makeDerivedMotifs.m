%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
params = defaultParams;
compositeReport = reportOnData('','',defaultParams,'verbose',true);
if ~iscell(compositeReport)
    compositeReport = {compositeReport};
end
nBirds = numel(compositeReport);
birdNeuronData = cell(nBirds, 1);

if ~exist('MUAinfo','var')
    MUAinfo = getSpikeMUAData;
end
%%
for ii = 1:nBirds
    allBirdSessions = compositeReport{ii};
    
    nSessions = numel(allBirdSessions);
    sessNeuronData = cell(nSessions,1);
    
    for jj = 1:nSessions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % setup
        sessionID = allBirdSessions(jj).sessionID;
        sessManifest = allBirdSessions(jj).manifest;
        sessSpikeFiles = allBirdSessions(jj).spikeFiles;
        if findInManifest(sessManifest, 'derivedMotifs');
            fprintf('derivedMotifs already constructed for session %s.\n', sessionID);
            continue;
        end
        % get the session/bird ID
        birdID = strtok(sessionID, '_');
        birdDataDir = [pwd filesep 'data' filesep birdID filesep];
        
        if ~findInManifest(sessManifest, 'metaStruct'),
            fprintf('Missing metadata, skipping session %s.\n',sessionID);
            continue;
        end
        if ~findInManifest(sessManifest, 'approvedSyllables'),
            fprintf('Missing approvedSyllables, skipping session %s.\n',sessionID);
            continue;
        end
        
        blankedSongStruct = loadFromManifest(sessManifest,'metaStruct');
        fs = 1/blankedSongStruct.interval;
        aS = loadFromManifest(sessManifest,'approvedSyllables');

        gapTime = params.derivedMotifMaxGap; % seconds, empirically agreed upon by JMA & JS, currently 0.3 (2/27/14).
        derivedMotifs = mergeGaps(sortBy(aS,'start'),gapTime);
        
        derivedMotifs = addPrePost(derivedMotifs, defaultParams, 'preroll', 50, 'postroll', 0); % adding premotor in
        [derivedMotifs.type]=deal('BOS');
        
        derivedMotifFile = [birdDataDir 'derivedMotifs-' sessionID '.mat'];        
        fprintf('Wrote derived motifs for session %s to %s.\n', sessionID, derivedMotifFile);
        save(derivedMotifFile,'derivedMotifs');
    end
end

