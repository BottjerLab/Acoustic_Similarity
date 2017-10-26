% unvalueAudiofile.m

% get the manifests
compositeReport = reportOnData;
% get rid of any files that don't have any spikes
%{
isNeuronBird = false(1,numel(comp));
for ii = 1:numel(comp)
    comp{ii}(cellfun('isempty',{comp{ii}.spikeFiles}))=[]; 
    isNeuronBird(ii) = ~isempty(comp{ii});
end
comp = comp(isNeuronBird);
%}
nBirds = numel(compositeReport);

%%
% put the useful info in 'Session Records' worksheet 
% Columns = 
% Bird ID, Session ID, Date of Recording, Age, #bouts, #motifs, #syllables, # core neurons, # shell neurons, Singing Time, Plastic/Subsong (Y/N)

for ii = 1:numel(compositeReport)
    allBirdSessions = compositeReport{ii};
    
    for jj = 1:numel(allBirdSessions)
        sessManifest = allBirdSessions(jj).manifest;
        
        sessionID = allBirdSessions(jj).sessionID;
        birdID = strtok(sessionID, '_');
        
        [sS, sFile] = loadFromManifest(sessManifest,'songStruct');
        if isempty(sS),
            fprintf('missing audio from session %s...\n', allBirdSessions(jj).sessionID);
            continue;
        end
        metaStruct = rmfield(sS,'values');
        [filePath, fileStem, fileExt] = fileparts(sFile);
        metaFile = [filePath filesep 'meta-' fileStem fileExt];
        save(metaFile, 'metaStruct')
        fprintf('Writing meta for session %s to file %s...\n', allBirdSessions(jj).sessionID, metaFile);       
    end
end
