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
        
        if numel(sessManifest) == 1 % we only have the audio data from spike,
            continue;
        end
        if isempty(allBirdSessions(jj).spikeFiles)
            continue;
        end
        
        sessionID = allBirdSessions(jj).sessionID;
        birdID = strtok(sessionID, '_');
        
        mM = loadFromManifest(sessManifest,'manualMotifs');
        
        motifNewTypes = {''};        
        % translate NaN values to strings
        if ~isempty(mM) 
            mTypes = {mM.type};
            mTypes(cellfun('isempty',mTypes))=[]; % get rid of empties
            for kk = 1:numel(mTypes)
                if isnumeric(mTypes{kk}), mTypes{kk} = num2str(mTypes{kk}); end;
            end
            motifTypes = unique(mTypes);

            if ~isempty(motifTypes)
                motifNewTypes = cell(1,2*numel(motifTypes));
                [motifNewTypes{1:2:end}] = motifTypes{:};
                [motifNewTypes{2:2:end-1}] = deal(', ');
                motifNewTypes{end} = '.';
            end
        end
        if ~isempty(motifNewTypes)
            fprintf('Session %s/%s has motif types %s...\n', birdID, sessionID, strcat(motifNewTypes{:}));
        else
            fprintf('Session %s/%s has empty motif types...\n');
        end
    end
end
