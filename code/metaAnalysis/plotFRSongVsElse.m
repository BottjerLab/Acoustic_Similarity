function neuronCollection = plotFRSongVsElse(session)

% goal: allow this to have the same argument list as reportOnData, so that
% we can always subselect by session according to some parameter (age, etc)

% TODO: REWRITE THIS FUNCTION TO FIT ITS NAME

%session = 'Lb277_3_20_1'; % TODO: batch preparation

% might be off here, test here
if nargin == 1, 
    if ~isempty(strfind(session,'_'))
        dataKeys = reportOnData(strtok(session, '_'), session);
    else
        dataKeys = reportOnData(session);
    end
else dataKeys = reportOnData;
end

neuronCollection = [];
for ii = 1:numel(dataKeys)
    birdAllRecords = dataKeys{ii};
    for jj = 1:numel(birdAllRecords)        
        manifest = birdAllRecords(jj).manifest;
        
        [~,keyEntry]=min(cellfun('length',{manifest.originalFile}));
        keyFile = manifest(keyEntry).originalFile;
        [~,thisSession] = fileparts(keyFile);
        % get motifs
        
        motifs = loadFromManifest(manifest, 'manualMotifs');
        if isempty(motifs), continue; end;
        
        spikes = loadSpikeData(birdAllRecords(jj).spikeFiles);
        if isempty(spikes), continue; end;
        
        songStruct = loadSpikeAudioStruct(keyFile, 'novalues');
        fprintf('Analyzing from manifest %s\n', keyFile);
        
        % assume BOS if unlabeled (yet)
        if all(cellfun('isempty',{motifs.type}))
            [motifs.type] = deal('BOS');
        end
        neuronData = spikingByCategory(songStruct, motifs, spikes);
        [neuronData.session] = deal(thisSession);
        neuronCollection = [neuronCollection; neuronData];
    end 
end
%{
    birdRecords = reportOnData(birdID, session); 
    manifest = birdRecords.manifest;
    % load just the syllables 
    
        
    % get spikes
    
    neuronData = spikingByCategory(songStruct, motifs, spikes);
    
    % plotting

    %title(['Syllable distribution, session ' session], 'Interpreter','none');
%}