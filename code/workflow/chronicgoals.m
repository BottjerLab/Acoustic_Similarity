%% Things to do for our Science/Nature paper

% Determine if bird's song is subsong or plastic song per day
    %% plot syllable length distributions
    session = 'Lb277_3_20_1'; % TODO: batch preparation
    plotDurationDistr(session); % FIXME: rewrote to take variable, not session string
    %%
    % can compare to distributions of syllables we can dig up from...
    %   existing songs we have from subsong and plastic song birds in...
    %   our colony (not from this experiment)   
    
% PER NEURON

%NB: subsong distribution of syllables is unimodal,
% plastic song is multimodal

% BUT this does not guarantee sequential arrangement of syllables
% For older birds in plastic song phase:
    % parse and label syllables
    % time-warp motifs
    % plot rasters/instan. FR
    % plot ISI distribution for singing vs not
    % determine FR differences for singing vs not
    % -> this can be done with spikingByCategory
    % determine FR per syllable (during how many syllables is it firing...
        % more than baseline
    % -> this can be done with spikingByCategory    
    % determine is ISI differs for different syllables
    %% run cross-correlation on spike trains for each motif rendition
    % we need the motifs and the spike trains
    
    session = 'Lb277_3_20_1'; % TODO: batch preparation
    birdID = strtok(session, '_');
    birdRecords = reportOnData(birdID, session); 
    manifest = birdRecords.manifest;
    
    varName = 'manualMotifs';
    correctRecord = manifest(strcmp(varName,{manifest.name}));
    lysis = load(correctRecord.originalFile, correctRecord.name);
    spikes = loadSpikeData(birdRecords.spikeFiles);
    
    [~, times] = countSpikes(lysis.(varName), spikes{1}, 'onset');
    %lengt
    binSize = 0.02; % in seconds
    
    %%
% For subsong birds:
    % parse syllables but do not attempt to label
    % determine FR differences for singing vs not
    session = 'Y231_12_11_030';
    plotFRSongVsElse(session);
    %
    %% plot ISI distribution for singing vs not
    session = 'Y231_12_11_030';
    plotISISong(session);
    %% have machine try to cluster syllables
    % plot FR vs the machine's clusters to see if neuron fires to...
        %specific kinds of syllables
    session = 'Y231_12_11_030'; 
    reportOnData('Y231', session);
    
%% For all birds:
    % Use parsed syllables but with no labels
    % Find similarity between a juvenile syllable and each tutor syllable
    % Plot FR vs "best" similarity score
    
    % Plot ISI distribution for top 25% similarity vs bottom 25% similarity
    % Determine FR for top 25% similarity scores vs bottom 25% similarity
    % Plot correlation between FR and similarity
    
    % Determine different neuron types

    % GROUP DATA (across neurons of the same category, e.g. excited separate...
    % from inhibited, core vs shell, across ages)
    % Plot correlation between FR and similarity
    % Plot FR for top 25% similarity vs bottom 25% 
    % Controls to be determined for similarity score
    