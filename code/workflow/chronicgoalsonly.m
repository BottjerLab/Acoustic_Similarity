%% Things to do for our Science/Nature paper

% Determine if bird's song is subsong or plastic song per day
    %% plot syllable length distributions
    session = 'Lb277_3_20_1'; % TODO: batch preparation
    birdID = strtok(session, '_');
    birdRecords = reportOnData(birdID, session); manifest = birdRecords.manifest;
    % load just the syllables 
    
    % TODO: conform data to standard set in RECORDINGDATA.txt
    varName = 'syllables';
    correctRecord = manifest(strcmp(varName,{manifest.name}));
    
    lysis = load(correctRecord.originalFile, correctRecord.name);
    
    syllableLengths = [lysis.(varName).stop] - [lysis.(varName).start];
    
    % plotting
    binWidth = 0.005;
    syllBins = 0:binWidth:0.15; binCenters = syllBins + binWidth/2;
    syllPDF = hist(syllableLengths,syllBins) / numel(syllableLengths);
    bar(binCenters, syllPDF, 1)
    title(['Syllable distribution, session ' session], 'Interpreter','none');
    %%
    % can compare to distributions of syllables we can dig up from...
    %   existing songs we have from subsong and plastic song birds in...
    %   our colony (not from this experiment)

    %% mini script to rename variables in a mat file, using load/save as struct
    srcVar = 'juvenileSylls';
    dstVar = 'syllables';
    foo = load(correctRecord(1).originalFile);
    foo.(dstVar) = foo.(srcVar); 
    foo = rmfield(foo,srcVar);
    save(correctRecord(1).originalFile,'-struct','foo');
    
% PER NEURON

% For older birds in plastic song phase:
    % parse and label syllables
    % time-warp motifs
    % plot rasters/instan. FR
    % plot ISI distribution for singing vs not
    % determine FR differences for singing vs not
    % determine FR per syllable (during how many syllables is it firing...
        % more than baseline
    % determine is ISI differs for different syllables
    %% run cross-correlation on spike trains for each motif rendition
    % we need the motifs and the spike trains
    
    %%
% For subsong birds:
    % parse syllables but do not attempt to label
    % plot ISI distribution for singing vs not
    % determine FR differences for singing vs nota
    % have machine try to cluster syllables
    % plot FR vs the machine's clusters to see if neuron fires to...
        %specific kinds of syllables
        
% For all birds:
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
    