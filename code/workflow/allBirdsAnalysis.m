% allBirdsAnalysis.m

profile off
% big bird style (fixed effects)
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
progressbar('Hi Jenny! I''m working hard!');
for ii = 1:nBirds
    sessionsForThisBird = compositeReport{ii};
    
    nSessions = numel(sessionsForThisBird);
    sessNeuronData = cell(nSessions,1);
    for jj = 1:nSessions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % setup
        tic
        sessManifest = sessionsForThisBird(jj).manifest;
        sessSpikeFiles = sessionsForThisBird(jj).spikeFiles;

        
        % get the session/bird ID
        sessionID = sessionsForThisBird(jj).sessionID;
        birdID = strtok(sessionID, '_');
        birdDataDir = [pwd filesep 'data' filesep birdID filesep];
        
        if numel(sessManifest) == 1 % we only have the audio data from spike,
            fprintf('Missing parsing data, skipping session %s.\n',sessionID);
            continue;
        end
        
        if isempty(sessionsForThisBird(jj).spikeFiles)
            fprintf('Missing stereotrode data, skipping session %s.\n',sessionID);
            continue;
        end
        
        % get metadata: age and plasticness for this session
        [thisAge, isThisPlastic] = getAgeOfSession(sessionID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load motifs (NB: should be in standard event struct, column format)
        mM = loadFromManifest(sessManifest,'derivedMotifs'); %includes 50ms premotor, from makeDerivedMotifs
        if isempty(mM) || ~any(strcmp('BOS', {mM.type}))
            fprintf('Missing BOS motifs, skipping session %s.\n',sessionID);
            continue;
        end
        
        aS = loadFromManifest(sessManifest,'approvedSyllables');
        if isempty(aS)
            warning('allBirdsAnalysis:missingApproved', 'Missing approved syllables, deferring to auto styllables in session %s.\n',sessionID);
            aS = loadFromManifest(sessManifest,'syllables');
            if isempty(aS)
                warning('allBirdsAnalysis:noSyllables','Missing syllables, continuing...');
                continue;
            end
        end
        
        % get first and last syllables within a motif
        
        [fAS, fFAS] = findFirstborn(mM, aS);
        [lAS, lLAS] = findLastborn (mM, aS);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load spikes
        [spikes,nNeuronsPerFile, isMUAUnit] = loadSpikeData(sessSpikeFiles, MUAinfo);
        
        % get metadata - core/shellness and MUA/SUA-ness of each neuron
        isCoreUnit = false(numel(spikes),1);
        
        cumIndex = [0 cumsum(nNeuronsPerFile)];
        for kk = 1:numel(sessSpikeFiles)
            unitsFromFile = (cumIndex(kk)+1):cumIndex(kk+1);
            
            % core/shell ness depends on the file folder
            isCoreUnit(unitsFromFile) = ~isempty(strfind(sessSpikeFiles{kk}, 'core'));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load the songstruct with only metadata, 'metaStruct' -
        blankedSongStruct = loadFromManifest(sessManifest,'metaStruct');
        
        if isempty(blankedSongStruct)
            fprintf('Missing song metadata, skipping session %s.\n',sessionID);
            continue;
        end
        
        % set up bird-specific parameters
        params = processArgs(defaultParams, 'fs', 1/blankedSongStruct.interval);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS: run for each section AND extract the
        % necessary summary statistics
        fprintf('Analyzing session %s/%s...', birdID, sessionID);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pre-loadable part
        % Compare FR during singing vs local baseline - constructBaseline
        % does automatic exclusions        
        baselinePeriodFile = [birdDataDir 'baselinePeriods-' sessionID '.mat'];
        if exist(baselinePeriodFile,'file')
            load(baselinePeriodFile);
        else % if this is the first time we look for silent baseline
            %just run it
            baselinePeriods = constructBaseline(mM, 'BOS', [-4, -2], [2, 4], params);
            [baselinePeriods.type] = deal('localBase'); %we might not use this anymore
            songStruct = loadFromManifest(sessManifest,'songStruct');
            baselinePeriodsS = constructBaselineSilence(songStruct, 2, aS, params);
            [baselinePeriodsS.type] = deal('silence');
            
            save(baselinePeriodFile,'baselinePeriods','baselinePeriodsS');
        end
        %%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % get response strengths comparing events to baselines        
        tmpNeuronData = getRS([mM; baselinePeriodsS], spikes, 'BOS', 'silence', ...
            't_p', 'meanRS', 'meanRSI', 'standRS', 'cvarFR', params, 'verbose',false); %compareTwoCatFiringRates
        
        preTime = (20:10:50); % pre-onset/offset times I might want in msec
        for pt = 1:length(preTime)
            preData = savingPreRS(baselinePeriodsS, spikes, params, aS, fAS, lAS, preTime(pt), tmpNeuronData);
            tmpNeuronData = preData;
        end
        neuronData = tmpNeuronData;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get event-related firing rate to syllables, onsets and offsets
        ERAOnsetWindow = [-0.2 0.2]; %s, the window in which to collect ERA hists
        onsetBinSize = 2e-3; % s
        %         binNumYouWant = find(ERAOnsetWindow(1):onsetBinSize:ERAOnsetWindow(2) >= -0.050,...
%             1);
%         binNumYouAlsoWant = find(ERAOnsetWindow(1):onsetBinSize:ERAOnsetWindow(2) >= 0,...
%             1);
        
        tempVar = eraFind(ERAOnsetWindow ,aS,'onset', onsetBinSize, spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).ERAOnsethist = tempVar(i,:)';
%             tempTemp = tempVar(i,:);
%             if mean(tempTemp(binNumYouWant:binNumYouAlsoWant)) > 0 %+ std(tempTemp(binNumYouWant:binNumYouAlsoWant)) %JMA added so wasn't set threshold
%                 neuronData(i).isERAPrePos = 1;
%             else
%                 neuronData(i).isERAPrePos = 0;
%             end
        end
        %PSTHs for onsets -working on it
                tempVar = eraFindFR(ERAOnsetWindow ,aS,'onset', onsetBinSize, spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).ERAOnsetFR = tempVar(i,:)';
        end
                       
        % Get event-related firing rate to syllable offsets
        tempVar = eraFind(ERAOnsetWindow,aS,'offset', onsetBinSize,spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).ERAOffsethist = tempVar(i,:)';
%             tempTemp = tempVar(i,:);
%             if mean(tempTemp(binNumYouWant:binNumYouAlsoWant)) < 0 %- std(tempTemp(binNumYouWant:binNumYouAlsoWant)) %JMA added so wasn't set threshold
%                 neuronData(i).isERAPostNeg = 1;
%             else
%                 neuronData(i).isERAPostNeg = 0;
%             end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get event-related firing rate to syllable onsets for first syllables in motif
        tempVar = eraFind(ERAOnsetWindow,fAS,'onset',onsetBinSize,spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).firstERAOnsethist = tempVar(i,:)';
%             if mean(tempTemp(binNumYouWant:binNumYouAlsoWant)) > 0 %+ std(tempTemp(binNumYouWant:binNumYouAlsoWant)) %JMA added so wasn't set threshold
%                 neuronData(i).isfERAPrePos = 1;
%             else
%                 neuronData(i).isfERAPrePos = 0;
%             end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get event-related firing rate to syllable onsets for last syllables in motif
        tempVar = eraFind(ERAOnsetWindow,lAS,'offset',onsetBinSize,spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).lastERAOffsethist = tempVar(i,:)';
%             if mean(tempTemp(binNumYouWant:binNumYouAlsoWant)) < 0 %- std(tempTemp(binNumYouWant:binNumYouAlsoWant)) %JMA added so wasn't set threshold
%                 neuronData(i).islERAPostNeg = 1;
%             else
%                 neuronData(i).islERAPostNeg = 0;
%             end
        end
        
        % get the normalized log ISI distributions for singing and baseline, Jenny wrote
        % this, ekkkkk!---Note I've embedded this in getRS now at least for singing, so can get
        % this info from there if I feel like cleaning this up
        nBins = 400;
        maxLog = 2;
        minLog = -3;
        binsUsed = logspace(minLog,maxLog,nBins); %time axis
        binsUsed = binsUsed';
       
        for kk = 1:numel(spikes)
            motifISIs = getISIs(spikes{kk}, mM);
            neuronData(kk).mISIsinging = mean(motifISIs);
            nonMotifISIs = getISIs(spikes{kk}, baselinePeriodsS);
            neuronData(kk).mISIbaseline = mean(nonMotifISIs);
            if isempty(motifISIs) || isempty(nonMotifISIs)
                warning ('allBirdsAnalysis:noSpikesinMotif',...
                    'there are no spikes in the motif intervals, so skipping this neuron')
                neuronData(kk).ISIsinging = NaN(size(binsUsed));
                neuronData(kk).ISIbaseline = NaN(size(binsUsed));
                neuronData(kk).burstFractionSinging = NaN;
                neuronData(kk).burstFractionBaseline = NaN;
                neuronData(kk).mISIsinging = NaN;
                neuronData(kk).mISIbaseline = NaN;
                continue
            end
            
            logMotifISIs    = log10(   motifISIs);
            logNonMotifISIs = log10(nonMotifISIs);
            
            H = ndhist(logMotifISIs'   , nBins, minLog, maxLog);
            I = ndhist(logNonMotifISIs', nBins, minLog, maxLog);
            
            neuronData(kk).ISIsinging = H;
            neuronData(kk).ISIbaseline = I;
            
            toBurst = find(binsUsed <= 0.01, 1,'last'); %hopefully in sec
            neuronData(kk).burstFractionSinging = sum(H(1:toBurst))/sum(H);
            neuronData(kk).burstFractionBaseline = sum(I(1:toBurst))/sum(I);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add metadata
        % label the metadata for each neuron (core/shell and MUA/SUA)
        foo = num2cell(isCoreUnit); [neuronData.isCore] = foo{:};
        foo = num2cell(isMUAUnit);  [neuronData.isMUA ] = foo{:};
        [neuronData.age] = deal(thisAge);
        [neuronData.isPlastic] = deal(isThisPlastic);
        
        saveFile = ['data' filesep birdID filesep 'neuronOnly-' sessionID '.mat'];
        save(saveFile,'neuronData');
        % stuff each session's data into a cell
        sessNeuronData{jj} = neuronData;
        [sessNeuronData{jj}.sessionID] = deal(sessionID);
        fprintf('done.\n');
        toc;
    end
    
    % stuff all the sessions for each bird into another cell
    if ~isempty(sessNeuronData) && any(~cellfun('isempty',sessNeuronData))
        birdNeuronData{ii} = vertcat(sessNeuronData{:});
    end
    progressbar(ii/nBirds);
end

% stuff all of the birds' data containing all sessions into one enormous
% structure.
%% aggregation
bigBirdNeuronData = vertcat(birdNeuronData{:});
%% flag making
isCore    = [bigBirdNeuronData.isCore   ];
isMUA     = [bigBirdNeuronData.isMUA    ];
isPlastic = [bigBirdNeuronData.isPlastic];

including = 'single'; %options to include all neurons, single unit only or multiunit only: 'all', 'single', 'multi'
if strcmp(including,'all')
    theCore  =  isCore;
    theShell = ~isCore;
end
if strcmp(including,'single')
    theCore  =  isCore & ~isMUA;
    theShell = ~isCore & ~isMUA;
end
if strcmp(including,'multi')
    theCore  =  isCore & isMUA;
    theShell = ~isCore & isMUA;
end

%% plotting
%for flagflip = 0:1 %0:1; to combine subsong and plastic
%     isSongStage = (isPlastic == flagflip); %to combine subsong and plastic
%     if flagflip == 0 %to combine subsong and plastic
%         stageName = 'subsong';% to combine subsong and plastic
%     else %to combine subsong and plastic
%         stageName = 'plastic'; %to combine subsong and plastic
%     end %to combine subsong and plastic
isSongStage = ~isMUA; %to combine subsong and plastic
stageName = 'combo';%to combine subsong and plastic
    
    theCoreFlag     =      theCore(isSongStage);
    theShellFlag    =     theShell(isSongStage);
    bigBirdStageData     = bigBirdNeuronData(isSongStage);
    
    isNeverSig = ...
        [bigBirdStageData.t_p]  > 0.05 & ...
        [bigBirdStageData.p_p_50]  > 0.05 & ...
        [bigBirdStageData.pp_p_50] > 0.05 & ...
        [bigBirdStageData.fp_p_50] > 0.05 & ...
        [bigBirdStageData.lp_p_50] > 0.05;

non1 = isNeverSig & theCoreFlag;
non2 = isNeverSig & theShellFlag;
nNon(1) = sum(non1);
nNon(2) = sum(non2);
nonProp(1) = nNon(1)/sum(theCoreFlag);
nonProp(2) = nNon(2)/sum(theShellFlag);

figure;
xSig = {'No Response';'Motif Exc'; 'Motif Inh';'50ms Pre-Syll Exc';...
    '50ms Pre-Syll Inh';'50ms Pre-Motif Exc'; '50ms Pre-Motif Inh'};
[sRS, NsRS, nRS, NnRS] = proporMake('t_p', 'meanRS', bigBirdStageData, theCoreFlag, theShellFlag); %neurons with significant response during motif
[sPR, NsPR, nPR, NnPR] = proporMake('p_p_50', 'preMeanRS_50', bigBirdStageData, theCoreFlag, theShellFlag); %neurons with significant pre-syllable response
[sFP, NsFP, nFP, NnFP] = proporMake('fp_p_50', 'preFMeanRS_50', bigBirdStageData, theCoreFlag, theShellFlag); %neurons with significant pre-motif response
handles_bars = bar(1:7,[nonProp;sRS;NsRS;sPR;NsPR;sFP;NsFP]);
set(gca, 'XTickLabel', xSig);
text(7,0.8,num2str([nNon;nRS;NnRS;nPR;NnPR;nFP;NnFP]));
set(handles_bars(1), 'FaceColor',[0.5 0.5 0.5]);
set(handles_bars(2), 'FaceColor',[1 0 0]);
ylim([0 1]); ylabel('Fraction of neurons');
title([stageName, ' total core n = ', num2str(sum(theCoreFlag)),' shell n = ', num2str(sum(theShellFlag))]);

    %% Firing rates baseline versus singing
    label = {'baseline'; 'singing (excited neurons)'; 'baseline'; 'singing (suppressed neurons)'};
    yLabel = 'Firing Rate (spikes/s) for song-excited neurons';
    fRPC = theCoreFlag & [bigBirdStageData.meanRS] > 0 & [bigBirdStageData.t_p] < 0.05;
    fRPS = theShellFlag & [bigBirdStageData.meanRS] > 0 & [bigBirdStageData.t_p] < 0.05;
    fRNC = theCoreFlag & [bigBirdStageData.meanRS] < 0 & [bigBirdStageData.t_p] < 0.05;
    fRNS = theShellFlag & [bigBirdStageData.meanRS] < 0 & [bigBirdStageData.t_p] < 0.05;
    fourSetsBar([bigBirdStageData(fRPC).FR_silence],[bigBirdStageData(fRPS).FR_silence],...  
        [bigBirdStageData(fRPC).FR_BOS],[bigBirdStageData(fRPS).FR_BOS],...
        [bigBirdStageData(fRNC).FR_silence],[bigBirdStageData(fRNS).FR_silence],...
        [bigBirdStageData(fRNC).FR_BOS],[bigBirdStageData(fRNS).FR_BOS], label, yLabel, stageName);
    
    %% Response Strength
    rSPC = theCoreFlag & [bigBirdStageData.meanRS] > 0 & [bigBirdStageData.t_p] < 0.05;
    rSPS = theShellFlag & [bigBirdStageData.meanRS] > 0 & [bigBirdStageData.t_p] < 0.05;
    rSNC = theCoreFlag & [bigBirdStageData.meanRS] < 0 & [bigBirdStageData.t_p] < 0.05;
    rSNS = theShellFlag & [bigBirdStageData.meanRS] < 0 & [bigBirdStageData.t_p] < 0.05;
    label = {'positive response strength'; 'negative response strength'};
    yLabel = 'Response Strength (spikes/s)';
    twoSetsBar([bigBirdStageData(rSPC).meanRS],[bigBirdStageData(rSPS).meanRS],...
        [bigBirdStageData(rSNC).meanRS],[bigBirdStageData(rSNS).meanRS], label, yLabel, stageName);
    
    %% Response Index
    rIPC = theCoreFlag & [bigBirdStageData.meanRSI] > 0 & [bigBirdStageData.t_p] < 0.05;
    rIPS = theShellFlag & [bigBirdStageData.meanRSI] > 0 & [bigBirdStageData.t_p] < 0.05;
    rINC = theCoreFlag & [bigBirdStageData.meanRSI] < 0 & [bigBirdStageData.t_p] < 0.05;
    rINS = theShellFlag & [bigBirdStageData.meanRS] < 0 & [bigBirdStageData.t_p] < 0.05;
    label = {'positive response index'; 'negative response index'};
    yLabel = 'Response Index';
    twoSetsBar([bigBirdStageData(rIPC).meanRSI],[bigBirdStageData(rIPS).meanRSI],...
        [bigBirdStageData(rINC).meanRSI],[bigBirdStageData(rINS).meanRSI], label, yLabel, stageName);
    %% Standardize Response
    sRPC = theCoreFlag  & [bigBirdStageData.standRS] > 0 & [bigBirdStageData.t_p] < 0.05;
    sRPS = theShellFlag & [bigBirdStageData.standRS] > 0 & [bigBirdStageData.t_p] < 0.05;
    sRNC = theCoreFlag  & [bigBirdStageData.standRS] < 0 & [bigBirdStageData.t_p] < 0.05;
    sRNS = theShellFlag & [bigBirdStageData.standRS] <= 0 & [bigBirdStageData.t_p] < 0.05;
    label = {'positive standardized response'; 'negative standardized response'};
    yLabel = 'Standardized Response (z-score)';
    twoSetsBar([bigBirdStageData(sRPC).standRS],[bigBirdStageData(sRPS).standRS],...
        [bigBirdStageData(sRNC).standRS],[bigBirdStageData(sRNS).standRS], label, yLabel, stageName);
    %% Coefficient of variation of FR
%     vPC = theCoreFlag  & [bigBirdStageData.meanRS] > 0 & [bigBirdStageData.t_p] < 0.05;
%     vPS = theShellFlag & [bigBirdStageData.meanRS] > 0 & [bigBirdStageData.t_p] < 0.05;
%     vNC = theCoreFlag  & [bigBirdStageData.meanRS] < 0 & [bigBirdStageData.t_p] < 0.05;
%     vNS = theShellFlag & [bigBirdStageData.meanRS] < 0 & [bigBirdStageData.t_p] < 0.05;
%     label = {'excited neurons'; 'suppressed neurons'};
%     yLabel = 'coefficent of variation of FR';
%     twoSetsBar([bigBirdStageData(vPC).cvarFR],[bigBirdStageData(vPS).cvarFR],...
%         [bigBirdStageData(vNC).cvarFR],[bigBirdStageData(vNS).cvarFR], label, yLabel, stageName);
%    
    %% Onset and Offset Responses
%     pRS = [bigBirdStageData.p_p_20] < 0.05;
%     ppRS = [bigBirdStageData.pp_p_20] < 0.05;
%      onOffsetPlot(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_20],[bigBirdStageData.postMeanRS_20],...
%         theCoreFlag, theShellFlag, bigBirdStageData,'all', stageName,'20 ms', pRS, ppRS); %all syllables
%         
%     pRS = [bigBirdStageData.p_p_30] < 0.05;
%     ppRS = [bigBirdStageData.pp_p_30] < 0.05; 
%     onOffsetPlot(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_30],[bigBirdStageData.postMeanRS_30],...
%         theCoreFlag, theShellFlag, bigBirdStageData,'all', stageName,'30 ms', pRS, ppRS); %all syllables
%     
%     pRS = [bigBirdStageData.p_p_40] < 0.05;
%     ppRS = [bigBirdStageData.pp_p_40] < 0.05;
%     onOffsetPlot(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_40],[bigBirdStageData.postMeanRS_40],...
%         theCoreFlag, theShellFlag, bigBirdStageData,'all', stageName,'40 ms', pRS, ppRS); %all syllables
    
    pRS = [bigBirdStageData.p_p_50] < 0.05; %onset
    ppRS = [bigBirdStageData.pp_p_50] < 0.05; %offset
    onOffsetPlot(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_50],[bigBirdStageData.postMeanRS_50],...
        theCoreFlag, theShellFlag, bigBirdStageData,'all', stageName,'50 ms', pRS, ppRS); %all syllables
    
    onOffsetPlot(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preFMeanRS_50],[bigBirdStageData.postLMeanRS_50],...
        theCoreFlag, theShellFlag, bigBirdStageData,'FL', stageName, '50 ms', pRS, ppRS); %first and last syllables in motifs
    %% plot all the onsets and offsets and look for trends
%     pRS = [bigBirdStageData.p_p_20] < 0.05;
%     ppRS = [bigBirdStageData.pp_p_20] < 0.05;
%    plotAllFRERA(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_20],[bigBirdStageData.postMeanRS_20],...
%         theCoreFlag, theShellFlag, bigBirdStageData, stageName,'20 ms', pRS, ppRS); %all syllables
%     
%     pRS = [bigBirdStageData.p_p_30] < 0.05;
%     ppRS = [bigBirdStageData.pp_p_30] < 0.05;
%     plotAllFRERA(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_30],[bigBirdStageData.postMeanRS_30],...
%         theCoreFlag, theShellFlag, bigBirdStageData, stageName,'30 ms', pRS, ppRS); %all syllables
%     
%     pRS = [bigBirdStageData.p_p_40] < 0.05;
%     ppRS = [bigBirdStageData.pp_p_40] < 0.05;
%     plotAllFRERA(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_40],[bigBirdStageData.postMeanRS_40],...
%         theCoreFlag, theShellFlag, bigBirdStageData, stageName,'40 ms', pRS, ppRS); %all syllables
    
    pRS = [bigBirdStageData.p_p_50] < 0.05;
    ppRS = [bigBirdStageData.pp_p_50] < 0.05;
    plotAllFRERA(ERAOnsetWindow,onsetBinSize,[bigBirdStageData.preMeanRS_50],[bigBirdStageData.postMeanRS_50],...
        theCoreFlag, theShellFlag, bigBirdStageData, stageName,'50 ms', pRS, ppRS); %all syllables
    
    
 %% compare ISI distribution for singing vs baseline   
   %{ 
    ISIbins = (0:0.03:1);
    
    coreMeanISI     = horzcat(bigBirdStageData( isCoreAndSUAFlag).mISIsinging );
    coreMeanNonISI  = horzcat(bigBirdStageData( isCoreAndSUAFlag).mISIbaseline);
    shellMeanISI    = horzcat(bigBirdStageData(isShellAndSUAFlag).mISIsinging );
    shellMeanNonISI = horzcat(bigBirdStageData(isShellAndSUAFlag).mISIbaseline);
    
    meanISICoreSinging  = histc(coreMeanISI,     ISIbins);
    meanISICoreNon      = histc(coreMeanNonISI,  ISIbins);
    meanISIShellSinging = histc(shellMeanISI,    ISIbins);
    meanISIShellNon     = histc(shellMeanNonISI, ISIbins);
    
    if any(isCoreAndSUAFlag)
        figure;
        subplot(211);
        
        bar(ISIbins,meanISICoreNon);
        hold on
        
        title(sprintf('mean ISIs during baseline in CORE for %s',stageName));
        hold off;
        subplot(212);
        bar(ISIbins,meanISICoreSinging);
        hold on
        
        title(sprintf('mean ISIs during singing in CORE for %s', stageName));
        hold off;
    end
    
    if any(~isCore)
        figure;
        subplot(211);
        bar(ISIbins,meanISIShellNon);
        hold on
        title(sprintf('mean ISIs during baseline in SHELL for %s', stageName));
        hold off;
        subplot(212);
        bar(ISIbins,meanISIShellSinging);
        hold on
        
        title(sprintf('mean ISIs during singing in SHELL for %s', stageName));
        hold off;
    end
        
    %}
    %% Burst Fraction Excited Neurons
    label = {'baseline'; ' excited singing'};
    yLabel = 'BurstFraction (ISI < 10ms)';
    bFCP = theCoreFlag & [bigBirdStageData.t_p] < 0.05 & [bigBirdStageData.meanRS] > 0;
    bFSP = theShellFlag & [bigBirdStageData.t_p] < 0.05 & [bigBirdStageData.meanRS] > 0;

    twoSetsBar([bigBirdStageData(bFCP).burstFractionBaseline],[bigBirdStageData(bFSP).burstFractionBaseline],...
        [bigBirdStageData(bFCP).burstFractionSinging],[bigBirdStageData(bFSP).burstFractionSinging], label, yLabel, stageName);
    
    forAnov = vertcat([bigBirdStageData(bFCP).burstFractionSinging]',[bigBirdStageData(bFCP).burstFractionBaseline]',...
        [bigBirdStageData(bFSP).burstFractionSinging]',[bigBirdStageData(bFSP).burstFractionBaseline]');
    g1 = zeros(length(forAnov),1);
    g1(1:length([bigBirdStageData(bFCP).burstFractionSinging]) + length([bigBirdStageData(bFCP).burstFractionBaseline])) = 1; %core = 1
    g2 = zeros(length(forAnov),1);
    g2(1:length([bigBirdStageData(bFCP).burstFractionSinging])) = 1; %singing = 1
    g2(length([bigBirdStageData(bFCP).burstFractionSinging]) + length([bigBirdStageData(bFCP).burstFractionBaseline]) + 1:length(forAnov) - length([bigBirdStageData(bFSP).burstFractionBaseline])) = 1;
    p = anovan(forAnov,{g1 g2},'model','interaction');
    
    %% Burst Fraction Suppressed Neurons
    label = {'baseline'; 'suppressed singing'};
    yLabel = 'BurstFraction (ISI < 10ms)';
    bFCN = theCoreFlag & [bigBirdStageData.t_p] < 0.05 & [bigBirdStageData.meanRS] < 0;
    bFSN = theShellFlag & [bigBirdStageData.t_p] < 0.05 & [bigBirdStageData.meanRS] < 0;

    twoSetsBar([bigBirdStageData(bFCN).burstFractionBaseline],[bigBirdStageData(bFSN).burstFractionBaseline],...
        [bigBirdStageData(bFCN).burstFractionSinging],[bigBirdStageData(bFSN).burstFractionSinging], label, yLabel, stageName);
    
    forAnov = vertcat([bigBirdStageData(bFCN).burstFractionSinging]',[bigBirdStageData(bFCN).burstFractionBaseline]',...
        [bigBirdStageData(bFSN).burstFractionSinging]',[bigBirdStageData(bFSN).burstFractionBaseline]');
    g1 = zeros(length(forAnov),1);
    g1(1:length([bigBirdStageData(bFCN).burstFractionSinging]) + length([bigBirdStageData(bFCN).burstFractionBaseline])) = 1; %core = 1
    g2 = zeros(length(forAnov),1);
    g2(1:length([bigBirdStageData(bFCN).burstFractionSinging])) = 1; %singing = 1
    g2(length([bigBirdStageData(bFCN).burstFractionSinging]) + length([bigBirdStageData(bFCN).burstFractionBaseline]) + 1:length(forAnov) - length([bigBirdStageData(bFSN).burstFractionBaseline])) = 1;
    p = anovan(forAnov,{g1 g2},'model','interaction');
   
%end %to combine subsong and plastic
%% ANOVA for standardized response across song stage
allStandRS = [bigBirdNeuronData(~isMUA).standRS];
allStandRS(isinf(allStandRS)) = NaN;
factorShell = ~isCore(~isMUA)';
factorPlastic = isPlastic(~isMUA)'; 
factorPos = (allStandRS>0)';

[anpvals, tabble,stats] = anovan(abs(allStandRS)', [factorShell factorPlastic factorPos], 'model','interaction');
fprintf(['p-values, main effect of core vs. shell (F(1,1) = %0.2f, p = %0.2g), ' ...
    'main effect of subsong vs. plastic, (F(1,1) = %0.2f, p = %0.2g), ' ... 
    'main effect of suppressed vs. excited, (F(1,1) = %0.2f, p = %0.2g)\n'], ...
    tabble{2,6}, anpvals(1), tabble{3,6}, anpvals(2), tabble{4,6}, anpvals(3));

%%
profile viewer
profile off

% hello, it's the bottom of the page