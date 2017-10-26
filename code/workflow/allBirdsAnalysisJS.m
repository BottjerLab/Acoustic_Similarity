% allBirdAnalysis.m

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
    allBirdSessions = compositeReport{ii};
    
    nSessions = numel(allBirdSessions);
    sessNeuronData = cell(nSessions,1);
    for jj = 1:nSessions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % setup
        sessManifest = allBirdSessions(jj).manifest;
        sessSpikeFiles = allBirdSessions(jj).spikeFiles;
        
        % get the session/bird ID
        sessionID = allBirdSessions(jj).sessionID;
        birdID = strtok(sessionID, '_');
        birdDataDir = [pwd filesep 'data' filesep birdID filesep];
        
        if numel(sessManifest) == 1 % we only have the audio data from spike,
            fprintf('Missing parsing data, skipping session %s.\n',sessionID);
            continue;
        end
        
        if isempty(allBirdSessions(jj).spikeFiles)
            fprintf('Missing stereotrode data, skipping session %s.\n',sessionID);
            continue;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load motifs (NB: should be in standard event struct, column format)
        mM = loadFromManifest(sessManifest,'derivedMotifs'); %going to change this to new motifs
        if isempty(mM)% || ~any(strcmp('BOS', {mM.type}))
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
        
        % get core/shellness and MUA/SUA-ness of each neuron
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
        fprintf('Analyzing session %s/%s...\n', birdID, sessionID);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pre-loadable part
        % Compare FR during singing vs local baseline - constructBaseline
        % does automatic exclusions
        myCodeIsPerfect = true;
        baselinePeriodFile = [birdDataDir 'baselinePeriods-' sessionID '.mat'];
        if myCodeIsPerfect && exist(baselinePeriodFile,'file')
            load(baselinePeriodFile);
        else % if this is the first time we look for silent baseline
            % just run it
            baselinePeriods = constructBaseline(mM, 'BOS', [-4, -2], [2, 4], params);
            [baselinePeriods.type] = deal('localBase');
            
            songStruct = loadFromManifest(sessManifest,'songStruct');
            baselinePeriodsS = constructBaselineSilence(songStruct, 2, aS, params);
            [baselinePeriodsS.type] = deal('silence');
            
            save(baselinePeriodFile,'baselinePeriods','baselinePeriodsS');
        end
        neuronData = getRS([mM; baselinePeriodsS], spikes, 'BOS', 'silence', params, 'verbose',false); %compareTwoCatFiringRates
        
        % Get event-related firing rate to syllable onsets
        ERAOnsetWindow = [-0.2 0.2]; %s, the window in which to collect ERA hists
        onsetBinSize = 2e-3; % s
        tempVar = eraFind(ERAOnsetWindow ,aS,'onset', onsetBinSize, spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).ERAOnsethist = tempVar(i,:)';
        end


        % Get event-related firing rate to syllable offsets
        ERAOffsetWindow = [-0.2 0.2]; %s, the window in which to collect ERA hists
        offsetBinSize = 2e-3; % s
        tempVar = eraFind(ERAOffsetWindow,aS,'offset', offsetBinSize,spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).ERAOffsethist = tempVar(i,:)';
        end

                
        % Get event-related firing rate to syllable onsets for first syllables in motif
        tempVar = eraFind([-0.2 0.2],fAS,'onset',onsetBinSize,spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).firstERAOnsethist = tempVar(i,:)';
        end
        
        % Get event-related firing rate to syllable onsets for last syllables in motif
        tempVar = eraFind([-0.2 0.2],lAS,'offset',offsetBinSize,spikes,baselinePeriodsS, params);
        for i = 1:size(tempVar,1)
            neuronData(i).lastERAOffsethist = tempVar(i,:)';
        end
        
        % get the normalized log ISI distributions for singing and baseline, Jenny wrote
        % this, ekkkkk!
        nBins = 400;
        maxLog = 2;
        minLog = -3;
        binsUsed = logspace(minLog,maxLog,nBins); %time axis
        binsUsed = binsUsed';
        
        for kk = 1:numel(spikes)
            motifISIs = getISIs(spikes{kk}, mM);
            nonMotifISIs = getISIs(spikes{kk}, baselinePeriodsS);
            if isempty(motifISIs) || isempty(nonMotifISIs)
                warning ('allBirdsAnalysis:noSpikesinMotif',...
                    'there are no spikes in the motif intervals, so skipping this neuron')
                neuronData(kk).ISIsinging = NaN(size(binsUsed));
                neuronData(kk).ISIbaseline = NaN(size(binsUsed));
                neuronData(kk).burstFractionSinging = NaN;
                neuronData(kk).burstFractionBaseline = NaN;
                continue
            end
            
            logMotifISIs    = log10(   motifISIs);
            logNonMotifISIs = log10(nonMotifISIs);
            
            H = ndhist(logMotifISIs', nBins, minLog, maxLog);
            I = ndhist(logNonMotifISIs', nBins, minLog, maxLog);
            
            neuronData(kk).ISIsinging = H;
            neuronData(kk).ISIbaseline = I;
            
            toBurst = find(binsUsed <= 0.01, 1,'last'); %hopefully in sec
            neuronData(kk).burstFractionSinging = sum(H(1:toBurst))/sum(H);
            neuronData(kk).burstFractionBaseline = sum(I(1:toBurst))/sum(I);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % label the metadata for each neuron (core/shell and MUA/SUA)
        foo = num2cell(isCoreUnit); [neuronData.isCore] = foo{:};
        foo = num2cell(isMUAUnit);  [neuronData.isMUA ] = foo{:};
        
        % stuff each session's data into a cell
        sessNeuronData{jj} = neuronData;
    end
    
    % stuff all the sessions for each bird into another cell
    if ~isempty(sessNeuronData) && any(~cellfun('isempty',sessNeuronData))
        birdNeuronData{ii} = vertcat(sessNeuronData{:});
    end
    progressbar(ii/nBirds);
end

% stuff all of the birds' data containing all sessions into one enormous
% structure.
bigBirdNeuronData = vertcat(birdNeuronData{:});
%% plotting/aggregating
isCore = [bigBirdNeuronData.isCore];
isMUA  = [bigBirdNeuronData.isMUA ];
bigBirdSUAData = bigBirdNeuronData(~isMUA);
isCoreSUAOnly = isCore(~isMUA);

ERAcumHistCore  = sum(horzcat(bigBirdSUAData( isCoreSUAOnly).ERAOnsethist ),2)/sum( isCoreSUAOnly);
ERAcumHistShell = sum(horzcat(bigBirdSUAData(~isCoreSUAOnly).ERAOnsethist ),2)/sum(~isCoreSUAOnly);
ERAcumHistCoreSem  = std(horzcat(bigBirdSUAData( isCoreSUAOnly).ERAOnsethist ),[],2)/sqrt(length([bigBirdSUAData( isCoreSUAOnly).ERAOnsethist])); %Jenny added
ERAcumHistShellSem = std(horzcat(bigBirdSUAData(~isCoreSUAOnly).ERAOnsethist ),[],2)/sqrt(length([bigBirdSUAData(~isCoreSUAOnly).ERAOnsethist])); %Jenny added

% local windowed baseline
subplot(211);
lWB = 0;
%lWB = mean([bigBirdSUAData.FR_localBase]) * onsetBinSize;
xx = (ERAOnsetWindow(1):onsetBinSize:ERAOnsetWindow(2)) * 1000;
xx = xx'; %Jenny added

%{
if ~any(isCore) %Jenny added
    ERAcumHistCore = zeros(length(xx),1);
end
if ~any(isShell) %Jenny added
    ERAcumHistShell = zeros(length(xx),1);
end
%}
% redone for transparency example
if any(isCore)
    plotLineTransparentSEM(xx(1:end-2), ERAcumHistCore(1:end-2), ERAcumHistCoreSem(1:end-2), [0.2 0.2 1]);
    hold on;
end
if any(~isCore)
    plotLineTransparentSEM(xx(1:end-2), ERAcumHistShell(1:end-2), ERAcumHistShellSem(1:end-2), [1 0.2 0.2]);
end
if any(isCore)
    hold off;
end


% plot(xx(1:end-2), ERAcumHistCore(1:end-2),'-b', ... %Jenny changed to -2
%      xx(1:end-2), ERAcumHistShell(1:end-2),'-r', 'LineWidth', 2);  %Jenny changed to -2
% hold on; plot(xx([1,end]), [lWB lWB], 'k--'); hold off;
% hold on; plot([0 0], ylim, 'k--'); hold off
% hold on; plot(xx(1:end-2),ERAcumHistCore(1:end-2) + ERAcumHistCoreStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),ERAcumHistCore(1:end-2) - ERAcumHistCoreStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),ERAcumHistShell(1:end-2) + ERAcumHistShellStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),ERAcumHistShell(1:end-2) - ERAcumHistShellStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold off; %Jenny added

xlabel('time relative to onset (ms)');


subplot(212);
ERAcumHistCore  = sum(horzcat(bigBirdSUAData( isCoreSUAOnly).ERAOffsethist ),2)/sum( isCoreSUAOnly);
ERAcumHistShell = sum(horzcat(bigBirdSUAData(~isCoreSUAOnly).ERAOffsethist ),2)/sum(~isCoreSUAOnly);
ERAcumHistCoreStd  = std(horzcat(bigBirdSUAData( isCoreSUAOnly).ERAOffsethist ),[],2)/sqrt(length([bigBirdSUAData( isCoreSUAOnly).ERAOffsethist])); %Jenny added
ERAcumHistShellStd = std(horzcat(bigBirdSUAData(~isCoreSUAOnly).ERAOffsethist ),[],2)/sqrt(length([bigBirdSUAData(~isCoreSUAOnly).ERAOffsethist])); %Jenny added


xx = (ERAOffsetWindow(1):offsetBinSize:ERAOffsetWindow(2)) * 1000;
xx = xx'; %Jenny added

%{
if isempty(ERAcumHistCore) %Jenny added
    ERAcumHistCore = zeros(length(xx),1);
end
if isempty(ERAcumHistShell) %Jenny added
    ERAcumHistShell = zeros(length(xx),1);
end
%}
% redone for transparency example
if any(isCore)
    plotLineTransparentSEM(xx(1:end-2), ERAcumHistCore(1:end-2), ERAcumHistCoreStd(1:end-2), [0.2 0.2 1]);
    hold on;
end
if any(~isCore)
    plotLineTransparentSEM(xx(1:end-2), ERAcumHistShell(1:end-2), ERAcumHistShellStd(1:end-2), [1 0.2 0.2]);
end
if any(isCore)
    hold off;
end

% plot(xx(1:end-2), ERAcumHistCore(1:end-2) , '-b',... %Jenny changed to -2
%      xx(1:end-2), ERAcumHistShell(1:end-2), 'r-', 'LineWidth', 2); %Jenny changed to -2
% hold on; plot(xx([1,end]), [lWB lWB], 'k--'); hold off;
% hold on; plot([0 0], ylim, 'k--'); hold off
% hold on; plot(xx(1:end-2),ERAcumHistCore(1:end-2) + ERAcumHistCoreStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),ERAcumHistCore(1:end-2) - ERAcumHistCoreStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),ERAcumHistShell(1:end-2) + ERAcumHistShellStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),ERAcumHistShell(1:end-2) - ERAcumHistShellStd(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold off; %Jenny added

xlabel('time relative to offset (ms)');
legend('core','shell');

%now for the first and last only
figure
fERAcumHistCore  = sum(horzcat(bigBirdSUAData( isCoreSUAOnly).firstERAOnsethist ),2)/sum( isCoreSUAOnly);
fERAcumHistShell = sum(horzcat(bigBirdSUAData(~isCoreSUAOnly).firstERAOnsethist ),2)/sum(~isCoreSUAOnly);
fERAcumHistCoreSem  = std(horzcat(bigBirdSUAData( isCoreSUAOnly).firstERAOnsethist ),[],2)/sqrt(length([bigBirdSUAData( isCoreSUAOnly).firstERAOnsethist])); %Jenny added
fERAcumHistShellSem = std(horzcat(bigBirdSUAData(~isCoreSUAOnly).firstERAOnsethist ),[],2)/sqrt(length([bigBirdSUAData(~isCoreSUAOnly).firstERAOnsethist])); %Jenny added

% local windowed baseline
subplot(211);
lWB = 0;
%lWB = mean([bigBirdSUAData.FR_localBase]) * onsetBinSize;
xx = (ERAOnsetWindow(1):onsetBinSize:ERAOnsetWindow(2)) * 1000;
xx = xx'; %Jenny added

%{
if isempty(fERAcumHistCore) %Jenny added
    fERAcumHistCore = zeros(length(xx),1);
end
if isempty(fERAcumHistShell) %Jenny added
    fERAcumHistShell = zeros(length(xx),1);
end
%}
% redone for transparency example
if any(isCore)
    plotLineTransparentSEM(xx(1:end-2), fERAcumHistCore(1:end-2), fERAcumHistCoreSem(1:end-2), [0.2 0.2 1]);
    hold on;
end
if any(~isCore)
    plotLineTransparentSEM(xx(1:end-2), fERAcumHistShell(1:end-2), fERAcumHistShellSem(1:end-2), [1 0.2 0.2]);
end
if any(isCore)
    hold off;
end


% plot(xx(1:end-2), fERAcumHistCore(1:end-2),'-b', ... %Jenny changed to -2
%      xx(1:end-2), fERAcumHistShell(1:end-2),'-r', 'LineWidth', 2);  %Jenny changed to -2
% hold on; plot(xx([1,end]), [lWB lWB], 'k--'); hold off;
% hold on; plot([0 0], ylim, 'k--'); hold off
% hold on; plot(xx(1:end-2),fERAcumHistCore(1:end-2) + fERAcumHistCoreSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),fERAcumHistCore(1:end-2) - fERAcumHistCoreSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),fERAcumHistShell(1:end-2) + fERAcumHistShellSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),fERAcumHistShell(1:end-2) - fERAcumHistShellSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold off; %Jenny added

xlabel('time relative to onset (ms)');
title('Relative to first and last syllables only');


subplot(212);
lERAcumHistCore  = sum(horzcat(bigBirdSUAData( isCoreSUAOnly).lastERAOffsethist ),2)/sum( isCoreSUAOnly);
lERAcumHistShell = sum(horzcat(bigBirdSUAData(~isCoreSUAOnly).lastERAOffsethist ),2)/sum(~isCoreSUAOnly);
lERAcumHistCoreSem  = std(horzcat(bigBirdSUAData( isCoreSUAOnly).lastERAOffsethist ),[],2)/sqrt(length([bigBirdSUAData( isCoreSUAOnly).lastERAOffsethist])); %Jenny added
lERAcumHistShellSem = std(horzcat(bigBirdSUAData(~isCoreSUAOnly).lastERAOffsethist ),[],2)/sqrt(length([bigBirdSUAData(~isCoreSUAOnly).lastERAOffsethist])); %Jenny added


xx = (ERAOffsetWindow(1):offsetBinSize:ERAOffsetWindow(2)) * 1000;
xx = xx'; %Jenny added

%{
if isempty(lERAcumHistCore) %Jenny added
    lERAcumHistCore = zeros(length(xx),1);
end
if isempty(lERAcumHistShell) %Jenny added
    lERAcumHistShell = zeros(length(xx),1);
end
%}
% redone for transparency example
if any(isCore)
    plotLineTransparentSEM(xx(1:end-2), lERAcumHistCore(1:end-2), lERAcumHistCoreSem(1:end-2), [0.2 0.2 1]);
    hold on;
end
if any(~isCore)
    plotLineTransparentSEM(xx(1:end-2), lERAcumHistShell(1:end-2), lERAcumHistShellSem(1:end-2), [1 0.2 0.2]);
end
if any(isCore)
    hold off;
end

% plot(xx(1:end-2), lERAcumHistCore(1:end-2) , '-b',... %Jenny changed to -2
%      xx(1:end-2), lERAcumHistShell(1:end-2), 'r-', 'LineWidth', 2); %Jenny changed to -2
% hold on; plot(xx([1,end]), [lWB lWB], 'k--'); hold off;
% hold on; plot([0 0], ylim, 'k--'); hold off
% hold on; plot(xx(1:end-2),lERAcumHistCore(1:end-2) + lERAcumHistCoreSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),lERAcumHistCore(1:end-2) - lERAcumHistCoreSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),lERAcumHistShell(1:end-2) + lERAcumHistShellSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold on; plot(xx(1:end-2),lERAcumHistShell(1:end-2) - lERAcumHistShellSem(1:end-2), 'Color', [0.7 0.7 0.7]); %Jenny added
% hold off; %Jenny added

xlabel('time relative to offset (ms)');
legend('core','shell');

%% compare ISI distribution for singing vs baseline, Jenny wrote this (reason why it is crazy looking)
%change to subplots
%  for mm = 1:length(bigBirdSUAData)
%      figure
%      subplot(211)
%      plot(binsUsed,bigBirdSUAData(mm).ISIbaseline,'b');
%      set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
%      xlabel('log time (seconds)')
%      ylabel('count')
%      if bigBirdSUAData(mm).isCore
%          thisName = 'CORE';
%      else
%          thisName = 'SHELL';
%      end
%      title(['ISI distribution for ',thisName, 'neuron'])
%      legend('Baseline')
%      subplot(212)
%      plot(binsUsed,bigBirdSUAData(mm).ISIsinging, 'k');
%      set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
%      xlabel('log time (seconds)')
%      ylabel('count')
%      legend('Singing')
%      hold off
%      saveCurrFigure(sprintf('ISI_%03d',mm));
%      close (gcf)
%  end

NonISICore      = nanmean(horzcat(bigBirdSUAData( isCoreSUAOnly).ISIbaseline), 2);
normNonISICore  = NonISICore/sum(NonISICore(:));
NonISIShell     = nanmean(horzcat(bigBirdSUAData(~isCoreSUAOnly).ISIbaseline), 2);
normNonISIShell = NonISIShell/sum(NonISIShell(:));
MotISICore      = nanmean(horzcat(bigBirdSUAData( isCoreSUAOnly).ISIsinging) , 2);
normMotISICore  = MotISICore/sum(MotISICore(:));
MotISIShell     = nanmean(horzcat(bigBirdSUAData(~isCoreSUAOnly).ISIsinging) , 2);
normMotISIShell = MotISIShell/sum(MotISIShell(:));
NonISICoreSem   = nanstd (horzcat(bigBirdSUAData( isCoreSUAOnly).ISIbaseline), [],2)./sqrt(sum(~isnan([bigBirdSUAData( isCoreSUAOnly).ISIbaseline]),2));
normNonISICoreSem = NonISICoreSem/sum(NonISICoreSem(:));
NonISIShellSem  = nanstd (horzcat(bigBirdSUAData(~isCoreSUAOnly).ISIbaseline), [],2)./sqrt(sum(~isnan([bigBirdSUAData(~isCoreSUAOnly).ISIbaseline]),2));
normNonISIShellSem = NonISIShellSem/sum(NonISIShellSem(:));
MotISICoreSem   = nanstd (horzcat(bigBirdSUAData( isCoreSUAOnly).ISIsinging) , [],2)./sqrt(sum(~isnan([bigBirdSUAData( isCoreSUAOnly).ISIsinging ]),2));
normMotISICoreSem = MotISICoreSem/sum(MotISICoreSem(:));
MotISIShellSem  = nanstd (horzcat(bigBirdSUAData(~isCoreSUAOnly).ISIsinging) , [],2)./sqrt(sum(~isnan([bigBirdSUAData(~isCoreSUAOnly).ISIsinging ]),2));
normMotISIShellSem = MotISIShellSem/sum(MotISIShellSem(:));

if any(isCore)
    figure;
    subplot(211);
    %plotLineTransparentSEM(binsUsed, NonISICore, NonISICoreSem, [0.8 0.8 0.8]);
    
plot(binsUsed, NonISICore, '-k'); %binsUsed comes from calcuation above
hold on;
plot(binsUsed, NonISICore + NonISICoreSem, 'Color', [0.7 0.7 0.7]);
plot(binsUsed, NonISICore - NonISICoreSem, 'Color', [0.7 0.7 0.7]);
    
    title(sprintf('ISIs during baseline in CORE'));
    set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
    hold off;
    subplot(212);
    plot(binsUsed, MotISICore, '-k');
    hold on;
    plot(binsUsed, MotISICore + MotISICoreSem, 'Color', [0.7 0.7 0.7]);
    plot(binsUsed, MotISICore - MotISICoreSem, 'Color', [0.7 0.7 0.7]);
    title(sprintf('ISIs during singing in CORE'));
    set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
    hold off;
end

if any(~isCore)
figure;
subplot(211);
plot(binsUsed, normNonISIShell, '-k'); %binsUsed comes from calcuation above
hold on;
plot(binsUsed, normNonISIShell + normNonISIShellSem, 'Color', [0.7 0.7 0.7]);
plot(binsUsed, normNonISIShell - normNonISIShellSem, 'Color', [0.7 0.7 0.7]);
title(sprintf('normalized ISIs during baseline in SHELL'));
set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
hold off;
subplot(212);
plot(binsUsed, normMotISIShell, '-k');
hold on;
plot(binsUsed, normMotISIShell + normMotISIShellSem, 'Color', [0.7 0.7 0.7]);
plot(binsUsed, normMotISIShell - normMotISIShellSem, 'Color', [0.7 0.7 0.7]);
title(sprintf('normalized ISIs during singing in SHELL'));
set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
hold off;
end

if any(isCore)
figure;
subplot(211);
plot(binsUsed, normNonISICore, '-k'); %binsUsed comes from calcuation above
hold on;
plot(binsUsed, normNonISICore + normNonISICoreSem, 'Color', [0.7 0.7 0.7]);
plot(binsUsed, normNonISICore - normNonISICoreSem, 'Color', [0.7 0.7 0.7]);
title(sprintf('normalized ISIs during baseline in CORE'));
set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
hold off;
subplot(212);
plot(binsUsed, normMotISICore, '-k');
hold on;
plot(binsUsed, normMotISICore + normMotISICoreSem, 'Color', [0.7 0.7 0.7]);
plot(binsUsed, normMotISICore - normMotISICoreSem, 'Color', [0.7 0.7 0.7]);
title(sprintf('normalized ISIs during singing in CORE'));
set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
hold off;
end

if any(~isCore)
figure;
subplot(211);
plot(binsUsed, normNonISIShell, '-k'); %binsUsed comes from calcuation above
hold on;
plot(binsUsed, normNonISIShell + normNonISIShellSem, 'Color', [0.7 0.7 0.7]);
plot(binsUsed, normNonISIShell - normNonISIShellSem, 'Color', [0.7 0.7 0.7]);
title(sprintf('ISIs during baseline in SHELL'));
set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
hold off;
subplot(212);
plot(binsUsed, normMotISIShell, '-k');
hold on;
plot(binsUsed, normMotISIShell + normMotISIShellSem, 'Color', [0.7 0.7 0.7]);
plot(binsUsed, normMotISIShell - normMotISIShellSem, 'Color', [0.7 0.7 0.7]);
title(sprintf('ISIs during singing in SHELL'));
set(gca, 'XScale', 'log', 'XLim', [10^-3 10^2]);
hold off;
end


%% anova!
%realBurstFractionCoreS  = bigBirdSUAData( isCoreSUAOnly).burstFractionSinging(~isnan([bigBirdSUAData( isCoreSUAOnly).burstFractionSinging]));
realBurstFractionCoreS  = [bigBirdSUAData( isCoreSUAOnly).burstFractionSinging];
realBurstFractionCoreS(isnan(realBurstFractionCoreS)) = [];
realBurstFractionCoreB  = [bigBirdSUAData( isCoreSUAOnly).burstFractionBaseline];
realBurstFractionCoreB(isnan(realBurstFractionCoreB)) = [];
realBurstFractionShellS = [bigBirdSUAData(~isCoreSUAOnly).burstFractionSinging];
realBurstFractionShellS(isnan(realBurstFractionShellS)) = [];
realBurstFractionShellB = [bigBirdSUAData(~isCoreSUAOnly).burstFractionBaseline];
realBurstFractionShellB(isnan(realBurstFractionShellB)) = [];

burstFractionCoreS  = mean(realBurstFractionCoreS);
burstFractionCoreB  = mean(realBurstFractionCoreB);
burstFractionShellS = mean(realBurstFractionShellS);
burstFractionShellB = mean(realBurstFractionShellB);

burstFractionCoreSsem  = std(realBurstFractionCoreS)/sqrt(length(realBurstFractionCoreS));
burstFractionCoreBsem  = std(realBurstFractionCoreB)/sqrt(length(realBurstFractionCoreB));
burstFractionShellSsem = std(realBurstFractionShellS)/sqrt(length(realBurstFractionShellS));
burstFractionShellBsem = std(realBurstFractionShellB)/sqrt(length(realBurstFractionShellB));

[coreT, coreP] = ttest2(realBurstFractionCoreS,realBurstFractionCoreB);
[shellT,shellP] = ttest2(realBurstFractionShellS,realBurstFractionShellB);
forAnov = vertcat(realBurstFractionCoreS',realBurstFractionCoreB',realBurstFractionShellS',realBurstFractionShellB');
g1 = zeros(length(forAnov),1);
g1(1:length(realBurstFractionCoreS) + length(realBurstFractionCoreB)) = 1; %core = 1
g2 = zeros(length(forAnov),1);
g2(1:length(realBurstFractionCoreS)) = 1; %singing = 1
g2(length(realBurstFractionCoreS) + length(realBurstFractionCoreB) + 1:length(forAnov) - length(realBurstFractionShellB)) = 1;
p = anovan(forAnov,{g1 g2})
tVals = [coreT,shellT];
figure;
barFraction = [burstFractionCoreS burstFractionShellS; burstFractionCoreB burstFractionShellB];
barError = [burstFractionCoreSsem burstFractionShellSsem; burstFractionCoreBsem burstFractionShellBsem];

plotBarError(barFraction,barError,[],tVals,[0.7 0.7 0.7; 1 0 0]); %using John's plotting
barLabel = {'singing', 'baseline'};
set(gca,'XTickLabel',barLabel);
title('Burst Fraction');
legend('CORE','SHELL');
ylabel('fraction of spikes with ISI <= 5ms');
%%
% determine FR characteristics surrounding syllable onset & offset

% characterize any neuron types based on those measures of FR during song
% and FR surrounding offsets and onsets

% determine if those "type" can be explained by song development/age

%% note: can be made into a plotting function
RSbins = -2:0.1:2;
RSI = [bigBirdSUAData.meanRSI];
hCoreRS = histc(RSI(isCoreSUAOnly), RSbins);
hShellRS = histc(RSI(~isCoreSUAOnly), RSbins);
figure;

subplot(211)
bar(RSbins, hCoreRS, 'b')
xlabel('normalized response strength'); ylabel('count');
xlim([-1 1])
title('SUA Only');
subplot(212)
bar(RSbins, hShellRS, 'r')
xlabel('normalized response strength'); ylabel('count');
xlim([-1 1])

profile viewer
profile off

% hello it's the bottom of the page