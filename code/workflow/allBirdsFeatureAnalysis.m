% allBirdsFeatureAnalysis.m
%profile on
compositeReport = reportOnData('','',defaultParams,'verbose',true);
if ~iscell(compositeReport)
    compositeReport = {compositeReport};
end
nBirds = numel(compositeReport);

birdUnitCorrData = cell(nBirds, 1);

if ~exist('MUAinfo','var')
    MUAinfo = getSpikeMUAData;
end
%%
progressbar;
birdID = cell(nBirds,1);
mdlSSE = []; glmSSE = [];
nNeurons = 351; % we know this from previous counts
for ii = 1:nBirds
    sessionsForThisBird = compositeReport{ii};
    nSessions = numel(sessionsForThisBird);
    sessNeuronData = cell(nSessions,1);
    
    firstSessionID = sessionsForThisBird(1).sessionID;
    birdID{ii} = strtok(firstSessionID, '_');
    birdDataDir = ['data' filesep birdID{ii} filesep];
    load([birdDataDir 'allSpecs-' birdID{ii} '.mat']);
    
    %% censor features:
    censorList = {'harmonicPitch', 'aperiodicity'}; % JMA added aperiodicity because we didn't do this for clustering?
    
    for jj = 1:numel(censorList)
        censorItem = censorList{jj};
        len = length(censorItem);
        fn = fieldnames(featureTable);
        toRemove = strncmp(censorItem, fn, len);
        featureTable = rmfield(featureTable, fn(toRemove));
    end
    %%
    for jj = 1:nSessions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % setup
        sessManifest = sessionsForThisBird(jj).manifest;
        sessSpikeFiles = sessionsForThisBird(jj).spikeFiles;
        
        % get the session/bird ID
        sessionID = sessionsForThisBird(jj).sessionID;
        
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
        % note: this is already loaded from allSpecs.
        %{
        mM = loadFromManifest(sessManifest,'derivedMotifs');
        if isempty(mM) || ~any(strcmp('BOS', {mM.type}))
            fprintf('Missing BOS motifs, skipping session %s.\n',sessionID);
            continue;
        end
        
        aS = loadFromManifest(sessManifest,'approvedSyllables');
        if isempty(aS)
            %warning('allBirdsAnalysis:missingApproved', 'Missing approved syllables, deferring to auto styllables in session %s.\n',sessionID);
            %aS = loadFromManifest(sessManifest,'syllables');
            %if isempty(aS)
                warning('allBirdsAnalysis:noSyllables','Missing syllables, continuing...');
                continue;
            %end
        end
        %}
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
        fprintf('Analyzing session %s/%s...\n', birdID{ii}, sessionID);
        
        matchingFileA = [birdDataDir sessionID '.mat'];
        matchingFileB = [birdDataDir filesep sessionID '.mat'];
        
        matchesFile = strcmpi(matchingFileA, {DRsylls.file}) | strcmpi(matchingFileB, {DRsylls.file});
        
        matchSylls = DRsylls(matchesFile);
        matchFeatures = featureTable(matchesFile);
        matchSpectra = spectra(matchesFile);
        matchSylls = addPrePost(matchSylls,[],'preroll', 50, 'postroll', 0);%adding in time lag for premotor-ness
        [matchSylls.type] = deal('syllable');
        
        % no point in regressing if we have fewer syllables than features
        % note: in reality we should have ~1000 syllables because we want
        % 10x the number of features
        if numel(matchSylls) < 3*numel(fieldnames(featureTable))
            warning('Not enough syllables to regress against features, skipping...');
            continue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % correlations to features are against response strength and use
        % silent baselines
        baselinePeriodsS = loadFromManifest(sessManifest, 'baselinePeriodsS');
        if isempty(baselinePeriodsS)
            error('can''t find baseline periods');
        end
        [baselinePeriodsS.file] = deal('whatever,dude'); % this doesn't matter
        [baselinePeriodsS.type] = deal('baseline');
        if size(baselinePeriodsS,1) == 1, baselinePeriodsS = baselinePeriodsS'; end
        if size(matchSylls,1) == 1, matchSylls = matchSylls'; end
        % get rid of extra fields from either structure
        fB = fieldnames(baselinePeriodsS);
        fM = fieldnames(matchSylls);
        
        baselinePeriodS = rmfield(baselinePeriodsS, setdiff(fB,fM));
        matchSylls = rmfield(matchSylls, setdiff(fM,fB));
        
        % get neural response to all syllables
        [dummyData, pairedFiring] = getRS([matchSylls; baselinePeriodsS], spikes, ...
            'syllable', 'baseline', ...
            'p_ttest', 'meanRS', 'meanRSI', 'standRS','cvarFR',...
            params, 'verbose',false);
        RS = cellfun(@(x) -diff(x), pairedFiring, 'UniformOutput', false);
        
        %correlateData = NaN(numel(spikes), numel(fieldnames(matchFeatures)));
        neuronFeatureData = struct([]);
        for ll = 1:numel(spikes)
            if isempty(RS{ll}), error('no response strengths...: ');  end
            %TODO: introduce different time delays?
            [stats, statFamily] = ...
                findSpikeCorrelates(matchSylls, matchFeatures, RS{ll});
            
            neuronFeatureData(ll).neuronIndex = ll;
            neuronFeatureData(ll).pResponsive = dummyData(ll).p_ttest;
            neuronFeatureData(ll).isExcited = (dummyData(ll).meanRS > 0);
            neuronFeatureData(ll).pCorr = stats.pCorr;
            neuronFeatureData(ll).tCorr = stats.tCorr;
            neuronFeatureData(ll).RSquare = stats.RSquare;
            neuronFeatureData(ll).pFull = stats.pAll; % p-value
            neuronFeatureData(ll).fFull = stats.fAll; %JMA added
            neuronFeatureData(ll).meanRS = dummyData(ll).meanRS;
            neuronFeatureData(ll).fPartial = stats.fPartial;
            neuronFeatureData(ll).pPartial = stats.pPartial;
            neuronFeatureData(ll).statFamily = statFamily;
            neuronFeatureData(ll).numSyll = numel(matchSylls);
            
            % if any of the family p values is low enough, plot the graph
            % between feature
            
            nTests = nNeurons * numel(statFamily);
            pMult = (1 - (1 - 0.05)^(1/nTests));
            if any([neuronFeatureData(ll).statFamily.P] < pMult)
                familyStats = neuronFeatureData(ll).statFamily;
                for mm = 1:numel(familyStats)
                    if (familyStats(mm).P > pMult), continue; end
                    
                    fldNames = fieldnames(matchFeatures);
                    familyName = familyStats(mm).feature;
                    lenN = length(familyName);
                    isFamily = strncmp(familyName, fldNames, lenN);
                    
                    % plotting all variables among the family against RS
                    %featureVals = cellfun(@(x) [matchFeatures.(x)], fldNames(isFamily), 'UniformOutput', false);
                    %featureVals = vertcat(featureVals{:});
                    %plotmatrix([featureVals', RS{ll}']);
                    
                    %plotting the most correlated variable among the family against RS
                    [~,bestFeatureIdx] = min(abs(neuronFeatureData(ll).pCorr(isFamily)));
                    find(isFamily); bestFeatureIdx = ans(bestFeatureIdx);
                    featureVals = [matchFeatures.(fldNames{bestFeatureIdx})];
                    hf = figure;
                    plot(featureVals, RS{ll},'k.');
                    xlabel(fldNames{bestFeatureIdx}, 'Interpreter','none');
                    ylabel('Response Strength');
                    title(sprintf('Session %s, Neuron %d, Family F test p = %0.2g, Individual p test = %0.2g', ...
                        sessionID, ll, familyStats(mm).P, neuronFeatureData(ll).pCorr(bestFeatureIdx)),...
                        'Interpreter','none');
                    saveCurrFigure(sprintf('figures/individualFeatureCorrelations/%s-n%d-f%s', ...
                        sessionID, ll, fldNames{bestFeatureIdx}));
                    close(hf)
                end
            end
            % get the top v bottom 25% for each feature
            featureNames = fieldnames(matchFeatures);
            % denominator for change in RS
            RS_syll = pairedFiring{ll}(1,:); RS_base = pairedFiring{ll}(2,:);
            covRS = cov(RS_syll, RS_base); if numel(covRS) > 1, covRS = covRS(2,1); end;
            
            
            neuronFeatureData(ll).diffRS = zeros(1, numel(featureNames));

            for mm = 1:numel(featureNames)
                featureCol = [matchFeatures.(featureNames{mm})];
                QRs = prctile(featureCol, [25 75]); %JMA changed from quantile(featureCol, [0.25 0.75]), but maybe it's the same
                lowRS  = mean(RS{ll}(featureCol <= QRs(1)));%JMA changed to <= because some quartiles are zero
                highRS = mean(RS{ll}(featureCol >= QRs(2))); %JMA added = to match the first quartile
                
                RSSE = sqrt(var(RS_syll) + var(RS_base) - 2 * covRS) / sqrt(numel(featureCol));
                neuronFeatureData(ll).diffRS(mm) = (lowRS - highRS) / RSSE; %I think maybe z-score should be calculated separately for low and high feature values. The SE isn't really of these values so shouldn't be denominator- JMA
               lowCV = std(RS_syll(featureCol <= QRs(1)))/mean(RS_syll (featureCol <= QRs(1))); %JMA added to compare CV for variance in FF; <= because some quartiles are zero
               highCV = std(RS_syll(featureCol >= QRs(2)))/mean(RS_syll (featureCol >= QRs(2)));%JMA added = to match the first quartile
               lowFR = mean(RS_syll (featureCol <= QRs(1)));
               highFR = mean(RS_syll (featureCol >= QRs(2)));
               lowVar = var(RS_syll (featureCol <= QRs(1)));
               highVar = var(RS_syll (featureCol >= QRs(2)));
               neuronFeatureData(ll).lowCV(mm) = lowCV;
               neuronFeatureData(ll).highCV(mm) = highCV;
               neuronFeatureData(ll).lowFR(mm) = lowFR;
               neuronFeatureData(ll).highFR(mm) = highFR;
               neuronFeatureData(ll).lowVar(mm) = lowVar;
               neuronFeatureData(ll).highVar(mm) = highVar;
            end

        end
        [neuronFeatureData.sessionID] = deal(sessionID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add metadata
        % label the metadata for each neuron (core/shell and MUA/SUA)
        foo = num2cell(isCoreUnit); [neuronFeatureData.isCore] = foo{:};
        foo = num2cell(isMUAUnit);  [neuronFeatureData.isMUA ] = foo{:};
        [neuronFeatureData.age] = deal(thisAge);
        [neuronFeatureData.isPlastic] = deal(isThisPlastic);
        % stuff each session's data into a cell
        sessNeuronData{jj} = neuronFeatureData';
        [sessNeuronData{jj}.sessionID] = deal(sessionID);
    end
    
    % stuff all the sessions for each bird into another cell
    if ~isempty(sessNeuronData) && any(~cellfun('isempty',sessNeuronData))
        birdUnitCorrData{ii} = vertcat(sessNeuronData{:});
    end
    progressbar(ii/nBirds);
end
%%
features = fieldnames(featureTable);
save('data/allUnitFeatureCorr.mat','birdUnitCorrData','features');
%% cache some category
allUnits = vertcat(birdUnitCorrData{:});
cc = [allUnits.statFamily];
statNames = {allUnits(1).statFamily.feature};
unitsAreCore = [allUnits.isCore];
unitsAreSUA  = ~[allUnits.isMUA];
unitsAreSig = [allUnits.pResponsive] < (0.05 * length(allUnits));

directP = reshape([cc(:).P], size(cc)); % rows are statistic families, columns are units
partialP = reshape([cc(:).Ppartial], size(cc));
fullP = reshape([fullcc(:).pAll], size(fullcc)); %JMA added

readableStatNames = strrep(statNames, 'Freq', ' frequency');
readableStatNames = strrep(readableStatNames, 'Goodness', ' goodness');
readableStatNames = strrep(readableStatNames, 'rawAM', 'power derivative');
readableStatNames = strrep(readableStatNames, 'totalPower', 'power');
readableStatNames = strrep(readableStatNames, 'Entropy', ' entropy');
readableStatNames = strrep(readableStatNames, 'mFD', 'maximum FM within frequency');
readableStatNames = strrep(readableStatNames, 'mTD', 'maximum AM within frequency');

%{
isCore_F = @(x) [x.isCore];
isMUA_F  = @(x) [x.isMUA ];
isPlastic_F = @(x) [x.isPlastic];

isCoreAndSUA_F = @(x) isCore_F(x) & ~isMUA_F(x);
%}
%%
colMap = jet(numel(statNames));
label = {'CORE','SHELL'};

for hh = 0:1
    figure(hh+1);
    unitsLocation = ~unitsAreCore == hh; %JMA added
    unitsInclude = unitsLocation & unitsAreSUA; %JMA added
    nNeurons = sum(unitsInclude); %JMA changed to unitsInclude
    pSidak = (1 - (1 - 0.05)^(1/(numel(cc(unitsInclude))))); %JMA changed to number of pvalue for single units and core or shell
    %     hlin = semilogx([0.05, 0.05], [1 nNeurons], '--','Color', [0.5 0.5 0.5]);
    %     set(hlin, 'HandleVisibility', 'off');
    %     hold on;
    hlin = semilogx([pSidak pSidak], [1 nNeurons], '--','Color', [0.3 0.3 0.3]);
    set(hlin, 'HandleVisibility', 'off');
    hold on; %JMA added
    for ii = 1:numel(statNames)
        thisFeature = statNames{ii};
        pValsFeature = partialP(ii, unitsInclude); %JMA changed to number of pvalue for single units and core or shell, change to partialP instead of directP to get partial values from full model
        
        semilogx(sort(pValsFeature), 1:nNeurons, 'Color', colMap(ii,:), 'LineWidth', 2);
    end

    xlim([1e-5 1]); ylim([0 110])
    hold off;
    set(gca,'XDir', 'reverse');
    if hh == 0
        legend(readableStatNames, 'Location','North');
    end
    xlabel(sprintf('linear correlation p-values of features to RS in %s',label{hh+1}));
    ylabel('Neuron rank');
    set(gca,'Box','off');
    set(gcf,'Color',[1 1 1]);
    saveCurrFigure(sprintf('figures/paper/featureCorrelations/pvalueCDFsFeatureFamily-%s.jpg',label{hh+1}));
end
%% plot the frequency of family feature significance for all neurons
pSidak = (1 - (1 - 0.05)^(1/(numel(cc))));
pctSigDirect  = sum(directP  < pSidak, 2) / size( directP, 2);
pctSigPartial = sum(partialP < pSidak, 2) / size(partialP, 2);
pctSigBoth = sum(directP < pSidak & partialP < pSidak, 2) / size(partialP, 2);
isDirectSig = directP < pSidak;
isPartialSig = partialP < pSidak;
isBothSig = isDirectSig & isPartialSig;

%{
figure;
subplot(311); bar(pctSigDirect ); ylim([0 1]); set(gca,'XTickLabel', statNames);
xlabel('Direct correlation');
subplot(312); bar(pctSigPartial); ylim([0 1]); set(gca,'XTickLabel', statNames);
xlabel('Partial correlation');
subplot(313); bar(pctSigBoth   ); ylim([0 1]); set(gca,'XTickLabel', statNames);
xlabel('Both correlation');
%}
%
% plot of each individual p-value
figure;
nFamilies = size(cc,1);
nUnits = size(cc,2);
jitteredXCoords = (1:numel(statNames))' * ones(1, nUnits) + 0.05*randn(size(cc));
coreCoord = logical(ones(nFamilies,1) * unitsAreCore); %JMA moved parentheses to correct
semilogy(jitteredXCoords(coreCoord(:)), directP(coreCoord(:)), '.','Color',[0.5 0.5 0.5]);
hold on;
semilogy(jitteredXCoords(~coreCoord(:)), directP(~coreCoord(:)), '.','Color',[1 0 0]);
plot(xlim, [pSidak pSidak], 'r-');
set(gca,'XTick', 1:numel(statNames), 'XTickLabel', statNames);
ylabel('family linear correlation p-value'); xlabel('feature family');

title('p-values of significant correlations from feature families to RS');
saveCurrFigure('figures\sigFeatCorrByBird\allFamilyCorrPValues.jpg');
%%
figure;
subplot(2,1,1);
coreSig  = sum(directP(:,  unitsAreCore) < pSidak, 2);
shellSig = sum(directP(:, ~unitsAreCore) < pSidak, 2);
plotBarError([coreSig shellSig],[],[],[],[0.5 0.5 0.5; 1 0 0]);
ylabel('Count of significant correlations to feature families');
set(gca, 'XTickLabel', statNames);
subplot(2,1,2);
corePartSig  = sum(partialP(:,  unitsAreCore) < pSidak, 2);
shellPartSig = sum(partialP(:, ~unitsAreCore) < pSidak, 2);
plotBarError([corePartSig shellPartSig],[],[],[],[0.5 0.5 0.5; 1 0 0]);
ylabel('Count of significant part correlations to feature families');
set(gca, 'XTickLabel', statNames);
saveCurrFigure('figures\sigFeatCorrByBird\allFamilyCorrCount.jpg');
% very few significant results
%% plot response strength differences for core vs. shell single-unit populations
% plotting only occurs if the RS differences are significant
allDiffRS = vertcat(allUnits.diffRS);
unitsAreCore = [allUnits.isCore];
unitsAreSUA  = ~[allUnits.isMUA];
unitsAreSig = [allUnits.pResponsive] < 0.05;
unitsAreExcited = [allUnits.isExcited];
nFeatures = length(features);

diffRSContrasts = ...
    initEmptyStructArray({'pCoreShell','pCore', 'pShell', ...
    'rCoreShell', 'rCore', 'rShell', 'feature'}, 3*nFeatures); % r for rank
diffRSContrasts = reshape(diffRSContrasts, nFeatures, 3);
pSidak = (1 - (1 - 0.05)^(1/nFeatures));
prefix = {'Inhibited', 'Excited', 'All'};
dRSFeatCoreMean = zeros(nFeatures,3);
dRSFeatShellMean = zeros(nFeatures,3);
dRSFeatCoreSEM = zeros(nFeatures,3);
dRSFeatShellSEM = zeros(nFeatures,3);
for ii = 1:nFeatures
    diffRSContrasts(ii).feature = features{ii};
    for jj = 0:2
        if jj == 2
            diffRSByFeature_Core{ii,jj+1}  = allDiffRS(unitsAreSUA & unitsAreSig &  unitsAreCore, ii);
            diffRSByFeature_Shell{ii,jj+1} = allDiffRS(unitsAreSUA & unitsAreSig & ~unitsAreCore, ii);
        else
            diffRSByFeature_Core{ii,jj+1}  = allDiffRS(unitsAreSUA & unitsAreSig &  unitsAreCore & unitsAreExcited == jj, ii);
            diffRSByFeature_Shell{ii,jj+1} = allDiffRS(unitsAreSUA & unitsAreSig & ~unitsAreCore & unitsAreExcited == jj, ii);
        end
        diffRSByFeature_Core{ii,jj+1}(isnan(diffRSByFeature_Core{ii,jj+1} )) = [];
        diffRSByFeature_Shell{ii,jj+1}(isnan(diffRSByFeature_Shell{ii,jj+1})) = [];
        
        if isempty(diffRSByFeature_Core{ii,jj+1}) || isempty(diffRSByFeature_Shell{ii,jj+1}),
            fprintf('Skipping feature %s...\n', features{ii});
            continue;
        end
        
        [diffRSContrasts(ii,jj+1).pCoreShell,~,stats] = ...
            ranksum(diffRSByFeature_Core{ii,jj+1}(~isnan(diffRSByFeature_Core{ii,jj+1} )), ...
            diffRSByFeature_Shell{ii,jj+1}(~isnan(diffRSByFeature_Shell{ii,jj+1})));
        diffRSContrasts(ii,jj+1).rCoreShell = stats.ranksum;
        [diffRSContrasts(ii,jj+1).pCore,~,stats] = signtest(diffRSByFeature_Core{ii,jj+1});
        diffRSContrasts(ii,jj+1).rCore = stats.sign;
        [diffRSContrasts(ii,jj+1).pShell,~,stats] = signtest(diffRSByFeature_Shell{ii,jj+1});
        diffRSContrasts(ii,jj+1).rShell = stats.sign;
        if diffRSContrasts(ii,jj+1).pCoreShell < pSidak
            fprintf('%s, single-unit & responsive: Core vs shell difference in dRS across %s, MW U: p = %0.3f\n',...
                prefix{jj+1}, features{ii}, diffRSContrasts(ii,jj+1).pCoreShell);
            plotFeatRS = true;
        end
        if diffRSContrasts(ii,jj+1).pCore < pSidak
            fprintf('%s, single-unit & responsive: Core non-zero dRS across %s, sign stat %d, sign test: p = %0.3g\n',...
                prefix{jj+1}, features{ii}, diffRSContrasts(ii,jj+1).rCore, diffRSContrasts(ii,jj+1).pCore);
            plotFeatRS = true;
        end
        if diffRSContrasts(ii,jj+1).pShell < pSidak
            fprintf('%s, single-unit & responsive: Shell non-zero dRS across %s, sign stat %d, sign test: p = %0.3g\n',...
                prefix{jj+1}, features{ii}, diffRSContrasts(ii,jj+1).rShell, diffRSContrasts(ii,jj+1).pShell);
            plotFeatRS = true;            
        end
        nCore  = numel(diffRSByFeature_Core{ii,jj+1});
        nShell = numel(diffRSByFeature_Shell{ii,jj+1});
        dRSFeatCoreMean(ii,jj+1) = mean(diffRSByFeature_Core{ii,jj+1}); 
        dRSFeatCoreSEM(ii,jj+1) = std(diffRSByFeature_Core{ii,jj+1})/sqrt(nCore - 1);
        dRSFeatShellMean(ii,jj+1) = mean(diffRSByFeature_Shell{ii,jj+1});
        dRSFeatShellSEM(ii,jj+1) = std(diffRSByFeature_Shell{ii,jj+1})/sqrt(nShell - 1);
    end
end
%%
nFeatures = length(features);
familyAssign = zeros(nFeatures,1);
zeroedFeatures = dRSFeatCoreMean(:,3) == 0 & dRSFeatShellMean(:, 3)  == 0;

for ii = 1:numel(statNames)
    familyAssign(strncmp(statNames{ii}, features, length(statNames{ii}))) = ii;
end
colMap = jet(numel(statNames));
dR = 0.7;
plot([-dR dR NaN 0 0 NaN -dR dR],[0 0 NaN -dR dR NaN -dR dR],'k-','HandleVisibility','off');
hold on;
pSig = 0.05;
pSidak = (1 - (1 - 0.05)^(1/nFeatures)); 
%{
isSig = diffRSContrasts(ii).pCore < pSig | ...
        diffRSContrasts(ii).pShell< pSig | ...
        diffRSContrasts(ii).pCoreShell < pSig;
isMultSig = diffRSContrasts(ii).pCore  < pSidak | ...
            diffRSContrasts(ii).pShell < pSidak | ...
            diffRSContrasts(ii).pCoreShell < pSidak;
fprintf('%d features significant after multiple comparisons...\n',sum(isMultSig));
%}
doCross = false;    
for ii = 1:numel(statNames)    
    %doPlot = (isSig & familyAssign == ii);
    doPlot = (familyAssign == ii & ~zeroedFeatures);
    %fprintf('%d features of family %s significant once...\n', sum(doPlot), statNames{ii});
    
%    if ~any(doPlot), continue; end    
    if doCross
%        hlin = ploterr(dRSFeatCoreMean(doPlot, 3), dRSFeatShellMean(doPlot, 3),...
%            dRSFeatCoreSEM(doPlot,3), dRSFeatShellSEM(doPlot, 3),'.');
        hlin = plot(dRSFeatCoreMean(doPlot, 3), dRSFeatShellMean(doPlot, 3), '.',...
            dRSFeatCoreMean(doPlot, 3) * [1 1]+ dRSFeatCoreSEM(doPlot, 3) * [-1 1], dRSFeatShellMean(doPlot, 3) * [1 1], '-', ...
            dRSFeatCoreMean(doPlot, 3) * [1 1], dRSFeatShellMean(doPlot, 3) * [1 1] + dRSFeatShellSEM(doPlot,3) * [-1 1], '-')
        
        set(hlin, 'Color', colMap(ii,:));
        set(hlin, 'MarkerEdgeColor', colMap(ii,:));
        set(hlin, 'MarkerFaceColor', colMap(ii,:));
        set(hlin([2 3]), 'HandleVisibility', 'off');
    else
        plot(dRSFeatCoreMean(familyAssign == ii, 3), ...
            dRSFeatShellMean(familyAssign == ii, 3), ...
            '.','Color',colMap(ii,:), 'MarkerSize', 14);
    end
    hold on;
end
xlim([-dR dR]); ylim([-dR dR]);
xlabel('zRS differences between low and high feature quartile in CORE');
ylabel('zRS differences between low and high feature quartile in SHELL');
hleg = legend(readableStatNames,'Location','Southeast');
set(hleg, 'FontSize', 10);
set(gca,'Box','off');
set(gcf,'Color', [1 1 1]);
plot(-0.3692,0.6130,'k+'); % the average zRS to the centrally chosen syllable for core, shell
hold off;

appString = [];
if doCross, appString = '-SEM'; end
saveCurrFigure(sprintf('figures/paper/featureCorrelations/syllableFeatureRSDiffs%s.jpg',appString));
%% plot the histogram of distribution of each neuron's dRS to a given feature
for ii = 1:nFeatures
    for jj = 1:3
        if ~isempty(diffRSByFeature_Core{ii,jj}) && ~isempty(diffRSByFeature_Shell{ii,jj}) %JMA added
        
        figure
        nBins = 100;
        binMin = max(-5, min([diffRSByFeature_Core{ii,jj}; diffRSByFeature_Shell{ii,jj}]));
        binMax = min( 5, max([diffRSByFeature_Core{ii,jj}; diffRSByFeature_Shell{ii,jj}]));
        bins = linspace(binMin, binMax, nBins);
        plotInterlaceBars(diffRSByFeature_Core{ii,jj}, diffRSByFeature_Shell{ii,jj}, bins);
        ytop = max(ylim);
        
        % todo maybe: use median/jackknife estimates of variance of median?
        hold on;
        plotHorzErrorBar(dRSFeatCoreMean(ii,jj), ytop, dRSFeatCoreSEM(ii,jj), [0.5 0.5 0.5]);
        plotHorzErrorBar(dRSFeatShellMean(ii,jj), ytop + 1, dRSFeatShellSEM(ii,jj), [1 0 0]);
        plot([0 0], ylim, 'k-');
        xlabel(sprintf('diff RS between low - high on %s',diffRSContrasts(ii).feature),...
            'Interpreter','none');
        ylabel('count');
        title(sprintf('%s units: Core/Shell MW U: p = %0.2g, Core p = %0.2g, Shell p = %0.2g',...
            prefix{jj}, diffRSContrasts(ii).pCoreShell, diffRSContrasts(ii).pCore, diffRSContrasts(ii).pShell))
        
        saveCurrFigure(sprintf('figures/featureRSDifferences/%sNeuronsRSDiff-%s.jpg', ...
           prefix{jj}, diffRSContrasts(ii).feature));
        hold off;
        close;
        else continue 
        end
       
    end
    % none of these are significant after multiple comparisons...
end

save('data\RSFeatureContrasts.mat', 'diffRSContrasts');
%% plotting for each bird

theFeatures = zeros(nBirds, length(features));
isCore    = [allUnits.isCore   ]; %JMA changed from birdUnitCorrData
isMUA     = [allUnits.isMUA    ];
isPlastic = [allUnits.isPlastic];

isCoreAndSUA  =  isCore & ~isMUA;
isShellAndSUA = ~isCore & ~isMUA;
isCoreAndMUA  =  isCore &  isMUA;
isShellAndMUA = ~isCore &  isMUA;

%% plot the percent of neurons whose FR correlates to different features-need to fix something in this section
for ii = 1:nBirds
    birdID{ii} = strtok(birdUnitCorrData{ii}(1).sessionID, '_');
end
label = {'CORE','SHELL'};
nFeatures = numel(features);
%nUnitsSigFeature
%pctUnitsSigMultFeature = zeros(nFeatures, nBirds);
pNormal = 0.05;
pSidak = (1 - (1 - pNormal)^(1/nFeatures));
nUnits = zeros(2, nBirds); % first dim is core/shell
nSigFeatures = cell(2, nBirds); % first dim is core/shell
nPartialSigFeatures = cell(2, nBirds); % first dim is core/shell
nSigUnits = zeros(nFeatures, nBirds, 2);% last  dim is core/shell
nSigPartialUnits = zeros(nFeatures, nBirds, 2);
for ii = 1:nBirds
    thisBird = birdUnitCorrData{ii}; % all neurons for a given bird
    isCore =  [thisBird.isCore];
    isSUA  = ~[thisBird.isMUA ];
    for jj = 0:1 % core/shell
        % the p values for each feature and each neuron
        pValFeatures = vertcat(thisBird(isCore == jj & isSUA).pCorr)';
        pPartialValFeatures = vertcat(thisBird(isCore == jj & isSUA).pPartial)';
        nSUAUnits = size(pValFeatures, 1);
        if nSUAUnits > 0
            
            whichSingleSig  = abs(       pValFeatures) <= pSidak;
            whichPartialSig = abs(pPartialValFeatures) <= pSidak;
            % fraction of units that show significant modulation by a given features,
            % for multiple comparisons
            
            % raw correlations uncorrected for independence
            nSigUnits(:,ii, jj+1) = sum(whichSingleSig ,2);
            nSigPartialUnits(:,ii, jj+1) = sum(whichPartialSig,2);
            
            nSigFeatures{jj+1, ii} = sum(whichSingleSig, 1);
            nPartialSigFeatures{jj+1, ii} = sum(whichPartialSig, 1);
            nUnits(jj+1, ii) = nSUAUnits;
        else
            fprintf('no %s SUA neurons for bird %s...\n', label{jj+1}, birdID{ii});
        end
    end
end
for jj = 1:2
    figure;
    nSigFeatsPerNeuron = horzcat(nSigFeatures{jj,:});
    maxSigFeats = max(nSigFeatsPerNeuron);
    hist(nSigFeatsPerNeuron, 0:maxSigFeats);
    xlabel(sprintf('# features significantly correlated to %s neurons', label{jj}));
    ylabel(sprintf('count of single unit %s neurons',label{jj}));
    saveCurrFigure(sprintf('figures/sigFeatCorrByBird/%s-directFRCorrFeatureCt.jpg', label{jj}));
    figure;
    nSigFeatsPerNeuron = horzcat(nPartialSigFeatures{jj,:});
    maxSigFeats = max(nSigFeatsPerNeuron);
    hist(nSigFeatsPerNeuron, 0:maxSigFeats);
    xlabel(sprintf('# features significantly part-correlated to %s neurons', label{jj}));
    ylabel(sprintf('count of single unit %s neurons',label{jj}));
    saveCurrFigure(sprintf('figures/sigFeatCorrByBird/%s-partialFRCorrFeatureCt.jpg', label{jj}));
    
    figure;
    barh(1:nFeatures, sum(nSigUnits(:,:,jj), 2));
    set(gca,'YTick', 1:nFeatures,'YTickLabel', features, 'xlim', [0 1], 'yLim', [0 length(features)]);
    xlabel('# neurons responsive to this feature');
    title(sprintf('single unit FR single correlations in %s to features,' ,label{jj}))
    saveCurrFigure(sprintf('figures/sigFeatCorrByBird/%s-directFRCorrNeuronCnt.jpg', label{jj}));
    
    figure;
    barh(1:nFeatures, sum(nSigPartialUnits(:,:,jj), 2));
    set(gca,'YTick', 1:nFeatures,'YTickLabel', features, 'xlim', [0 1], 'yLim', [0 length(features)]);
    xlabel('# neurons responsive to this feature');
    title(sprintf('single unit FR partial correlations in %s to features',label{jj}));
    saveCurrFigure(sprintf('figures/sigFeatCorrByBird/%s-partialFRCorrNeuronCnt.jpg', label{jj}));
    
    
end
%% plotting across birds
acrossFeatures = mean(theFeatures,1);
for ff = (1:10:120);
    figure;
    bar(1:10,acrossFeatures(ff:ff+9));
    set (gca,'xTickLabel', features(ff:ff+9), 'xLim', [0 10]);
    title('Mean proportion of  core neurons correlated to features');
end
%%

%% find example neuron
%{
mysteryField = 'totalPower_release';
assert(any(strcmpi(mysteryField,fieldnames(featureTable))));
mFIdx = find(strcmpi(mysteryField,fieldnames(featureTable)));

for ii = 1:numel(birdNeuronData) % loop through the birds
    % to find the best correlated neuron in each one
    nSessions = numel(birdNeuronData{ii});
        
    corrVals = horzcat(birdNeuronData{ii}(:).correlates);
    [bestCorrelationVal, bestNeuronNum] = min(abs(corrVals(mFIdx,:)));
    
    % match the correct neuron to the correct session (because all the
    % neurons embedded in sections are smushed into one bird
    
    birdDataDir = [pwd filesep 'data' filesep birdID{ii} filesep]; %Jenny added because old def. was in a diff loop
    load([birdDataDir 'allSpecs-' birdID{ii} '.mat']); %Jenny added because old def. was in diff loop

    % converts best neuron number (in the whole bird) to best neuron number
    % (within its session)
    justBird = horzcat(birdNeuronData{ii}(:));
    bestSessID = justBird(bestNeuronNum).sessionID;
    bestNeuronNum = find(bestNeuronNum == find(strcmp(bestSessID, {justBird.sessionID})));
    
    % find the syllables that match that session
    matchingFileA = [birdDataDir bestSessID '.mat'];
    matchingFileB = [birdDataDir filesep bestSessID '.mat'];
    matchesFile = strcmpi(matchingFileA, {DRsylls.file}) | strcmpi(matchingFileB, {DRsylls.file});
    mysteryFeatures = [featureTable(matchesFile).(mysteryField)];

    % get the firing rates for that special neuron
    corrSylls = DRsylls(matchesFile);
    corrSpikes = loadSpikeData(compositeReport{ii}(strcmp(bestSessID, {compositeReport{ii}.sessionID})).spikeFiles);
    specialNeuronSpikes = corrSpikes{bestNeuronNum};
    
    bestRates = countSpikes(corrSylls, specialNeuronSpikes); %Jenny added the extra ~,
    
    figure;
        plot(bestRates, mysteryFeatures, 'r.');
    xlabel('Firing Rate (Hz)')
    ylabel(mysteryField);
    title(birdID{ii});
    hold off
end
%}
%% plot all neurons
%{
mysteryField = 'wienerEntropy_max';
assert(any(strcmpi(mysteryField,fieldnames(featureTable))));
mFIdx = find(strcmpi(mysteryField,fieldnames(featureTable)));

for ii = 1:numel(birdUnitCorrData) % loop through the birds
    % to find the best correlated neuron in each one
    nSessions = numel(birdUnitCorrData{ii});
        
    corrVals = horzcat(birdUnitCorrData{ii}(:).correlates);
    [bestCorrelationVal, bestNeuronNum] = min(abs(corrVals(mFIdx,:)));
    
    % match the correct neuron to the correct session (because all the
    % neurons embedded in sections are smushed into one bird
    
    birdDataDir = [pwd filesep 'data' filesep birdID{ii} filesep];
    load([birdDataDir 'allSpecs-' birdID{ii} '.mat']);

    % converts best neuron number (in the whole bird) to best neuron number
    % (within its session)
    justBird = horzcat(birdUnitCorrData{ii}(:));
    %run through all the sessions
    bestSessID = justBird(bestNeuronNum).sessionID;
    bestNeuronNum = find(bestNeuronNum == find(strcmp(bestSessID, {justBird.sessionID})));
    
    % find the syllables that match that session
    matchingFileA = [birdDataDir bestSessID '.mat'];
    matchingFileB = [birdDataDir filesep bestSessID '.mat'];
    matchesFile = strcmpi(matchingFileA, {DRsylls.file}) | strcmpi(matchingFileB, {DRsylls.file});
    mysteryFeatures = [featureTable(matchesFile).(mysteryField)];

    % get the firing rates for that special neuron
    corrSylls = DRsylls(matchesFile);
    corrSpikes = loadSpikeData(compositeReport{ii}(strcmp(bestSessID, {compositeReport{ii}.sessionID})).spikeFiles);
    specialNeuronSpikes = corrSpikes{bestNeuronNum};
    
    bestRates = countSpikes(corrSylls, specialNeuronSpikes);
    
    figure;
    plot(bestRates, mysteryFeatures, 'r.');
    xlabel('Firing Rate (Hz)')
    ylabel(mysteryField);
    title(birdID{ii});
    hold off
end
%}
%%
profile off
profile viewer