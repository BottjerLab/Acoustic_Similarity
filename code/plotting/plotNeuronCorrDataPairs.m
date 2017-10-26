function plotNeuronCorrData(allNeuronCorrData, params, varargin)
if nargin < 2 || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});

% plot difference between firing rates for near tutor/far from tutor
% and also p-values for correlations between neurons and firing rates

% this data is compiled in correlateDistanceToFiring
if nargin < 1 || isempty(allNeuronCorrData)
    load('data/allNeuronCorrelations.mat');
end

%% here we load the cluster quality
% get birds and ages first
sessionIDs = {allNeuronCorrData.sessionID};
birdIDs = strtok(sessionIDs, '_');
[uSessions, ~, rIdxSession] = unique(sessionIDs);  % index through ages can go back to sessions
uAges = getAgeOfSession(uSessions);
sessionAges = zeros(size(sessionIDs));
for ii = 1:numel(uAges)
    sessionAges(rIdxSession == ii) = uAges(ii);
end
[sessionQ   , allSubj]  = getClusterQuality(birdIDs, sessionAges, [allNeuronCorrData.syllID]);
[sessionObjQ, allObj ]  = getClusterQuality(birdIDs, sessionAges, [allNeuronCorrData.syllID], true);
foo = num2cell(sessionQ   ); [allNeuronCorrData.clusterQ   ] = foo{:};
foo = num2cell(sessionObjQ); [allNeuronCorrData.clusterObjQ] = foo{:};
foo = allSubj'; allSubj = [allSubj(:)];
foo = allObj' ; allObj  = [allObj(:) ];

missingData = isnan(allSubj) | isnan(allObj);
qualityFit = polyfit(allSubj(~missingData), allObj(~missingData), 1);
plot(1:5, polyval(qualityFit,1:5),'r-');
hold on;
%boxplot(sessionObjQ', sessionQ', 'notch', 'on');
plot(allSubj, allObj, 'k.');

[rQualCorr, pQualCorr] = corrcoef(allSubj(~missingData), allObj(~missingData));
legend(sprintf('r^2 = %0.3f, p = %0.3g', rQualCorr(2,1), pQualCorr(2,1)));
xlim([0.5 5.5])
xlabel('Subjective Cluster Quality');
ylabel('Davies-Bouldin Index');
title('Subjective vs. objective cluster quality correlations');
if params.saveplot
    saveCurrFigure('figures\A_keeper\objSubjClusterQuality.jpg');
end
%% flags
isCore        = [allNeuronCorrData.isCore];
isMUA         = [allNeuronCorrData.isMUA];
isPlastic     = [allNeuronCorrData.isPlastic];
isSignificant = [allNeuronCorrData.sigResponse];
nSylls        = [allNeuronCorrData.nSylls];
%isSignificant = true(1,numel(allNeuronCorrData));
isExcited     = [allNeuronCorrData.isExcited];
%%
isSubjGood    = [allNeuronCorrData.clusterQ]  < 1.5; % < 2.5
isObjGood     = [allNeuronCorrData.clusterObjQ] < 0.8; % < 1
%%
% criteria for cluster inclusion
% must have at least three syllables per quartile
isPresel = ~isMUA & nSylls >= 40 & isObjGood; %JMA
%isPresel = ~isMUA & isObjGood; %~isMUA & isObjGood
%%isPresel = isSignificant; % for mostresponsive neurons in todo0410.m
%isPresel = nSylls >= 12;

% subplot rows / columns
nR = 3; nC = 2;
% measures
distanceTypes = {'tutor', 'intra', 'inter', 'consensus', 'central','humanMatch'}';
distanceDescriptions = {'closest tutor', 'cluster center', 'normed center', ...
    'closest tutor to cluster consensus', 'closest tutor to cluster center', 'expert-designated tutor'};
dFieldsRS  = [strcat(distanceTypes, '_nearMeanRS') strcat(distanceTypes, '_farMeanRS')];
dFieldsSEM = [strcat(distanceTypes, '_nearMeanSEM') strcat(distanceTypes, '_farMeanSEM')];
eiTitle = {'Significantly inhibited single unit-syllable pairs', ...
    'Significantly excited single unit-syllable pairs', ...
    'All significant single unit-syllable pairs'};
xlabels = strcat({'RS for near - far to '}, distanceDescriptions);
filsuff = {'inh','exc','all'};
%%
%{
for hh = 1:3 % inhibited, excited, all
    for ii = 1:numel(distanceTypes) % six of them
        figure;
        diffTutorMeanRS = [allNeuronCorrData.(dFieldsRS{ii,1})] - [allNeuronCorrData.(dFieldsRS{ii,2})];
        nearMeanRS  = [allNeuronCorrData.(dFieldsRS{ii,1})];
         farMeanRS  = [allNeuronCorrData.(dFieldsRS{ii,2})];
        nearMeanSEM = [allNeuronCorrData.(dFieldsRS{ii,1})];
         farMeanSEM = [allNeuronCorrData.(dFieldsRS{ii,2})];
        if hh < 3
            selHereCore  = isPresel & isExcited == hh-1 &  isCore;
            selHereShell = isPresel & isExcited == hh-1 & ~isCore;
             coreDiffRS = diffTutorMeanRS(selHereCore);
            shellDiffRS = diffTutorMeanRS(selHereShell);
            
             coreSubDiffRS = diffTutorMeanRS(selHereCore  & ~isPlastic);
            shellSubDiffRS = diffTutorMeanRS(selHereShell & ~isPlastic);
            
             corePlastDiffRS = diffTutorMeanRS(selHereCore & isPlastic);
            shellPlastDiffRS = diffTutorMeanRS(selHereShell & isPlastic);
        else
            coreDiffRS  = diffTutorMeanRS(isPresel &  isCore);
            shellDiffRS = diffTutorMeanRS(isPresel & ~isCore);

             coreSubDiffRS = diffTutorMeanRS(isPresel &  isCore & ~isPlastic);
            shellSubDiffRS = diffTutorMeanRS(isPresel & ~isCore & ~isPlastic);
            
             corePlastDiffRS = diffTutorMeanRS(isPresel &  isCore & isPlastic);
            shellPlastDiffRS = diffTutorMeanRS(isPresel & ~isCore & isPlastic);
        end
        pc = signrank( coreDiffRS);
        ps = signrank(shellDiffRS);
        pMannU = ranksum(coreDiffRS(~isnan(coreDiffRS)), shellDiffRS(~isnan(shellDiffRS)));
        fprintf(['%s-%s:\n\tsign-rank test p-value for core: %0.3f' ...
                 '      \n\tsign-rank test p-value for shell: %0.3f',...
                 '      \n\tMann-Whitney U test p-value for core v shell: %0.3f\n'],...
            eiTitle{hh},xlabels{ii},pc,ps,pMannU);

        % clunky way just to get the top histogram value
        RSdiffBins = -10:0.2:10;
        plotInterlaceBars(coreDiffRS, shellDiffRS, RSdiffBins);
        ytop = ylim * [0 1]';
        
        % todo: plot significance on graph
        hold on;
        plotSEMBar(      coreDiffRS, ytop  , [0.5 0.5 0.5]);
        plotSEMBar(   coreSubDiffRS, ytop+1, [0.5 0.5 0.5]);
        plotSEMBar( corePlastDiffRS, ytop+2, [0.5 0.5 0.5]);
        plotSEMBar(     shellDiffRS, ytop+3, [  1   0   0]);
        plotSEMBar(  shellSubDiffRS, ytop+4, [  1   0   0]);
        plotSEMBar(shellPlastDiffRS, ytop+5, [  1   0   0]);
        plot([0 0], ylim, 'k--');
        hold off;
        % redo y axis labels
        yt = get(gca,'YTick');
        yt = [yt(yt < ytop) ytop:ytop+5];
        ytl = cellfun(@(x) sprintf('%d',x),num2cell(yt),'UniformOutput',false);
        ytl(end-5:end) = {'Core','Core/Subsong','Core/Plastic','Shell','Shell/Subsong','Shell/Plastic'};
        set(gca,'YTick',yt,'YTickLabel',ytl);
                
        % figure formatting
        xlabel(xlabels{ii});
        ylabel('Count');
        xlim([min(RSdiffBins) max(RSdiffBins)]);
        set(gca,'Box','off');
        set(gca, 'FontSize', 14);
        set(get(gca,'XLabel'),'FontSize', 14);
        set(get(gca,'YLabel'),'FontSize', 14);
        set(get(gca,'Title' ),'FontSize', 14);
        title(eiTitle{hh});
        set(gcf,'Color',[1 1 1]);
        
        if params.saveplot
            saveCurrFigure(sprintf('figures/distanceCorrelations/RSdiffs-SUA-%s-%s.jpg', distanceTypes{ii}, filsuff{hh}));
        end
    end
end
%}

%% JMA added this section to compare neurons with compiled cluster data

for aa = 1: length(allNeuronCorrData)
    allNeuronCorrData(aa).intra_DistanceAll = allNeuronCorrData(aa).intra_DistanceAll';%these were in columns
    allNeuronCorrData(aa).inter_DistanceAll = allNeuronCorrData(aa).inter_DistanceAll';
    allNeuronCorrData(aa).burstFraction = allNeuronCorrData(aa).burstFraction';
    allNeuronCorrData(aa).quality = allNeuronCorrData(aa).clusterObjQ;
end
ccc = allNeuronCorrData(1);%trying to get unique syllable iterations to compare similarity to tutor with song stage goodness-of-fit coefficient
for aaa = 2: length(allNeuronCorrData)
    ddd = strcmp(allNeuronCorrData(aaa).sessionID,{ccc.sessionID});
    if ~any(ddd)
        ccc = [ccc;allNeuronCorrData(aaa)];%I know, I know, I can't pre-allocate
    else
        eee = ccc(ddd);
        bbb = ~ismember(allNeuronCorrData(aaa).syllID,[eee.syllID]);
        if bbb
            ccc = [ccc;allNeuronCorrData(aaa)];
        end
    end
end
load('goodnessOfFitCo.mat')
for aa = 1: length(ccc)
   tt = zeros(ccc(aa).nSylls,1);%needed to get syllable type for every syllable
   tt(:,1) = deal(ccc(aa).syllID);
   ccc(aa).syllID = tt';
   for bb = 1: length(goodness.fit)
   sess(bb) = strcmp(ccc(aa).sessionID,goodness.sessionID(bb));
   end
   gof = cell2mat(goodness.fit(sess));
   ttt = zeros(ccc(aa).nSylls,1);
   ttt(:,1) = deal(gof);
   ccc(aa).indGoF = ttt';
end

usablePairs = allNeuronCorrData(isPresel);
% usableClusterSessions = {usablePairs.sessionID};
% [uUCSessions, ~, ~] = unique(usableClusterSessions);
% corrByNeuron = struct([]);
% for mm = 1: length(uUCSessions)
%     isCurrentSession = strcmp(uUCSessions(mm),{usablePairs.sessionID});
%     currentSessionPairs = usablePairs(isCurrentSession);
% %     [neuronsHere, ~, ~] = unique([currentSessionPairs.unitNum]);
% %     for nn = 1: length(neuronsHere)
% %         isCurrentNeuron = [currentSessionPairs.unitNum] == neuronsHere(nn);
% %         currentNeuronPairs = currentSessionPairs(isCurrentNeuron);
% %         compiledNeuron.isCore = currentNeuronPairs(1).isCore;
% %         compiledNeuron.isMUA = currentNeuronPairs(1).isMUA;
% %         compiledNeuron.isPlastic = currentNeuronPairs(1).isPlastic;
% %         compiledNeuron.sessionID = currentNeuronPairs(1).sessionID;
% %         compiledNeuron.unitNum = currentNeuronPairs(1).unitNum;
% %         compiledNeuron.nSylls = sum([currentNeuronPairs.nSylls]);
% %         compiledNeuron.sigSyll = sum([currentNeuronPairs.sigResponse]);
% %         compiledNeuron.RSAll = horzcat([currentNeuronPairs.RSAll]);
% %         compiledNeuron.FRSyll = horzcat([currentNeuronPairs.FRSyll]);
% %         compiledNeuron.FRBase = horzcat([currentNeuronPairs.FRBase]);
% %         compiledNeuron.burstFraction = horzcat([currentNeuronPairs.burstFraction]);
% %         compiledNeuron.tutor_DistanceAll = horzcat([currentNeuronPairs.tutor_DistanceAll]);
% %         compiledNeuron.consensus_DistanceAll = horzcat([currentNeuronPairs.consensus_DistanceAll]);
% %         compiledNeuron.central_DistanceAll = horzcat([currentNeuronPairs.central_DistanceAll]);
% %         compiledNeuron.intra_DistanceAll = horzcat([currentNeuronPairs.intra_DistanceAll]);
% %         compiledNeuron.inter_DistanceAll = horzcat([currentNeuronPairs.inter_DistanceAll]);
% %         compiledNeuron.quality = horzcat([currentNeuronPairs.clusterObjQ]);
% %         compiledNeuron.syllID = horzcat([currentNeuronPairs.syllID]);
% %         numClassSyll = length(unique(compiledNeuron.syllID));
% %         compiledNeuron.classSyll = deal(numClassSyll);
%         corrByNeuron = [corrByNeuron; compiledNeuron];
% %     end
% end
% isenoughSyll = [corrByNeuron.nSylls] > 39; %want at least 10 syllables in each quartile
% usableNeuron = [corrByNeuron.sigSyll] > 0 & isenoughSyll; %neuron has to respond to at least one syllable cluster (but maybe shouldn't do this)
% usableNeuron = isenoughSyll; 
% corrByNeuron = corrByNeuron(usableNeuron);
corrByNeuron = usablePairs;

%correlation of distance to response strength
for bb = 1: length(corrByNeuron)
    yDist = corrByNeuron(bb).tutor_DistanceAll'; %can try other distances
    xRS = corrByNeuron(bb).RSAll'; %response strength, does it make sense to use firing rate?
[linfit, ~,~,~, fitStats] = regress(yDist, [ones(numel(xRS),1) xRS]); %checked with corrcoef and gives same p value
%if params.plot
% figure
% plot(xRS, yDist, 'k.', 'HandleVisibility', 'off');
% hold on;
% plot(xRS, linfit(1) + xRS * linfit(2), '--','Color',[1 0 0]);
% legend(sprintf('r^2 = %0.3g, F = %0.3g, p = %0.3g\n',...
%    fitStats(1), fitStats(2),fitStats(3)));
% xlabel('Response Strength'); ylabel('Matched distance');
%end
corrByNeuron(bb).linfit = linfit;
corrByNeuron(bb).fitStats = fitStats;
end

isCore2 = [corrByNeuron.isCore];
CorrPRS = zeros(length(corrByNeuron),1);
for cc = 1: length(corrByNeuron)
    CorrPRS(cc) = corrByNeuron(cc).fitStats(3);
end
isCorrel = CorrPRS < 0.05;
cSC = isCore2 & isCorrel';
cSS = ~isCore2 & isCorrel';
fprintf('Number neurons with significant correlation of RS in core %s out of %s core neurons \n', num2str(sum(cSC)), num2str(sum(isCore2)))
fprintf('Number neurons with significant correlation of RS in shell %s out of %s shell neurons \n', num2str(sum(cSS)), num2str(sum(~isCore2)))

%correlation of distance to burst fraction
for bb = 1: length(corrByNeuron)
    yDist = corrByNeuron(bb).tutor_DistanceAll'; 
    xRS = corrByNeuron(bb).burstFraction';
[linfit, ~,~,~, fitStats] = regress(yDist, [ones(numel(xRS),1) xRS]); %checked with corrcoef and gives same p value
corrByNeuron(bb).linfitBF = linfit;
corrByNeuron(bb).fitStatsBF = fitStats;
end
isCore2 = [corrByNeuron.isCore];
CorrP = zeros(length(corrByNeuron),1);
for cc = 1: length(corrByNeuron)
    CorrP(cc) = corrByNeuron(cc).fitStatsBF(3);
end
isCorrel = CorrP < 0.05;
cSC = isCore2 & isCorrel';
cSS = ~isCore2 & isCorrel';
fprintf('Number neurons with significant correlation of BF in core %s out of %s core neurons \n', num2str(sum(cSC)), num2str(sum(isCore2)))
fprintf('Number neurons with significant correlation of BF in shell %s out of %s shell neurons \n', num2str(sum(cSS)), num2str(sum(~isCore2)))

%compare population response, standardized response to top and bottom 50%
%similarity to tutor song
for dd = 1: length(corrByNeuron)
    dists = corrByNeuron(dd).tutor_DistanceAll; %tutor_DistanceAll
    iqDists = prctile(dists, 50);
    near_medianFR = corrByNeuron(dd).FRSyll(dists < iqDists); nMedian = numel(near_medianFR);
    far_medianFR  = corrByNeuron(dd).FRSyll(dists > iqDists); fMedian = numel( far_medianFR);
    near_medianFRBase = corrByNeuron(dd).FRBase(dists < iqDists);
    far_medianFRBase  = corrByNeuron(dd).FRBase(dists > iqDists);
    near_medianBF = corrByNeuron(dd).burstFraction(dists < iqDists);
    far_medianBF = corrByNeuron(dd).burstFraction(dists > iqDists);
    near_medianDist = corrByNeuron(dd).tutor_DistanceAll(dists < iqDists);
    far_medianDist = corrByNeuron(dd).tutor_DistanceAll(dists > iqDists);
    allDist = corrByNeuron(dd).tutor_DistanceAll;
%     near_quartileQual = corrByNeuron(dd).quality(dists < iqDists(1));
%     far_quartileQual = corrByNeuron(dd).quality(dists > iqDists(2));
%     near_quartileType = corrByNeuron(dd).syllID(dists < iqDists(1));
%     far_quartileType = corrByNeuron(dd).syllID(dists > iqDists(2));
    farMedian = nanmean(far_medianDist);
    nearMedian = nanmean(near_medianDist);
    meanDist = nanmean(allDist);
%     farQual = nanmean(far_quartileQual);
%     nearQual = nanmean(near_quartileQual);
%     farTypes = length(unique(far_quartileType));
%     nearTypes = length(unique(near_quartileType));
    nq_meanFR = nanmean(near_medianFR); nq_varFR = nanvar(near_medianFR);
    fq_meanFR = nanmean( far_medianFR); fq_varFR = nanvar( far_medianFR);
    nq_meanBF = nanmean(near_medianBF);
    fq_meanBF = nanmean(far_medianBF);
    nq_CV = nanstd(near_medianFR)/nq_meanFR;
    fq_CV = nanstd(far_medianFR)/fq_meanFR;
    nq_meanFRBase = nanmean(near_medianFRBase); nq_varFRBase = nanvar(near_medianFRBase);
    fq_meanFRBase = nanmean( far_medianFRBase); fq_varFRBase = nanvar( far_medianFRBase);
    [~, pVal, ~, tValStruct] = ttest2(near_medianFR- near_medianFRBase, far_medianFR- far_medianFRBase);
    tVal = tValStruct.tstat;
    fcovar = nancov(far_medianFR,far_medianFRBase);if numel(fcovar) > 1, fcovar = fcovar(2,1); end;
    ncovar = nancov(near_medianFR,near_medianFRBase);if numel(ncovar) > 1, ncovar = ncovar(2,1); end;
    farZDenom = sqrt(fq_varFR + fq_varFRBase - 2*fcovar);
    nearZDenom = sqrt(nq_varFR + nq_varFRBase - 2*ncovar);
    farZ = ((fq_meanFR - fq_meanFRBase)* sqrt(fMedian))/farZDenom;
    nearZ = ((nq_meanFR - nq_meanFRBase)* sqrt(nMedian))/nearZDenom;
    corrByNeuron(dd).medianPValue = pVal;
    corrByNeuron(dd).farZ = farZ;
    corrByNeuron(dd).nearZ = nearZ;
    corrByNeuron(dd).farCV = fq_CV;
    corrByNeuron(dd).nearCV = nq_CV;
    corrByNeuron(dd).nearBF = nq_meanBF;
    corrByNeuron(dd).farBF = fq_meanBF;
    corrByNeuron(dd).farQuart = farMedian;
    corrByNeuron(dd).nearQuart = nearMedian;
    corrByNeuron(dd).meanDist = meanDist;
    corrByNeuron(dd).farRS = (fq_meanFR - fq_meanFRBase);
    corrByNeuron(dd).nearRS = (nq_meanFR - nq_meanFRBase);
%     corrByNeuron(dd).farQual = farQual;
%     corrByNeuron(dd).nearQual = nearQual;
%     corrByNeuron(dd).farTypes = farTypes;
%     corrByNeuron(dd).nearTypes = nearTypes;
%     corrByNeuron(dd).farTypeID = {unique(far_quartileType)};
%     corrByNeuron(dd).nearTypeID = {unique(near_quartileType)};
end

isQS = [corrByNeuron.medianPValue] < 0.05;
SQC = isCore2 & isQS;
SQS = ~isCore2 & isQS;
fprintf('Number neurons with significant difference in RS to near vs far in core %s out of %s core neurons \n', num2str(sum(SQC)), num2str(sum(isCore2)))
fprintf('Number neurons with significant difference in RS to near vs far in shell %s out of %s shell neurons \n', num2str(sum(SQS)), num2str(sum(~isCore2)))

meanFarZC = nanmean([corrByNeuron(isCore2).farZ]);
meanFarZS = nanmean([corrByNeuron(~isCore2).farZ]);
SEMFarC = nanstd([corrByNeuron(isCore2).farZ])/sqrt(length(corrByNeuron(isCore2)));
SEMFarS = nanstd([corrByNeuron(~isCore2).farZ])/sqrt(length(corrByNeuron(~isCore2)));
meanNearZC = nanmean([corrByNeuron(isCore2).nearZ]);
meanNearZS = nanmean([corrByNeuron(~isCore2).nearZ]);
SEMNearC = nanstd([corrByNeuron(isCore2).nearZ])/sqrt(length(corrByNeuron(isCore2)));
SEMFarS = nanstd([corrByNeuron(~isCore2).nearZ])/sqrt(length(corrByNeuron(~isCore2)));


% normalized RS values--not using this part
dNormedRSFields = strcat(distanceTypes, '_dRSnorm');
for hh = 1:3 % inhibited, excited, all
    for ii = 1:numel(distanceTypes) % six of them
        figure;
        dRSNormed = [allNeuronCorrData.(dNormedRSFields{ii})];
        if hh < 3
            coreDiffRSNorm = dRSNormed(isPresel & isExcited == hh-1 &  isCore);
            shellDiffRSNorm = dRSNormed(isPresel & isExcited == hh-1 & ~isCore);
            
            coreSubDiffRSNorm = dRSNormed(isPresel & isExcited == hh-1 &  isCore & ~isPlastic);
            shellSubDiffRSNorm = dRSNormed(isPresel & isExcited == hh-1 & ~isCore & ~isPlastic);
            
            corePlastDiffRSNorm = dRSNormed(isPresel & isExcited == hh-1 &  isCore & isPlastic);
            shellPlastDiffRSNorm = dRSNormed(isPresel & isExcited == hh-1 & ~isCore & isPlastic);
        else
            coreDiffRSNorm  = dRSNormed(isPresel &  isCore);
            shellDiffRSNorm = dRSNormed(isPresel & ~isCore);
            
            coreSubDiffRSNorm = dRSNormed(isPresel &  isCore & ~isPlastic);
            shellSubDiffRSNorm = dRSNormed(isPresel & ~isCore & ~isPlastic);
            
            corePlastDiffRSNorm = dRSNormed(isPresel &  isCore & isPlastic);
            shellPlastDiffRSNorm = dRSNormed(isPresel & ~isCore & isPlastic);
        end
        fprintf('%s-normed %s: ', eiTitle{hh},xlabels{ii});
        if ~(all(isnan(   coreSubDiffRSNorm)) || ...
                all(isnan( corePlastDiffRSNorm)) || ...
                all(isnan(  shellSubDiffRSNorm)) || ...
                all(isnan(shellPlastDiffRSNorm)))
            
            % significance tests: two-way anova, permuted
            % this function is not consistent with matlab's anovan, so it
            % won't be used until we can see why the inconsistency's there
            %[stats, df, pvals] = statcond(...
            %    {noNaN( coreSubDiffRSNorm), noNaN( corePlastDiffRSNorm); ...
            %     noNaN(shellSubDiffRSNorm), noNaN(shellPlastDiffRSNorm)}, ...
            %'mode','param');
            % test against anova
            
            xx = [coreSubDiffRSNorm corePlastDiffRSNorm shellSubDiffRSNorm shellPlastDiffRSNorm]';
            grps = [zeros(size(coreSubDiffRSNorm)) zeros(size(corePlastDiffRSNorm)) ...
                ones(size(shellSubDiffRSNorm)) ones(size(shellPlastDiffRSNorm)); ...
                zeros(size(coreSubDiffRSNorm)) ones(size(corePlastDiffRSNorm)) ...
                zeros(size(shellSubDiffRSNorm)) ones(size(shellPlastDiffRSNorm))]';
            if isreal(xx) %JMA added
                [pAnova, tAnova] = anovan(xx,grps,'model','interaction','display', 'off');
                fprintf(['\n\tANOVA (fixed model), 2-way: effect of core/shell, p = %0.2f, '...
                    'effect of subsong/plastic, p = %0.2f, interaction, p = %0.2f'], ...
                    pAnova(1), pAnova(2), pAnova(3))
            else
                fprintf('Skipping two-way permutation ANOVA');
            end
        end
        % significance tests: post-hoc, core vs shell
        if ~isempty(coreDiffRSNorm) && ~isempty(shellDiffRSNorm)
            pc = signrank( coreDiffRSNorm);
            ps = signrank(shellDiffRSNorm);
            pMannU = ranksum(coreDiffRSNorm(~isnan(coreDiffRSNorm)), shellDiffRSNorm(~isnan(shellDiffRSNorm)));
            % no subsong/plastic significance tests
            
            fprintf(['\n\tsign-rank test p-value for core: %0.3f' ...
                '\n\tsign-rank test p-value for shell: %0.3f',...
                '\n\tMann-Whitney U test p-value for core v shell: %0.3f\n'],...
                pc,ps,pMannU);
        end
        
        % set bins for histogram
        RSdiffBins = -10:0.5:10;
        plotInterlaceBars(coreDiffRSNorm, shellDiffRSNorm, RSdiffBins);
        ytop = ylim * [0 1]';
        
        % todo: plot significance on graph
        hold on;
        plotSEMBar(      coreDiffRSNorm, ytop  , [0.5 0.5 0.5]);
        plotSEMBar(     shellDiffRSNorm, ytop+1, [  1   0   0]);
        plotSEMBar(   coreSubDiffRSNorm, ytop+2, [0.5 0.5 0.5]);
        plotSEMBar(  shellSubDiffRSNorm, ytop+3, [  1   0   0]);
        plotSEMBar( corePlastDiffRSNorm, ytop+4, [0.5 0.5 0.5]);
        plotSEMBar(shellPlastDiffRSNorm, ytop+5, [  1   0   0]);
        plot([0 0], ylim, 'k--');
        hold off;
        % redo y axis labels
        yt = get(gca,'YTick');
        yt = [yt(yt < ytop) ytop:ytop+5];
        ytl = cellfun(@(x) sprintf('%d',x),num2cell(yt),'UniformOutput',false);
        ytl(end-5:end) = {'Core','Shell','Core/Subsong','Shell/Subsong','Core/Plastic','Shell/Plastic'};
        set(gca,'YTick',yt,'YTickLabel',ytl);
        
        % figure formatting
        xlabel(['normalized ' xlabels{ii}]);
        ylabel('Count');
        legend(sprintf('CORE:  n = %d', numel( coreDiffRSNorm(~isnan( coreDiffRSNorm)))),...
            sprintf('SHELL: n = %d', numel(shellDiffRSNorm(~isnan(shellDiffRSNorm)))));
        xlim([min(RSdiffBins) max(RSdiffBins)]);
        set(gca,'Box','off');
        set(gca, 'FontSize', 14);
        set(get(gca,'XLabel'),'FontSize', 14);
        set(get(gca,'YLabel'),'FontSize', 14);
        set(get(gca,'Title' ),'FontSize', 14);
        title(sprintf('%s, core/shell diff p = %0.2g, core from zero p = %0.2g, shell from zero p = %0.2g', eiTitle{hh}, pMannU, pc, ps));
        set(gcf,'Color',[1 1 1]);
        
        %mean(coreDiffRSNorm)
        %mean(shellDiffRSNorm)
        %pause;
        if params.saveplot
            imFile = sprintf('figures/paper/distanceCorrelations-subjectiveScoreFilter/normedRSdiffs-SUA-%s-%s.pdf', distanceTypes{ii}, filsuff{hh});
            fprintf('Writing image to %s', imFile);
            scrsz = get(0,'ScreenSize');
            set(gcf, 'Position', [1 1 scrsz(3) scrsz(4)]);
    
            export_fig(imFile); 

%            saveCurrFigure(sprintf('figures/A_keeper/mostResponsive/normedRSdiffs-SUA-%s-%s.jpg', distanceTypes{ii}, filsuff{hh}));
        end
    end
end
%{
fprintf('\n\np-values of FR correlation to distance types');
pFields = strcat(distanceTypes, 'Distance_p');
R2Fields = strcat(distanceTypes, 'DistanceR2');
xlabels = strcat({'Linear trend p-values of FR to '}, distanceDescriptions);
   
for hh = 1:3 % inhibited, excited, all
    figure;
    for ii = 1:numel(distanceTypes)
        subplot(nR,nC,ii)
        corrPVals = [allNeuronCorrData.(pFields{ii})];
        corrR2Vals = [allNeuronCorrData.(R2Fields{ii})];
        if hh < 3
             coreCPVs = corrPVals(isPresel & isExcited == hh-1 &  isCore);
            shellCPVs = corrPVals(isPresel & isExcited == hh-1 & ~isCore);
        
            % get % variance explained
            [ coreR2M,  coreR2SEM] = meanSEM(corrR2Vals(isPresel & isExcited == hh-1 &  isCore));
            [shellR2M, shellR2SEM] = meanSEM(corrR2Vals(isPresel & isExcited == hh-1 & ~isCore));
             %coreSubCPVs = corrPVals(isPresel & isExcited == hh-1 &  isCore & ~isPlastic);
            %shellSubCPVs = corrPVals(isPresel & isExcited == hh-1 & ~isCore & ~isPlastic);
            
             %corePlastCPVs = corrPVals(isPresel & isExcited == hh-1 &  isCore & isPlastic);
            %shellPlastCPVs = corrPVals(isPresel & isExcited == hh-1 & ~isCore & isPlastic);
        else
             coreCPVs = corrPVals(isPresel &  isCore);
            shellCPVs = corrPVals(isPresel & ~isCore);

            [ coreR2M,  coreR2SEM] = meanSEM(corrR2Vals(isPresel &  isCore));
            [shellR2M, shellR2SEM] = meanSEM(corrR2Vals(isPresel & ~isCore));

             %coreSubCPVs = corrPVals(isPresel &  isCore & ~isPlastic);
            %shellSubCPVs = corrPVals(isPresel & ~isCore & ~isPlastic);
            
             %corePlastCPVs = corrPVals(isPresel &  isCore & isPlastic);
            %shellPlastCPVs = corrPVals(isPresel & ~isCore & isPlastic);
        end
        
        % test the difference between core and shell
        pMannU = ranksum(coreCPVs(~isnan(coreCPVs)), shellCPVs(~isnan(shellCPVs)));
        fprintf('%s - %s: Mann-Whitney U p-value for core v shell: %0.3f\n',...
            eiTitle{hh}, xlabels{ii},pMannU);
        fprintf(['\tCore variance explained: %0.2f +/- %0.2f, '...
            'shell variance explained: %0.2f +/- %0.2f\n'], ...
        coreR2M, coreR2SEM, shellR2M, shellR2SEM);
        pBins = logspace(-3,0,30);
        plotInterlaceBars(coreCPVs, shellCPVs, pBins);
        
        legend(sprintf('CORE:  n = %d', numel( coreCPVs(~isnan( coreCPVs)))),...
               sprintf('SHELL: n = %d', numel(shellCPVs(~isnan(shellCPVs)))));
        xlabel(xlabels{ii});
        ylabel('Count');
        xlim([0 1])
        
        % figure formatting
        set(gca, 'XTick', [0.05 0.1 0.2 0.4 0.6 0.8]);
        set(gca, 'Box', 'off');
        set(gca, 'FontSize', 12);
%        xticklabel_rotate([],45); % this messes up the figure subplots
        set(get(gca,'XLabel'),'FontSize', 14);
        set(get(gca,'YLabel'),'FontSize', 14);
        set(get(gca,'Title' ),'FontSize', 14);

        % todo: plot the significance markers
        ytop = ylim * [0 1]'; ylim([0 ytop+2]);
        hold on;
        plotSEMBar( coreCPVs, ytop-0.2, [0.5 0.5 0.5]);
        plotSEMBar(shellCPVs, ytop-0.1, [1 0 0]);
        hold off;
    end
    subplot(nR,nC,1);
    title(eiTitle{hh});
    set(gcf,'Color',[1 1 1]);
    if params.saveplot
        saveCurrFigure(sprintf('figures/distanceCorrelations/neuronDistanceCorr-SUA-%s.jpg', filsuff{hh}));
    end

end
%}
end

function [m, v] = meanSEM(set1)
m = nanmean(set1);
if numel(set1) >= 2
    v = nanstd(set1 )/sqrt(numel(set1)-1);
else
    v = 0;
end
end

% plot the error bars w/ SEM on top
function plotSEMBar(set, y, col)
[m,v] = meanSEM(set);
plotHorzErrorBar(m,y,v,col);
end

function x = noNaN(x)
x(isnan(x)) = [];
end
