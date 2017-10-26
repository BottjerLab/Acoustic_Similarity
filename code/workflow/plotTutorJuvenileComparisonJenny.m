function plotTutorJuvenileComparisonJenny(sessionID, neuronNum, clusterNum, isCore, params, varargin)
% plot juvenile syllables in a cluster that are similar to/far from a tutor
% and show mean/sem of their differences
if nargin < 5 || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});

birdID = strtok(sessionID,'_');
% lots of preloading
%% load syllables and labels for a given session
rep = reportOnData(birdID, sessionID);
fprintf('Loading syllables for session %s...\n', sessionID)
juvSylls = loadFromManifest(rep.manifest, 'approvedSyllables');
juvLabels = loadFromManifest(rep.manifest, 'acceptedLabels');

fprintf('Total syllables IDed / total: (%d/%d)\n', sum(~isnan(juvLabels)), numel(juvLabels));
%% load the particular neuron within that session
spikes = loadSpikeData(rep.spikeFiles);
neuronSpikes = spikes{neuronNum};

%% load tutor syllables
tutorFil = ['data' filesep birdID filesep 'tutor-' birdID '.mat'];
if exist(tutorFil,'file') ~= 2
    error('compareTutorJuvenile:tutorFileNotFound',...
        'Tutor file %s for bird %s not found...',tutorFil, birdID);
end
load(tutorFil);

%% load distances to tutor syllables
centMatch = loadFromManifest(rep.manifest,'centMatch');
distToCentral  = loadFromManifest(rep.manifest,'distToCentral');
distClusterInter = loadFromManifest(rep.manifest,'interClusterDists');
% get the first consensus tutor syllable to the cluster in question
isCluster = (juvLabels == clusterNum);
assert(any(isCluster));
tutorMatchToClust = centMatch(isCluster);
% does this cluster have all the same consensus syllable matched?
assert(all(tutorMatchToClust == tutorMatchToClust(1)));
tutorNum = tutorMatchToClust(1);

tutorExample = tutorSylls(find(...
    strcmp({tutorSylls.type},char('a'-1+tutorNum)),...
    1));
%% get the longest syllable length out of the cluster and syllables to make a common scale
tutorLen = tutorExample.stop - tutorExample.start;
juvLens = [juvSylls(isCluster).stop] - [juvSylls(isCluster).start];
if size(juvLens, 2) == 1, juvLens = juvLens'; end;
maxLen = max([tutorLen juvLens]);

%% get the distances for each cluster
juvCluster = juvSylls(isCluster);
[~,DRfile] = findInManifest(rep.manifest,'sAudio');
[juvCluster.file] = deal(DRfile);
num2cell(distToCentral(isCluster));  [juvCluster.distance] = ans{:};
num2cell(distClusterInter(isCluster)); [juvCluster.interClusterDist] = ans{:};

%% load baseline periods and get the closest ones
baselinePeriodsS = loadFromManifest(rep.manifest, 'baselinePeriodsS');

for ii = 1:numel(juvCluster)
    distToFront = abs([baselinePeriodsS.stop]- juvCluster(ii).start);
    distToBack = abs([baselinePeriodsS.start]- juvCluster(ii).stop);
    dist = min([distToFront; distToBack], [], 1);
    [~, idxMatchedBaseline] = min(dist);
    matchedBaseline(ii) = baselinePeriodsS(idxMatchedBaseline);
end
%% plot raw spikes
%{
[~,relSpikeTimes] = countSpikes(juvCluster, neuronSpikes);
subplot(211);
plotRaster(juvCluster, relSpikeTimes);
ylabel('juvenile syllables')
subplot(212);
[~, relSpikeTimes] = countSpikes(matchedBaseline, neuronSpikes);
plotRaster(matchedBaseline, relSpikeTimes);
ylabel('matched baseline syllables')
%}
%% get the firing data
nSData = loadFromManifest(rep.manifest,'neuronSyllableData');
nSEntry = nSData(neuronNum, clusterNum); 
if numel(nSEntry.syllIndex) ~= numel(juvCluster)
    fprintf('Number of syllables in neuron data: %d, in cluster, %d\n', numel(nSEntry.syllIndex), numel(juvCluster));
    return;
end
RSs = -diff(nSEntry.rawRates,1,1);

RS_syll = nSEntry.rawRates(1,:); RS_base = nSEntry.rawRates(2,:);
covRS = cov(RS_syll, RS_base); if numel(covRS) > 1, covRS = covRS(2,1); end;
dRS_SE = sqrt(var(RS_syll) + var(RS_base) - 2 * covRS) / sqrt(numel(nSEntry.syllIndex));
zRSs = (RSs - mean(RSs)) / dRS_SE;
num2cell(RSs); [juvCluster.bcFR] = ans{:};
num2cell(zRSs); [juvCluster.zRS] = ans{:};

%% order the syllables by similarity
% these are sorted clusters
quartiles = prctile(distToCentral(isCluster), [25 75]);
nearCluster = sortBy(juvCluster([juvCluster.distance] < quartiles(1)), 'distance');
farCluster  = sortBy(juvCluster([juvCluster.distance] > quartiles(2)), 'distance', 'descend');

%% plot in the same window style as checkStereotypy, one row at a time, w/ common scale
nR = 15;
assert(mod(nR,2)==1);

fMaxX = 0.4;
fMaxY = 1/nR;
plotParams = processArgs(defaultParams, 'preroll', 5, 'postroll', 5, ...
    'dgram.minContrast', 3e-11);
% get the first tutor syllable of that type
figBox = [0 (nR-1)/nR fMaxX*(tutorLen/maxLen) fMaxY];
subplot('Position', figBox);

if isCore, textLabel = 'CORE'; else, textLabel = 'SHELL'; end
% plotting the spectrogram of the tutor clip in a figure box
plotDGFromRegion(tutorStruct, tutorExample, plotParams);
set(gca,'XTick',[], 'YTick', [], 'Box', 'off');
text(1,0,['central tutor syllable, ' textLabel], 'Units','normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
for jj = 1:min(nR/2,numel(nearCluster))
    % plotting the spectrogram of a near syllable in a figure box
    nearLen = nearCluster(jj).stop - nearCluster(jj).start;
    figBox = [0 (nR-1-jj)/nR fMaxX*(nearLen/maxLen) fMaxY];
    subplot('Position', figBox);
    plotDGFromRegion([], nearCluster(jj), plotParams);
    
    % text label of distance and baseline-corrected firing rate
    text(1,0,sprintf('Dist: %0.3f, zRS: %0.3f', nearCluster(jj).distance, nearCluster(jj).zRS), ...
        'Units','normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    set(gca,'XTick',[], 'YTick', [], 'Box', 'off', 'Color', [1 0 0]);
    drawnow;
end
for jj = 1:min(nR/2,numel(farCluster))
    % plotting the spectrogram of a near syllable in a figure box
    farLen = farCluster(jj).stop - farCluster(jj).start;
    figBox = [0 ((nR-1)/2-jj)/nR fMaxX*(farLen/maxLen) fMaxY];
    subplot('Position', figBox);
    plotDGFromRegion([], farCluster(jj), plotParams);
    
    % text label of distance and baseline-corrected firing rate
    text(1,0,sprintf('Dist: %0.3f, zRS: %0.3f', farCluster(jj).distance, farCluster(jj).zRS), ...
        'Units','normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    set(gca,'XTick',[], 'YTick', [], 'Box', 'off', 'Color', [1 0 0]);
    drawnow;
end

% plotting firing data
subplot(1,2,2); 

nearMean = mean([nearCluster.zRS]); nearSEM = std([nearCluster.zRS])/sqrt(numel(nearCluster)-1);
 farMean = mean([ farCluster.zRS]);  farSEM = std([ farCluster.zRS])/sqrt(numel( farCluster)-1);
 
dRSNearFar = (nearMean - farMean) / dRS_SE;

[~,pT] = ttest2([nearCluster.zRS], [farCluster.zRS]);
errorbar([1 2], [nearMean farMean], [nearSEM farSEM], 'k.');
ylabel('z-standardized response strength'); set(gca,'XTick',[1 2], 'XTickLabel', {'Near', 'Far'});
xlabel('Syllable distance to central tutor');
set(gca,'Box','off');
title(sprintf('%s, RS Difference = %0.2f (p = %0.2f)', textLabel, dRSNearFar, pT) ,'FontSize',14);
%iCDq = prctile([juvCluster.interClusterDist],[25 75]);

%% draw line to distinguish tutor syllables from juvenile syllables
set(gcf,'Name',sprintf('%s, tutor %d, cluster %d', sessionID, tutorNum, clusterNum),'Units','normalized','Position',[0 0 1 1]);
annotation('line', [0 fMaxX], (1-1/nR)*[1 1],'LineWidth', 2, 'Color', [0.2 0.8 0.3])
% draw line to distinguish near syllables from far syllables
annotation('line', [0 fMaxX], (1/2-1/(2*nR))*[1 1], 'LineWidth', 2, 'Color', [0.8 0.2 0.2])

if params.saveplot
    imFile = ['figures\paper\individualNeuronClusterCorrelations\tutorJuvExample' '-' sessionID '-neuron' ...
        num2str(neuronNum) '-cluster' num2str(clusterNum) '.pdf'];
    fprintf('Saving figure to file %s...\n',imFile);
    export_fig(imFile);
end
end