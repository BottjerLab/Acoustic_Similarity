function plotTutorJuvenileComparison(sessionID, neuronNum, clusterNum, params, varargin)
% plot juvenile syllables in a cluster that are similar to/far from a tutor
% and somehow show their firing rates
if nargin < 4 || isempty(params)
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
%% load tutor syllables
tutorFil = ['data' filesep birdID filesep 'tutor-' birdID '.mat'];
if exist(tutorFil,'file') ~= 2
    error('compareTutorJuvenile:tutorFileNotFound',...
        'Tutor file %s for bird %s not found...',tutorFil, birdID);
end
load(tutorFil);

%% load distances to tutor syllables
consMatch = loadFromManifest(rep.manifest,'centMatch');
distToConsensus  = loadFromManifest(rep.manifest,'distToCentral');
distClusterInter = loadFromManifest(rep.manifest,'interClusterDists');
% get the first consensus tutor syllable to the cluster in question
isCluster = (juvLabels == clusterNum);
assert(any(isCluster));
tutorMatchToClust = consMatch(isCluster);
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

%% order the syllables by similarity
juvCluster = juvSylls(isCluster);
[~,DRfile] = findInManifest(rep.manifest,'sAudio');
[juvCluster.file] = deal(DRfile);
quartiles = prctile(distToConsensus(isCluster), [25 75]);
num2cell(distToConsensus(isCluster));  [juvCluster.distance] = ans{:};
num2cell(distClusterInter(isCluster)); [juvCluster.interClusterDist] = ans{:};

nearCluster = sortBy(juvCluster([juvCluster.distance] < quartiles(1)), 'distance');
farCluster  = sortBy(juvCluster([juvCluster.distance] > quartiles(2)), 'distance', 'descend');

%%
nSData = loadFromManifest(rep.manifest,'neuronSyllableData');
nSEntry = nSData(neuronNum, clusterNum); 
%FIXME: why are these not the same? for instance,  Gy217 clust #1/5/6
if numel(nSEntry.syllIndex) ~= numel(juvCluster)
    fprintf('Number of syllables in neuron data: %d, in cluster, %d\n', numel(nSEntry.syllIndex), numel(juvCluster));
    return;
end
%% plot in the same window style as checkStereotypy, one row at a time, w/ common scale
nR = 15;
assert(mod(nR,2)==1);

fMaxX = 1/3;
fMaxY = 1/nR;
plotParams = processArgs(defaultParams, 'preroll', 5, 'postroll', 5, ...
    'dgram.minContrast', 3e-11);
% get the first tutor syllable of that type
figBox = [0 (nR-1)/nR fMaxX*(tutorLen/maxLen) fMaxY];
subplot('Position', figBox);
% plotting the spectrogram of a clip in a figure box
plotDGFromRegion(tutorStruct, tutorExample, plotParams);
set(gca,'XTick',[], 'YTick', [], 'Box', 'off');
text(1,0,'consensus tutor syllable', 'Units','normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
for jj = 1:min(nR/2,numel(nearCluster))
    nearLen = nearCluster(jj).stop - nearCluster(jj).start;
    figBox = [0 (nR-1-jj)/nR fMaxX*(nearLen/maxLen) fMaxY];
    subplot('Position', figBox);
    plotDGFromRegion([], nearCluster(jj), plotParams);
    
    %%
    text(1,0,num2str(nearCluster(jj).distance, '%0.3f'), 'Units','normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    set(gca,'XTick',[], 'YTick', [], 'Box', 'off', 'Color', [1 0 0]);
    drawnow;
end
for jj = 1:min(nR/2,numel(farCluster))
    farLen = farCluster(jj).stop - farCluster(jj).start;
    figBox = [0 ((nR-1)/2-jj)/nR fMaxX*(farLen/maxLen) fMaxY];
    subplot('Position', figBox);
    plotDGFromRegion([], farCluster(jj), plotParams);
    
    %%
    text(1,0,num2str(farCluster(jj).distance, '%0.3f'), 'Units','normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    set(gca,'XTick',[], 'YTick', [], 'Box', 'off', 'Color', [1 0 0]);
    drawnow;
end
% firing data

RSs = -diff(nSEntry.rawRates,1,1);
subplot(2,2,2); plot(RSs, [juvCluster.distance], 'k.');
hold on; plot(xlim, quartiles(1) * [1 1], 'g--', xlim, quartiles(2) * [1 1], 'r--'); hold off;
xlabel('syllable - baseline firing (Hz)'); ylabel('distance to consensus closest tutor');
set(gca,'Box','off');
title('Tutor similarity','FontSize',14);
iCDq = prctile([juvCluster.interClusterDist],[25 75]);
subplot(2,2,4); plot(RSs, [juvCluster.interClusterDist], 'k.');
hold on; plot(xlim, iCDq(1) * [1 1], 'g--', xlim, iCDq(2) * [1 1], 'r--'); hold off;
xlabel('syllable - baseline firing (Hz)'); ylabel('distance within cluster');
title('Cluster variability','FontSize',14);

%% draw line to distinguish tutor syllables from juvenile syllables
set(gcf,'Name',sprintf('%s, tutor %d, cluster %d', sessionID, tutorNum, clusterNum),'Units','normalized','Position',[0 0 1 1]);
annotation('line', [0 fMaxX], (1-1/nR)*[1 1],'LineWidth', 2, 'Color', [0.2 0.8 0.3])
% draw line to distinguish near syllables from far syllables
annotation('line', [0 fMaxX], (1/2-1/(2*nR))*[1 1], 'LineWidth', 2, 'Color', [0.8 0.2 0.2])

if params.saveplot
    imFile = ['figures' filesep 'individualNeuronClusterCorrelations\tutorJuvExample' '-' sessionID '-neuron' ...
        num2str(neuronNum) '-cluster' num2str(clusterNum) '.jpg'];
    fprintf('Saving figure to file %s...\n',imFile);
    export_fig(imFile);
end
end