function recalcClustersJenny(birdID)
% nb: as of 2/20 all these recalc functions are obsoleted
% these were experiments to determine what cluster values were best%%
birdID = 'R204';
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];

prefix = 'altClustDataAge-';
filNames = dir([clusterDir prefix '*.mat']); filNames = {filNames.name};
filNames = {filNames{1:end-1}};

clustSessions = strrep(strrep(filNames, '.mat',''), prefix,''); % take off the prefix and the '.mat' suffix;

fprintf('Loading syllable lengths for bird %s...\n', birdID);    
load([dataDir 'allSpecs-' birdID '.mat'], 'DRsylls');

%% normalize cluster local distances by length
DRages = [DRsylls.age];
DRlens = [DRsylls.stop] - [DRsylls.start];
%%
for ii = 1:numel(filNames)
    thisAge = str2double(strtok(clustSessions{ii},'-'));
    
    fprintf('Loading clustering data for session %s...\n', clustSessions{ii});
    load([clusterDir filNames{ii}]);

    seld = find(DRages == thisAge);
    nSeld = numel(seld);
    seldLens = DRlens(seld);
    NC2 = nchoosek(nSeld,2);
    
    avg = ones(nSeld,1) * seldLens + (ones(nSeld,1) * seldLens)' / 2; 
    avg(1:nSeld+1:end) = 0;
    seldMeanLen = squareform(avg);
    % normalize cluster local distances by length
    distMats.warpedLocalMean = distMats.warpedLocal ./ seldMeanLen;
    
    % re-evaluate empirical distributions
    nEmpD = min(NC2,2e4);
    
    fprintf('Realculating warped local mean empirical scores...\n');
    
    fld = 'warpedLocalMean';
    [empMats.(fld), empDistrs.(fld)] = makeSelfEmpirical(distMats.(fld), nEmpD);
    
    % base clustering on
    % recalculated cosimilarity
    %%
    distMats.localCosim  = pdist(squareform(empMats.warpedLocalMean), 'correlation');
    distMats.globalCosim = pdist(squareform(empMats.global), 'correlation');
    
    % do the clustering - the easiest part
    nClusters = 4:10;
    pairLinks = linkage(distMats.localCosim,'complete');
    localCosimIdxs = cluster(pairLinks,'maxclust',nClusters);
    
    pairLinks = linkage(distMats.globalCosim,'complete');
    globalCosimIdxs = cluster(pairLinks,'maxclust',nClusters);
    
    % now do some plotting
    params = defaultParams;
    clusterFile = [clusterDir 'expClusters-' clustSessions{ii} '-' ];
    
    seldSylls = DRsylls(seld);
    for jj = 1:nClusters(end)
        fig = figure(jj);
        fprintf('Plotting mean-adjusted clusters for syllable %d...\n',jj);
        
        mosaicDRSpec(seldSylls(localCosimIdxs(:,end)==jj), params, ...
            'dgram.minContrast', 3e-11, ...
            'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5);
        set(fig, 'Name', sprintf('Syllable #%d, Testing',jj));
        figFilName = [clusterFile 'localOnly-c' num2str(jj) '.jpg'];
        saveCurrFigure(figFilName);
        close(fig);
        
        fig = figure(jj);
        fprintf('Plotting fused empirical clusters for syllable %d...\n',jj);
        
        mosaicDRSpec(seldSylls(globalCosimIdxs(:,end)==jj), params, ...
            'dgram.minContrast', 3e-11, ...
            'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5);
        set(fig, 'Name', sprintf('Syllable #%d, Testing',jj));
        figFilName = [clusterFile 'globalOnly-c' num2str(jj) '.jpg'];
        saveCurrFigure(figFilName);
        
        close(fig);
    end
end
% optional: examine proportional weighting of global/local