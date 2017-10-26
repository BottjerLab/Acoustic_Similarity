function recalcClustersBestFeature(birdID)
% nb: as of 2/20 all these recalc functions are obsoleted
% these were experiments to determine what cluster values were best

clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];
timeFlag = ['T-' datestr(clock, 'mm_dd_HH_MM')];
timeClusterDir = [clusterDir timeFlag filesep];


fullFiles = uigetfile([clusterDir 'altClustDataAge-*.mat'],'Select the file to recalculate','MultiSelect','on');
if ~iscell(fullFiles), fullFiles = {fullFiles}; end
% grab a particular session's spectra
fprintf('Loading spectra for bird %s...\n', birdID);
load([dataDir 'allSpecs-' birdID '.mat']);
mkdir(clusterDir, timeFlag);

%% normalize cluster local distances by length
for hh = 1:numel(fullFiles)
    filName = fullFiles{hh};
    
    clustSession = strrep(filName(1:end-5), 'altClustDataAge-',''); % take off the prefix and the '.mat' suffix;
    fprintf('Re-analyzing session %s...\n', clustSession);
    
    distMats = []; empMats = []; empDistrs = []; load([clusterDir filName]);
    
    seldSylls = DRsylls(seld);
    %%
    nSeld = numel(seldSylls);
    seldLens = [seldSylls.stop] - [seldSylls.start];
    
    % find the mean length for each pair of syllables (not the most eff. but
    % whatev)
    
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
    
    % get the list of features
    fn = fieldnames(featureTable);
    fnRoots = strtok(fn,'_');

    seldFeatureTable = cellfun(@(x) [featureTable(seld).(x)]', fn', 'UniformOutput',false);
    seldFeatureTable = [seldFeatureTable{:}];
    zNormedTable = zscore(seldFeatureTable);
    
    % pick only a few good features
    approvedFeatureRoots = {'length', 'wienerEntropy', 'fundamentalFreq', 'FM'};
    isApproved = false(numel(fn),1);
    for ii = 1:numel(approvedFeatureRoots)
        iRoot = approvedFeatureRoots{ii};
        isApproved = isApproved | strcmp(iRoot, fnRoots);
    end
       
    fprintf('Calculating global dissimilarity scores for cultivated list...\n' );
    
    glFld = 'approvedGlobal';
    distMats.(glFld) = pdist(zNormedTable(:,isApproved));
    [empMats.(glFld), empDistrs.(glFld)] = makeSelfEmpirical(distMats.(glFld), nEmpD);
    distMats.(['cosim_' (glFld)]) = pdist(squareform(empMats.(glFld)), 'correlation');
    
    % recalculate clustering
    frac = 0.5; % even
    fusedHalfVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.approvedGlobal));
    frac = 2/3;
    fusedThirdVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.approvedGlobal));
    frac = 1/3;
    fusedLowThirdVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.approvedGlobal));
    
    % recalculated cosimilarity
    distMats.adjustedCosimEven = pdist(squareform(fusedHalfVals), 'correlation');
    distMats.adjustedCosimThird = pdist(squareform(fusedThirdVals), 'correlation');
    distMats.adjustedCosimLessThird = pdist(squareform(fusedLowThirdVals), 'correlation');
    %%
    
    % do the clustering - the easiest part
    nClusters = 4:25;    
    cosimTypes = {'adjustedCosimEven', 'adjustedCosimThird', 'adjustedCosimLessThird'};
        
    for ii = 1:numel(cosimTypes)
        pairLinks = linkage(distMats.(cosimTypes{ii}),'complete');        
        clusterIdxs.(cosimTypes{ii}) = cluster(pairLinks,'maxclust',nClusters);
    end
    
    %% write the clusters to reflect all the different kind of clusterings
    save([timeClusterDir filName],'empMats','empDistrs','distMats', 'clusterIdxs');
end
% optional: examine proportional weighting of global/local