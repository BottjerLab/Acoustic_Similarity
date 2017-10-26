function recalcClustersLengthWeight(birdID)
% nb: as of 2/20 all these recalc functions are obsoleted
% these were experiments to determine what cluster values were best
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];
timeFlag = ['T-' datestr(clock, 'mm_dd_HH_MM')];
timeClusterDir = [clusterDir timeFlag filesep];


fullFiles = uigetfile([clusterDir 'altClustDataAge-*.mat'],'Select the file to recalculate','MultiSelect','on');
if ~iscell(fullFiles), fullFiles = {fullFiles}; end

mkdir(clusterDir, timeFlag);
% grab a particular session's spectra
fprintf('Loading spectra for bird %s...\n', birdID);
load([dataDir 'allSpecs-' birdID '.mat']);

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
    dLens = abs(ones(nSeld,1) * seldLens - (ones(nSeld,1) * seldLens)');
    dLens = (dLens / max(dLens(:))); % normalization
    
    %dSim = 1 - dLens; dSim = dSim - diag(diag(dSim));
    avg(1:nSeld+1:end) = 0;
    seldMeanLen = squareform(avg);
    
    % normalize cluster local distances by length
    distMats.warpedLocalMean = distMats.warpedLocal ./ seldMeanLen;
    
    % re-evaluate empirical distributions
    nEmpD = min(NC2,2e4);
    
    fprintf('Realculating warped local mean empirical scores...\n');
    
    fld = 'warpedLocalMean';
    [empMats.(fld), empDistrs.(fld)] = makeSelfEmpirical(distMats.(fld), nEmpD);
    
    % independent cosimilarity measures <- note: these are actually
    % dissimilarity/distance
    distMats.localCosim  = pdist(squareform(empMats.warpedLocalMean), 'correlation');
    distMats.globalCosim = pdist(squareform(empMats.global), 'correlation');
    distMats.lengthCosim = squareform(dLens);

    %%
    % optional: examine length infusion
    fld = 'lenAdded';
    ratios = 0:0.05:0.5;
    specFlds = cell(1,numel(ratios));
    for ii = 1:numel(ratios)
        lenInfusRatio = ratios(ii);
        otherRatio = (1.0 - lenInfusRatio)/2;
        specFlds{ii} = [fld '_' num2str(lenInfusRatio)];
        specFlds{ii} = strrep(specFlds{ii}, '.','');
          distMats.(specFlds{ii}) = exp((log(empMats.warpedLocalMean) + log(empMats.global)) * otherRatio + ...
                              log(distMats.lengthCosim) * lenInfusRatio);                          
    end
    
    % do the clustering - the easiest part
    nClusters = 4:25;
    cosimTypes = ['localCosim', 'globalCosim', specFlds];
    %fileLabels = {'localOnly', 'globalOnly', specFlds};
        
    for ii = 1:numel(cosimTypes)
        pairLinks = linkage(distMats.(cosimTypes{ii}),'complete');        
        clusterIdxs.(cosimTypes{ii}) = cluster(pairLinks,'maxclust',nClusters);
    end
    
    % write the clusters to reflect all the different kind of clusterings
    save([timeClusterDir filName],'empMats','empDistrs','distMats', 'clusterIdxs');
end