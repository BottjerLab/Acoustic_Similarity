function clustResultFiles = recalcClustersInnerLenWeight(birdID, fullFiles)
% nb: as of 2/20 all these recalc functions are obsoleted
% these were experiments to determine what cluster values were best
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];
timeFlag = ['T-' datestr(clock, 'mm_dd_HH_MM')];
timeClusterDir = [clusterDir timeFlag filesep];

if nargin < 2
    fullFiles = uigetfile([clusterDir 'altClustDataAge-*.mat'],'Select the file to recalculate','MultiSelect','on');
else
    fullFiles = strrep(fullFiles, clusterDir, '');
end
if ~iscell(fullFiles), fullFiles = {fullFiles}; end

% grab a particular session's spectra
fprintf('Loading spectra for bird %s...\n', birdID);
load([dataDir 'allSpecs-' birdID '.mat']);
mkdir(clusterDir, timeFlag);

clustResultFiles = cell(1,numel(fullFiles)); 
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
    
    % reweight length and the other features
    weightedFeatureRoots = {'length', 'wienerEntropy', 'fundamentalFreq', 'FM'};
    weights = [3 2 2 2];
    isApproved = false(numel(fn),1);
    zWeightTable = zNormedTable;
    for ii = 1:numel(weights)
        thisIndex = strcmp(weightedFeatureRoots{ii}, fnRoots);
        zWeightTable(:,thisIndex) = weights(ii) * zWeightTable(:,thisIndex);
    end
       
    fprintf('Calculating global dissimilarity scores w/o and with reweighting...\n' );
    
    glFld = 'global';
    distMats.(glFld) = pdist(zNormedTable(:,isApproved));
    [empMats.(glFld), empDistrs.(glFld)] = makeSelfEmpirical(distMats.(glFld), nEmpD);
    distMats.(['cosim_' (glFld)]) = pdist(squareform(empMats.(glFld)), 'correlation');

    glFld = 'weightedGlobal'; % note: i don't think these work
    distMats.(glFld) = pdist(zWeightTable(:,isApproved));
    [empMats.(glFld), empDistrs.(glFld)] = makeSelfEmpirical(distMats.(glFld), nEmpD);
    distMats.(['cosim_' (glFld)]) = pdist(squareform(empMats.(glFld)), 'correlation');
    
    % recalculate clustering
    ratios = [1/4 0.5 3/4];
    suffixes = {'LessQuart','Even','Quart'};
    for ii = 1:numel(ratios)
        % recalculated cosimilarity
        frac = ratios(ii); 
        
        fusedVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.global));
        distMats.(['adjustedCosim' suffixes{ii}]) = pdist(squareform(fusedVals), 'correlation');
        
        fusedVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.weightedGlobal));
        distMats.(['weightedCosim' suffixes{ii}]) = pdist(squareform(fusedVals), 'correlation');
    end
    %% do the clustering - the easiest part
    nClusters = 4:25;    
    cosimTypes = [strcat('adjustedCosim', suffixes), strcat('weightedCosim', suffixes)];
        
    for ii = 1:numel(cosimTypes)
        pairLinks = linkage(distMats.(cosimTypes{ii}),'complete');        
        clusterIdxs.(cosimTypes{ii}) = cluster(pairLinks,'maxclust',nClusters);
    end
    
    %% write the clusters to reflect all the different kind of clusterings
    clustResultFiles{hh} = [timeClusterDir filName];
    save(clustResultFiles{hh}, 'empMats','empDistrs','distMats', 'clusterIdxs');
end
% optional: examine proportional weighting of global/local