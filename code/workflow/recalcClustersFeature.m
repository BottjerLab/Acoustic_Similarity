function recalcClustersFeature(birdID)

clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];
timeFlag = ['T-' datestr(clock, 'mm_dd_HH_MM')];
timeClusterDir = [clusterDir timeFlag filesep];

mkdir(clusterDir, timeFlag);

fullFiles = uigetfile([clusterDir 'altClustDataAge-*.mat'],'Select the file to recalculate','MultiSelect','on');
if ~iscell(fullFiles), fullFiles = {fullFiles}; end
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
    lenField = find(strcmp(fn,'length')); % remove the length field
    
    [pres, posts] = strtok(fn,'_');
    posts = strrep(posts,'_','');
    
    seldFeatureTable = cellfun(@(x) [featureTable(seld).(x)]', fn', 'UniformOutput',false);
    seldFeatureTable = [seldFeatureTable{:}];
    zNormedTable = zscore(seldFeatureTable);
    [fnRoots, ~, fnRootTable] = unique(pres );
    [fnReds,  ~, fnRedTable ] = unique(posts); 

    % make exception for length field
%{
    lenIdx = fnRootTable(lenField);
    fnRootTable(lenField) = NaN;
    fnRootTable(fnRootTable > lenIdx) = fnRootTable(fnRootTable > lenIdx) - 1;
    fnRoots(lenIdx) = [];
%}
    % length field has no reduction, so it can't be included
    lenIdx = fnRedTable(lenField);
    fnRedTable(lenField) = NaN;
    fnRedTable(fnRedTable > lenIdx) = fnRedTable(fnRedTable > lenIdx) - 1;
    fnReds(lenIdx) = [];
    clear lenIdx
    
    for ii = 1:numel(fnRoots)
        fprintf('Calculating global dissimilarity scores among root feature %s...\n', fnRoots{ii});
        
        glFld = ['global_' fnRoots{ii}];
        seldFeatures = (ii == fnRootTable);% seldFeatures(lenField) = true;
        
        distMats.(glFld) = pdist(zNormedTable(:,seldFeatures));
        [empMats.(glFld), empDistrs.(glFld)] = makeSelfEmpirical(distMats.(glFld), nEmpD);
        distMats.(['cosim_' fnRoots{ii}]) = pdist(squareform(empMats.(glFld)), 'correlation');        
        %keyboard
        %localAugScore = sqrt(empMats.(glFld) .* empMats.warpedLocalMean);
        %rawScore = empMats.(glFld);        
    end
    
    for ii = 1:numel(fnReds)
        fprintf('Calculating global dissimilarity scores among reduction %s...\n', fnReds{ii});
        
        glFld = ['global_' fnReds{ii}];
        seldFeatures = (ii == fnRedTable); % seldFeatures(lenField) = true;
        distMats.(glFld) = pdist(zNormedTable(:,seldFeatures));
        [empMats.(glFld), empDistrs.(glFld)] = makeSelfEmpirical(distMats.(glFld), nEmpD);
        
        %localAugScore = sqrt(empMats.(glFld) .* empMats.warpedLocalMean);
        distMats.(['cosim_' fnReds{ii}]) = pdist(squareform(empMats.(glFld)), 'correlation');        
    end
    %%
    % independent cosimilarity measures
    distMats.localCosim  = pdist(squareform(empMats.warpedLocalMean), 'correlation');
    distMats.globalCosim = pdist(squareform(empMats.global), 'correlation');
    
    % recalculate clustering
    frac = 0.5; % even
    fusedHalfVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.global));
    frac = 2/3;
    fusedThirdVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.global));
    frac = 1/3;
    fusedLowThirdVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * log(empMats.global));
    % recalculated cosimilarity
    distMats.adjustedCosimEven = pdist(squareform(fusedHalfVals), 'correlation');
    distMats.adjustedCosimThird = pdist(squareform(fusedThirdVals), 'correlation');
    distMats.adjustedCosimLessThird = pdist(squareform(fusedLowThirdVals), 'correlation');
    %%
    % optional: examine empirical value of joint multiple
%    fld = 'fused';
%    distMats.(fld) = sqrt(distMats.warpedLocalMean .* distMats.global);
%    [empMats.(fld), empDistrs.(fld)] = makeSelfEmpirical(distMats.(fld), nEmpD);
    %distMats.fusedEmpCosim = pdist(squareform(distMats.fused), 'correlation');
    
    % do the clustering - the easiest part
    nClusters = 4:25;    
    cosimTypes = cat(2,{'cosim', 'adjustedCosimEven', 'adjustedCosimThird', 'adjustedCosimLessThird', ...
        'localCosim', 'globalCosim'}, ...
        strcat('cosim_',fnRoots)', strcat('cosim_',fnReds)');
    fileLabels = cat(2,'original', 'meanAdj', 'meanAdjMoreLocal', 'meanAdjLessLocal',...
        'localOnly', 'globalOnly',...
        strcat(fnRoots, 'Only')', strcat(fnReds, 'Only')');
        
    for ii = 1:numel(cosimTypes)
        pairLinks = linkage(distMats.(cosimTypes{ii}),'complete');        
        clusterIdxs.(cosimTypes{ii}) = cluster(pairLinks,'maxclust',nClusters);
    end
    
    %% write the clusters to reflect all the different kind of clusterings
    save([timeClusterDir filName],'empMats','empDistrs','distMats', 'clusterIdxs');
end
% optional: examine proportional weighting of global/local