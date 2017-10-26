function recalcClusters(birdID)
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];

fullFiles = dir([clusterDir 'altClustDataAge-*.mat']);
fullFiles = {fullFiles.name};
%[filName, pathName] = uigetfile([clusterDir 'altClustDataAge-*.mat'],'Select the file to recalculate');
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
    
    % independent cosimilarity measures
    distMats.localCosim  = pdist(squareform(empMats.warpedLocalMean), 'correlation');
    distMats.globalCosim = pdist(squareform(empMats.global), 'correlation');
    
    % recalculate clustering
    frac = 0.5; % even
    fusedVals = exp(frac * log(empMats.warpedLocalMean) + (1 - frac) * (empMats.global));
    
    % recalculated cosimilarity
    distMats.adjustedCosim = pdist(squareform(fusedVals), 'correlation');
    
    % optional: examine empirical value of joint multiple
    fld = 'fused';
    distMats.(fld) = distMats.warpedLocalMean .* distMats.global;
    [empMats.(fld), empDistrs.(fld)] = makeSelfEmpirical(distMats.(fld), nEmpD);
    
    distMats.fusedEmpCosim = pdist(squareform(distMats.fused), 'correlation');
    
    % do the clustering - the easiest part
    nClusters = 4:25;
    
    cosimTypes = {'cosim', 'adjustedCosim', 'localCosim', 'globalCosim', 'fused'};
    fileLabels = {'original', 'meanAdj', 'localOnly', 'globalOnly', 'fusedEmp'};
    logPrefix = [clusterDir filesep 'clustIdxLog-' clustSession '-'];
    for ii = 1:numel(cosimTypes)
        pairLinks = linkage(distMats.(cosimTypes{ii}),'complete');
        specificIdxs = cluster(pairLinks,'maxclust',nClusters);
        clustIdxs.(cosimTypes{ii}) = specificIdxs;
        
        diary([logPrefix fileLabels{ii} '.txt']);
        fprintf('Cross-tabulation of syllables for %s cosimilarity clustering...\n', fileLabels{ii});
        nClusterings = numel(nClusters);
        nMaxClust = max(nClusters);
        splitNodes = NaN(1,nClusterings);
        parentRef = NaN(nMaxClust, nClusterings);
        branches = zeros(2,nClusterings);
        for jj = 1:nClusterings-1
            ct = crosstab(specificIdxs(:,jj), specificIdxs(:,jj+1));
            splitNodes(jj) = find(sum(ct > 0, 2)==2,1);
            branches(:, jj) = find(ct(splitNodes(jj),:) > 0);
            for kk = 1:size(ct,2)
                parentRef(kk,jj+1) = find(ct(:, kk) > 0, 1);
            end
        end
        
        splitNodes
        parentRef
        diary off
    end
    
    %%%%%%%% now do some plotting
    params = defaultParams;
    clusterFile = [clusterDir filesep 'expClusters-' clustSession '-'];
    for ii = 1:nClusters(end)
        for jj = 1:numel(cosimTypes)
            fig = figure(ii);
            fprintf('Plotting %s clusters for syllable %d...\n', cosimTypes{jj}, ii);
            
            mosaicDRSpec(seldSylls(clustIdxs.(cosimTypes{jj})(:,end)==ii), params, ...
                'dgram.minContrast', 3e-11, ...
                'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5);
            set(fig, 'Name', sprintf('Syllable #%d, %s', ii, cosimTypes{jj}));
            figFilName = [clusterFile fileLabels{jj} '-c' num2str(ii) '.jpg'];
            saveCurrFigure(figFilName);
            close(fig);
        end
    end
end
% optional: examine proportional weighting of global/local