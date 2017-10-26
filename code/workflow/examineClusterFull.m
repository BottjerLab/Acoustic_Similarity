%function examineClusterFull(birdID, sessFilter)
function clustQuality = examineClusterFull(birdID, ages)
% get the objective cluster quality for the given ages of a bird

% TODO: test rewrite - make objective cluster score
% get the files
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];

clustQuality = cell(1,numel(ages));
for ii = 1:numel(ages) 
    thisAge = ages(ii);
    [fil, filExist] = getLatestFile([clusterDir 'altClustDataAge-' ...
        num2str(thisAge) '*.mat']);
    
    if ~filExist
        error('examineClusterFull:missingDistances',...
            'Missing distances file for bird %s, age %d', birdID, thisAge);
    end
    load(fil, 'distMats');
    
    labelsFil = sprintf('%sacceptedLabels-%s-age%d.mat', dataDir, birdID, thisAge);
    if ~exist(labelsFil,2)
        error('examineClusterFull:missingClusters',...
            'Missing clusterfile for bird %s, age %d', birdID, thisAge);        
    end
    syllLabels = loadAcceptedLabels(birdID, thisAge);
    
    % assign the types
    num2cell(syllLabels); [thisAgeSylls.type] = ans{:}; %#ok<NOANS>
        
    [~,clustQuality{ii}] = clusterQuality(distMats.cosim, syllLabels);
end
    %{
    % cross tabulate / tree out types
    fprintf('Cross-tabulation of syllables...\n');
    nClusterings = size(syllLabels,2);
    nMaxClust = max(syllLabels(:,end));
    splitNodes = zeros(1,nClusterings);
    parentRef = NaN(nMaxClust, nClusterings);
    for jj = 1:nClusterings-1
        ct = crosstab(syllLabels(:,jj), syllLabels(:,jj+1));
        splitNodes(jj) = find(sum(ct > 0, 2)==2,1);
        for kk = 1:size(ct,2)
            parentRef(kk,jj+1) = find(ct(:, kk) > 0, 1);
        end
    end

    nSylls = sum(isThisAge);
    mConform = zeros(nSylls, nClusterings);
    vConform = zeros(nSylls, nClusterings);
    clustQuality = NaN(nMaxClust, nClusterings);
    clustMConform = NaN(nMaxClust, nClusterings);
    clustVConform = NaN(nMaxClust, nClusterings);
    
    for jj = 1:nClusterings
        thisIdxs = syllLabels(:,jj);
    
        [mConform(:,jj), vConform] = conformity(distMats.cosim, thisIdxs);
        nClusters = max(thisIdxs);
        [~,clustQuality(1:nClusters, jj)] = clusterQuality(distMats.cosim, thisIdxs);
        for kk = 1:nClusters
            clustMConform(kk,jj) = mean(mConform(thisIdxs == kk));
            clustVConform(kk,jj) = mean(vConform(thisIdxs == kk));
        end
    end
    
    parentRef
    clustQuality
    clustMConform
    clustVConform
    %{
    figureDir = [pwd filesep 'figures' filesep 'cluster-' birdID filesep clustSession{ii} filesep];
    mkdir([pwd filesep 'figures' filesep 'cluster-' birdID filesep], clustSession{ii});    
    nClusts = max([thisAgeSylls.type]);
    for jj = 1:nClusts
        trainClust = thisAgeSylls([thisAgeSylls.type]==jj);
        if ~isempty(trainClust), 
            figure
            hf = mosaicDRSpec(trainClust, [], 'dgram.minContrast', 1e-10, 'maxMosaicLength', 5.5, 'noroll');
            set(hf,'Name', sprintf('Trained cluster %d, age %d, bird %s', jj, currAge, birdID));
            saveCurrFigure(sprintf('%s%s_a%d_c%d-train.jpg', figureDir, birdID, currAge, jj));
            close(hf);
        end
    end
    %}
end
%{
% get all cluster files that match the pattern for 
files = dir([clusterDir sessFilter]); files = {files.name}';
fPats = {'altClustDataAge-'};

isAssigned = strncmp(fPats{1}, files, length(fPats{1}));
clustSession = strrep(strrep(files(isAssigned), fPats{1}, ''),'.mat','');

fprintf('Loading library for bird %s...\n',birdID);
load([dataDir 'allSpecs-' birdID])
for ii = 1:numel(clustSession)
    if ~exist([clusterDir fPats{1} clustSession{ii} '.mat'], 'file') %&& ...
         %exist([clusterDir fPats{2} clustSession{ii}],'file'))
        continue;
    end
    fprintf('Loading session %s...\n', clustSession{ii});
    clusterIdxs = []; load([clusterDir fPats{1} clustSession{ii}]);
    %load([clusterDir fPats{2} clustSession{ii}]);
    
    currAge = str2double(clustSession{ii}(1:2));
    
    isThisAge = [DRsylls.age]==currAge;
    thisAgeSylls = DRsylls(isThisAge);

    % assign the types
    num2cell(clusterIdxs(:,end)); [thisAgeSylls.type] = ans{:}; %#ok<NOANS>
        
    % cross tabulate / tree out types
    fprintf('Cross-tabulation of syllables...\n');
    nClusterings = size(clusterIdxs,2);
    nMaxClust = max(clusterIdxs(:,end));
    splitNodes = zeros(1,nClusterings);
    parentRef = NaN(nMaxClust, nClusterings);
    for jj = 1:nClusterings-1
        ct = crosstab(clusterIdxs(:,jj), clusterIdxs(:,jj+1));
        splitNodes(jj) = find(sum(ct > 0, 2)==2,1);
        for kk = 1:size(ct,2)
            parentRef(kk,jj+1) = find(ct(:, kk) > 0, 1);
        end
    end

    nSylls = sum(isThisAge);
    mConform = zeros(nSylls, nClusterings);
    vConform = zeros(nSylls, nClusterings);
    clustQuality = NaN(nMaxClust, nClusterings);
    clustMConform = NaN(nMaxClust, nClusterings);
    clustVConform = NaN(nMaxClust, nClusterings);
    
    for jj = 1:nClusterings
        thisIdxs = clusterIdxs(:,jj);
    
        [mConform(:,jj), vConform] = conformity(distMats.cosim, thisIdxs);
        nClusters = max(thisIdxs);
        [~,clustQuality(1:nClusters, jj)] = clusterQuality(distMats.cosim, thisIdxs);
        for kk = 1:nClusters
            clustMConform(kk,jj) = mean(mConform(thisIdxs == kk));
            clustVConform(kk,jj) = mean(vConform(thisIdxs == kk));
        end
    end
    
    parentRef
    clustQuality
    clustMConform
    clustVConform
    figureDir = [pwd filesep 'figures' filesep 'cluster-' birdID filesep clustSession{ii} filesep];
    mkdir([pwd filesep 'figures' filesep 'cluster-' birdID filesep], clustSession{ii});
    
    nClusts = max([thisAgeSylls.type]);
    for jj = 1:nClusts
        trainClust = thisAgeSylls([thisAgeSylls.type]==jj);
        if ~isempty(trainClust), 
            figure
            hf = mosaicDRSpec(trainClust, [], 'dgram.minContrast', 1e-10, 'maxMosaicLength', 5.5, 'noroll');
            set(hf,'Name', sprintf('Trained cluster %d, age %d, bird %s', jj, currAge, birdID));
            saveCurrFigure(sprintf('%s%s_a%d_c%d-train.jpg', figureDir, birdID, currAge, jj));
            close(hf);
        end
    end
end
%}