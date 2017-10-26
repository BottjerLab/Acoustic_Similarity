%% step 1: gather regions (in this case, syllables)
birdID = 'Lb277';
dataPath = [pwd filesep 'data' filesep birdID filesep];

catalog = reportOnData(birdID, '', [], 'verbose',false);
%%
features = cell(1,numel(catalog));
sessionHasSpikeData = ~cellfun('isempty',{catalog.spikeFiles});
% features is a horz cell vector of feature struct arrays 
% these are 'reduced features', i.e. features that capture a statistical
% reduction of the entire timeseries for that 
% within struct arrays, each syllable is represented once
%%
DRsylls = cell(1,numel(catalog));
for ii = 1:numel(catalog)
    iSession = catalog(ii);
    hasData = any(findInManifest(iSession.manifest, {'bosSyllables','manualSyllables','syllables'}));
    if ~hasData, continue; end;
 
    if ~sessionHasSpikeData(ii), continue; end;
    
    % sources = {'bosSyllables', 'manualSyllables', 'syllables'};
    sources = {'approvedSyllables'}; 
    jj = 0;
    while jj + 1 <= numel(sources) && isempty(DRsylls{ii})
        jj = jj + 1;
        DRsylls{ii} = loadFromManifest(iSession.manifest, sources{jj});
    end
    loadedSource = sources{jj};
    
    if isempty(DRsylls{ii}) % no variable, no plot
        warning('findInManifest fell through or value is empty, couldn''t load...');
        continue;
    end
    
    % add metadata for syllables
    keyFile = [dataPath iSession.sessionID '.mat'];
    assert(2==exist(keyFile,'file'));
    [DRsylls{ii}.file] = deal(keyFile);
    [DRsylls{ii}.age] = deal(getAgeOfSession(iSession.sessionID));
    
    % put in column form
    if size(DRsylls{ii},1) == 1 && size(DRsylls{ii},2) > 1, DRsylls{ii} = DRsylls{ii}'; end
    
    fprintf('Loaded session %s (%d/%d), source %s (# = %d)...\n', iSession.sessionID, ii, numel(catalog), loadedSource, numel(DRsylls{ii}));
    %{
    noiseMask = loadFromManifest(iSession.manifest, 'noiseMask');
    % step 2: collect relevant features for each region
    params = processArgs(defaultParams, birdID, 'preroll',1,'postroll',1, ...
        'doFilterNoise', true, 'noiseFilter', noiseMask, ...
        'fs', 1/getfield(loadFromManifest(iSession.manifest,'metaStruct'),'interval'));
    
    features{ii} = getFeatures([], DRsylls{ii},params);
    %}
end

%% step 3: collect spiking data for each section - times and rates
spikes = cell(1,numel(catalog)); % each spike 
% spikes is a horz cell vector of cell arrays spike times within a segmented syllable
% within each cell, the inner cell array represents each syllable crossed with
% each neuron 
for ii = 1:numel(catalog)
    iSession = catalog(ii);
    fprintf('Loading spikes from session %s (%d/%d)...\n',iSession.sessionID, ii, numel(catalog));
    
    spikeTrains = loadSpikeData(iSession.spikeFiles); % we don't look at if they're MUA or not...
    % skip if there aren't any spikes...
    if isempty(spikeTrains), continue; end;
        
    % NOTE: for all subsequent analyses here, we treat all units in LMAN as equivalent and
    % from the same 'big neuron'
    bigNeuronSpikeTrain = vertcat(spikeTrains{:});
    [~,spikes{ii}] = countSpikes(DRsylls{ii}, bigNeuronSpikeTrain);    
end

%% step 3b: combine syllables across sessions
allFeatures = [features{sessionHasSpikeData}];
allSpikes = [spikes{:}];
allDRsylls = vertcat(DRsylls{sessionHasSpikeData});
%% step 3c: save regions associated with these features
featureSpikeFile = [dataPath 'featureSpike-' birdID '.mat'];

save(featureSpikeFile,'allFeatures','allSpikes','allDRsylls');
%uisave({'allFeatures','allSpikes','allDRsylls'},featureSpikeFile);
%% step 3C: load regions associated with these features
featureSpikeFile = [dataPath 'featureSpike-' birdID '.mat'];

%save(featureSpikeFile,'allFeatures','allSpikes','allDRsylls');
load(featureSpikeFile);

%% step 4a: correlate features to spiking activity
nSpikes = cellfun('length',allSpikes);
[pCorr, tCorr] = findCorrelates(allFeatures, nSpikes ./ [allFeatures.length]);

% find significance level for number of tests
alpha = 0.05;
pBonf = 1-(1-alpha).^(1/numel(fieldnames(allFeatures)));
tBonf = norminv(1-pBonf);

tCorr(isnan(tCorr)) = 0;
reducedStatTable(fieldnames(allFeatures), tCorr);
title(sprintf('signed t-values of spiking correlations (# sylls = %d), sig @ \\pm%0.2f',numel(allDRsylls),tBonf));

% based on jma's preference in colors
coralCol = [255 127 80]/255;
sfGreenCol = [127 255 212]/255;
jmacm = [interpColormap(coralCol, [0.25,0.25,0.25], 128); ...
    flipud(interpColormap(sfGreenCol, [0.25,0.25,0.25], 128))];
colormap(jmacm);

%% step 4b: PCA syllable features
% represent each syllable now as an anonymous list of features in
% allFeaturesTable, which has one row per syllable and one column per
% feature

fn = fieldnames(allFeatures);
allFeaturesTable = cellfun(@(x) [allFeatures.(x)]', fn', 'UniformOutput',false);
allFeaturesTable = [allFeaturesTable{:}];
%%
% the actual pca
[pc,score,latent,tsquare] = princomp(zscore(allFeaturesTable));
% find 99.5% of the variance captured
thresh = 0.995; 
nSuffDims = find(cumsum(latent) / sum(latent) > thresh,1);
fprintf('Reduced from %d -> %d dims...\n',numel(fn), nSuffDims);
reducedFeatures = score(:,1:nSuffDims);
%% plotting some pca analysis
figure; 
subplot(211); 

semilogy(1:numel(latent),latent,'k-'); 
ylabel('comp. eigenvector magnitudes');
xlabel('# components');
subplot(212); 
plot(1:numel(latent),cumsum(latent)/sum(latent),'k-');
ylabel('fraction of explained variance')
xlabel('# components');
%% first try on clustering: simple k-means on features
nObs = numel(allFeatures);

nClusters = 8;
[kClustIdx, clustCenters, distToAssignCluster, distToAllClusters] = kmeans(reducedFeatures,nClusters,...
    'distance', 'sqEuclidean');

reductionCoeff = mean(distToAllClusters(sub2ind([nObs,nClusters],1:nObs,kClustIdx'))./sum(distToAllClusters,2)');
figure; 
cols = jet(nClusters);
for ii = 1:nClusters
    isClust = (kClustIdx == ii);
    plot3(reducedFeatures(isClust,1), reducedFeatures(isClust,2), reducedFeatures(isClust,3), '.', 'Color', cols(ii,:));
    hold on;
end
hold off;
fprintf('Ratio of distances minimized = %0.3f\n', reductionCoeff);

fprintf('Quit me here!');
keyboard
%% step 4.5: subselect syllables for further clustering
% pick a random sample of 1000 syllables
nForCluster = min(ceil(numel(allDRsylls)/2),1000);
trainIdxs = randperm(numel(allDRsylls)); trainIdxs = trainIdxs(1:nForCluster);
clustSylls = allDRsylls(trainIdxs);

fieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'};
% store the feature-based spectra for all of them
params.fine.features = {'wienerEntropy','deriv','harmonicPitch','fundamentalFreq'};

progressbar(sprintf('Calculating features for sample (# = %d)',nForCluster));
for ii = 1:nForCluster
    tmpSpec = getMTSpectrumStats(getClipAndProcess([],clustSylls(ii),params), params.fine);
    for jj = 1:numel(fieldsToKeep)
        spec(ii).(fieldsToKeep{jj}) = tmpSpec.(fieldsToKeep{jj});
    end
    progressbar(ii/nForCluster);
end
%% calculate local distances, fixed-time version, no
%entry = find(strcmp('aperiodicity',{params.featureCatalog.}
yn = input('Would you like to load pre-loaded unwarped similarity from the file, y/n [no by default]? ','s');
if strncmpi(yn,'y',1)
    load([dataPath 'localSim-' birdID '.mat'], 'distM');
else
    yn = input('Would you like to load pre-loaded unwarped similarity from the base workspace, y/n [yes by default]? ','s');
    if strncmpi(yn,'n',1),
        distM = zeros(1,ii*(ii-1)/2);
    end
end

innerIdx = 0;
progressbar('Unwarped Distance Calcs');
for ii = 1:nForCluster-1
    for jj = ii+1:nForCluster
        innerIdx = innerIdx + 1;
        distM(innerIdx) = min(standardDistance(spec(ii), spec(jj)));
        % not sure why this is needed
        if nForCluster-ii==1, continue; end;
        progressbar(innerIdx/nchoosek(nForCluster,2));
    end
end
progressbar(1);
save([dataPath 'localSimSyllsOnly-' birdID '.mat'],'clustSylls');
save([dataPath 'localSim-' birdID '.mat'],'distM');
%% calculate local distances, TIME WARPED version
%entry = find(strcmp('aperiodicity',{params.featureCatalog.}
tic
progressbar('Saves','Time Warped Distance Calcs');
twDistM = NaN(1,nchoosek(nForCluster,2));

innerIdx = 0;
for ii = 1:nForCluster-1
    for jj = ii+1:nForCluster
        innerIdx = innerIdx + 1;

        twDistM(innerIdx) = timeWarpedDistance(spec(ii), spec(jj));
        if ii==nForCluster-1, continue; end;
        
        progressbar([],innerIdx/nchoosek(nForCluster,2));
        
        if rem(innerIdx, floor(sqrt(nchoosek(nForCluster,2)))) == 0
            save([dataPath 'localSimTW-' birdID '.mat'],'clustSylls','twDistM');
            progressbar(floor(innerIdx/floor(sqrt(nchoosek(nForCluster,2)))) / ...
                floor(nchoosek(nForCluster,2)/floor(sqrt(nchoosek(nForCluster,2)))))
        end
    end
end

progressbar(1);
toc

save([dataPath 'localSimTW-' birdID '.mat'],'clustSylls','twDistM');
%% step 5: measure global distances within pairs of syllables

%seldFeaturesTable = allFeaturesTable;
%clustSylls = allDRsylls(trainIdxs);
seldFeaturesTable = allFeaturesTable(trainIdxs,:);

% start with unnormalized table of features
% step 1: normalize to z-scores
fprintf('Calculating global dissimilarity scores...\n');
zNormedFeatures = zscore(seldFeaturesTable);
seldGlobalDistVec = pdist(zNormedFeatures);

globalDistVec = pdist(zscore(allFeaturesTable));

%% get non-parametric (probability-rank) ordering of local similarity
fprintf('Calculating empirical local dissimilarity scores');
[~,rord] = sort(twDistM);
twDistPvals(rord) = [1:numel(twDistM)] / numel(twDistM);
%% get non-parametric (probability-rank) ordering of global similarity

nBoot = 10000;
bootstrapSubseld = randi(numel(globalDistVec), nBoot);
[~,mini] = min(globalDistVec); [~, maxi] = max(globalDistVec);
bootstrapSubseld(1:2) = [mini maxi];

fprintf('Constructing empirical cdf for global dissimilarity from whole-session sample...\n');
[empSimDistrY,empSimDistrX] = ecdf(globalDistVec(trainIdxs));
empSimDistrX(1) = 0;
globalDistPVals = interp1(empSimDistrX, empSimDistrY,seldGlobalDistVec);
globalDistPVals(isnan(globalDistPVals)) = 1; % maximum dissimalirity
%% construct co-similarity as fusion of local and global p-values
fprintf('Calculating co-dissimilarity (correlation of dissimilarities, which is a similarity score)...\n');
fusedPVals = sqrt(twDistPvals .* globalDistPVals);
coSimVector = pdist(squareform(fusedPVals), 'correlation');

%%
subplot(2,2,1);
imagetext(triu(squareform(globalDistPVals)), false, [], get(gca,'Position'), false);
set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Global distances')
subplot(2,2,2);
imagetext(triu(squareform(twDistPvals)), false, [], get(gca,'Position'), false);
set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Local time-warped distances');
subplot(2,2,3);
imagetext(triu(squareform(fusedPVals)), false, [], get(gca,'Position'), false);
set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Fused distances');
subplot(2,2,4);
imagetext(triu(squareform(coSimVector)), false, [], get(gca,'Position'), false);
set(gca,'XTickLabel',[],'YTickLabel',[]);
title('Co-similarity distances');

%% save
batchSimFile = [pwd filesep 'data' filesep birdID filesep 'batchSimilarity-' birdID '.mat'];

uisave({'clustSylls','twDistM','twDistPvals','seldGlobalDistVec','globalDistPVals','coSimVector'},batchSimFile);
%% step 5.b: cluster on co-similarity (note: this step is fast)
nClusters = 4:10;
pairLinks = linkage(coSimVector,'complete');
clustIdxs = cluster(pairLinks,'maxclust',nClusters);
%% visualize: clusters
typedDRsylls = allDRsylls(trainIdxs); 
num2cell(clustIdxs(:,7)); [typedDRsylls.type] = ans{:};

for ii = 1:max(nClusters)
    figure(ii);
    fprintf('Plotting syllable %d...\n',ii);
    theseSylls = typedDRsylls([typedDRsylls.type]==ii);
    cumLens = cumsum([theseSylls.stop]-[theseSylls.start]);
    
    % get some user input
    %timeLimit = input(sprintf('Cluster #%d has %0.2fs of syllables, how many (s) would you like to show? (default is 5s.)> ',...
    %    ii,sum([theseSylls.stop] - [theseSylls.start])));
    %if isempty(timeLimit), timeLimit = 5.0; end;
    timeLimit = 5.0; %seconds
    
    seld = find([cumLens Inf] >= timeLimit, 1, 'first') - 1; % get 3.8 s worth of syllables
    
    mosaicDRSpec(theseSylls, params, 'dgram.minContrast', 1e-10, ...
        'preroll', 3, 'postroll', 3);
    title(sprintf('Syllable #%d',ii)); 
    %pause;
    %{
    exportPath = [dataPath 'wav' filesep sprintf('cluster%02d',ii) filesep sprintf('%s-c%02d', birdID, ii)];
    fprintf('Writing waves to %s...\n',exportPath);  
    writeDRWAVClips(theseSylls, exportPath);
    pause;
    %}
end
%% check ages of different clusters
[uFiles,~,uFilesrIdx] = unique({typedDRsylls.file});
uSessions = strrep(strrep(uFiles, dataPath, ''),'.mat','');
uAges = getAgeOfSession(uSessions);
num2cell(uAges(uFilesrIdx)); [typedDRsylls.age] = ans{:};

%% now: assign test syllables to clusters and see how well they do

