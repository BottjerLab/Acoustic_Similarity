function assignClustersToTutor(birdID, thisAge)
% pulls up the figures for the tutor and the different clusters, and 
% takes user input to decide which tutor syllable matches those clusters
dataDir  = ['data' filesep birdID filesep];
% open viewer for the tutor 
imshow(['figures' filesep 'tutor-' birdID '.jpg']);
set(gca,'Position', [0 0 1 1]);
%% 
report = reportOnData(birdID,'',defaultParams,'verbose',false);
ageSylls = loadAgeSylls(birdID, thisAge); 

% get the session ID for each syllable
sForSylls = cell(1,numel(ageSylls));
for ii = 1:numel(ageSylls)
    [~,sForSylls{ii}] = fileparts(ageSylls(ii).file);
end
uSessions = unique(sForSylls);
fprintf('Bird %s, age %d\n', birdID, thisAge);

%% get the appropriate type of cluster labels
clusterIdxs = loadAcceptedLabels(birdID, thisAge);
fprintf('Total syllables IDed / total: (%d/%d)\n', sum(~isnan(clusterIdxs)), numel(clusterIdxs));

nTypes = max(clusterIdxs);
typeMatch = cell(1,nTypes);
hf = figure;
for kk = 1:nTypes % loop through syllable types    
    figure(hf)
    sType = kk;
    thisCluster = ageSylls(clusterIdxs == sType);
    
    fprintf('Type %d, # = %d\n', kk, numel(thisCluster));
    
    mosaicDRSpec(thisCluster, defaultParams, ...
        'dgram.minContrast', 1e-11, 'doFilterNoise', false,...
        'preroll', 5, 'postroll', 5, 'maxMosaicLength', 2.5);
   typeMatch{kk} = input('Which cluster?', 's');    
end
%% 
save([dataDir 'matchToTutor-age' num2str(thisAge) '.mat'],'typeMatch')