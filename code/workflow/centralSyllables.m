% get centrality / peripherality
birdID = 'Lb189';

clustDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir  = [pwd filesep 'data' filesep            birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
%% make all neurons psths for individual syllables for R204/age 53, R247, Lb277, giving coreness
% neuron loop - first loop over sessions, then gather spikes, then gather
report = reportOnData(birdID,'',defaultParams,'verbose',false);   

for hh = 1:numel(uAges)
    thisAge = uAges(hh);
    ageSylls = DRsylls([DRsylls.age]==thisAge);
        
    fprintf('Bird %s, age %d\n', birdID, thisAge);
    
    % get the last type of cluster labels
    acceptedLabelFile = [dataDir 'acceptedLabels-' birdID '-age' num2str(thisAge) '.mat'];
    if ~exist(acceptedLabelFile, 'file')
        fprintf('No acceptedLabels file, # sylls = %d\n', numel(ageSylls));
        continue;
    end
    clusterIdxs = []; load(acceptedLabelFile);
    clTypes = fieldnames(clusterIdxs); 
    if numel(clTypes) == 3, prefType = clTypes{2};
    else prefType = clTypes{1};
    end
    clusterIdxs = clusterIdxs.(prefType);
    fprintf('Total syllables IDed / total: (%d/%d)\n', sum(~isnan(clusterIdxs)), numel(clusterIdxs));

    % get the cluster matrices
    glob = [clustDir 'mfccClustDataAge-' num2str(thisAge) '*.mat'];
    fils = dir(glob); 
    if isempty(fils), error('No cluster file found'); end
    fils = {fils.name};
    clustFile = [clustDir fils{1}];
    load(clustFile, 'distMats')
    clusterMat = distMats.cosim;
    centralOrder = findMostCentral(clusterMat, clusterIdxs);
    
end
