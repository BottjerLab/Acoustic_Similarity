function saveStereotypy(birdID, thisAge)
% plots stereotypical and not stereotypical syllables
clustDir = ['data' filesep 'cluster-' birdID filesep];
dataDir  = ['data' filesep            birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
if ~any(thisAge == uAges), error('saveStereotypy:AgeNotRepresented','Age %d not represented', thisAge); end

plotParams = processArgs(defaultParams, 'preroll', 5, 'postroll', 5, ...
    'dgram.minContrast', 3e-11);
%% neuron loop - first loop over sessions, then gather spikes, then gather
report = reportOnData(birdID,'',defaultParams,'verbose',false);
ageSylls = DRsylls([DRsylls.age]==thisAge);

% get the session ID for each syllable
sForSylls = cell(1,numel(ageSylls));
for ii = 1:numel(ageSylls)
    [~,sForSylls{ii}] = fileparts(ageSylls(ii).file);
end
[uSessions, ~, rIdx] = unique(sForSylls);
fprintf('Bird %s, age %d\n', birdID, thisAge);

% get the cluster labels
try
    clusterIdxs = loadAcceptedLabels(birdID, thisAge);
catch err
    if strcmp('loadAcceptedLabels:fileNotFound', err.identifier) 
        clusterIdxs = [];
    else
        rethrow(err);        
    end
    [juvSylls.type] = deal(NaN);
end

fprintf('Total syllables IDed / total: (%d/%d)\n', sum(~isnan(clusterIdxs)), numel(clusterIdxs));
if sum(~isnan(clusterIdxs)) == 0, fprintf('No IDed syllables, exiting...\n'); return; end;

% get the cluster distances matrices for ordering
fils = dir([clustDir 'altClustDataAge-' num2str(thisAge) '*.mat']);
if isempty(fils), error('saveStereotypy:noClusterFile','No cluster file found'); end
clustFile = [clustDir fils(1).name]; load(clustFile, 'distMats');

clusterMat = distMats.cosim;
[intraCDists, interCDists] = clusterDistances(clusterMat, clusterIdxs);
fprintf('Cluster distances found.\n');

for ii = 1:numel(uSessions)
    intraClusterDists = intraCDists(rIdx==ii);
    interClusterDists = interCDists(rIdx==ii);
    saveFile = ([dataDir 'labelStereotypy-' uSessions{ii} '.mat']);
    fprintf('Saving file %s...\n', saveFile);
    save(saveFile, 'intraClusterDists', 'interClusterDists');    
end
