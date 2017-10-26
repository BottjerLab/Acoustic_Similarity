function recalcAllClusters
profile on
dataDir = [pwd filesep 'data'];

% find all the bird IDs from directory
files = dir(dataDir);
birdIDs = {files([files.isdir]).name};
isID = ~cellfun('isempty',regexp(birdIDs, '^[A-Z][a-z]?\d{1,3}')); %could make tighter by using cap
birdIDs = birdIDs(isID);

for ii = 1:numel(birdIDs)   
    clusterDir = [pwd filesep 'data' filesep 'cluster-' birdIDs{ii} filesep];
    fullFiles = dir([clusterDir 'altClustDataAge-*.mat']);
    fullFiles = strcat(clusterDir, {fullFiles.name}');
    clustResultFiles = recalcClustersInnerLenWeight(birdIDs{ii},fullFiles);
    
    drawClusters(birdIDs{ii}, clustResultFiles, 'plotall');
    
end
profile viewer