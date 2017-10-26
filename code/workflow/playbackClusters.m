function playbackClusters(birdID)

% note: obsoleted by drawClusters
if nargin < 1
    birdID = 'Lb277';
end

load(['data\' birdID '\allSpecs-' birdID '.mat']);
%dataDir = [pwd filesep 'data' filesep birdID filesep];
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
%%

[allFiles dirClustFile] = uigetfile([clusterDir 'altClustDataAge-*.mat'], 'MultiSelect','on');
if ~iscell(allFiles), allFiles = {allFiles}; end
%%
%{
fullFiles = dir([clusterDir 'altClustDataAge-*.mat']);
fullFiles = {fullFiles.name};
% select one of the files
%}
for hh = 1:numel(allFiles)
    sessionID = strrep(allFiles{hh}(1:end-5),'altClustDataAge-','');
    seldAge = str2double(strtok(sessionID, '-'));
    load([dirClustFile allFiles{hh}]);
    
    seldSylls = DRsylls([DRsylls.age] == seldAge);
    
    
    if exist('clusterIdxs','var')==1
        nClusters = 4:max(clusterIdxs.cosim(:));
    else
        nClusters = 4:13;
        cosimTypes = fieldnames(distMats);
        typesToRemove = {'unwarpedLocal','global','warpedLocal'};
        cosimTypes = setxor(cosimTypes, typesToRemove);
        
        fprintf('Recovering clustering...\n');
        for ii = 1:numel(cosimTypes)
            pairLinks = linkage(distMats.(cosimTypes{ii}),'complete');
            specificIdxs = cluster(pairLinks,'maxclust',nClusters);
            clusterIdxs.(cosimTypes{ii}) = specificIdxs;
        end
    end
    nToSave = 12;
    idxsIn = nToSave == nClusters;
    
    clustTypes = fieldnames(clusterIdxs);
    seldTypes = clustTypes(listdlg('ListString', clustTypes));
    for ii = 1:nToSave
        for jj = 1:numel(seldTypes)
            subDir = sprintf('wav-%s-%02d',seldTypes{jj}, ii);
            mkdir(dirClustFile, subDir);
            clustWavDir = [dirClustFile subDir filesep];
            
            fprintf('Writing wavs for cluster %d, type %s into folder [%s]...\n', ii, seldTypes{jj}, clustWavDir);
            
            theseTypes = seldSylls(clusterIdxs.(cosimTypes{jj})(:,idxsIn)==ii);
            wavFile = [clustWavDir 'example-'];
            writeBunchWAVClip([],theseTypes(kk), wavFile);            
        end
    end
end