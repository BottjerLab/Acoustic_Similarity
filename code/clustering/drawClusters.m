function drawClusters(birdID, allFiles, reportType)
%
% note: this file is to be obsoleted.  use drawClustersUpdated instead for
% now
% allFiles refers to the altClustDataAge files within the 'cluster' folder.
% the argument may be a string glob, i.e. a pattern matching string with
% asterisks
%
% reportType can be 'plotall', 'plot', and 'wav'
% 'plotall' chooses all the clusterings, automatic and manual, otherwise
% the choice is given to the user which clusterings to use
% 'plot' plots the syllables, up to 5 seconds at a time, and saves the mosaics in 
% the 'data/cluster-' folder
% 'wav' saves the syllables as 'wavs', and saves them in a 'wav' file.  

if nargin < 1
    birdID = 'Lb277';
end

%dataDir = [pwd filesep 'data' filesep birdID filesep];
clusterDir = ['data' filesep 'cluster-' birdID filesep];

if nargin < 2 || isempty(allFiles)
    [allFiles dirClustFile] = uigetfile([clusterDir 'altClustDataAge-*.mat'], 'MultiSelect','on');
else
    dirClustFile = [fileparts(allFiles{1}) filesep];
    allFiles = strrep(allFiles, dirClustFile,'');
end

if ~iscell(allFiles), allFiles = {allFiles}; end

if nargin < 3
    reportType = 'plot'; %or 'wav'
end;
load(['data\' birdID '\allSpecs-' birdID '.mat']);
%%
%%
for hh = 1:numel(allFiles)
    clustSession = strrep(allFiles{hh}(1:end-5),'altClustDataAge-','');
    seldAge = str2double(strtok(clustSession, '-'));
    fprintf('Loading cluster session %s...\n',clustSession);
    load([dirClustFile allFiles{hh}]);
    
    seldSylls = DRsylls([DRsylls.age] == seldAge);
    
    if exist('clusterIdxs','var')==1 && isfield(clusterIdxs,'cosim')        
        nClusters = 4:max(clusterIdxs.cosim(:));
    else
        nClusters = 4:13;
        cosimTypes = fieldnames(distMats);
        isCosim = ~cellfun(@isempty, strfind(cosimTypes, 'cosim')) | ~cellfun(@isempty, strfind(cosimTypes, 'Cosim')) ;
        cosimTypes = cosimTypes(isCosim);
        
        fprintf('Recovering clustering...\n');
        for ii = 1:numel(cosimTypes)
            pairLinks = linkage(distMats.(cosimTypes{ii}),'complete');
            specificIdxs = cluster(pairLinks,'maxclust',nClusters);
            clusterIdxs.(cosimTypes{ii}) = specificIdxs;
        end
    end
    
    nToAnalyze = 12;
    idxsIn = (nToAnalyze == nClusters);
    clusterFile = [dirClustFile filesep 'expClusters-' clustSession '-'];
    
    params = defaultParams;
    
    if hh == 1 
        if strcmpi('plotall', reportType)
            seldTypes = fieldnames(clusterIdxs);
        else
            clustTypes = fieldnames(clusterIdxs);
            seldTypes = clustTypes(listdlg('ListString', clustTypes));
        end
    end
    
    if any(strncmpi('plot',reportType,4))
        for ii = 1:numel(seldTypes) % different kinds of clustering type
            for jj = 1:nToAnalyze % number of clusters
                fig = figure(jj);
                fprintf('Plotting %s clusters for syllable %d...\n', seldTypes{ii}, jj);
                
                theseSylls = seldSylls(clusterIdxs.(seldTypes{ii})(:,idxsIn)==jj);
                mosaicDRSpec(theseSylls, params, ...
                    'dgram.minContrast', 1e-11, 'doFilterNoise', false,...
                    'preroll', 5, 'postroll', 5, 'maxMosaicLength', 4.5);
                set(fig, 'Name', sprintf('Syllable #%d, %s', jj, seldTypes{ii}));
                figFilName = [clusterFile seldTypes{ii} '-c' num2str(jj) '.jpg'];
                saveCurrFigure(figFilName);
                close(fig);
            end
        end
    end
    if any(strcmpi('wav',reportType))
        for ii = 1:nToAnalyze
            for jj = 1:numel(seldTypes)
                subDir = sprintf('wav-%s-%s-%02d',clustSession, seldTypes{jj}, ii);
                mkdir(dirClustFile, subDir);
                clustWavDir = [dirClustFile subDir filesep];
                
                fprintf('Writing wavs for cluster %d, type %s into folder [%s]...\n', ii, seldTypes{jj}, clustWavDir);
                
                theseSylls = seldSylls(clusterIdxs.(seldTypes{jj})(:,idxsIn)==ii);
                wavFile = [clustWavDir 'example-'];
                writeBunchWAVClip([],theseSylls, wavFile);
            end
        end
    end
end