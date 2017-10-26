function fullClusterByDay(birdID, startAge)
% we are looking for the following:
% DRsylls
% featureTable
% spectra
%params = processArgs(defaultParams, 'warpingCost', 0.7);
%%
params = processArgs(defaultParams, 'warpingCost', 1.0);
%birdID = 'R247';
cachedResFile = ['data' filesep birdID filesep 'allSpecs-' birdID '.mat'];
necessaryVars = {'DRsylls', 'featureTable' 'spectra'};
containedVars = who('-file',cachedResFile);
hasValue = false(1,numel(necessaryVars));
for ii = 1:numel(necessaryVars),
    hasValue(ii) = any(strcmp(necessaryVars{ii}, containedVars));
end
if ~hasValue(1)
    error('DRsylls not found...');
end
if ~all(hasValue(2:3)),
    warning('Spectra values not cached...');
end
fprintf('Loading cached data from %s...\n', cachedResFile);
load(cachedResFile, necessaryVars{:})
%%
ages = unique([DRsylls.age]);

% limit for test purposes/brevity
isThisAge = ([DRsylls.age]==startAge);

clustFolder = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
timeFlag = datestr(clock, 'mm_dd_HH_MM');

diary([clustFolder 'diary-' timeFlag '-age' num2str(startAge) '.txt']);
nTrain = sum(isThisAge);
seld = find(isThisAge);

t1 = clock;
%if numel(seld) < 5000
    [clusterIdxs, empMats, distMats, empDistrs] = ...
        DRcluster(DRsylls(seld), featureTable(seld), spectra(seld), params);
%else
%    fprintf('Running disk-write version of clustering...\n');
%    [clusterIdxs, empMats, distMats, empDistrs] = ...
%        DRcluster(DRsylls(seld), featureTable(seld), spectra(seld), params);    
%end
%%
fprintf('Time for total clustering: %0.2fs\n',etime(clock, t1));

clusterFile = [clustFolder 'altClustDataAge-' num2str(startAge) '-' timeFlag];
fprintf('Writing clustering to file %s...\n', clusterFile);
% make some figures and save work
% note: seld indexes into the cached syllables in the allSpecs file
save([clusterFile '.mat'], 'clusterIdxs', 'empMats', 'distMats', 'empDistrs', 'seld');

%%    
%{
    typedDRsylls = DRsylls(seld);
    takenIdxs = clusterIdxs(:,end);
    nTypes = max(takenIdxs);
    
    for jj = 1:nTypes
        fig = figure(jj);
        theseSylls = typedDRsylls(takenIdxs==jj);
        fprintf('Plotting trained clusters for syllable %d (#=%d)...\n',jj, sum(takenIdxs==jj));
        
        mosaicDRSpec(theseSylls, params, 'dgram.minContrast', 1e-10, ...
            'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5);
        set(fig, 'Name', sprintf('Syllable #%d, Full',jj));
        figFileName = [clusterFile '-a' num2str(startAge) '-c' num2str(jj) '-train.jpg'];
        fprintf('Saving figure to %s...\n',figFileName);
        saveCurrFigure(figFileName);
        close(fig);
    end
    %}
    %%
    diary off
