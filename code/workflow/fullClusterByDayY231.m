% we are looking for the following:
% DRsylls
% featureTable
% spectra
%params = processArgs(defaultParams, 'warpingCost', 0.7);
%%
params = processArgs(defaultParams, 'warpingCost', 1.0);
birdID = 'Y231';
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
    warning('Spectra Values not cached...');
end
fprintf('Loading cached data from %s...\n', cachedResFile);
load(cachedResFile, necessaryVars{:})
%%
ages = unique([DRsylls.age]);

% note: Y231 day 2 (age 49) is stored in one of the tmpDists from
% 9/13-9/16, plz restore
% day 3 is also stored as tmpDists-09-16-11-29
% day 4 is done
% day 5 is interrupted from 9-21
startAge = 5;
for ii = startAge:numel(ages)
    % limit for test purposes/brevity
    thisAge = ages(ii);
    isThisAge = ([DRsylls.age]==thisAge);
    
    clustFolder = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
    timeFlag = datestr(clock, 'mm_dd_HH_MM');
    
    diary([clustFolder 'diary-' timeFlag '-age' num2str(thisAge) '.txt']);
    nTrain = sum(isThisAge);
    seld = find(isThisAge);

    t1 = clock;
    [clusterIdxs, empMats, distMats, empDistrs] = DRcluster(DRsylls(seld), featureTable(seld), spectra(seld), params);
    fprintf('Time for total clustering: %0.2fs\n',etime(clock, t1));
    
    
    clusterFile = [clustFolder 'altClustDataAge-' num2str(ages(ii)) '-' timeFlag];
    
    % make some figures and save work
    % note: seld indexes into the cached syllables in the allSpecs file
    save([clusterFile '.mat'], 'clusterIdxs', 'empMats', 'distMats', 'empDistrs', 'seld');
    
%%    
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
        figFileName = [clusterFile '-a' num2str(thisAge) '-c' num2str(jj) '-train.jpg'];
        fprintf('Saving figure to %s...\n',figFileName);
        saveCurrFigure(figFileName);
        close(fig);
    end
    %%
    diary off
end