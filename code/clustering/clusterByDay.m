% we are looking for the following:
% DRsylls
% featureTable
% spectra
params = defaultParams;

birdID = 'Lb277';
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

ages = unique([DRsylls.age]);

for ii = 1:numel(ages)
    % limit for test purposes/brevity
    thisAge = ages(ii);
    isThisAge = ([DRsylls.age]==thisAge);
    
    nTrain = 800; % JS, 082913: 500 seems to be not enough, also random sampling may not be great
    seld = zeros(1,nTrain);
    if nTrain > sum(isThisAge)
        nTrain = sum(isThisAge);
        seld = find(isThisAge);
    else
        % boring stratified random sampling code
        [uSessions, ~, rSessionIdx] = unique({DRsylls(isThisAge).file});
        nSessions = numel(uSessions);
        nPerSession = histc(rSessionIdx, 1:nSessions);
        
        % divvy up the quotas
        quotas = floor(nTrain / nSessions) * ones(1, nSessions);
        if rem(nTrain, nSessions) > 0
            randExcess = randperm(nSessions, rem(nTrain, nSessions));
            quotas(randExcess) = quotas(randExcess)+1;
        end
        
        % while at least one session can't meet quota, then shuffle the remainder
        % evenly and continue
        while any(quotas > nPerSession)
            overQuota = (quotas > nPerSession);
            excess = sum(quotas(overQuota) - nPerSession(overQuota));
            quotas(overQuota) = nPerSession(overQuota);
            
            isFull = (quotas == nPerSession);
            quotas(~isFull) = quotas(~isFull) + evenRDistr(excess, sum(~isFull));
        end
        bounds = [0 cumsum(quotas)];
        
        % get seld indices
        for jj = 1:nSessions
            inSession = find(rSessionIdx==jj);
            seld((bounds(jj)+1):bounds(jj+1)) = inSession(randperm(nPerSession(jj), quotas(jj)));
        end
        
        agesIdxs = find(isThisAge);
        seld(seld==0)=[];
        seld = agesIdxs(seld); % move index from smaller (thisAge) to larger (DRsylls) scope
    end
    t1 = clock;
    [clusterIdxs, empMats, distMats, empDistrs] = DRcluster(DRsylls(seld), featureTable(seld), spectra(seld));
    fprintf('Time for total clustering: %0.2fs\n',etime(clock, t1));
    
    clustFolder = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
    timeFlag = datestr(clock, 'mm_dd_HH_MM');
   
    clusterFile = [clustFolder 'clustDataAge-' num2str(ages(ii)) '-' timeFlag];
    
    % make some figures and save work
    % note: seld indexes into the cached syllables in the allSpecs file 
    save([clusterFile '.mat'], 'clusterIdxs', 'empMats', 'distMats', 'empDistrs', 'seld'); 
    
    typedDRsylls = DRsylls(seld);
    takenIdxs = clusterIdxs(:,end);
    nTypes = max(takenIdxs);
    
    for jj = 1:nTypes
        fig = figure(jj);
        fprintf('Plotting trained clusters for syllable %d...\n',jj);
        theseSylls = typedDRsylls(takenIdxs==jj);
        
        mosaicDRSpec(theseSylls, params, 'dgram.minContrast', 1e-10, ...
            'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5);
        set(fig, 'Name', sprintf('Syllable #%d, Training',jj));
        saveCurrFigure([clusterFile '-a' num2str(thisAge) '-c' num2str(jj) '-train.jpg']);
        close(fig);
    end
    
    tmpLblSylls = DRsylls; [tmpLblSylls.type] = deal(NaN);
    num2cell(takenIdxs); [tmpLblSylls(seld).type] = deal(ans{:}); %#ok<NOANS>
    
    % pick the most central exemplars from each category
    nEx = 5;
    exemplars = NaN(1, nEx * nTypes);
    for jj = 1:nTypes
        isJType = (takenIdxs == jj); jType = find(isJType);
        mostCentralList = sort(sum(empMats(isJType, isJType)));
        mostCentral = jType(mostCentralList(1:min(end,nEx)));        
        exemplars((nEx*jj-nEx)+(1:numel(mostCentral))) = mostCentral;
    end
    exemplars(isnan(exemplars))=[];
    
    isTestSeld = isThisAge & isnan([tmpLblSylls.type]); % pick all unlabeled syllables of this age
    isTestSeld(seld(exemplars)) = true; % pick up all the exemplars

    seldGrid = false(nTrain); seldGrid(exemplars, exemplars) = true; 
    seldGrid(1:nTrain+1:end) = false; % set diagonal to 0
    seldVec = squareform(seldGrid);
    
    fn = fieldnames(empMats);
    for jj = 1:numel(fn)
        testEmpMatrices.(fn{jj}) = empMats.(fn{jj})(seldVec);
    end
    
    [assignedSylls, cosims] = DRassignCluster(...
        tmpLblSylls(isTestSeld),testEmpMatrices, empDistrs, ...
        featureTable(isTestSeld), spectra(isTestSeld));
    
    assignedFile = [clustFolder 'assignedClustAge-' num2str(ages(ii)) '-' timeFlag];
    fprintf('Assignment complete! Saving to %s.\n', assignedFile);
    
    % make some figures and save work
    save([assignedFile '.mat'], 'assignedSylls', 'cosims');

    for jj = 1:nTypes
        fig = figure(jj);
        fprintf('Plotting trained clusters for syllable %d...\n',jj);
        theseSylls = assignedSylls([assignedSylls.type]==jj);
        
        mosaicDRSpec(theseSylls, params, 'dgram.minContrast', 1e-10, ...
            'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5);
        set(fig, 'Name', sprintf('Syllable #%d, Testing',jj));
        saveCurrFigure([clusterFile '-a' num2str(thisAge) '-c' num2str(jj) '-test.jpg']);
        close(fig);
    end
end