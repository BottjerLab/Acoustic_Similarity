function examineCluster(birdID, sessFilter)

% get the files
if nargin == 0, birdID = 'R204'; end
if nargin <  2, sessFilter = '*.mat'; end    
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];
files = dir([clusterDir sessFilter]);
files = {files.name}'; files(1:2)=[]; % get rid of '.' and '..'
fPats = {'assignedClustAge-', 'clustDataAge-'};

isAssigned = strncmp(fPats{1}, files, length(fPats{1}));

clustSession = strrep(files(isAssigned), fPats{1}, '');

load([dataDir 'allSpecs-' birdID])
for ii = 1:numel(clustSession)
    if ~(exist([clusterDir fPats{1} clustSession{ii}], 'file') && ...
         exist([clusterDir fPats{2} clustSession{ii}],'file'))
        continue;
    end
    load([clusterDir fPats{1} clustSession{ii}]);
    load([clusterDir fPats{2} clustSession{ii}]);
    
    currAge = str2double(clustSession{ii}(1:2));
    %seld refers to the training set
    
    isThisAge = [DRsylls.age]==currAge;
    thisAgeSylls = DRsylls(isThisAge);
    if numel(seld) == numel(thisAgeSylls) % take the training set
        num2cell(clusterIdxs(:,end)); [thisAgeSylls.type] = ans{:}; %#ok<NOANS>
        [thisAgeSylls.inTrain] = deal(true);
    else
        isAllTrain = false(1,numel(DRsylls));
        isAllTrain(seld) = true;
        isThisAgeTrain = isAllTrain(isThisAge);
        
        nTest = sum(~isThisAgeTrain);
        % fuse the training and test set
        num2cell(clusterIdxs(:,end)); [thisAgeSylls(isThisAgeTrain).type] = ans{:}; %#ok<NOANS>
        thisAgeSylls = [thisAgeSylls(isThisAgeTrain) assignedSylls];
        num2cell([true(1,numel(seld)) false(1, nTest)]); [thisAgeSylls.inTrain] = ans{:}; %#ok<NOANS>        <- error here
    end
    fprintf('Cross-tabulation of syllables...\n');
    crosstab([thisAgeSylls.inTrain], [thisAgeSylls.type])

    figureDir = [pwd filesep 'figures' filesep clustSession filesep];
    mkdir([pwd filesep 'figures' filesep], clustSession);
    
    nClusts = max([thisAgeSylls.type]);
    for jj = 1:nClusts
        trainClusts = thisAgeSylls([thisAgeSylls.type]==jj & [thisAgeSylls.inTrain]);
        if ~isempty(trainClusts), 
            figure
            hf = mosaicDRSpec(trainClusts, [], 'dgram.minContrast', 1e-10, 'maxMosaicLength', 5.5, 'noroll');
            set(hf,'Name', sprintf('Trained cluster %d, age %d, bird %s', jj, currAge, birdID));
            saveCurrFigure(sprintf('%s%s_a%d_c%d-train.jpg', figureDir, birdID, currAge, jj));
            close(hf);
        end
        
        testClusts = thisAgeSylls([thisAgeSylls.type]==jj & ~[thisAgeSylls.inTrain]);
        if ~isempty(testClusts), 
            figure
            hf = mosaicDRSpec(testClusts, [], 'dgram.minContrast', 1e-10, 'maxMosaicLength', 5.5, 'noroll');
            set(hf,'Name', sprintf('Tested cluster %d, age %d, bird %s', jj, currAge, birdID));
            saveCurrFigure(sprintf('%s%s_a%d_c%d-test.jpg', figureDir, birdID, currAge, jj));
            close(hf);
        end
    end
end