function [acceptedLabels, augmentedLabels] = browseAndAccept(birdID)
% interactive part of clustering post-machine step
% several prompts based on the output results from recalcClusters* family
%%
close all;
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir    = [pwd filesep 'data' filesep            birdID filesep];

plotParams = processArgs(defaultParams,...
                            'dgram.minContrast', 1e-11, 'doFilterNoise', false,...
                            'preroll', 3, 'postroll', 3);
        
% load the subset of syllables
fprintf('Loading spectra for bird %s...', birdID);
DRsylls = []; featureTable = []; spectra = [];

load([dataDir 'allSpecs-' birdID '.mat'], 'DRsylls','featureTable');
fprintf(' done loading.\n');

% choose the cluster file
[filName, pathName] = uigetfile([clusterDir '*.mat'],'Select the file to grab clusters from');
fprintf('Loading cluster file...');
%%
% parse the age, subselect the syllables to save memory
[~,foooo] = strtok(filName,'-');
thisAge = sscanf(foooo(2:end), '%d'); thisAge = thisAge(1); % make me less hacky.......... use a regexp
clustSession = strrep(filName(1:end-5), 'altClustDataAge-','');
%%
isAge = ([DRsylls.age] == thisAge);
DRsylls = DRsylls(isAge);
featureTable = featureTable(isAge);

%% load cluster file
clusterIdxs = [];
load([pathName filesep filName], 'clusterIdxs');
fprintf('Now loading distance matrix for additional refinement... ');
load([pathName filesep filName], 'distMats');
fprintf('done loading.\n');

%%
% pick the kind of clustering
if isstruct(clusterIdxs)
    fn = fieldnames(clusterIdxs);
    seldClusterType = fn{listdlg('ListString', fn, 'SelectionMode', 'single','Name', 'Which type of clustering?')};
    if isempty(seldClusterType), error('Did not select option, buggin'' out.'); end
%    cIdxs = clusterIdxs.(seldClusterType);
else
%    cIdxs = clusterIdxs;
    seldClusterType = 'cosim';
end
%maxes = max(cIdxs,[],1);
%%
dists = squareform(distMats.(seldClusterType)); clear distMats;

%% pick which level of clustering you want:
nDefClusts = 12;
maxes = 4:30;
%%
linktree = linkage(dists, 'complete');
cIdxs = cluster(linktree,'maxclust',maxes);
%%
nClusts = str2double(inputdlg(...
    sprintf('How many clusters are you looking for? (%d-%d)',min(maxes),max(maxes)),...
    'Number of clusters', 1, {num2str(nDefClusts)}));

thisSetPtr = find(maxes==nClusts);

thisSetIdxs = cIdxs(:,thisSetPtr);
counts = hist(thisSetIdxs,1:nClusts);

%% loop through the different clusters in order of decreasing size
[~, sortPerm] = sort(counts, 'descend');
acceptedLabels = -ones(size(thisSetIdxs));

isFused = false(1,nClusts);
%%
grpCtr = 1;
for ii = 1:nClusts
    if isFused(sortPerm(ii)), continue; end;
    %% inspect cluster
    clf    
    mosaicDRSpec(DRsylls(thisSetIdxs == sortPerm(ii)), plotParams, 'maxMosaicLength', 4.5);
    set(gcf,'Name',sprintf('Cluster #%d (%d/%d)',sortPerm(ii),ii,nClusts));    
    
    retry = true;
    while retry
        retry = false;
        switch nm_questdlg('How to treat this cluster?', ...
                sprintf('Cluster %d, # = %d', sortPerm(ii), sum(thisSetIdxs == sortPerm(ii))),...
                'Accept', 'Reject', 'Edit', 'Accept')
            case 'Accept'
                acceptedLabels(thisSetIdxs == sortPerm(ii)) = grpCtr;
                grpCtr = grpCtr + 1;
            case 'Reject'
                acceptedLabels(thisSetIdxs == sortPerm(ii)) = -1; % sentinel for rejection
                continue;
            case 'Edit'
                % todo: report metrics of conformity/similarity?
                [newLabelGroups, retry] = editSub(thisSetIdxs, sortPerm(ii));
                % implement labeling groups
                for jj = 1:numel(newLabelGroups)
                    acceptedLabels(newLabelGroups{jj}) = grpCtr;
                    grpCtr = grpCtr + 1;
                end
        end
    end
end

acceptedLabels(acceptedLabels == -1) = NaN;

% save to an acceptedLabels file
clusterIdxs = struct('accepted', acceptedLabels);
saveFileName = [pwd filesep 'data' filesep birdID filesep 'acceptedLabels-' birdID '-age' num2str(thisAge) '.mat'];
fprintf('Saving cluster identifies to %s...', saveFileName)
save(saveFileName, 'clusterIdxs');

%% step 1.5: encourage review and merging of labels
% goal: make more representative syllables
regroupedLabels = mergeClustersByHand(DRsylls, dists, acceptedLabels);

%% second pass: try to match labels to establish groups
% find matches by looking at the distance
% loop through the unlabeled syllables
% look for matches within the range of acceptable distances
augmentedLabels = semiLabelSyllables(DRsylls, dists, regroupedLabels);

%% save to an acceptedLabels file
clusterIdxs = struct('accepted', acceptedLabels, 'regrouped', regroupedLabels, 'augmented', augmentedLabels);
saveFileName = [pwd filesep 'data' filesep birdID filesep 'acceptedLabels-' birdID '-age' num2str(thisAge) '.mat'];
fprintf('Saving cluster identifies to %s...', saveFileName);
save(saveFileName, 'clusterIdxs');
fprintf('done. hooray. \n');
%%
    function [newLabelGrps, tryAgain] = editSub(labels, oldLabel)
        nThisType = sum(labels == oldLabel);
        ttl = sprintf('Cluster %d, # = %d', oldLabel, nThisType);
        
        openFigures = [];
        
        tryAgain = false;
        newLabelGrps = {};
        switch nm_questdlg('How to edit this cluster?', ttl, 'Merge', 'Split', 'Sort by Ear', 'Merge')            
            case 'Merge'
                % follow the tree backwards
                mergedCluster = findMerge(thisSetPtr, oldLabel, cIdxs);
                
                % now presentIden and clustIden should fuse
                fprintf('Fusing candidate...');
                
                figureName = sprintf('%sexpClusters-%s-%s-c%d.jpg', pathName, clustSession, ...
                    seldClusterType, mergedCluster);

                if exist(figureName,'file') == 2;
                    fOp = figure; openFigures = [openFigures fOp];
                    imshow(figureName);
                    set(gcf,'Name',[figureName ' - candidate cluster for merging']);
                else
                    
                end
                switch nm_questdlg('Fuse?',ttl,'Yes','No, reject all', 'Retry', 'Yes')
                    case 'Yes'
                        isFused(oldLabel) = true;
                        isFused(mergedCluster) = true;
                        newLabelGrps = [find(labels == oldLabel | labels == mergedCluster) newLabelGrps];
                    case 'Accept original'
                        newLabelGrps = [find(labels == oldLabel) newLabelGrps];                        
                    case 'No, reject'
                        % do nothing
                    case 'Retry'
                        tryAgain = true;
                end
            case 'Split'
                % how do we split?
                splitFeatOptions = ['tree'; fieldnames(featureTable)];
                % follow the tree forwards
                
                repeatSplit = true;
                while repeatSplit % loop around try-catchtry for clustering robustness
                    repeatSplit = false;
                    [splitOpts, hitOk] = listdlg('ListString', splitFeatOptions, ...
                        'Name', 'Which features to split?');
                    
                    if ~hitOk
                        tryAgain = true; 
                        closeAll(openFigures);
                        return; 
                    end
                    
                    if any(splitOpts == 1)
                        [typeA, typeB] = findSplit(thisSetPtr, oldLabel, cIdxs);                        
                    else
                        % try a two-component mixture of gaussian clusters
                        nFeats = numel(splitOpts);
                        stats = zeros(nThisType, nFeats);
                        omitFeature = false(1, nFeats);
                        varStat = zeros(1,nFeats);
                        for kk = 1:nFeats
                            stats(:,kk) = [featureTable(labels == oldLabel).(splitFeatOptions{splitOpts(kk)})];                            
                        end
                        varStat = var(stats);
                        omitFeature = (varStat < eps(max(varStat))*nThisType);
                        
                        if all(omitFeature), warning('All features are constant...'); end;
                        if any(omitFeature)
                            omitStr = '';
                            omitList = splitOpts(omitFeature);
                            
                            % create the join string
                            for kk = 1:numel(omitList),
                                omitStr = strcat(omitStr , splitFeatOptions(omitList(kk)));
                                if kk < nFeats, omitStr = strcat(omitStr,', '); end
                            end
                            omitStr = omitStr{1};
                            fprintf('Removing features %s...\n', omitStr);
                            
                            stats(:,omitFeature) = [];
                            splitOpts(omitFeature) = [];
                            nFeats = numel(splitOpts);
                        end
                        
                        try
                            fitObj = gmdistribution.fit(stats,2,'Regularize',0.0001);
                            newIdxs = cluster(fitObj,stats);
                            
                            % do some plotting
                            fTab = figure('Name','Clustergram');
                            openFigures = [openFigures fTab];                            
                            if size(stats,2) == 1
                                nBins = 40; bins = zeros(1,nBins);                                
                                bins(2:end-1) = linspace(min(stats),max(stats),nBins-2);
                                bins(end) = 2 * bins(end-1) - bins(end-2);
                                bins(1) = 2 * bins(2) - bins(3);
                                bar(bins,histc(stats, bins),1); 
                                xlim(bins([1 end]));
                                hold on;
                                plot(bins, pdf(fitObj,bins'),'r-', 'LineWidth', 2);
                                hold off;
                                xlabel(splitFeatOptions{splitOpts}, 'Interpreter','none'); ylabel('Count');
                            else
                                % pick the two most diagnostic features
                                dprime = zeros(1,nFeats);
                                for kk = 1:size(stats,2)
                                    dMu = diff(fitObj.mu(:,kk));
                                    sumSigma = sqrt(sum(fitObj.Sigma(kk,kk,:)));
                                    dprime(kk) = dMu/sumSigma;
                                end
                                [~,bestDims] = sort(dprime, 'descend'); bestDims = bestDims(1:2);
                                plot(stats(:,bestDims(1)), stats(:,bestDims(2)),'k.');
                                xlabel(splitFeatOptions(splitOpts(bestDims(1))), 'Interpreter','none');
                                ylabel(splitFeatOptions(splitOpts(bestDims(2))), 'Interpreter','none');
                            end
                            
                            % split done
                            subset = find(labels == oldLabel);
                            typeA = subset(newIdxs == 1);
                            typeB = subset(newIdxs == 2);
                        catch err
                            repeatSplit = questdlg(['Error in clustering: [', err.message, ']; try again?'],...
                                'Clustering Error', ...
                                'Yes','No','Yes');
                            repeatSplit = strcmp('Yes', repeatSplit);                            
                        end                        
                        
                    end
                        
                end                
                % give option to view/sing/approve blind
                switch nm_questdlg(sprintf('How to review the split (A = %d,B = %d)?', ...
                        numel(typeA), numel(typeB)), ...
                        ttl,'View','Listen','Continue w/o review', 'View')
                    case 'View'                                        
                        totalLenA = sum([DRsylls(typeA).stop] - [DRsylls(typeA).start]);
                        totalLenB = sum([DRsylls(typeB).stop] - [DRsylls(typeB).start]);
                        
                        defaultVal = min(totalLenA, 4.0);
                        tPrev = nm_inputdlg(...
                            sprintf('Number of seconds to preview (max %.1fs): ', totalLenA),...
                            'Cluster A', 1, {num2str(defaultVal)});
                        
                        if isempty(tPrev), tryAgain = true; return; 
                        else tPrev = str2double(tPrev); end

                        fprintf('Plotting split cluster A (# = %d)...\n', numel(typeA));                                                                
                        figA = figure;
                        openFigures = [openFigures figA];
                        
                        mosaicDRSpec(DRsylls(typeA), plotParams, 'maxMosaicLength', tPrev);
                        set(gcf,'Name',sprintf('Cluster A (# = %d)',numel(typeA)));

                        defaultVal = min(totalLenB, tPrev);
                        tPrev = nm_inputdlg(...
                            sprintf('Number of seconds to preview (max %.1fs): ', totalLenB),...
                            'Cluster B', 1, {num2str(defaultVal)});
                        
                        if isempty(tPrev), tryAgain = true; return; 
                        else tPrev = str2double(tPrev); end

                        fprintf('Plotting split cluster B (# = %d)...\n', numel(typeB));
                        figB = figure;
                        openFigures = [openFigures figB];
                        
                        mosaicDRSpec(DRsylls(typeB), plotParams, 'maxMosaicLength', tPrev);
                        set(gcf,'Name',sprintf('Cluster B (# = %d)',numel(typeB)));                        
                    case 'Listen'
                        % simply play the sounds: it's faster to load but
                        % slower to be complete
                        fprintf('Splitting candidate... branch A audio: \n');
                        markRegions([],DRsylls(typeA));
                        fprintf('Splitting candidate... branch B audio: \n');
                        markRegions([],DRsylls(typeB));
                
                    case 'Accept w/o review'
                end
                
                choices = {'Both','Branch A', 'Branch B','None'};
                [respStr, isAccept] = nm_listdlg('Name', 'Which branches to accept (none is an option)?', 'ListString',choices, ...
                    'SelectionMode', 'Single','OKString','Accept',...
                    'CancelString','Retry');
                
                if ~isAccept, 
                    tryAgain = true;
                    closeAll(openFigures);
                    return;
                end;
                respStr = choices{respStr};
                switch respStr
                    case 'Both'
                        % will flatten these later
                        newLabelGrps = {typeA, typeB};
                    case 'Branch A'
                        newLabelGrps = {typeA};
                    case 'Branch B'
                        newLabelGrps = {typeB};                        
                    case 'None'
                        newLabelGrps= {};
                end
            case 'Sort by Ear'
                fprintf('Split into no more than 10..., based on audio: \n');
                newIdxs = multiMark([], DRsylls(labels == oldLabel));
                if any(isnan(newIdxs))
                    if strcmp('Start over',questdlg(['Hand labeling is stopped early, accept or start over?'],...
                                'Labeling stopped', ...
                                'Accept','Start over','Start over'))
                        tryAgain = true;
                        closeAll(openFigures);
                        return
                    end
                end
                [~,~,rIdxs] = unique(newIdxs);
                
                newLabelGrps = cell(1,numel(rIdxs));
                for kk = 1:numel(rIdxs)
                    newLabelGrps{kk} = find(rIdxs == kk);
                end
        end
        closeAll(openFigures);
    end
end

function closeAll(figH)
        for kk = 1:numel(figH)
            close(figH(kk));
        end
end
function [otherCluster, mergeDepth] = findMerge(grpPtr, origClust, clusterLabels)
isMerged = false;

origPtr = grpPtr;
clustTrav = origClust;
while ~isMerged && grpPtr > 1
    ct = crosstab(clusterLabels(:,grpPtr), clusterLabels(:, grpPtr-1));
    % the row that contains the current cluster - does it
    % contain another cluster?
    dstClusters = find(ct(clustTrav,:)>0);
    if sum(ct(:,dstClusters) > 0) > 1
        otherCluster = sum(find(ct(:,dstClusters))) - clustTrav;
        isMerged = true;
    else
        clustTrav = dstClusters;
        grpPtr = grpPtr - 1;
    end
end
mergeDepth = grpPtr;

% follow the new cluster back to the current division? (we
% don't have to do this, it might split again)
for jj = grpPtr:origPtr-1
    ct = crosstab(clusterLabels(:,jj), clusterLabels(:,jj+1));
    [~, otherCluster] = max(ct(otherCluster,:));
end
end

function [clusterInds, splitInds, splitDepth] = findSplit(grpPtr, origClust, clusterLabels)
hasSplit = false;
clustTrav = origClust;

nClusterings = size(clusterLabels, 2);
while ~hasSplit && grpPtr < nClusterings
    ct = crosstab(clusterLabels(:,grpPtr), clusterLabels(:, grpPtr+1));
    % the row that contains the current cluster - does it
    % contain another cluster?
    dstClusters = find(ct(clustTrav,:) > 0);
    
    if numel(dstClusters) == 2
        clustTrav = dstClusters(1);
        otherCluster = dstClusters(2);
        hasSplit = true;
    else
        clustTrav = dstClusters;
    end
    grpPtr = grpPtr + 1;
end
if hasSplit
    clusterInds = find(clusterLabels(:,grpPtr) == clustTrav);
    splitInds = find(clusterLabels(:,grpPtr) == otherCluster);
    splitDepth = grpPtr;
else
    clusterInds = find(clusterLabels(:,grpPtr) == clustTrav);
    splitInds = [];
    splitDepth = nClusterings; % excepted
end
end