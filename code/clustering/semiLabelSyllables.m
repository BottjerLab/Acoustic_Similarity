function improvedLabels = semiLabelSyllables(syllSet, dists, currLabels)
% currLabels should have some clusters and some NaNs
if isvector(dists), dists = squareform(dists); end;

passOpts = {'manual','median', 'kNN'};
plotParams = processArgs(defaultParams,...
                            'dgram.minContrast', 1e-11, 'doFilterNoise', false,...
                            'preroll', 3, 'postroll', 3);
        

secondPass = nm_listdlg('ListString', passOpts, 'Name','Choose selection mode','SelectionMode', 'single');
secondPass = passOpts{secondPass};
improvedLabels = currLabels;

nLabels = nanmax(currLabels);
isUnlabeled = isnan(currLabels);

unlabbed = find(isUnlabeled);
nUn = numel(unlabbed);

if nUn == 0
    fprintf('No unlabeled syllables, sweet..');
    improvedLabels = currLabels;
    return;
end

elgDists = zeros(nUn,nLabels);
likenessScore = NaN(nUn, nLabels);
switch secondPass
    case 'median'
        % better would be do to one vs one for all pairs and vote that way, but
        % it'd likely be the same
        [diams, cMeans, cStds] = clusterDiameter(dists, currLabels, 1:nLabels);
        for ii = 1:nUn
            for jj = 1:nLabels
                elgDists(ii,jj) = median(dists(unlabbed(ii), currLabels == jj));
            end
            % inclusion would be if these distances are less than the mean? try
            zLikeScore = (cMeans - elgDists(ii,:) ) ./ cStds;
            likenessScore(ii,:) = normcdf(zLikeScore);
        end
        
        likenessScore = normcdf((ones(nUn,1) * cMeans - elgDists) ./ (ones(nUn,1) * cStds));
        % we need these three fields
        [bestLikeness, bestMatch] = max(likenessScore,[],2);
        confidenceScore =  bestLikeness ./ sum(likenessScore, 2);
    case 'kNN'
        [diams, cMeans, cStds] = clusterDiameter(dists, currLabels, 1:nLabels);
        
        % standardize distances to known quantities
        knownDists = dists(isUnlabeled, ~isUnlabeled);
        extantLabels = currLabels(~isUnlabeled);
        indivMeans = cMeans(extantLabels); if size(indivMeans, 2) == 1, indivMeans = indivMeans'; end
        indivStds = cStds(extantLabels); if size(indivStds, 2) == 1, indivStds = indivStds'; end
        knownDists = (knownDists - ones(nUn,1) * indivMeans) ./ (ones(nUn,1) * indivStds); % can be negative
        
        defaultVal = 5;
        kNeigh = nm_inputdlg(...
            'Value of k (roughly 1-50, number of neigbhors):' ,...
            'kNN parameters', 1, {num2str(defaultVal)});
        
        if isempty(kNeigh), tryAgain = true; return;
        else kNeigh = str2double(kNeigh); end
        
        bestMatch = zeros(1,nUn);
        bestLikeness = zeros(1,nUn);
        confidenceScore = zeros(1,nUn);
        % now take the first kk entries
        for ii = 1:nUn
            [lowScores, bestOrder] = sort(knownDists(ii,:));
            closestNeighbors = extantLabels(bestOrder(1:kNeigh));
            matchCounts = histc(closestNeighbors, 1:nLabels);
            [nMostMatch, bestMatch(ii)] = max(matchCounts);
            if sum(matchCounts == nMostMatch) > 1                                
                candidates = find(matchCounts == nMostMatch);
                % pick the one that appears the first? 
                firstAppear = zeros(1,numel(candidates));
                for jj = 1:numel(candidates)
                    firstAppear(jj) = find(closestNeighbors == candidates(jj),1);
                end
                [~,ord] = min(firstAppear); 
                bestMatch(ii) = candidates(ord);
                % or the one with
                % smaller sum of ranks?                
            end
            
            % get confidence via 1 - prod(ratios)
            matchPos = find(closestNeighbors == bestMatch(ii), kNeigh); 
            matchScore = lowScores(matchPos);
            
            nonMatchPos = find(closestNeighbors ~= bestMatch(ii));
            
            for jj = 1:nMostMatch
                npos = find(nonMatchPos > matchPos(jj), 1);
                if isempty(npos), 
                    nonMatchScore(jj) = lowScores(kNeigh+1);
                else
                    nonMatchScore(jj) = lowScores(npos);                    
                end
            end
            bestLikeness(ii)    = 1 / (1 + exp(matchScore(1)));
            confidenceScore(ii) = 1 / (1 + exp(sum(matchScore) - sum(nonMatchScore)));
        end
end

%% interactive to set thresholds
if ~strcmp(secondPass,'manual')
    cols = jet(nLabels);
    legStr = cell(1,nLabels);
    for ii = 1:nLabels
        if ~any(bestMatch == ii)
            continue;
        end
        legStr{ii} = sprintf('Cluster %d, # = %d', ii, sum(bestMatch==ii));
        subplot(211);
        h=cdfplot(confidenceScore(bestMatch==ii));
        set(h,'Color', cols(ii,:));
        xlabel('Confidence');
        hold on;
        if ii == nLabels, hold off; end;
        
        subplot(212);
        h=cdfplot(bestLikeness(bestMatch==ii));
        set(h,'Color', cols(ii,:));
        xlabel('Raw score');
        hold on;
        if ii == nLabels, hold off; end;
    end
    legStr(cellfun(@isempty,legStr)) = [];
    legend(legStr);
    thresholds = nm_inputdlg({'What threshold for confidence?','What threshold for raw score?','AND or OR to accept'}, ...
        'Threshold for automatic', 1,{'0.6','0.5','or'});
    
    confThresh = str2double(thresholds{1});
    rawThresh  = str2double(thresholds{2});
    threshOp = thresholds{3};
end
%% try looking through


if strcmp(secondPass,'manual')
    %{
    % crib from syllables
    [~,birdSession] = fileparts(syllSet(1).file);
    birdID = strtok(birdSession ,'_');
    thisAge = getAgeOfSession(birdSession);
    tmpDir = [pwd filesep 'tmp-' birdID '-age' num2str(thisAge) filesep];
    mostCentral = findMostCentral(dists, currLabels); % cell array
    for ii = 1:nLabels
        % create temporary figures
        nEx = 5;
        bestFive = syllSet(mostCentral{ii}(1:nEx));
        mosaicDRSpec(bestFive, plotParams, 'doFilterNoise', true);
        
        filName = [tmpDir 'clusterEx-' num2str(ii,'%02d') '.jpg'];
        saveCurrFigure(filName);
        clf
    end
    system(['explorer ' tmpDir]);
    %}
end
% loop through each syllable, but shunt  & make recommendations
if strcmp(secondPass,'manual')
    figPerm = figure('MenuBar', 'none', 'ToolBar','none');
end
for ii = 1:nUn
    % check the scores
    
    % assign w/o check
    if ~strcmp(secondPass, 'manual') 
        if (strcmpi('and',threshOp) && (confidenceScore(ii) >= confThresh && bestLikeness(ii) >= rawThresh)) || ...
                (strcmpi('or' ,threshOp) && (confidenceScore(ii) >= confThresh || bestLikeness(ii) >= rawThresh))
            improvedLabels(ii) = bestMatch(ii);
            continue;
        end
    end
    
    % plot the clip
    if strcmp(secondPass, 'manual')
        figure(figPerm);
        subplot('Position',[0 0 1 1]);
        [cl,fs] = getClipAndProcess([],syllSet(unlabbed(ii)), ...
            defaultParams, 'doFilterNoise', false);
        
        fineSpecParams = getfield(defaultParams,'fine');
        fineSpecParams.features = {'deriv'};
        fineSpecParams.fs = fs;
        spec = getMTSpectrumStats(cl, fineSpecParams);
        plotDerivGram(spec, defaultParams,'dgram.minContrast', 1e-9);
        set(gca,'YTickLabel',[],'YTick',[]);
        set(figPerm,'Name',sprintf('Candidate %d/%d for cluster #%d',jj,nUn,ii));
        str2double(nm_inputdlg(...
            sprintf('Which cluster should we assign this to? (1-%d, NaN for none)',nLabels),...
            sprintf('Approvals, %d/%d',ii,nUn), 1, {num2str(bestMatch(ii))}));
        
        clf(figPerm);
    end
end
if strcmp(secondPass,'manual')
    delete([tmpDir '*.jpg']); deleteStat = rmdir(tmpDir);
end
end