function drawClustersUpdated(birdID, age, params, varargin)

if nargin < 3 || isempty(params)
    params = defaultParams;
end
params = processArgs(params, varargin{:});

dataDir = [pwd filesep 'data' filesep birdID filesep];
% clusterDir = ['data' filesep 'cluster-' birdID filesep];

% get sessions for a certain age
rep = reportOnData(birdID);
birdSessions = {rep.sessionID};

sessAges = getAgeOfSession(birdSessions);
rep = rep(age == sessAges);

% get syllables for a certain age by concatenating across different groups
ageSylls = [];
for ii = 1:numel(rep)
    mani = rep(ii).manifest;
    [~,fil] = findInManifest(mani, 'sAudio');
   
    theseSylls = loadFromManifest(mani, 'approvedSyllables');

    % minimum length used in DRcluster would be smaller than the
    % spectrogram window
    % remove syllables that are too short
    minLen = params.fine.windowSize / 1000;
    isTooShort = minLen > [theseSylls.stop] - [theseSylls.start];
    theseSylls(isTooShort) = [];
    
    [theseSylls.file] = deal(fil);
    
    % get the labels
    labels = loadFromManifest(mani, 'acceptedLabels');
    if isempty(labels)
        fprintf('No acceptedLabels for session %s, not retrieving...\n', rep(ii).sessionID);
        continue;
    end
    
    % apply labels
    num2cell(labels); [theseSylls.type] = ans{:};    
    ageSylls = [ageSylls theseSylls];
end

% plot all syllables within each cluster as a mosaic
labels = [ageSylls.type];
nClusts = nanmax([ageSylls.type]);
clustLabels = unique(labels(~isnan(labels)));
if any(isnan(labels)), clustLabels(end+1) = NaN; end
    
mkdir('figures/',sprintf('%s-age%d', birdID, age));
redoneLabels = labels;
for ii = 1:nClusts
    clustSylls = ageSylls(labels == ii);
    if isempty(clustSylls), continue; end;
    
    clustSylls = addPrePost(clustSylls, [], 'preroll', 5, 'postroll',5);
    
    maxMosaicLength = 6.5;
    [hfs, newLabels] = plotInPages(clustSylls, maxMosaicLength, ii, ...
        clustLabels, birdID, age, params);
    close(hfs)
    
    % updating syllable labels
    if params.manuallyUpdateClusters              
        redoneLabels(labels == ii) = newLabels;
        redoneFile = sprintf('%sredoneLabels-%s-age%d.mat', dataDir, birdID, age, params);
        fprintf('Intermediate saving to %s...\n', redoneFile);
        save(redoneFile, 'redoneLabels');    
    end
end
unIDedSylls = ageSylls(isnan([ageSylls.type]));
unIDedSylls = addPrePost(unIDedSylls, [], 'preroll', 5, 'postroll',5);

maxMosaicLength = 6.5;
[hfs, newLabels] = plotInPages(unIDedSylls, maxMosaicLength, NaN, clustLabels, birdID, age, params);
% updating syllable labels
if params.manuallyUpdateClusters
    redoneLabels(isnan(labels)) = newLabels;
    redoneFile = sprintf('%sredoneLabels-%s-age%d.mat', dataDir, birdID, age);
    fprintf('Final saving to %s...\n', redoneFile);
    save(redoneFile, 'redoneLabels');
end

close(hfs)
end

function [hfs, newLabels] = plotInPages(sylls, pageLength, syllNum, clustLabels, birdID, age, params)
% pageLength in seconds, max amount to plot, accounting for wasted margins on the mosaic
% break down clustSylls into groups no longer than maxTime
lens = [sylls.stop] - [sylls.start];
grpStarts = 1;
grpEnds = numel(sylls);
while sum(lens(grpStarts(end):end)) > pageLength
    newGrpStart = find(cumsum(lens) > pageLength, 1);
    grpStarts(end+1) = newGrpStart;
    grpEnds   = [grpEnds(1:end-1) (newGrpStart-1) grpEnds(end)];
    lens(1:newGrpStart-1) = 0; % zero out the lengths that are already accounted for
end

isClosed = false(1,numel(grpStarts));
hfs = zeros(1,numel(grpStarts));
newLabels = ones(1,numel(sylls)) * syllNum;
for jj = 1:numel(hfs)
    figure;
    % get the figure handle and image handles
    [hfs(jj) hIms] = mosaicDRSpec(sylls(grpStarts(jj):grpEnds(jj)), [],...
        'dgram.minContrast', 1e-11, 'doFilterNoise', false,...
        'noroll', 'maxMosaicLength', Inf);
    set(hfs(jj), 'Name', sprintf('%s, syllable #%d, instances #%d-%d', ...
        birdID, syllNum, grpStarts(jj), grpEnds(jj)));
    
    % interactive portion %
    if params.manuallyUpdateClusters
        % get the axes objects
        hRowAxes = get(hIms,'Parent');
        % de-cell
        hRowAxes = [hRowAxes{:}];
        for kk = 1:numel(hRowAxes)
            hThisRowAxes = hRowAxes(kk);
            hThisIm = hIms(kk);
            sepH = findall(hThisRowAxes, 'Type','line','Color',[1 1 1]);
            nSeps = numel(sepH);
            
            % boundaries of the regions in figure x-coordinates
            xSeps = zeros(1,nSeps + 2);
            % ends of the separators
            xl = xlim(hThisRowAxes);
            xSeps(1:2) = xl;
            for ll = 1:nSeps
                xl = get(sepH(ll), 'XData');
                xSeps(ll+2) = xl(1); % first x-coordinate of the line is x position
            end
            xSeps = sort(xSeps);
            
            % make red text labels
            axes(hThisRowAxes);
            for ll = 1:nSeps+1 %each segment
                hlabel = createTextLabel(mean(xSeps(ll:ll+1)), ...
                    num2str(syllNum));
                set(hlabel, 'UserData', ll);
            end
            % save bounds and current labels in figure data structure
            set(hThisIm,'UserData', ...
                struct('bounds', xSeps, ...
                'currLabels', syllNum * ones(1,nSeps + 1)));
            
            hcmenu = uicontextmenu;
            for ll = 1:numel(clustLabels)
                uimenu(hcmenu,'Label',num2str(clustLabels(ll)),...
                    'Callback', @(hEntry, data) changeLabels(...
                    clustLabels(ll), hThisIm)); % both input variables are dummy variables
            end
            uimenu(hcmenu,'Label','Done',...
                'Callback', @(hEntry, data) ...
                stopLabeling(gcf)); %the data input dummy variables
            set(hThisIm,'uicontextmenu', hcmenu);
        end
        
        % set the close request function and a cleanup in case of
        % ctrl-c - this needs to be fixed so that c cleans up on ctrl-c
        % instead of on request to delete hfs(jj)
        set(hfs(jj),'CloseRequestFcn',@(hfig,event) set(hfig, 'Tag','Done'));
        c = onCleanup(@() resetCloseFcn(hfs(jj)));
        waitfor(hfs(jj),'Tag','Done')
        
        % change the labels by grabbing userdata from the images
        ptr = grpStarts(jj);
        for kk = 1:numel(hIms)
            foo = getfield(get(hIms(kk),'UserData'),'currLabels');
            newLabels(ptr:ptr+numel(foo)-1) = foo;
            ptr = ptr + numel(foo);
        end
        % todo: delete the labels
        assert(ptr == grpEnds(jj)+1);
        
        nChanged = sum(newLabels~=syllNum);
        if isnan(syllNum)
            nChanged = sum(isnan(newLabels));
        end
        fprintf('\tcluster #%d: %d labels changed...\n', syllNum, nChanged);
    end
    % end interactive portion %
    
    if params.saveplot
        fName = sprintf('figures/%s-age%d/%s-age%d-clust%02d-i%04d-%04d.jpg',...
            birdID,age, birdID,age, syllNum,grpStarts(jj), grpEnds(jj));
        saveCurrFigure(fName);
    end
    isClosed(jj) = true;
    close(hfs(jj));
end
hfs(isClosed) = [];
end

function resetCloseFcn(hfig)
    fprintf('Calling reset fxn\n');
    if isvalid(hfig), 
        set(hfig,'CloseRequestFcn','closereq'); 
    end
end

function changeLabels(seldLabel, hIm)
    % seldLabel is the name of the label
    % hIm is the image handle that contains the track of labels in userData
    
    % get the current point of the click in normalized coordinates
    hCurrAxes = get(hIm,'Parent');
    hCurrFig = get(hCurrAxes, 'Parent');
    currPt = get(hCurrFig, 'CurrentPoint');
    figDims = get(hCurrFig, 'Position');
    currPt = currPt ./ figDims([3 4]);
    
    % bounds of array
    userData = get(hIm,'UserData');
    xb = userData.bounds;

    % get the segment which was selected
    xpos = currPt(1);
    segNum = find(xb(1:end-1) < xpos & xb(2:end) >= xpos, 1);
    
    % reassign the label from the segment which was clicked    
    userData.currLabels(segNum) = seldLabel;
    set(hIm,'UserData', userData);    
    
    % rewrite the label
    hCurrLabel = findall(hCurrAxes, 'Type', 'text', 'UserData', segNum);
    set(hCurrLabel, 'String', num2str(seldLabel));
end

function stopLabeling(gCurrFig)
    % set the done tag
    set(gCurrFig,'Tag','Done');
end
