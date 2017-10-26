function eventSpikeTimes = plotAlignedPSTH(events, spikeTimes, subEvents, labelsOfInterest)
% function plotAlignedPSTH(events, spikeTimes) plot time-warped rasters/PSTH
%
%
% spikeTimes is expected to be a raw array of spike times

if nargin < 3
    subEvents = initEvents(1);
end
if nargin < 4
    labelsOfInterest = {};
elseif ~iscell(labelsOfInterest)
    labelsOfInterest = {labelsOfInterest};
end

[counts, eventSpikeTimes] = countSpikes(events, spikeTimes,'onset');

syllableLabels = {subEvents.type};
syllableLengths = [subEvents.stop] - [subEvents.start];
[uLabels, foo, rLabelIdx] = unique(syllableLabels);

[parentage, eventsInContext] = findParent(events,subEvents);

nAlign = numel(labelsOfInterest);
rAlignIdx = NaN(size(rLabelIdx));
doWarping = (nargin == 4);
if doWarping
    % set index of event types according to the align list    
    for ii = 1:nAlign
        labelsIndex = find(strcmp(uLabels,labelsOfInterest{ii}));
        if ~isempty(labelsIndex)
            rAlignIdx(rLabelIdx == labelsIndex) = ii;
        end
    end
    
    % get the average onset and length of each syllable
    avgLengths = zeros(1,nAlign);
    for ii = 1:nAlign
        % get average length
        isThisLabel = (rAlignIdx == ii);
        avgLengths(ii) = mean(syllableLengths(isThisLabel));
        avgStarts(ii) = mean([eventsInContext(isThisLabel).start]);
    end
    
    warpedEventsInContext = eventsInContext;
    for ii = 1:numel(events)
        isChild = (parentage == ii);
        childEvs = eventsInContext(isChild);
        controlPoints = NaN(2*nAlign, 2); % first column is sample, second is standard;
        controlSet = false(1,2*nAlign);
        for jj = 1:nAlign
            % in case of multiply labeled syllables, just pick the first
            % one for alignment purposes
            matchAlign = find(strcmp({childEvs.type},labelsOfInterest{jj}),1);
            if ~isempty(matchAlign)
                controlPoints(2*jj-1:2*jj,:) = ...
                    [childEvs(matchAlign).start avgStarts(jj); ...
                    childEvs(matchAlign).stop  avgStarts(jj) + avgLengths(jj)];
                controlSet(2*jj-1) = true;
                controlSet(2*jj) = true;
            end
        end
        controlPoints(~controlSet,:) = []; % get rid of nan points
        
        % now use control points to linearly interpolate times
        adjEventSpikeTimes{ii} = interpLinearPoints(eventSpikeTimes{ii},...
            controlPoints(:,2), controlPoints(:,1));
        
        % for plotting purposes: warp event boundaries themselves
        warpedStarts = interpLinearPoints([childEvs.start],...
            controlPoints(:,2), controlPoints(:,1));
        warpedStops  = interpLinearPoints([childEvs.stop],...
            controlPoints(:,2), controlPoints(:,1));
        
        foo = num2cell(warpedStarts); 
        [warpedEventsInContext(isChild).start] = foo{:};
        foo = num2cell(warpedStops); 
        [warpedEventsInContext(isChild).stop]  = foo{:};        
    end
    
end

% prepare the bins - all important binning parameter
binSize = 1e-3; % in seconds

% concatenate the spiketimes
if doWarping
    eventSpikeTimes = adjEventSpikeTimes;    
end
xcoords = vertcat(eventSpikeTimes{:});
maxTime = max(xcoords);

% plot raster
subplot(3,1,1:2)
plotRaster;
title(sprintf('Raster for song'));

% plot average PSTH
subplot(3,1,3)
title(sprintf('PSTH'));
plotPeriHist

    function plotRaster
        % option
        plotColors = true;
        
        % this code tells us how to sort out the y-coordinates so that
        % each 'spike' lands on the correct row
        stepsUp = zeros(size(xcoords));
        stepsUp(cumsum(counts(1:end-1))+1) = 1; stepsUp(1) = 1;
        ycoords = cumsum(stepsUp);
        
        % actually plot the spikes
        plot(xcoords, ycoords, 'ks','MarkerSize',2,'MarkerFaceColor','k');
        xlim([0 maxTime]);
        ylim([0.5 numel(events) + 0.5])
        
        % optional: plot marks behind
        for kk = 1:numel(events)
            isKthChild = (parentage == kk);
            subEvs = eventsInContext(isKthChild);
            if doWarping
                subEvs = warpedEventsInContext(isKthChild);
            end
            
            % assign the new colors only to the labels that correspond to alignment ones
            if plotColors
                typeColors = autumn(nAlign); % color assignment for each label
                evTypes = rAlignIdx(isKthChild);
                markColors = 0.5*ones(numel(subEvs),3);
                markColors(~isnan(evTypes),:) = typeColors(evTypes(~isnan(evTypes)),:);
            else
                markColors = [0.5 0.5 0.5];
            end
            plotAreaMarks(subEvs, markColors, false, [-0.5 0.5] + kk);
        end
    end

    function plotPeriHist
        histBar = histc(xcoords, 0:binSize:maxTime);
        bar(0:binSize:maxTime, histBar, 'histc');
        xlim([0 maxTime]);
    end
end

function yy = interpLinearPoints(xx,y1,x1)
% interpolate linearly between each pair of points & interpolate 
assert(numel(x1) == numel(y1))
if numel(x1) == 0, yy = xx; return; end;

[x1, idxs] = sort(x1);
y1 = y1(idxs);
intervalPtr = 1;

yy = zeros(size(xx));
for kk = 1:numel(xx)
    intervalPtr = find(xx(kk) <= x1,1);
    if isempty(intervalPtr)
        yy(kk) = xx(kk) - x1(end) + y1(end); % no warping after last control point
    elseif intervalPtr == 1
        yy(kk) = xx(kk) - x1(1) + y1(1); % no warping before first control point
    else
        yy(kk) = interp1(x1(intervalPtr-1:intervalPtr),...
            y1(intervalPtr-1:intervalPtr),xx(kk));
    end
end
end
