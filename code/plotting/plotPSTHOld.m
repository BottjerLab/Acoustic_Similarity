function eventSpikeTimes = plotPSTHOld(events, spikeTimes, subEvents, labelsOfInterest, prePostParams)
% function plotPSTH(events, spikeTimes) plot time-warped rasters/PSTH
% TODO: revise/refactor this - too specific
% 
% This function plots the raster and PSTH for each event.  Spike times
% should be input uncorrected, that is, raw times with respect to the same 
% timeframe as that of the events.
%
% plotPSTH(events, spikeTimes, subEvents) plots the subEvents as gray
% areas.  Tip: Use events as an expanded version of subEvents to plot a PSTH
% with baseline.
%
% plotPSTH(events, spikeTimes, subEvents, labelsOfInterest) plots a separate
% raster/PSTH with time warping for each event of type label of interest.  Times are
% warped so that events of each label are warped.
%
% labelsOfInterest is a cell array of strings for the syllables that should
% be aligned.  everything will be time warped to those syllables. (These
% should be inorder)
%
% spikeTimes is expected to be a column vector of spike times.
%
% if warping is requested (by adding labels)
if nargin < 3
    subEvents = initEvents;
end
if nargin < 4
    labelsOfInterest = {};
elseif ~iscell(labelsOfInterest)
    labelsOfInterest = {labelsOfInterest};
end

% find which subevents belong to which events
% note: if the events and subevents are in one to one correspondence, (i.e.
% one is the pre/post rolled version of another, then the parentage is
% probably just one-to-one
if numel(events) == numel(subEvents)
    parentage = 1:numel(subEvents);
    for ii = 1:numel(subEvents)
        eventsInContext(ii) = adjustTimeStamps(subEvents(ii), ...
            -events(ii).start);
    end
else
    if nargin < 5
        [parentage, eventsInContext] = findParent(events,subEvents);
    else %nargin == 5
        parentage = findParent(events,subEvents);
        events = addPrePost(events, prePostParams);
        % take out orphan events :[
        subEvents(isnan(parentage)) = [];
        parentage(isnan(parentage)) = [];
        for ii = 1:numel(subEvents)
            eventsInContext(ii) = adjustTimeStamps(subEvents(ii), ...
                    -events(parentage(ii)).start);
        end
    end

end

% take out orphan events :[
subEvents(isnan(parentage)) = [];
eventsInContext(isnan(parentage)) = [];
parentage(isnan(parentage)) = [];

% equalize outside events so that spikes are evenly counted
maxEvLength = max([events.stop] - [events.start]);

% now extend/alter the events so that they all have the same postroll /
% preroll
% if ~isempty(events)
%     num2cell(maxEvLength + [events.start]); [events.stop] = ans{:};
%     [events.length] = deal(maxEvLength);
% end
[counts, eventSpikeTimes] = countSpikes(events, spikeTimes,'onset');

syllableLabels = {subEvents.type};
syllableLengths = [subEvents.stop] - [subEvents.start];
[uLabels, foo, rLabelIdx] = unique(syllableLabels);

nAlign = numel(labelsOfInterest);
rAlignIdx = NaN(size(rLabelIdx));

% should we warp? (according to whether or not user gave syllables list)
doWarping = (nargin >= 4);
if doWarping
    % set index of event types according to the align list    
    for ii = 1:nAlign
        labelsIndex = find(strcmp(uLabels,labelsOfInterest{ii}));
        if ~isempty(labelsIndex)
            rAlignIdx(rLabelIdx == labelsIndex) = ii;
        end
    end
    
    % get the average onset and length of each syllable
    avgEnds = zeros(1,nAlign);
    avgStarts = zeros(1,nAlign);
    for ii = 1:nAlign
        % get average length
        isThisLabel = (rAlignIdx == ii);
        avgEnds(ii) = mean([eventsInContext(isThisLabel).stop]);
        avgStarts(ii) = mean([eventsInContext(isThisLabel).start]);
    end
    
    % repair cases where average events overlap for some strange reason
    % (needs to be more reasonable)
    if any(avgEnds(1:end-1) > avgStarts(2:end))
        warning('distorted average event profile, repairing...');
        defaultGap = 0.005;
        overlaps = find(avgEnds(1:end-1) > avgStarts(2:end));
        for ll = 1:numel(overlaps)
            gapJump = avgEnds(ll) - avgStarts(ll+1) + defaultGap;
            avgEnds((ll+1):end) = avgEnds((ll+1):end) + gapJump;
            avgStarts((ll+1):end) = avgStarts((ll+1):end) + gapJump;
        end
    end
    avgLengths = avgEnds - avgStarts;

    warpedEventsInContext = eventsInContext;
    for ii = 1:numel(events)
        isChild = (parentage == ii);
        childEvs = eventsInContext(isChild);
        controlPoints = NaN(2*nAlign, 2); % first column is sample, second is standard;
        controlSet = false(1,2*nAlign);
        for jj = 1:nAlign
            % NB: in case of multiply labeled syllables, just pick the first
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

% concatenate the context spiketimes
if doWarping
    eventSpikeTimes = adjEventSpikeTimes;    
end

% these values are needed by both plotting subroutines
xcoords = vertcat(eventSpikeTimes{:});
maxTime = max(max(xcoords), maxEvLength);

% plot raster
subplot(3,1,1:2)
plotRasterInner;
title(sprintf('Raster for song'));

% plot average PSTH
subplot(3,1,3)
title(sprintf('PSTH'));
plotPeriHist
% for title positioning
subplot(3,1,1:2)

%%%%%%%%%%%%%%%%%%%%%%%% plotting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotRasterInner
        plotColors = true;
        
        % this code tells us how to sort out the y-coordinates so that
        % each 'spike' lands on the correct row
        ycoords = [];
        for ii = 1:numel(counts)
            ycoords = [ycoords; ii * ones(counts(ii), 1)]; 
        end
        
        % actually plot the spikes
        plot(xcoords, ycoords, 'ks','MarkerSize',2,'MarkerFaceColor','k');
        xlim([0 maxTime]);
        ylim([0.5 numel(events) + 0.5])
        
        % optional: plot gray/colored subEvent backgrounds behind
        for kk = 1:numel(events)
            %isKthChild = (parentage == kk);
            
            % subEvs = eventsInContext(isKthChild);
            
            % plot all gray/colored 
            isKthChild = isWithinEvent(subEvents, events(kk));
            subEvs = adjustTimeStamps(subEvents(isKthChild), -events(kk).start);
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
        % counts per bin
        spikesPerBin = histc(xcoords, 0:binSize:maxTime);
        spikeRate = spikesPerBin / (binSize * numel(events));
        
        % smooth bins with gaussian curve:
        tt = 0:binSize:maxTime;
        smoothedSpikeCount = zeros(size(tt));
        smoothWidth = 4e-3;
        for ii = 1:numel(tt)
            smoothedSpikeCount(ii) = 1/sqrt(2*pi) * sum(exp(-((xcoords-tt(ii))/smoothWidth).^2));
        end
        
        smoothedRate = smoothedSpikeCount / (binSize * numel(events));
        
        % plot the raw PSTH
        %bar(0:binSize:maxTime, spikeRate, 'histc');
        %hold on
        
        % plot the smoothed stuff
        plot(tt,smoothedRate, 'r-');
        xlim([0 maxTime]);
        xlabel('Time(s)');
        ylabel('Mean firing rate (Hz)');
        %hold off;
    end
end

function yy = interpLinearPoints(xx,y1,x1)
% interpolate linearly between each pair of points & interpolate, without
% warping before and after the control points
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
