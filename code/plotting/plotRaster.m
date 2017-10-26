function xcoords = plotRaster(events, spikeTimes, alignmentOpt, customAlign)
% function plotRaster(events, spikeTimes) plot time-warped rasters/PSTH
% TODO: revise/refactor this - too specific
%
% This function plots the raster over a single event.  Spike times
% should be in a corresponding cell array (spikeTimes) 
% and absolute (as provided by countSpikes).  
%
% If events has a 'children' field, plots children as gray areas.  Children
% should NOT be timed relative to their parents.  
%
% alignmentOpt determines how the events are aligned.  
% It can be 'onset', 'offset', or 'custom', and is 'onset' by default.  

if nargin < 3
    alignmentOpt = 'onset';
end
assert(iscell(spikeTimes) && numel(events) == numel(spikeTimes));
if strcmp(alignmentOpt, 'custom')
    assert(nargin == 4 && numel(events) == numel(customAlign));
end

% these values are needed by both plotting subroutines
maxEvLength = max([events.stop] - [events.start]);

N = numel(events);
counts = cellfun(@(x) length(x), spikeTimes);

csum = [0 cumsum(counts)];
xcoords = zeros(1,csum(end));
ycoords = zeros(1,csum(end));

% set up warping to children
if strcmp(alignmentOpt, 'warping')
    % gather up each of the children as control points - assumes everyone
    % has the same children... (no checking for now)
    controlPoints = cell(1,N);
    zeroOffControlPoints = cell(1,N);
    for ii = 1:N
        cntrl = [0 events(ii).stop - events(ii).start];
        for jj = 1:numel(events(ii).children)
            cntrl = [cntrl (events(ii).children(jj).start - events(ii).start) ...
                (events(ii).children(jj).stop - events(ii).start)];
        end
        zeroOffControlPoints{ii} = sort(cntrl);
        controlPoints{ii} = zeroOffControlPoints{ii} + events(ii).start;
    end
    allPts = vertcat(zeroOffControlPoints{:});
    registerPts = mean(allPts, 1);
end

for ii = 1:N
    if counts(ii) == 0, continue; end;
    switch alignmentOpt
        case 'onset'
            xcoords((csum(ii)+1):csum(ii+1)) = spikeTimes{ii} - events(ii).start;
        case 'offset'
            xcoords((csum(ii)+1):csum(ii+1)) = spikeTimes{ii} - events(ii).stop + maxEvLength;
        case 'custom'
            xcoords((csum(ii)+1):csum(ii+1)) = spikeTimes{ii} - customAlign(ii);
        case 'warping'
            xcoords((csum(ii)+1):csum(ii+1)) = ...
                interp1(controlPoints{ii}, registerPts, spikeTimes{ii},...
                'linear', 'extrap');
    end
    ycoords((csum(ii)+1):csum(ii+1)) = ii; 
end

%%%%%%%%%%%%%%%%%%% plot raster
% this code tells us how to sort out the y-coordinates so that
% each 'spike' lands on the correct row
% actually plot the spikes as black squares
for ii = 1:numel(xcoords)
    plot([1 1] * xcoords(ii), [1 1] * ycoords(ii) + [-0.3 0.3], 'k', 'MarkerSize',3);
    hold on;
end
hold off;
% correct sizing
if strcmp(alignmentOpt, 'warping')
    xlim([0 max(registerPts)]);
elseif ~isempty(xcoords)
    maxTime = max(max(xcoords), maxEvLength);
    xlim([0 maxTime]);
else 
    maxTime = maxEvLength;
end
ylim([0.5 numel(events) + 0.5])

% optional: plot gray/colored subEvent backgrounds behind
%{
if isfield(events, 'children')    
    for ii = 1:numel(events)
        assert(isEvent(events(ii).children) || isempty(events(ii)));
        offset = 0;
        toDraw = events(ii).children;
        switch alignmentOpt
            case 'onset'
                offset = events(ii).start;
            case 'offset'
                offset = events(ii).stop - maxEvLength;
            case 'custom'
                offset = customAlign(ii);
            case 'warping'
                toDraw = eventFromTimes(registerPts(2:2:end-1), registerPts(3:2:end-1),10000); %fs is arbitrary
        end
        subEvs = adjustTimeStamps(toDraw, -offset);

        markColors = [0.5 0.5 0.5];
        plotAreaMarks(subEvs, markColors, false, [-0.5 0.5] + ii);
    end
end
%}