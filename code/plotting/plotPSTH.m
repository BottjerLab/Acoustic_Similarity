function plotPSTH(events, spikeTimes, alignmentOpt, customAlign, params)
% function plotPSTH(events, spikeTimes) plot PSTHs
%
% This function plots the PSTH within all events.  Spike times
% should be in a corresponding cell array (spikeTimes) 
% and absolute (as provided by countSpikes).  
%

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
aggrTimes = zeros(1,csum(end));

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
            aggrTimes((csum(ii)+1):csum(ii+1)) = spikeTimes{ii} - events(ii).start;
        case 'offset'
            aggrTimes((csum(ii)+1):csum(ii+1)) = spikeTimes{ii} - events(ii).stop + maxEvLength;
        case 'custom'
            aggrTimes((csum(ii)+1):csum(ii+1)) = spikeTimes{ii} - customAlign(ii);
        case 'warping'
            aggrTimes((csum(ii)+1):csum(ii+1)) = ...
                interp1(controlPoints{ii}, registerPts, spikeTimes{ii},...
                'linear', 'extrap');
    end
end


% smooth bins with gaussian curve:
binSize = 0.01; %seconds
sigma = 0.010;

% calculate values
if strcmp(alignmentOpt, 'warping')
    maxTime = max(registerPts);
elseif ~isempty(aggrTimes)
    maxTime = max(max(aggrTimes), maxEvLength);
else 
    maxTime = maxEvLength;
end
tt = 0:binSize:maxTime;
spikesPerBin = histc(aggrTimes, tt);
%spikeRate = spikesPerBin / (binSize * numel(events));
%smoothedRate = spikeRate;
%smoothedRate = smoothSignal(spikeRate, sigma/binSize, 'SD');

% plot the raw PSTH
%bar(0:binSize:maxTime, spikeRate, 'histc');
%hold on

% plot the smoothed and raw rate
%bar(tt,smoothedRate, 1, 'k'); % be a bar
bar(tt, spikesPerBin, 1, 'k');
xlim([0 maxTime]);
xlabel('Time(s)');
ylabel('Mean firing rate (Hz)');
%hold off;
