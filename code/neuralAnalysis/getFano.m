function [fano, timeBins] = getFano(events, eventSpikeTimes)
% count variance / mean (= fano) in each window, stepping according to
% window size
% eventSpikeTimes is extracted as if from countSpikes, but with onset
stepSize = 0.002; %s, 2 ms bins
windowSize = 0.03; %s, 30 ms window

% get relative spike times
%[~, eventSpikeTimes] = countSpikes(events, spikeTrain, 'onset');

% set up bins for each event
evLengths = [events.stop] - [events.start];
minL = min(evLengths);

% bin spike trains and smooth
% note: this code truncates the spike trains according to the event stops,
% which may lead to distortion/bias
nEvents = numel(events);

sth = zeros(ceil(max(evLengths)/stepSize), nEvents);
boxSize = floor(windowSize / stepSize);
for ii = 1:nEvents
    bins = 0:stepSize:evLengths(ii);
    foo = histc(eventSpikeTimes{ii}, bins);

    % 'predictive' boxcar smooth each spike train
    foo = conv(foo,ones(boxSize,1));
    
    foo(1:boxSize-1) = [];
    sth(1:numel(foo),ii) = foo;    
end

% at the end, sth contains all the spikes in the window from
% currTime->currTime + windowSize

timeBins = 0:stepSize:max(evLengths);
fano = var(sth,0,2)./mean(sth, 2);
fano(isnan(fano)) = 0;

end