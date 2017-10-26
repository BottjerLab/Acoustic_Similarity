function [basalRate, meanSpikeRate, rateSEM, tp] = findSpikeRates(events, spikeTimes)
%function findSpikeCorrelates(events, statistics, spikeTimes)
% events is a struct array of events
% statistics is a struct array of scalars, one for each event
% spikeTimes is a vector
% test if spiking is significantly higher/lower during the event than
% during other times

maxTime = max(spikeTimes); % an approximation
spikeCounts = countSpikes(events, spikeTimes);
evLengths = [events.stop] - [events.start];
spikeRates = spikeCounts ./ evLengths;
basalRate = (numel(spikeTimes) - sum(spikeCounts)) / (maxTime - sum(evLengths));
rateSEM = std(spikeRates)/sqrt(numel(events)-1);
[h,tp,ci,stats]=ttest(spikeRates - basalRate);

fprintf('Basal rate %0.2f Hz, event rate %0.2f +/- %0.2f Hz, (p = %0.3f)\n', ...
    basalRate, mean(spikeRates), rateSEM, tp);

meanSpikeRate = mean(spikeRates);