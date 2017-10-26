function [histOut, baseIFR] = eraFindPSTH( ERAWindow, events, type, binSize, spikes, baseline, params )
%get PSTH in relation to event onsets or offsets
%   Jenny wrote this 12/29/14

%ERAWindow is time in seconds from event of interest at time 0, e.g. [-0.2 0.2]
%Event is type of event, such as aS for approved syllables, etc.
%binSize is in seconds
%type is a string 'onset' or 'offset'

if strcmp(type,'onset')
    sEnd = 'start';
end
if strcmp(type,'offset')
    sEnd = 'stop';
end

ERAPeriods = eventFromTimes([events.(sEnd)]+ERAWindow(1), [events.(sEnd)]+ERAWindow(2), params.fs);
nUnits = numel(spikes);
bins = ERAWindow(1):binSize:ERAWindow(2);
histOut = zeros(length(spikes),length(bins(1:end-3))); %had to do end-3 to match
baseIFR = zeros(1,nUnits);

for kk = 1:nUnits
    
    [~,spikeTimes] = countSpikes(ERAPeriods, spikes{kk},'onset');
    spikeTimes = vertcat(spikeTimes{:}) + ERAWindow(1);
    
    % get total histogram for each neuron
    tmpHist = histc(spikeTimes,bins); %/numel(events)
    %tmpHist = tempHist/binSize;
    
    [~,spikeTimes] = countSpikes(baseline, spikes{kk},'onset'); %get baseline FR
    spikeTimes = vertcat(spikeTimes{:});
    tmpHist2 = histc(spikeTimes,bins)/ numel(baseline);
    tmpHist2 = tmpHist2/binSize;
    baseIFR(kk) = mean(tmpHist2);
    
    histOut(kk,:) = tmpHist(1:end-3);
    
end

end

