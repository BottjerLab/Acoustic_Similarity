function [histOut, baseIFR] = eraFind( ERAWindow, events, type, binSize, spikes, baseline, params )
%get firing rates in relation to event onsets or offsets


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
    
    % get total histogram for each neuron, and smooth
    tmpHist = histc(spikeTimes,bins) / numel(events);
    tmpHist = tmpHist/binSize;
    
    [~,spikeTimes] = countSpikes(baseline, spikes{kk},'onset'); %get baseline FR
    spikeTimes = vertcat(spikeTimes{:});
    tmpHist2 = histc(spikeTimes,bins)/ numel(baseline);
    tmpHist2 = tmpHist2/binSize;
    baseIFR(kk) = mean(tmpHist2);
    
    %if we want to subtract baseline
    %tmpHist = tmpHist - baseIFR(kk);
    
    %if we want to mean subtract to center on zero like Goldberg & Fee
    tmpHist = tmpHist(1:end-3); %to remove weird stuff at the end
    tmpHist = tmpHist - mean(tmpHist); 
    
    %maybe we want to normalize the response strength?
    
    histOut(kk,:) = smoothSignal(tmpHist, 15, 'SD'); % smoothing window size in bins, 15 binSize-width smoothing
%     if size(histOut(kk,:), 2) > 1, % force column vector output
%         histOut(kk,:) = histOut(kk,:)';
%     end

end

end

