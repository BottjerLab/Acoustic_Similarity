function [averageCC, ccList, shuffledCC, shuffledList] = getCC(events, eventSpikeTimes)
% note: NaNs arise out of the list when there are no spikes in an event
binSize = 0.001; %s, 1 ms bins
smoothKernel = 0.01; %s, 10 ms smoothing kernel

% get relative spike times
%[~, eventSpikeTimes] = countSpikes(events, spikeTrain, 'onset');

% set up bins for each event
evLengths = [events.stop] - [events.start];
nEvs = numel(evLengths);

% bin spike trains and smooth
% note: this code truncates the spike trains according to the event stops,
% which may lead to distortion/bias
[ccList, averageCC] = innerCC(eventSpikeTimes, evLengths);

% shufflling?
if nargout >= 3
    shuffledEventTimes = cell(size(eventSpikeTimes));
    for ii = 1:numel(events);
        if isempty(eventSpikeTimes{ii}), continue; end
        % shuffle at least 100 ms to 5-- ms in either direction
        randOffset = randIn([min(0.1,ceil(evLengths(ii))) 0.5]) * sign(rand - 0.5);
        shuffledEventTimes{ii} = rem(eventSpikeTimes{ii} + randOffset, evLengths(ii));
    end
    [shuffledList, shuffledCC] = innerCC(shuffledEventTimes, evLengths);
end

% plot comparison histograms
% bins for the value of the cross correlations
ccBins = -0.1:0.02:1;
truCol = [0.9 0.2 0.2];
hCCList = hist(ccList,ccBins);
barh(ccBins, hCCList, 'FaceColor', truCol);
hold on;
plot(1, averageCC, '<', 'MarkerFaceColor', truCol);
hold off;
if nargout >= 3
    backCol = [0.5 0.5 0.5];
    hCCShuff = hist(shuffledList,ccBins);
    hold on; 
    barh(ccBins, -hCCShuff, 'FaceColor', backCol);
    plot(-1, shuffledCC, '>', 'MarkerFaceColor', backCol * 0.5, 'MarkerEdgeColor', backCol * 0.5);
    hold off;
end
ylim(ccBins([1 end]));
drawnow;

    function ans = randIn(interval)
        rand(1) * diff(interval) + interval(1);
    end

    function [ccList, meanCC] = innerCC(spikeTimes, lens)
        maxlen = ceil(max(lens)/binSize);
        
        % bin spike trains and smooth
        % note: this code truncates the spike trains according to the event stops,
        % which may lead to distortion/bias
        
        % note: if there is no 
        N = numel(lens);
        trains = zeros(maxlen, N);
        for jj = 1:N
            bb = 0:binSize:lens(jj);
            if isempty(spikeTimes{jj}), continue; end
            
            trains(1:numel(bb),jj) = histc(spikeTimes{jj}, bb);
            
            % smooth each spike train
            trains(:,jj) = smoothSignal(trains(:,jj), smoothKernel/binSize, 'sd');
            
            % subtract mean firing rate
            trains(:,jj) = trains(:,jj) - mean(trains(:,jj));
        end
        
        % calculate zero-lag cross correlation across all the spike trains
        ccList = corrcoef(trains);
        
         % trick to get the upper half only in column form 
        ccList(1:N+1:end) = 0;
        ccList = squareform(ccList);

        meanCC = nanmean(ccList);
    end
end