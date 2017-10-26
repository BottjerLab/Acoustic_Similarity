function [counts, times, rates] = countSpikes(evs, spikeTimes, opt)

% function [counts,times] = countSpikes(evs, spikeTimes, opt)
% Counts the spikes within each event in the evs struct array.
% Returns counts, the number of spikes within each event, 
% times, a cell array of the spike times within each event, and
% rates, the firing rate in Hz within each event

% opt can take different string values and changes the times value to be
% relative to which part of the event

if nargin < 3, opt = 'absolute'; end

nEvs = numel(evs);
spikeTimes = sort(spikeTimes);
counts = zeros(1, nEvs);
times = cell(1, nEvs);
nSpikes = numel(spikeTimes);

starts = vertcat(evs.start);
stops  = vertcat(evs.stop );
flags  = [zeros(1,nSpikes) 1:nEvs 1:nEvs];
tPoses = [spikeTimes; starts; stops];

[sTposes, poseSortIdx] = sort(tPoses);
sTflags = flags(poseSortIdx);

isInterestPt = (sTflags ~= 0);
for ii = 1:nEvs
    iPoles = find(sTflags==ii);
    iRange = (iPoles(1)+1):(iPoles(2)-1);
    iRange(isInterestPt(iRange)) = [];
    times{ii} = sTposes(iRange);    

    switch lower(opt)
        case 'onset'
            times{ii} = times{ii} - evs(ii).start;
        case 'offset'
            times{ii} = times{ii} - evs(ii).stop;
        case 'absolute'
            % do nothing
        otherwise
            error('BadOption', 'opt must be ''onset'', ''offset'', or ''absolute'''); %#ok<CTPCT>
    end
    counts(ii) = numel(times{ii});
end

rates = counts ./ ([evs.stop]-[evs.start]);

