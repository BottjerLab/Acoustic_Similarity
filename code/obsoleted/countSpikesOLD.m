function [counts, times, rates] = countSpikes(evs, spikeTimes, opt)

% function [counts,times] = countSpikes(evs, spikeTimes, opt)
% Counts the spikes within each event in the evs struct array.
% Returns counts, the number of spikes within each event, 
% times, a cell array of the spike times within each event, and
% rates, the firing rate in Hz within each event

% opt can take different string values and changes the times value to be
% relative to which part of the event

if nargin < 3, opt = 'absolute'; end

counts = zeros(1, numel(evs));
times  =  cell(1, numel(evs));
rates = zeros(1, numel(evs));
if isempty(spikeTimes), return; end;
spikeTimes = sort(spikeTimes);

for ii = 1:numel(evs)
    % turns out this is already pretty darn fast - must be some internal
    % JIT optimization
    times{ii} = [spikeTimes(spikeTimes >= evs(ii).start & ...
        spikeTimes <= evs(ii).stop)];
    switch lower(opt)
        case 'onset'
            times{ii} = times{ii} - evs(ii).start;
        case 'offset'
            times{ii} = times{ii} - evs(ii).stop;
        case 'absolute'
            % do nothing
        otherwise
            error('BadOption', 'opt must be ''onset'', ''offset'', or ''absolute''');
    end
    counts(ii) = numel(times{ii});
end

rates = counts ./ ([evs.stop]-[evs.start]);
end