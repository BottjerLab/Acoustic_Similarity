spikesIn = spikes;
spikesIn{end+1} = vertcat(spikesIn{:});

for ii = 1:numel(spikesIn)
    allspikesIn = spikesIn{ii};
% get song firing rates
spikesInSong = countSpikes(juvenileSylls,allspikesIn);
songLength = [juvenileSylls.stop] - [juvenileSylls.start];
ratesInSong = spikesInSong ./ songLength;

% get non song firing rates
totalTime = (2125-800);
nspikesInNotInSong = numel(allspikesIn) - sum(spikesInSong);
nonSongTime = totalTime - sum(songLength); % seconds
baseRate = nspikesInNotInSong / nonSongTime;

[h,p,ci,stats]=ttest(ratesInSong - baseRate);
fprintf('Neuron #%d: baseline rate %0.2f Hz, song rate %0.2f Hz, p value for diff %0.3f\n', ...
    ii, baseRate, mean(ratesInSong), p); 
end
% find control periods
%{
nControl = 100;
controlEvents = initEvents(nControl);
for ii = 1:nControl
    randStart = rand * totalTime;
    if 
    controlEvents(ii).start = ;
end
%}