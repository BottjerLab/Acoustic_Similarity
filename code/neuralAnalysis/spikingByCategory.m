function neuronData = spikingByCategory(songStruct, events, spikes)
% spikes is a cell array of different neurons
% get non song firing rates

nNeurons = numel(spikes);

syllLengths = [events.stop] - [events.start];
for ii = 1:nNeurons
    spikeCounts{ii} = countSpikes(events, spikes{ii},'onset');
    spikeRates{ii} = spikeCounts{ii} ./ syllLengths; 
end
totalSpikeRates = sum(cellfun(@sum,spikeRates));

totalTime = songStruct.interval * songStruct.length;
fs = 1/songStruct.interval;

% construct local baseline %%%%
% get local baseline rate [-3, -2] seconds before, that do not
% overlap with these events
% todo: make this customizable
localBaseStarts = [max(0, [events.start] - 4)];%; min(songStruct.length/fs, [ROIs.stop] + 2)];
localBaseStops =  [max(0, [events.start] - 2)];%; min(songStruct.length/fs, [ROIs.stop] + 3)];
localBaseEvents = eventFromTimes(localBaseStarts, localBaseStops, fs);
localBaseNotSong = findUniqueAreas(localBaseEvents, events);

fullDur = sum([localBaseEvents.stop] - [localBaseEvents.start]);
notSongDur = sum([localBaseNotSong.stop] - [localBaseNotSong.start]);
%isSongDur = sum([localBaseInSong.stop] - [localBaseInSong.start]);

% temporary: add premotor events
premotorEvents = eventFromTimes([events.start]-0.1, [events.start],fs);
[premotorEvents.type] = deal('premotor');
events = [events premotorEvents];

% get event types
evTypes = {events.type};
[uTypes, ~, rIdxType] = unique(evTypes);


fprintf('Baseline samples %0.2fs, %0.2f%% unoccluded by other events...\n', ...
    fullDur, notSongDur/fullDur * 100);
if notSongDur == 0
    warning('No legit local baseline.');
    keyboard
end
%fprintf('Neuron | Global | Local | (p)   | Adult | (p)   | Juvenile | (p)  \n');

spikesInEvent = cell(nNeurons,numel(uTypes));
ratesInEvent = cell(nNeurons,numel(uTypes));
spikesInValidRate = cell(nNeurons,1);
ratesInValidRate = cell(nNeurons,1);
mR = zeros(nNeurons, numel(uTypes));

for ii = 1:nNeurons   
    % get local baseline rate
    [spikesInValidBase{ii},~, ratesInValidBase{ii}] = countSpikes(localBaseNotSong,spikes{ii});
    localValidBaseRate(ii) = sum(spikesInValidBase{ii}) / notSongDur; 

    % get event firing rates
    for jj = 1:numel(uTypes)
        theseEvents = events(rIdxType == jj);
        [spikesInEvent{ii,jj}, ~, ratesInEvent{ii,jj}] = countSpikes(theseEvents,spikes{ii});
    end

   % get global (non-event) baseline rate
    baseRate(ii) = (numel(spikes{ii}) - sum([spikesInEvent{ii,:}])) / ...
        (totalTime - sum([events.stop]-[events.start]));
    
    [hbase{ii},pbase{ii},cibase{ii},statsbase{ii}]=ttest(ratesInValidBase{ii} - baseRate(ii));
    for jj = 1:numel(uTypes)
        [h(ii,jj),p(ii,jj),ci{ii,jj},stats(ii,jj)]=ttest(ratesInEvent{ii,jj} - localValidBaseRate(ii));
    end
    
    % report
    fprintf('Neuron %d: (Local) baseline rate = %0.3f Hz\n',ii,localValidBaseRate(ii));
    for jj = 1:numel(uTypes)
        mR(ii,jj)=mean(ratesInEvent{ii,jj});
%         fprintf('Syllable %s (n = %03d): Rate = %0.3f Hz (p = %0.3f)%s\n', ...
%             uTypes{jj}, sum(rIdxType==jj), mR(ii,jj), p(ii,jj),'*'*(p(ii,jj)<0.05 && mR(ii,jj) ~= 0) );
    end
%      fprintf('%6d | %6.3f | %4.3f | %4.3f | %4.3f | %4.3f | %8.3f | %4.3f \n', ...
%          ii, baseRate, ...
%          localValidBaseRate, p,  ...
%          mean(ratesInAdultSong),p2, ...
%          mean(ratesInJuvieSong),p3);
end

dataTable = [localValidBaseRate' mR p];
fieldNames = ['Local baseline', strcat('FR ', uTypes), strcat('p ', uTypes)];

hTable = uitable('Data', dataTable, 'ColumnName',fieldNames);
hExtent = get(hTable,'Extent'); %hExtent(4) = 800;
set(hTable, 'Position',hExtent);

kosherFieldNames = strrep(fieldNames, ' ', '_');
neuronData = initEmptyStructArray(kosherFieldNames, nNeurons);
for ii = 1:numel(kosherFieldNames)
    foo = num2cell(dataTable(:,ii));
    [neuronData.(kosherFieldNames{ii})] = foo{:};
end
% TODO: what to send back?
