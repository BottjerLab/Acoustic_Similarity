function neuronData = spikingByCategoryTwoTail(songStruct, events, spikes)
% spikes is a cell array of different neurons
% get non song firing rates

nNeurons = numel(spikes);

syllLengths = [events.stop] - [events.start];
for ii = 1:nNeurons
% count the spikes within each event
    spikeCounts{ii} = countSpikes(events, spikes{ii},'onset');
    spikeRates{ii} = spikeCounts{ii} ./ syllLengths; 
end

% some more fundamental numbers
totalSpikeRates = sum(cellfun(@sum,spikeRates));
totalTime = songStruct.interval * songStruct.length;
fs = 1/songStruct.interval;

% construct local baseline %%%%
localBaseStarts = [max(0, [events.start] - 4)];
localBaseStops =  [max(0, [events.start] - 2)];

localBaseEvents = eventFromTimes(localBaseStarts, localBaseStops, fs);
localBaseNotSong = findUnion(findIntersection(localBaseEvents, events));

fullDur = sum([localBaseEvents.stop] - [localBaseEvents.start]);
notSongDur = sum([localBaseNotSong.stop] - [localBaseNotSong.start]);

%isSongDur = sum([localBaseInSong.stop] - [localBaseInSong.start]);

% temporary: add premotor events
premotorEvents = eventFromTimes([events.start]-0.1, [events.start],fs);
[premotorEvents.type] = deal('premotor');
events = [events premotorEvents];

% get event types
evTypes = {events.type};
[uEvTypes, ~, rIdxType] = unique(evTypes);

fprintf('Baseline samples %0.2fs, %0.2f%% not intersecting with other events...\n', ...
    fullDur, notSongDur/fullDur * 100);

spikesInEvent = cell(nNeurons,numel(uEvTypes));
ratesInEvent = cell(nNeurons,numel(uEvTypes));
spikesInValidRate = cell(nNeurons,1);
ratesInValidRate = cell(nNeurons,1);
mR = zeros(nNeurons, numel(uEvTypes));

for ii = 1:nNeurons   
    % get local baseline rate
    [spikesInValidBase{ii},~, ratesInValidBase{ii}] = countSpikes(localBaseNotSong,spikes{ii});
    localValidBaseRate(ii) = sum(spikesInValidBase{ii}) / notSongDur; 

    % get event firing rates
    for jj = 1:numel(uEvTypes)
        theseEvents = events(rIdxType == jj);
        [spikesInEvent{ii,jj}, ~, ratesInEvent{ii,jj}] = countSpikes(theseEvents,spikes{ii});
    end

   % get global (non-event) baseline rate
%    baseRate(ii) = (numel(spikes{ii}) - sum([spikesInEvent{ii,:}])) / ...
%        (totalTime - sum([events.stop]-[events.start]));
    
%    [hbase{ii},pbase{ii},cibase{ii},statsbase{ii}]=ttest(ratesInValidBase{ii} - baseRate(ii));

% compare each event type to the local baseline
    for jj = 1:numel(uEvTypes)
%        [h(ii,jj),p(ii,jj),ci{ii,jj},stats(ii,jj)]=ttest(ratesInEvent{ii,jj} - localValidBaseRate(ii));
        [h(ii,jj),p(ii,jj),ci{ii,jj},stats(ii,jj)]=ttest2(ratesInEvent{ii,jj}, ratesInValidBase{ii});
    end
    
    % report
    fprintf('Neuron %d: (Local) baseline rate = %0.3f Hz\n',ii,localValidBaseRate(ii));
    for jj = 1:numel(uEvTypes)
        mR(ii,jj)=mean(ratesInEvent{ii,jj});
end

% package data for output argument
dataTable = [localValidBaseRate' mR p];
fieldNames = ['Local baseline', strcat('FR ', uEvTypes), strcat('p ', uEvTypes)];

% matlab tables are the worst
%hTable = uitable('Data', dataTable, 'ColumnName',fieldNames);
%hExtent = get(hTable,'Extent'); %hExtent(4) = 800;
%set(hTable, 'Position',hExtent);

kosherFieldNames = strrep(fieldNames, ' ', '_');
neuronData = initEmptyStructArray(kosherFieldNames, nNeurons);
for ii = 1:numel(kosherFieldNames)
    foo = num2cell(dataTable(:,ii));
    [neuronData.(kosherFieldNames{ii})] = foo{:};
end
