function neuronData = compareTwoCatFiringRates(events, spikes, keyType1, keyType2, params, varargin)
% this compares the firing rates between events of key type 1 and events of
% key type 2

if nargin < 5
    params = defaultParams;
end
if nargin > 5
    params = processArgs(params, varargin{:});
end
% songstruct does not need values
% spikes is a cell array of different neurons

% to do: make multi version w/o variable-embedded indices
nNeurons = numel(spikes);

%syllLengths = [events.stop] - [events.start]; Jenny commented out
isKey1 = strcmp(keyType1, {events.type});
isKey2 = strcmp(keyType2, {events.type});

meanRS = zeros(1,nNeurons);
for ii = 1:nNeurons
    % count the spikes within each event
    [~,~, spikeRates1{ii}] = countSpikes(events(isKey1), spikes{ii});
    [~,~, spikeRates2{ii}] = countSpikes(events(isKey2), spikes{ii});
    % to do: add SEMs to output

    % get the rates and p-values for each neuron 
    % normality testing? we don't need no stinking normality testing
    % (unless our Ns are low);
    
    [~,pbase(ii)]=ttest2(spikeRates1{ii}, spikeRates2{ii});

%     % normalized response strength index, Jenny commented out
%     RSI(ii) = (mean(spikeRates1{ii}) - mean(spikeRates2{ii}))/(mean(spikeRates1{ii}) + mean(spikeRates2{ii}));
%     % report
%     if params.verbose
%         fprintf('Neuron %d: ''%s'' rate = %0.3f Hz, ''%s'' rate = %0.3f Hz, p = %0.2f, RSI = %0.2f \n', ...
%             ii, keyType1, mean(spikeRates1{ii}), keyType2, mean(spikeRates2{ii}), pbase(ii), RSI(ii));
%     end
      %Response strength (FR during event 1 - FR during baseline/event2), using average of nearest two baseline events as
      %baseline, Jenny wrote
      key1 = events(isKey1);
      key2 = events(isKey2);
      RS = zeros(1,length(key1));
      for jj = 1:length(key1) 
          startList = abs(vertcat(key2.start) - key1(jj).start); %baselines closest to event start
          stopList =  abs(vertcat(key2.start) - key1(jj).stop); %baselines closest to event stop
          [bothList, bothInd] = sort([startList; stopList]);
          bothInd(bothInd > length(startList)) = bothInd(bothInd > length(startList)) - length(startList);
          base1 = bothInd(1);
          if bothInd(2) == bothInd(1) %don't use the same baseline event twice
              base2 = bothInd(3);
          else
              base2 = bothInd(2);
          end
          matchBase = (spikeRates2{ii}(base1) + spikeRates2{ii}(base2)) / 2;
          RS(jj) = spikeRates1{ii}(jj) - matchBase;
      end
      meanRS(ii) = mean(RS);    
end

%%% plating
% package data for output argument
dataTable =   [cellfun(@mean,spikeRates1); cellfun(@mean, spikeRates2); pbase;  meanRS];
fieldNames = {['FR_' keyType1],           ['FR_' keyType2],             't_p', 'meanRS'};

kosherFieldNames = strrep(fieldNames, ' ', '_');
neuronData = initEmptyStructArray(kosherFieldNames, nNeurons);
for ii = 1:numel(kosherFieldNames)
    foo = num2cell(dataTable(ii,:));
    [neuronData.(kosherFieldNames{ii})] = foo{:};
end
