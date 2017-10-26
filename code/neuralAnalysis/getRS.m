function [neuronData, pairedFiring] = getRS(events, spikes, keyType1, keyType2, ptype, mRSType, mRSIType, standType, cvarFRType, params, varargin)
% this compares the firing rates between events of key type 1 and events of
% key type 2
% ptype is the field name of the probability (by t-test) that the syllables
% are different
% mRStype is the field name of the response strength (difference in firing)
% mRSItype is the field name of the response strength normalized by the sum
% standType is the field name of the response strength index normalized by
% the standard deviation
%
% warning: this method only takes string types
%
% pairedFiring is a cell array of the firing rate of each event and its paired baseline
if nargin <= 9 || isempty(params)
    params = defaultParams;
end
if nargin > 9
    params = processArgs(params, varargin{:});
end
% songstruct does not need values
% spikes is a cell array of different neurons

% to do: make multi version w/o variable-embedded indices
nNeurons = numel(spikes);

%syllLengths = [events.stop] - [events.start]; Jenny commented out
isKey1 = strcmp(keyType1, {events.type});
isKey2 = strcmp(keyType2, {events.type});

% if there are no intervals that are labeled with the correct key
if sum(isKey1) < 1 || sum(isKey2) < 1
    fieldNames = {['FR_' keyType1],           ['FR_' keyType2],            ptype, mRSType, mRSIType, standType, cvarFRType};
    kosherFieldNames = strrep(fieldNames, ' ', '_');
    neuronData = initEmptyStructArray(kosherFieldNames, nNeurons);
    for ii = 1:numel(kosherFieldNames)
        [neuronData.(kosherFieldNames{ii})] = deal(NaN);
    end
    pairedFiring = cell(1,nNeurons);
    return
end
meanRS = NaN(1,nNeurons);
meanRSI = NaN(1,nNeurons);
standRS = NaN(1,nNeurons);
cvarFR = NaN(1,nNeurons);
burstF = cell(1, nNeurons); %JMA
pbase = NaN(1,nNeurons);
matchBase = cell(1,nNeurons);
pairedFiring = cell(1,nNeurons);
for ii = 1:nNeurons
    % count the spikes within each event
    [~,~, spikeRates1{ii}] = countSpikes(events(isKey1), spikes{ii});
    [~,~, spikeRates2{ii}] = countSpikes(events(isKey2), spikes{ii});
    % to do: add SEMs to output        
   
    if all(spikeRates1{ii} == 0) && all(spikeRates2{ii} == 0), 
        % no spikes in either type of event
        % the default for each paired firing vector
        pairedFiring{ii} = zeros(2,sum(isKey1));
        continue;   
    end;    
    % we can't fill out the fields if there's no firing ever - there should
    % be a more specific place for this
    
    %[~,pbase(ii)]=ttest2(spikeRates1{ii}, spikeRates2{ii});

      %Response strength (FR during event 1 - FR during baseline/event2), using average of nearest two baseline events as
      %baseline, Jenny wrote
      key1 = events(isKey1);
      key2 = events(isKey2);
      RS = NaN(1,length(key1));
      BF = NaN(1,length(key1));
      RSI = zeros(1,length(key1));
      FR = zeros(1,length(key1));
      matchBase{ii} = zeros(1,length(key1)); %the events of baseline/event 2 that match the original event
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
          matchBase{ii}(jj) = (spikeRates2{ii}(base1) + spikeRates2{ii}(base2)) / 2;
          
          
          RS(jj) = spikeRates1{ii}(jj) - matchBase{ii}(jj);
          FR(jj) = spikeRates1{ii}(jj);
          RSI(jj) = (spikeRates1{ii}(jj) - matchBase{ii}(jj)) / (spikeRates1{ii}(jj) + matchBase{ii}(jj));
          if isnan(RSI(jj)), RSI(jj) = 0; end
      end
      pairedFiring{ii} = [spikeRates1{ii}; matchBase{ii}];
      % t-test to distinguish spiking from baseline in events
      [~,pbase(ii)]=ttest(spikeRates1{ii}, matchBase{ii});
      
      % response strengths - raw (RS), normalized to sum of response (RSI), and
      % standardized to variance of response (standRS)
      meanRS(ii) = mean(RS);
      meanRSI(ii) = mean(RSI);
      cvarFR(ii) = std(FR)/mean(FR); %coeffient of variation
      burstF{ii} = BF;
      
      % get the covariance of the spikes
      covarSB = cov(spikeRates1{ii},matchBase{ii});     
      % if no spikes ever occur in either baseline regions or ROI, 
      % then the covariance matrix will be scalar
      if numel(covarSB) == 1, 
          covarSB = 0; 
      else
          covarSB = covarSB(2,1);
      end
      
      %(singing - baseline) events in that order
      standRS(ii) = sqrt(numel(key1)) * (mean(spikeRates1{ii}) - mean(matchBase{ii}))/sqrt(var(spikeRates1{ii}) + ... 
          var(matchBase{ii}) - 2 * covarSB); %Jenny wrote this, standardized response
      
      % the standardize response could be 0 or infinite if the variance is
      % poorly behaved
      if isinf(standRS(ii)) || isnan(standRS(ii)) 
          standRS(ii) = 0.0;
      end
      
     % report
     if params.verbose
         fprintf('Neuron %d: ''%s'' rate = %0.3f Hz, matched ''%s'' rate = %0.3f Hz, p = %0.2f, RSI = %0.2f, standRS = %0.2f \n', ...
             ii, keyType1, mean(spikeRates1{ii}), keyType2, mean(matchBase{ii}), pbase(ii), meanRSI(ii), standRS(ii));
     end
end

%%% plating
% package data for output argument

% NB: the p-value here is unsigned
dataTable =   [cellfun(@mean,spikeRates1); cellfun(@mean, matchBase);   pbase;  meanRS; meanRSI; standRS; cvarFR];
fieldNames = {['FR_' keyType1],           ['FR_' keyType2],             ptype, mRSType, mRSIType, standType, cvarFRType}; 
kosherFieldNames = strrep(fieldNames, ' ', '_');
neuronData = initEmptyStructArray(kosherFieldNames, nNeurons);
for ii = 1:numel(kosherFieldNames)
    foo = num2cell(dataTable(ii,:));
    [neuronData.(kosherFieldNames{ii})] = foo{:};
end