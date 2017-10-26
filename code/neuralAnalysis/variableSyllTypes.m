function allClusterRecords = variableSyllTypes(birdID, birdAge)
% get the variability/inter-cluster distance for each labeled syllable of a
% given bird/age
clustDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir  = [pwd filesep 'data' filesep            birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
if ~any(uAges==birdAge)
    error('variableSyllTypes:ageNotFound', 'Age %d not found for bird %s', birdAge, birdID)
end

%% neuron loop - first loop over sessions, then gather spikes, then gather
report = reportOnData(birdID,'',defaultParams,'verbose',false);

%load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');
ageSylls = DRsylls([DRsylls.age]==birdAge);

% get the session ID for each syllable
sForSylls = cell(1,numel(ageSylls));
for ii = 1:numel(ageSylls)
    [~,sForSylls{ii}] = fileparts(ageSylls(ii).file);
end
uSessions = unique(sForSylls);

fprintf('Bird %s, age %d\n', birdID, birdAge);

% get the last type of cluster labels
acceptedLabelFile = [dataDir 'acceptedLabels-' birdID '-age' num2str(birdAge) '.mat'];
if ~exist(acceptedLabelFile, 'file')
    error('variableSyllTypes:labelsNotFound', 'No acceptedLabels file, # sylls = %d\n', numel(ageSylls));    
end
clusterIdxs = []; load(acceptedLabelFile);
clTypes = fieldnames(clusterIdxs);
if numel(clTypes) == 3,
    prefType = clTypes{2};
else
    prefType = clTypes{1};
end
clusterIdxs = clusterIdxs.(prefType);
fprintf('Total syllables IDed / total: (%d/%d)\n', sum(~isnan(clusterIdxs)), numel(clusterIdxs));

%%
% get the cluster matrices for ordering
fils = dir([clustDir 'altClustDataAge-' num2str(birdAge) '*.mat']);
if isempty(fils), error('variableSyllTypes:clusterDistssNotFound','No cluster file found'); end
clustFile = [clustDir fils(1).name];
load(clustFile, 'distMats');
clusterMat = distMats.cosim;

% these values are celled by syllable
[centralOrder, innerDists] = findMostCentral(clusterMat, clusterIdxs);

% these values are not celled by syllable - over all syllables
[intraClusterDists, interClusterDists] = clusterDistances(clusterMat, clusterIdxs);

fprintf('Central-peripheral order found.\n');

% retrieve multiunit data
fprintf('Retrieving MUA data...\n');
MUAinfo = getSpikeMUAData;

allClusterRecords = [];
neuronIndex = 1;
for ii = 1:numel(uSessions) % loop through sessions
    thisSession = uSessions{ii};
    % grab the neurons and metadata for each session
    matchingRecord = report(strcmp({report.sessionID},thisSession));
    matchingSpikeFiles = matchingRecord.spikeFiles;
    [spikes, nNeuronsPerFile, isMUA] = loadSpikeData(matchingSpikeFiles, MUAinfo);
    nNeurons = numel(spikes);
    isCoreUnit = false(nNeurons,1);
    
    % get the core/shell identity of each neuron
    cumIndex = [0 cumsum(nNeuronsPerFile)];
    for jj = 1:numel(matchingSpikeFiles)
        unitsFromFile = (cumIndex(jj)+1):cumIndex(jj+1);
        % core/shell ness depends on the file folder
        isCoreUnit(unitsFromFile) = ~isempty(strfind(matchingSpikeFiles{jj}, 'core'));
    end
    
    % get subsong/plastic ness of session
    [~, isThisPlastic] = getAgeOfSession(thisSession);
    
    % get syllables specific to this session
    isThisSession = strcmp(sForSylls, thisSession);
    syllables = ageSylls(isThisSession);
    clusterNum = clusterIdxs(isThisSession);
    
    % transfer the syllable labels to the age-specific syllables
    num2cell(clusterNum); [ageSylls(isThisSession).type] = ans{:};
    
    minSyllInstances = 20;
    
    % load baseline
    baselinePeriods = [];
    baselinePeriodFile = [dataDir 'baselinePeriods-' thisSession '.mat'];
    load(baselinePeriodFile);
    [baselinePeriods.type] = deal('silence');
    
    for jj = 1:nNeurons
        thisSpikeTrain = spikes{jj}; 
        nTypes = max(clusterNum); %
        
        % prescreen for minimal number of syllables
        included = false(1,nTypes);
        for kk = 1:nTypes
            included(kk) = (sum(clusterNum == kk) >= minSyllInstances);
        end
        
        acceptedTypes = find(included);
        nTypes = numel(acceptedTypes); % number of syllables
        nC = nTypes;
        
        preMotorT = 50; %ms
        postRollT = 0;
        eachClusterRecord = initEmptyStructArray({'neuronIndex', 'syllableIndex', 'avgRate','semRate','rates', ...
            'centralRate', 'peripheralRate', 'isMUA','isCore', 'isPlastic',...
            'meanVariability', 'medianVariability', ...
            ... %'allNormRS', 'allProb', ...
            'centralNormRS', 'centralProb', 'periphNormRS', 'periphProb'},nTypes);
        for kk = 1:nTypes
            % loop through syllable types
            sType = acceptedTypes(kk);
            
            eachClusterRecord(kk).neuronIndex   = neuronIndex;
            eachClusterRecord(kk).syllableIndex = sType;
            fprintf('Plotting neuron #%d/%d, syllable #%d...\n',jj,numel(spikes),sType);
            
            % sort the cluster only within this session
            isInSession = isThisSession(centralOrder{sType});            
            orderInSession = centralOrder{sType}(isInSession);
            thisCluster = ageSylls(orderInSession)';
            
            theseIntraDistances = intraClusterDists(orderInSession);
            theseInterDistances = interClusterDists(orderInSession);
            nInCluster = numel(thisCluster);            
            num2cell((1:nInCluster)/nInCluster); [thisCluster.rank] = ans{:};
            
            thisClusterWithLag = addPrePost(thisCluster, defaultParams, ...
                'preroll', preMotorT, 'postroll', postRollT);
            
            [~,spikesInEvent, eachClusterRecord(kk).rates] = countSpikes(thisClusterWithLag, thisSpikeTrain);
            
            eachClusterRecord(kk).rankedness = [thisCluster.rank]; %of all syllables
            eachClusterRecord(kk).avgRate = mean(eachClusterRecord(kk).rates);
            if numel(thisCluster) > 1
                eachClusterRecord(kk).semRate = std(eachClusterRecord(kk).rates)/sqrt(numel(thisCluster)-1);
            end
            isCentral = [thisCluster.rank] < 0.25;
            isPeripheral = [thisCluster.rank] > 0.75;
            eachClusterRecord(kk).centralRate    = mean([eachClusterRecord(kk).rates(isCentral)]);
            eachClusterRecord(kk).peripheralRate = mean([eachClusterRecord(kk).rates(isPeripheral)]);
            eachClusterRecord(kk).meanVariability   = mean  (innerDists{sType});
            eachClusterRecord(kk).medianVariability = median(innerDists{sType});
            eachClusterRecord(kk).intraDists = theseIntraDistances;
            eachClusterRecord(kk).interDists = theseInterDistances;
            
            [thisClusterWithLag(isCentral).type]    = deal('central');
            [thisClusterWithLag(isPeripheral).type] = deal('peripheral');
            
            stdClustLag = thisClusterWithLag;
            fn = fieldnames(thisClusterWithLag);
            for ll = 1:numel(fn)
                if ~isfield(baselinePeriods, fn{ll})
                    stdClustLag = rmfield(stdClustLag, fn{ll});
                end
            end
            
            centralStat = getRS([baselinePeriods; stdClustLag], {thisSpikeTrain}, 'central', 'silence',...
                'p', 'RS', 'RSI', 'normRS');
            periphStat  = getRS([baselinePeriods; stdClustLag], {thisSpikeTrain}, 'peripheral', 'silence',...
                'p', 'RS', 'RSI', 'normRS');
            %[stdClustLag.type] = deal('syllable'); % now test 
            %allSyllsStat = getRS([baselinePeriods; stdClustLag], {thisSpikeTrain}, 'syllable', 'silence',...
            %    'p', 'RS', 'RSI', 'normRS');
            %eachClusterRecord(kk).allNormRS    = allSyllsStat.normRS;            
            eachClusterRecord(kk).centralNormRS =  centralStat.normRS;
            eachClusterRecord(kk).periphNormRS  =   periphStat.normRS;
            % probability that firing rate is significantly above baseline for these groups
            %eachClusterRecord(kk).allProb       = allSyllsStat.p * sign(allSyllsStat.RS);                        
            eachClusterRecord(kk).centralProb   =  centralStat.p * sign( centralStat.RS);
            eachClusterRecord(kk).periphProb    =   periphStat.p * sign(  periphStat.RS);
        end 
        % the attributes of each neuron/session pair
        [eachClusterRecord.isCore]    = deal(isCoreUnit(jj));
        [eachClusterRecord.isMUA]     = deal(isMUA(jj)); 
        [eachClusterRecord.isPlastic] = deal(logical(isThisPlastic));
        [eachClusterRecord.sessionID] = deal(thisSession);
        allClusterRecords = [allClusterRecords; eachClusterRecord];
        neuronIndex = neuronIndex + 1;
    end
end
clear clusterIdxs
