function allSyllables = stereotypedSylls(birdID, birdAge)

clustDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir  = [pwd filesep 'data' filesep            birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
if ~any(uAges==birdAge)
    error('Age %d not found for bird %s', birdAge, birdID)
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
    error('No acceptedLabels file, # sylls = %d\n', numel(ageSylls));    
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
allSyllables = [];
for ii = 1:numel(uSessions) % loop through sessions
    % grab the neurons and metadata for each session
    matchingRecord = report(strcmp({report.sessionID},uSessions{ii}));
    matchingSpikeFiles = matchingRecord.spikeFiles;
    [spikes, nNeuronsPerFile, isMUA] = loadSpikeData(matchingSpikeFiles);
    nNeurons = numel(spikes);
    isCoreUnit = false(nNeurons,1);
    
    % get the core/shell identity of each neuron
    cumIndex = [0 cumsum(nNeuronsPerFile)];
    for jj = 1:numel(matchingSpikeFiles)
        unitsFromFile = (cumIndex(jj)+1):cumIndex(jj+1);
        % core/shell ness depends on the file folder
        isCoreUnit(unitsFromFile) = ~isempty(strfind(matchingSpikeFiles{jj}, 'core'));
    end
    
    % get the cluster matrices for ordering
    fils = dir([clustDir 'mfccClustDataAge-' num2str(birdAge) '*.mat']); 
    if isempty(fils), error('No cluster file found'); end   
    clustFile = [clustDir fils(1).name];
    load(clustFile, 'distMats');
    clusterMat = distMats.cosim;
    centralOrder = findMostCentral(clusterMat, clusterIdxs);
    fprintf('Central-peripheral order found.\n');
    % get subsong/plastic ness of session
    [~, isThisPlastic] = getAgeOfSession(uSessions{ii});
    
    % get syllables specific to this session
    isThisSession = strcmp(sForSylls, uSessions{ii});
    syllables = ageSylls(isThisSession);
    clusterNum = clusterIdxs(isThisSession);
    
    % transfer the syllable labels to the age-specific syllables
    num2cell(clusterNum); [ageSylls(isThisSession).type] = ans{:};
    
    % assign cluster types in string form to this subset of syllables
    % for a given region
    % NB: is this redundant?
    %strrep(cellstr(strcat('cluster_', num2str(clusterNum))),' ' ,''); [syllables.type] = ans{:};
    
    nR = 6; % number of rows for figure
    minSyllInstances = 20;
    cols = jet(nNeurons);
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
        postRollT = 10; %ms
        eachSyllable = initEmptyStructArray({'avgRate','semRate','rates', ...
            'centralities', 'centralRate', 'peripheralRate', 'isCore', 'isPlastic'},nTypes);
        for kk = 1:nTypes
            % loop through syllable types
            sType = acceptedTypes(kk);
            fprintf('Plotting neuron #%d/%d, syllable #%d...\n',jj,numel(spikes),sType);
            
            % sort the cluster only within this session
            isInSession = isThisSession(centralOrder{sType});            
            orderInSession = centralOrder{sType}(isInSession);
            thisCluster = ageSylls(orderInSession);
            nInCluster = numel(thisCluster);            
            num2cell((1:nInCluster)/nInCluster); [thisCluster.rank] = ans{:};
            
            thisClusterParent = addPrePost(thisCluster, defaultParams, ...
                'preroll', preMotorT, 'postroll', postRollT);
            [~,spikesInEvent, eachSyllable(kk).rates] = countSpikes(thisClusterParent, thisSpikeTrain);
            eachSyllable(kk).centralities = [thisCluster.rank];
            eachSyllable(kk).avgRate = mean(eachSyllable(kk).rates);
            if numel(thisCluster) > 1
                eachSyllable(kk).semRate = std(eachSyllable(kk).rates)/sqrt(numel(thisCluster)-1);
            end
            eachSyllable(kk).centralRate = mean([eachSyllable(kk).rates([thisCluster.rank] < 0.5)]);
            eachSyllable(kk).peripheralRate = mean([eachSyllable(kk).rates([thisCluster.rank] > 0.5)]);
            % plot central vs peripheral
            %{
            figure('Color',[1 1 1]); % whiten figure
            plot(1:2, [eachSyllable(kk).centralRate eachSyllable(kk).peripheralRate], '-', 'Color', cols(jj,:));
            xlabel('Syllable centeredness');
            xlim([0.5 2.5])
            set(gca,'XTick', [1 2]);
            set(gca,'XTickLabel',{'Central','Peripheral'});
            ylabel('Firing rate (Hz)');
            hold on;
            %}
        end 
        [eachSyllable.isCore] = deal(isCoreUnit(jj));
        [eachSyllable.isPlastic] = deal(isThisPlastic);
        allSyllables = [allSyllables; eachSyllable];
    end
end
clear clusterIdxs
