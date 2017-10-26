birdID = 'Lb277';
compositeReport = reportOnData(birdID,'',defaultParams,'verbose',true);
if ~iscell(compositeReport)
    compositeReport = {compositeReport};
end
nBirds = numel(compositeReport);
birdNeuronData = cell(nBirds, 1);

if ~exist('MUAinfo','var')
    MUAinfo = getSpikeMUAData;
end
%%
birdData = struct([]);
for hh = 1:nBirds
    sessionsForThisBird = compositeReport{hh};
    birdID = strtok(compositeReport{hh}(1).sessionID, '_');
    sessionData = reportOnData(birdID, '', [],'verbose',false);
    sessions = {sessionData.sessionID};
    nSessions = numel(sessions);
    ages = zeros(1,numel(sessions));
    for ii = 1:nSessions
        ages(ii) = getAgeOfSession(sessions{ii});
    end
    %uAges = unique(ages(~isnan(ages))); 
    uAges = 53;
    
    % FIX ME PLEASE JOHN
    %uAges(uAges == 57) = []; %for R247
    
    realAllNeurons = struct([]);
    for gg = 1: numel(uAges)
        
    pickAnAge = uAges(gg);
    
    birdDir    = [pwd filesep 'data' filesep birdID filesep];
    clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
    
    % syllable feature information
    load([birdDir 'allSpecs-' birdID '.mat']);
    
    %%
    
    %clusterIdxs = loadAcceptedLabels(birdID, pickAnAge)
    % if we have acceptedLabels, then use those
    
    acceptedLabelFile = [birdDir 'acceptedLabels-' birdID '-age' num2str(pickAnAge) '.mat'];
    if exist(acceptedLabelFile, 'file')
        load(acceptedLabelFile, 'clusterIdxs');
        if isfield(clusterIdxs,'augmented')
            clusterIdxs = clusterIdxs.augmented; %.augmented is best, accepted is ok
        else
            clusterIdxs = clusterIdxs.accepted;
        end
    else
        % if we don't have acceptedLabels, then look for machine-done labels
        % load the cluster file for a given day - a little bit of file
        % directory/string matching
        % pick the latest modified directory
        contents = dir(clusterDir);
        batchDirs = contents([contents.isdir]);
        
        [~, perm] = sort([batchDirs.datenum],'descend');
        batchDirs = {batchDirs(perm).name};
        batchDirs(strncmp('.',batchDirs,1)) = [];
        
        foundClusterData = false;
        ii = 0;
        
        while ~foundClusterData
            ii = ii + 1;
            if ii > numel(batchDirs), error('Couldn''t find any cluster index files, =('); end;
            
            wholeDir = [clusterDir batchDirs{ii} filesep];
            searchGlobStr = [wholeDir 'altClustDataAge-' num2str(pickAnAge) '*.mat'];
            contents = dir(searchGlobStr);
            if isempty(contents), continue; end;
            
            % check to see if it has clusterIdxs
            contents = {contents.name};
            for jj = 1:numel(contents)
                if any(strcmp('clusterIdxs', who('-file', [wholeDir contents{jj}]))) % we found some
                    clusterFile = [wholeDir contents{jj}];
                    clusterIdxs = []; % keep matlab happy
                    load(clusterFile, 'clusterIdxs');
                    if isfield(clusterIdxs,'adjustedCosimEven')
                        clusterIdxs = clusterIdxs.adjustedCosimEven(:,3); % pick a few clusters only                        
                        foundClusterData = true;
                        break;
                    else
                        continue;
                    end
                end
            end
        end
    end
    
    %%
    DRsylls = DRsylls([DRsylls.age]==pickAnAge);
    featureTable = featureTable ([DRsylls.age]==pickAnAge);
    spectra = spectra([DRsylls.age]==pickAnAge);
        
    sForSylls = cell(1,numel(DRsylls));
    for ii = 1:numel(DRsylls)
        [~,sForSylls{ii}] = fileparts(DRsylls(ii).file);
    end
    %%
    uSessions = unique(sForSylls);
    neuronData = struct([]);
    
    meanRS = zeros(max(clusterIdxs), numel(uSessions));
    
    for ii = 1:numel(uSessions) % loop through sessions
        matchingRecord = sessionData(strcmp(uSessions{ii}, sessions));
        matchingSpikeFiles = matchingRecord.spikeFiles;
        [spikes, nNeuronsPerFile, isMUA] = loadSpikeData(matchingSpikeFiles);
        isCoreUnit = false(numel(spikes),1);
        
        cumIndex = [0 cumsum(nNeuronsPerFile)];
        for jj = 1:numel(matchingSpikeFiles)
            unitsFromFile = (cumIndex(jj)+1):cumIndex(jj+1);
            
            % core/shell ness depends on the file folder
            isCoreUnit(unitsFromFile) = ~isempty(strfind(matchingSpikeFiles{jj}, 'core'));
        end
        
        [~, isThisPlastic] = getAgeOfSession(uSessions{ii});
        
        isThisSession = strcmp(sForSylls, uSessions{ii}); 
        syllables = DRsylls(isThisSession); 
        clusterNum = clusterIdxs(isThisSession); 
        
        % transfer the labels to the age-specific syllables
        num2cell(clusterNum); [DRsylls(isThisSession).type] = ans{:};
        
        % assign cluster types in string form to this subset of syllables
        % for a given region
        strrep(cellstr(strcat('cluster_', num2str(clusterNum))),' ' ,''); [syllables.type] = ans{:};

        % preallocate struct array
        finalNames = {'FR_silence'};
        for jj = 1:max(clusterNum)
            if any(clusterNum == jj)
                finalNames = [finalNames strcat({'FR_cluster_', 't_p_cluster_', 'standRS_cluster_'}, num2str(jj))];
            end
        end
        
        allNeuronsData = initEmptyStructArray(finalNames, numel(spikes));
        meanRS = zeros(max(clusterNum),numel(uSessions));
                
        for jj = 1:max(clusterNum) % loop through syllable types
            thisCluster = syllables(clusterNum== jj);
            if isempty(thisCluster), continue; end;
            if length(thisCluster) == 1
                continue; end;
            
            baselinePeriodFile = [birdDir 'baselinePeriods-' uSessions{ii} '.mat'];
            load(baselinePeriodFile);
            
            [baselinePeriodsS.age] = deal(pickAnAge);
            [baselinePeriodsS.file] = deal('whatever,dude');
            [baselinePeriodsS.type] = deal('silence');
            syllsWPre = addPrePost(thisCluster', params, 'preroll', 50, 'postroll', 0); %adding 50ms premotor
            
            neuronData{jj} = getRS([syllsWPre; baselinePeriodsS], spikes, ... %this put all neurons into neuronData for each syllable type, rewrite so each neuron gets responses to all syllable types
                ['cluster_' num2str(jj)], 'silence', ...
                't_p', 'meanRS', 'meanRSI', 'standRS',...
                params, 'verbose',false);  
            
            for kk = 1:numel(spikes)
                if jj == 1
                    [allNeuronsData(kk).FR_silence] = neuronData{jj}(kk).FR_silence;
                end
                
                newField = ['FR_cluster_' num2str(jj)];
                allNeuronsData(kk).(newField) = neuronData{jj}(kk).(newField);
                newField = ['standRS_cluster_' num2str(jj)];
                allNeuronsData(kk).(newField) = neuronData{jj}(kk).standRS;
                newField = ['t_p_cluster_' num2str(jj)];
                allNeuronsData(kk).(newField) = neuronData{jj}(kk).t_p;
         
                %maybe also store a psth for each syllable type?
            end
            meanRS(jj,ii) = nanmean([neuronData{jj}.standRS]);
        end %end syllable type loop
        %plot RS's for each neuron for each syllable type for this
        %session

        figure; nR = 2; nC = 4;
        figureCtr = 1;            
        for kk = 1: numel(spikes)
            responses = NaN(max(clusterNum),1);
            for jj = 1: max(clusterNum)
                thisField = ['standRS_cluster_' num2str(jj)];
                thatField = ['FR_cluster_' num2str(jj)];
                if isfield(allNeuronsData(kk), thisField) && ...
                        ~isempty(getfield(allNeuronsData(kk),thisField))
                    responses(jj) = getfield(allNeuronsData(kk),thisField);
                end
                if isfield(allNeuronsData(kk), thatField) && ...
                        ~isempty(getfield(allNeuronsData(kk),thatField))
                    firingRate(jj) = getfield(allNeuronsData(kk),thatField);
                end
            end
            subplot(nR,nC, figureCtr)
            bar(1:max(clusterNum),responses);
            title([birdID ' ' num2str(pickAnAge) ' dph, neuron ' num2str(kk)], 'Interpreter', 'none');
            figureCtr = figureCtr + 1;
            if figureCtr > nR * nC, figure; figureCtr = 1; end
            
            allNeuronsData(kk).activityFrac = (1- (sum(firingRate)/numel(firingRate))^2/sum(firingRate.^2/numel(firingRate)))/(1 - (1/numel(firingRate))); %no selectivity = 0, max selectivity = 1
        end
        
        foo = num2cell(isCoreUnit); [allNeuronsData.isCore] = foo{:};
        foo = num2cell(isMUA);  [allNeuronsData.isMUA ] = foo{:};
        [allNeuronsData.age] = deal(uAges(gg));
        [allNeuronsData.isPlastic] = deal(isThisPlastic);
        %%
        acrossNeuron = struct('activityFrac',cell(1,numel(allNeuronsData)));
        acrossNeuron = copyPartialStruct(allNeuronsData, acrossNeuron, ...
            {'activityFrac', 'isCore', 'isPlastic', 'isMUA'});
        realAllNeurons = [realAllNeurons acrossNeuron];
        %%
    end %end session loop
    
    
    
    %     figure;
    %     bar(1:max(clusterNum),mean(meanRS,2));
%     title([birdID ' ' num2str(pickAnAge) ' dph'], 'Interpreter', 'none');
    
    %% let me look at a syllable type    
%     figure;
%     plotParams = processArgs(defaultParams,...
%                             'dgram.minContrast', 1e-11, 'doFilterNoise', false,...
%                             'preroll', 200, 'postroll', 200);
%     mosaicDRSpec(DRsylls([DRsylls.type] == 4), plotParams, 'maxMosaicLength', 2.5);
    end
    
    isCore    = [realAllNeurons.isCore];
    isMUA     = [realAllNeurons.isMUA];
    isPlastic = [realAllNeurons.isPlastic];
    isCoreAndSUA  =  isCore & ~isMUA;
    isShellAndSUA = ~isCore & ~isMUA;
    
    figure;
    bins = 0:.1:1;
    y = histc([realAllNeurons(isCoreAndSUA).activityFrac],bins);
    bar(bins,y);
    title('CORE');
    
    figure;
    bins = 0:.1:1;
    y = histc([realAllNeurons(isShellAndSUA).activityFrac],bins);
    bar(bins,y);
    title('SHELL');
end

% Do neurons fire to one syllable type?
% I want to look at PSTHs of neurons firing to each syllable type
% Where are the significant peaks and troughs in the PSTHs

% Do neurons fire more when there is comparatively more variability in a
% Syllable compared to the other syllables of its type?
% Do neurons fire more similarly temporally for more similar syllables?

% Do the neurons fire randomly during a syllable type or in a time-aligned
%manner? cross-correlations might suck i know




