function multiPSTHselect(birdID, birdAge)

clustDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir  = [pwd filesep 'data' filesep            birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
if ~any(uAges==birdAge)
    error('Age not found for bird %s', birdID)
end

%% make all neurons psths for individual syllables for R204/age 53, R247, Lb277, giving coreness
% neuron loop - first loop over sessions, then gather spikes, then gather
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
for ii = 1:numel(uSessions) % loop through sessions
    % grab the neurons and metadata for each session
    matchingRecord = report(strcmp({report.sessionID},uSessions{ii}));
    matchingSpikeFiles = matchingRecord.spikeFiles;
    [spikes, nNeuronsPerFile, isMUA] = loadSpikeData(matchingSpikeFiles);
    isCoreUnit = false(numel(spikes),1);
    
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
    strrep(cellstr(strcat('cluster_', num2str(clusterNum))),' ' ,''); [syllables.type] = ans{:};
    
    nR = 6; % number of rows for figure
    minSyllInstances = 50;
    for jj = 1:numel(spikes)
        figure('Color',[1 1 1]); % whiten figure
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
        perSyllable = initEmptyStructArray({'avgRate','semRate','rates','hspec','maxYFR'},nTypes);
        for kk = 1:nTypes
            % loop through syllable types
            sType = acceptedTypes(kk);
            fprintf('Plotting neuron #%d/%d, syllable #%d...\n',jj,numel(spikes),sType);
            thisCluster = syllables(centralOrder{sType});
            
            % toss all but the first N syllables to equalize counts
            %{
            if numel(thisCluster) > minSyllInstances
                thisCluster(minSyllInstances+1:end) = [];
            end
            %}
            thisClusterParent = addPrePost(thisCluster, defaultParams, ...
                'preroll', preMotorT, 'postroll', postRollT);
            [~,spikesInEvent, perSyllable(kk).rates] = countSpikes(thisClusterParent, thisSpikeTrain);
            
            perSyllable(kk).avgRate = mean(perSyllable(kk).rates);
            if numel(thisCluster) > 1
                perSyllable(kk).semRate = std(perSyllable(kk).rates)/sqrt(numel(thisCluster)-1);
            end
            
            % example syllable
            dieRoll = randi(numel(thisCluster));
            exEvent = thisCluster(dieRoll);
            [cl,fs] = getClipAndProcess([], exEvent, defaultParams, 'noroll');
            fP = getfield(defaultParams,'best'); fP.fs = fs;
            exSpec = getMTSpectrumStats(cl, fP);
            
            % set up children events
            for ll = 1:numel(thisCluster)
                thisClusterParent(ll).children = thisCluster(ll);
            end
            
            % plot rasters in middle figure
            subplot(nR, nC, (kk + nC):nC:(kk + (nR-2)* nC))
            plotRaster(thisClusterParent, spikesInEvent, 'warping');
            set(gca,'Box','off','FontName','MyriadPro-Regular');
            if kk > 1, ylabel(''); else ylabel('Iteration number (central first)'); end
            figPos = get(gca,'Position');
            xlims = xlim;
            absPos = figPos(1) + figPos(3) * ([preMotorT/1000 (xlims(2) - postRollT/1000)] - xlims(1)) / diff(xlims);
            
            % add zero to xTicks
            set(gca,'XTick', 0:preMotorT/1000:max(xlim));
            set(gca,'XTickLabel', num2str(str2num(get(gca,'XTickLabel')) * 1000 - preMotorT));
            set(gca,'TickDir','out');
            
            % plot the PSTH for a syllable
            subplot(nR, nC, kk + (nR-1) * nC)
            plotPSTH(thisClusterParent, spikesInEvent, 'warping');
            title(sprintf('FR = %0.3g \\pm %0.3g', perSyllable(kk).avgRate, perSyllable(kk).semRate));
            set(gca,'Box','off','FontName','MyriadPro-Regular');
            perSyllable(kk).maxYFR = max(ylim);
            if kk > 1, ylabel(''); else ylabel('Count'), end
            xlabel('Time (ms)');
            
            % add zero to xTicks
            set(gca,'XTick', 0:preMotorT/1000:max(xlim));
            set(gca,'XTickLabel', num2str(str2num(get(gca,'XTickLabel')) * 1000 - preMotorT));
            set(gca,'TickDir','out');
            title(sprintf('FR = %0.3g \\pm %0.3g', perSyllable(kk).avgRate, ...
                perSyllable(kk).semRate));
            % get the x axis coords:
            
            % find any syllable - todo: align syllable size to absolute
            % position of shaded area
            % how? find absolute figure position of 0.5 to end-0.1 in
            % PSTH/raster
            perSyllable(kk).hspec = subplot('Position', [absPos(1) 0.80 diff(absPos) 0.12]);
            
            plotSpectrogram(exSpec);
            set(gca,'XTick',[],'YTick', []);
            set(gca,'Box','off','FontName','MyriadPro-Regular');
            xlabel(''); ylabel('');
        end %end syllable type loop
        
        % adjusting tick marks
        for kk = 1:nTypes
            subplot(nR, nC, kk + (nR-1) * nC)
            ylim([0 max([perSyllable.maxYFR])]);
            subplot(nR, nC, (kk + nC):nC:(kk + (nR-2)* nC))
        end
        
        % finishing - calculate 1-factor anova on syllable type
        syllId = cell(1,nTypes);
        for kk = 1:nTypes
            syllId{kk} = ones(size(perSyllable(kk).rates)) * kk;
        end
        [anpval, anstats, antable] = anova1([perSyllable.rates]', [syllId{:}]', 'off');
        
        % finishing - title
        axes(perSyllable(1).hspec);
        avgRate = [perSyllable.avgRate];
        aF = (1- (sum(avgRate)/numel(avgRate))^2/sum(avgRate.^2/numel(avgRate)))/(1 - (1/numel(avgRate)));
        if isCoreUnit(jj), neurtypeID = ' CORE neuron';
        else neurtypeID = 'SHELL neuron'; end;
        title(sprintf(' %s, %s, session %s, %d dph, SI = %0.2g, 1-factor anova F(1,196) = %0.3f, p = %0.3f', ...
            neurtypeID, birdID, uSessions{ii}, birdAge, aF, anstats{2,5}, anpval), ...
            'Interpreter','none', 'HorizontalAlignment', 'left');
        
        % font niceness
        set(findall(gcf,'type','text'),'FontSize', 12, 'FontName', 'MyriadPro-Regular')
        fprintf('Neuron #%d/%d done...\n',jj,numel(spikes));
    end
end
clear clusterIdxs
pause;
