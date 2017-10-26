function checkSyllableLabels(birdID)

clustDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir  = [pwd filesep 'data' filesep            birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
%% make all neurons psths for individual syllables for R204/age 53, R247, Lb277, giving coreness
% neuron loop - first loop over sessions, then gather spikes, then gather
report = reportOnData(birdID,'',defaultParams,'verbose',false);   
for hh = 1:numel(uAges)
    thisAge = uAges(hh);
    %load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');
    ageSylls = DRsylls([DRsylls.age]==thisAge);
        
    % get the session ID for each syllable
    sForSylls = cell(1,numel(ageSylls));
    for ii = 1:numel(ageSylls)
        [~,sForSylls{ii}] = fileparts(ageSylls(ii).file);
    end
    uSessions = unique(sForSylls);
    fprintf('Bird %s, age %d\n', birdID, thisAge);
    
    % get the last type of cluster labels
    acceptedLabelFile = [dataDir 'acceptedLabels-' birdID '-age' num2str(thisAge) '.mat'];
    if ~exist(acceptedLabelFile, 'file')
        fprintf('No acceptedLabels file, # sylls = %d\n', numel(ageSylls));
        continue;
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
    
    % assume that clusterIdxs all have the same structure for now - not safe
%{
    while ~exist('clusterIdxs', 'var') || isstruct(clusterIdxs)        
        fprintf('Import me the syllable labels for bird %s, age %d, and place in clusterIdxs as a numeric vector, any way you can... ', ...
            birdID, thisAge);
        keyboard
    end
%}
    nR = 18; % number of examples / number of rows
    nTypes = max(clusterIdxs);
    for kk = 1:nTypes
        % loop through syllable types
        sType = kk;
        thisCluster = ageSylls(clusterIdxs == sType);
        fprintf('Type %d, # = %d\n', kk, numel(thisCluster));
        nEx = min(numel(thisCluster),nR);
        for ll = 1:nEx
            exEvent = thisCluster(ll);
            [cl,fs] = getClipAndProcess([], exEvent, defaultParams, ...
                'preroll', 5, 'postroll', 5);
            fP = getfield(defaultParams, 'best'); fP.fs = fs;
            exSpec = getMTSpectrumStats(cl, fP);
            
            figBox = [(kk-1)/nTypes (nR-ll)/nR 1/nTypes 1/nR];
            hspec(kk) = subplot('Position', figBox);
            plotDerivGram(exSpec, [], 'dgram.minContrast', 3e-12);
            set(gca,'XTick', [], 'YTick', [], 'Box', 'off');
            xlabel(''); ylabel('');
            
            % plot 10 ms scale bar
            xx = xlim; yy = ylim; 
            ybar = yy * [0.1 0.9]'; 
            xbar = xx * [0.1 0.9]'; 
            hold on;
            plot([(xbar-0.01) xbar], [ybar ybar], 'w-', 'LineWidth', 3);
            hold off;
            drawnow
        end
    end %end syllable type loop
    set(gcf,'Name',sprintf('%s, Age %d', birdID, thisAge));
    pause;
end