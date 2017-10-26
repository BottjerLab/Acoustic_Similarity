function checkStereotypy(birdID, thisAge)
% plots stereotypical and not stereotypical syllables
clustDir = ['data' filesep 'cluster-' birdID filesep];
dataDir  = ['data' filesep            birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
if ~any(thisAge == uAges), error('checkStereotypy:AgeNotRepresented','Age %d not represented', thisAge); end

plotParams = processArgs(defaultParams, 'preroll', 5, 'postroll', 5, ...
    'dgram.minContrast', 3e-11);

report = reportOnData(birdID,'',defaultParams,'verbose',false);
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
    error('checkStereotypy:noLabels','No acceptedLabels file, # sylls = %d\n', numel(ageSylls));
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

% get the cluster distances matrices for ordering
fils = dir([clustDir 'altClustDataAge-' num2str(thisAge) '*.mat']);
if isempty(fils), error('checkStereotypy:noClusterFile','No cluster file found'); end
clustFile = [clustDir fils(1).name]; load(clustFile, 'distMats');

clusterMat = distMats.cosim;
[centralOrder, ~, innerClusterDists] = findMostCentral(clusterMat, clusterIdxs);
fprintf('Central-peripheral order found.\n');

nR = 18; % number of examples / number of rows
nTypes = max(clusterIdxs);
for kk = 1:nTypes
    % loop through syllable types
    sType = kk;
    % thisCluster = ageSylls(clusterIdxs == sType);
    thisOrderedCluster = ageSylls(centralOrder{sType});
    thisOrderedDists = sort(innerClusterDists{sType});
    fprintf('Type %d, # = %d\n', kk, numel(thisOrderedCluster));
    nEx = min(numel(thisOrderedCluster),nR);
    
    for ll = 1:nEx
        if ll <= nEx/2
            exEvent = thisOrderedCluster(ll);
            innerDist = thisOrderedDists(ll);
        else
            exEvent = thisOrderedCluster(ll + end - nEx);
            innerDist = thisOrderedDists(ll + end - nEx);
        end
%        [cl,fs] = getClipAndProcess([], exEvent, defaultParams, ...
%            'preroll', 5, 'postroll', 5);
%        fP = getfield(defaultParams, 'best'); fP.fs = fs;
%        exSpec = getMTSpectrumStats(cl, fP);
        
        figBox = [(kk-1)/nTypes (nR-ll)/nR 1/nTypes 1/nR];
        hspec(kk) = subplot('Position', figBox);
        plotDGFromRegion([], exEvent, plotParams);
%        plotDerivGram(exSpec, [], 'dgram.minContrast', 3e-12);
        set(gca,'XTick', [], 'YTick', [], 'Box', 'off');
        xlabel(''); ylabel('');
        
        text(1,0,num2str(innerDist, '%0.3f'), 'Units','normalized', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
            'Color', [1 1 1]);
        
        drawnow
    end
end %end syllable type loop
set(gcf,'Name',sprintf('%s, Age %d', birdID, thisAge));
annotation('line', [0 1], [0.5 0.5],'LineWidth', 2, 'Color', [0.8 0.1 0.1])