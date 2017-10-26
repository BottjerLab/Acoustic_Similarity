birdID = 'Lb277';

clustDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep];
dataDir = [pwd filesep 'data' filesep birdID filesep];

load([dataDir 'allSpecs-' birdID '.mat'],'DRsylls');

uAges = unique([DRsylls.age]);
wantedAge = 59; % adhoc for now

load([clustDir 'T-09_23_09_10' filesep 'altClustDataAge-59-09_13_17_11.mat'])

load([dataDir 'acceptedLabels-' birdID '-age59.mat']);

%DRsylls = DRsylls([DRsylls.age] == wantedAge);

%distMats = distMats.adjustedCosimEven;
%%
params = defaultParams;
sOI = 2;
centralOrders = findMostCentral(distMats, clusterIdxs.accepted);% assemble a clip
% concatenate representative clips:

nEx = 5;
lens = zeros(1,nEx);
cl = cell(1,nEx);
for ii = 1:nEx
    rollLen = 3;    
    roi = addPrePost(DRsylls(centralOrders{sOI}(ii)), params, 'postroll', rollLen, 'preroll', rollLen);
    lens(ii) = roi.stop - roi.start;
    [cl{ii}, fs] = getClipAndProcess([], roi, params, 'noroll');
end
%%
concatCl = vertcat(cl{:});
fP = getfield(defaultParams, 'fine'); fP.fs = fs;
fP.features = [fP.features 'harmonicPitch'];

spec = getMTSpectrumStats(concatCl, fP);
hax = plotAllFigures(spec, [], params, 'optGraphs', {'deriv', 'AM', 'wienerEntropy', 'pitchGoodness', 'fundamentalFreq'});

barsx = cumsum(lens);
barsx = [barsx; barsx; NaN(size(barsx))];
%%
for ii = 1:numel(hax)
    axes(hax(ii));
    hold on;
    yy = ylim;
    barsy = barsx; barsy(1:2,:) = ylim' * ones(1,numel(lens));
    if ii ~= 1
        set(get(gca,'Children'),'LineWidth', 2);
        plot(barsx(:), barsy(:), 'k-' ,'LineWidth', 3);
    else        
        set(gca, 'XTick', [], 'YTick', []);
    end
    if ii == numel(hax)
        xlabel('Time (s)');
    end
    hold off;
    set(gca,'Box','off');
    set(gca,'FontName', 'MyriadPro-Regular');
    
end    
set(gcf,'Color', [1 1 1]);
set(findall(gcf,'type','text'),'FontSize', 12, 'FontName', 'MyriadPro-Regular')

%% 
export_fig([pwd filesep 'figures\JMAthesis\sampleRepeatedSyllableStatistics.pdf'])

