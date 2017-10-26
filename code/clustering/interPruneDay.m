function interPruneDay(birdID, age)
%% load day's worth of syllables
dataDir = 'C:\Users\Johnny\Documents\LMAN\data\';
load([dataDir birdID filesep 'allSpecs-' birdID '.mat']);
isAge = ([DRsylls.age] == age);
%featureTable = featureTable(isAge);
clear spectra; clear featureTable;

unlabeledSylls = DRsylls(isAge);
[unlabeledSylls.type] = deal('');

%% load appropriate cluster file

clustDir = [dataDir 'cluster-' birdID filesep]
appropFileGlob = [clustDir '*-' num2str(age) '-*.mat'];
fils = dir(appropFileGlob);
fils = sortBy(fils, 'datenum', 'descend');
clustFil = [clustDir fils(1).name];
load(clustFil);

%% recalculate the clusters - this is important for browseClusters should we choose to use it
nClusts = 6:16;
cIdxs = zeros(numel(unlabeledSylls),numel(nClusts));
linkArr = linkage(distMats.cosim,'complete');
for ii = 1:numel(nClusts);
    cIdxs(:,ii) = cluster(linkArr, 'maxclust', nClusts(ii));
end

%% run multiMark w/ new interactive set
close all

params = defaultParams;
params.doNoiseFilter = false;
params.dgram.minContrast = 3e-11;
%%
len = num2cell([unlabeledSylls.stop] - [unlabeledSylls.start]);
[unlabeledSylls.length] = len{:};
%
close all
[unlabeledSylls, rInd] = sortBy(unlabeledSylls,'length');
%%
dMat = squareform(distMats.cosim);

%%
% isolate different types of sounds
cementedSylls = initEmptyStructArray(fieldnames(unlabeledSylls),0); 
contLoop = true;
dSubMat = dMat(rInd, rInd);
while contLoop && numel(unlabeledSylls) > 0
    fprintf('Isolate given sounds...\n');
    isSoundType = markForSim([], unlabeledSylls, dSubMat, params);
    
    soundLabel = input('What kind of sound is this? ', 's');
    [unlabeledSylls(isSoundType).type] = deal(soundLabel);
    cementedSylls = [cementedSylls unlabeledSylls(isSoundType)];
    
    unlabeledSylls = unlabeledSylls(~isSoundType);  
    dSubMat = dSubMat(~isSoundType, ~isSoundType);
     
    contLoop = ~strcmp('n', input('Continue isolating sounds [Y/n]? ', 's'));
end

%% save work
cementedSylls = sortBy([cementedSylls unlabeledSylls], 'start');
cemFil = [dataDir birdID filesep 'cementedSylls-' birdID '-' num2str(age) '.mat'];
uisave('cementedSylls', cemFil);