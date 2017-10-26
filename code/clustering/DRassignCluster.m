function [assignedSylls, cosims] = DRassignCluster(syllables, empMatrices, empDistrs, featureTable, spectra)
% precondition: syllables contains a mix of labeled and unlabled syllables
% the empirical distribution and empirical distances correspond in order to the labeled syllables
% also, featureTable and spectra must be correctly formatted
% featureTable has to have cell of same features, with empty cells for
% uncalculated features
% spectra must have struct of arrays of the same size

featuresCached = (nargin >= 4) & numel(syllables) == numel(featureTable);
specsCached    = (nargin >= 5) & numel(syllables) == numel(spectra);

% todo: make robust to charcter input
N = numel(syllables);
isLabeled = ~isempty([syllables.type]) & ~isnan([syllables.type]);  

unlabIdx = find(~isLabeled);
  labIdx = find( isLabeled);

nUnlabeled = sum(~isLabeled);
  nLabeled = sum( isLabeled);

NLC2 = nchoosek(nLabeled,2);

% precond: the empirical distribution and empirical distances correspond in order to the labeled syllables
assert(all(NLC2 == cellfun(@(x) numel(empMatrices.(x)), fieldnames(empMatrices))));

% preparation for calculating spectra
fieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'};
% store the feature-based spectra for all of them
params.fine.features = {'wienerEntropy','deriv','harmonicPitch','fundamentalFreq'};

% get sampling rate
[filePath, fileStem] = fileparts(syllables(1).file);
metaFile = [filePath filesep 'meta-' fileStem];
metaStruct = []; load(metaFile);
params.fine.fs = 1/metaStruct.interval;

if nargin < 4
    featureTable = cell(1,N);
end
if nargin < 5
    spectra = initEmptyStructArray(fieldsToKeep, N);
end
    
% calculate spectra

if ~specsCached || ~featuresCached    
    missingFeatures = cellfun('isempty', featureTable);
    missingSpecs    = cellfun('isempty', {spectra.times});
    missingData = missingFeatures | missingSpecs;
    
    [remadeFeats, remadeSpecs] = getFeatures([],syllables(missingData), params);
    spectra(missingData) = remadeSpecs;
    featureTable(remadeFeats) = remadeFeats;
    featureTable = [featureTable{:}];
end

%% calculate local distances, TIME WARPED version
tic
progressbar('Unlabeled time-warped','Labeled');
dists.warpedLocal = zeros(nUnlabeled, nLabeled);
for ii = 1:nUnlabeled
    for jj = 1:nLabeled        
        dists.warpedLocal(ii,jj) = timeWarpedDistance(...
            spectra(unlabIdx(ii)), spectra(labIdx(jj)));        
        progressbar((ii-1)/nUnlabeled,jj/nLabeled);        
    end
end            
progressbar(1,1);
tt=toc;
fprintf('Time warping took %0.2f s...\n', tt);

%% step 5: measure global distances within pairs of syllables
% convert features from struct array to 2D array
fn = fieldnames(featureTable);
featureTable = cellfun(@(x) [featureTable.(x)]', fn', 'UniformOutput',false);
featureTable = [featureTable{:}];

% rows represent features, columns represent syllables
fprintf('Calculating global dissimilarity scores...\n');
zNormedFeatures = zscore(featureTable);
dists.global = zeros(nUnlabeled, nLabeled);
progressbar('Unlabeled global','Labeled');
for ii = 1:nUnlabeled
    for jj = 1:nLabeled
        dV = zNormedFeatures(unlabIdx(ii),:) - zNormedFeatures(labIdx(jj),:);
        dists.global(ii,jj) = sqrt(dot(dV,dV));
        progressbar((ii-1)/nUnlabeled,jj/nLabeled);
    end
end
progressbar(1,1);

%% calculate empirical distributions

% lookup on eCDF
emp.warpedLocal = interp1(empDistrs.warpedLocal(2,:), empDistrs.warpedLocal(1,:), dists.warpedLocal);
emp.global      = interp1(empDistrs.global     (2,:), empDistrs.global     (1,:), dists.global     );
emp.fused       = sqrt(emp.warpedLocal .* emp.global); 

%% lookup co-dissimilarity to the labeled vectors
cosims = zeros(nUnlabeled, nLabeled);
fusedLabeled = squareform(empMatrices.warpedLocal .* empMatrices.global);

% voting for membership is simple max, not sum
assignedSylls = syllables(~isLabeled);
labTypes = [syllables(labIdx).type];
nTypes = max(labTypes);
for ii = 1:nUnlabeled
    for jj = 1:nLabeled
        foo = corrcoef(emp.fused(ii,:), fusedLabeled(jj,:));
        cosims(ii,jj) = foo(2,1); % cosims is the final similarity metric
    end
    [~, foo] = sort(cosims(ii,:));
    ordSim(foo) = nLabeled:-1:1;
    
    cumScore = NaN(1, nTypes);
	ordScore = NaN(1, nTypes);    
    for jj = 1:nTypes
        cumScore(jj) = mean(cosims(ii,labTypes == jj));
        ordScore(jj) = mean(ordSim(labTypes == jj));
    end
    [~, maxOrdType ] = max(cumScore);
    [~, maxPropType] = max(ordScore);
    [~, maxID]    = max(cosims(ii,:));
    
    assignedSylls(ii).type = labTypes(maxID);
    assignedSylls(ii).propVoteType = labTypes(maxPropType);
    assignedSylls(ii).ordVoteType  = labTypes(maxOrdType);    
end

end