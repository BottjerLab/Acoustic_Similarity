function [clustIdxs, empMatrices, distMatrices, empDistrs] = DRcluster(songStruct,DRsylls,noiseMask,featureTable, spectra, params, varargin)

timeFlag = ['T-' datestr(clock, 'mm_dd_HH_MM')];
   
featuresCached = (nargin >= 2);
specsCached    = (nargin >= 3);
if nargin < 4
    params = defaultParams;
end
params = processArgs(params, varargin{:});

% remove any syllables that are too short
isTooShort = (params.fine.windowSize / 1000 > [DRsylls.stop] - [DRsylls.start]);
DRsylls(isTooShort) = [];
fprintf('Removing %d syllables that are too short...\n', sum(isTooShort));

N = numel(DRsylls);% = min(ceil(numel(allDRsylls)/2),1000);
NC2 = nchoosek(N,2);

nEmpD = min(NC2,2e4);

distMatrices = struct(...
    'warpedLocal'  , zeros(1,NC2),...
    'global'       , zeros(1,NC2));
empMatrices = struct(...
    'warpedLocal'  , zeros(1,NC2),...
    'global'       , zeros(1,NC2));
empDistrs = struct(...,
    'warpedLocal', zeros(2,nEmpD),...
    'global'     , zeros(2,nEmpD));

fieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'};
% store the feature-based spectra for all of them
params.fine.features = {'wienerEntropy','deriv','harmonicPitch','fundamentalFreq'};

% get sampling rate
params.fine.fs = 1/songStruct.interval;

% calculate spectra
if ~specsCached
    spectra = initEmptyStructArray(fieldsToKeep, N);
    if ~featuresCached
        featureTable = cell(1,N);
    end
    
    progressbar(sprintf('Calculating spectra & features for regions (# = %d)',N));
    for ii = 1:N
        cl = getClipAndProcess([],DRsylls(ii), params, 'noroll','doFilterNoise',true,'noiseFilter', noiseMask);
        tmpSpec = getMTSpectrumStats(cl, params.fine);
        for jj = 1:numel(fieldsToKeep)
            spectra(ii).(fieldsToKeep{jj}) = tmpSpec.(fieldsToKeep{jj});
        end
        if ~featuresCached
            featureTable{ii} = extractFeatures(tmpSpec);
        end
        progressbar(ii/N);
    end
end
if ~featuresCached
    featureTable = [featureTable{:}];
    save(['tmpFeatures-' timeFlag],'DRsylls','spectra','featureTable');
end

% convert features from struct array to 2D array
fn = fieldnames(featureTable);
featureTable = cellfun(@(x) [featureTable.(x)]', fn', 'UniformOutput',false);
featureTable = [featureTable{:}];
%% calculate local distances, fixed-time version, no
%{
innerIdx = 0;
progressbar('Unwarped Distance Calcs');
for ii = 1:N-1
    for jj = ii+1:N
        innerIdx = innerIdx + 1;
        distMatrices.unwarpedLocal(innerIdx) = min(standardDistance(spectra(ii), spectra(jj), params));
        progressbar(innerIdx/nchoosek(N,2));
    end
end
save(['tmpDists-' timeFlag],'DRsylls','distMatrices');
%}
%% calculate local distances, TIME WARPED version
tic
innerIdx = 0;
progressbar('Saves','Time Warped Distance Calcs');
for ii = 1:N-1
    iLen = DRsylls(ii).stop - DRsylls(ii).start;
    for jj = ii+1:N
        jLen = DRsylls(jj).stop - DRsylls(jj).start;
        innerIdx = innerIdx + 1;       
        distMatrices.warpedLocal(innerIdx) = ...
            timeWarpedDistance(spectra(ii), spectra(jj), params) / mean([iLen, jLen]);
        progressbar([],innerIdx/nchoosek(N,2));
        
        if rem(innerIdx, floor(sqrt(nchoosek(N,2)))) == 0
            %save(['tmpDists-' timeFlag],'DRsylls','distMatrices');
            progressbar(floor(innerIdx/floor(sqrt(nchoosek(N,2)))) / ...
                floor(nchoosek(N,2)/floor(sqrt(nchoosek(N,2)))))
        end
        
    end
end
% save(['tmpDists-' timeFlag],'DRsylls','distMatrices');
            
progressbar(1);
tt=toc;
fprintf('Time warping took %0.2f s...\n', tt);
%save([dataPath 'localSimTW-' birdID '.mat'],'clustSylls','twDistM');
%% step 5: measure global distances within pairs of syllables

%seldFeaturesTable = allFeaturesTable;
%clustSylls = allDRsylls(trainIdxs);
%featureTable = allFeaturesTable(trainIdxs,:);

% start with unnormalized table of features
% step 1: normalize to z-scores
fprintf('Calculating global dissimilarity scores...\n');
zNormedFeatures = zscore(featureTable);
distMatrices.global = pdist(zNormedFeatures);
% save(['tmpDists-' timeFlag],'DRsylls','distMatrices');

%% get non-parametric (probability-rank) ordering of similarity scores
fprintf('Calculating empirical scores...\n');
scoreTypes = fieldnames(distMatrices);
for ii = 1:numel(scoreTypes)
    fld = scoreTypes{ii};
    arr = distMatrices.(fld);
    [sArr,rord] = sort(arr);
    empMatrices.(fld)(rord) = [1:numel(arr)] / numel(arr);
        
    xx = linspace(0,1,nEmpD);
    if nEmpD == numel(arr)
        yy = sArr;
    else
        yy = interp1(linspace(0,1,numel(arr)), sArr, xx);
    end
    % prepare for interp1 by removing redundant entries
    redun = [diff(yy)==0 false];
    if any(redun)
        xx = xx(~redun); % might be better to take a mean of the p-values instead of the max (as this implies)
        yy = yy(~redun);
    end
    empDistrs.(fld) = zeros(2,numel(xx));
    empDistrs.(fld)(1,:) = xx;
    empDistrs.(fld)(2,:) = yy;
%     save(['tmpEmp-'  timeFlag],'DRsylls','empMatrices', 'empDistrs');
end
%% construct co-similarity as fusion of local and global p-values
fprintf('Calculating co-dissimilarity (correlation of dissimilarities, which is a similarity score)...\n');
fusedPVals = sqrt(empMatrices.warpedLocal .* empMatrices.global);
distMatrices.cosim = pdist(squareform(fusedPVals), 'correlation');

% do the clustering - the easiest part
nClusters = 4:25;
pairLinks = linkage(distMatrices.cosim,'complete');
clustIdxs = cluster(pairLinks,'maxclust',nClusters);

% undo sorting step
scoreTypes = fieldnames(distMatrices);
for ii = 1:numel(scoreTypes)
    distMatrices.(scoreTypes{ii}) = squareform(unsort2D(squareform(distMatrices.(scoreTypes{ii})), ...
        ones(1,length(squareform(distMatrices.(scoreTypes{ii}))))));
    if isfield(empMatrices, scoreTypes{ii})
         empMatrices.(scoreTypes{ii}) = squareform(unsort2D(squareform(empMatrices.(scoreTypes{ii})), ...
             ones(1,length(squareform(distMatrices.(scoreTypes{ii}))))));
    end
end

function mat = unsort2D(mat, sI)    
    [coordsi, coordsj] = meshgrid(sI,sI);
    mat(sub2ind(size(mat),coordsi, coordsj)) = mat;
end
end
