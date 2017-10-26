% load acceptedLabels, DRsylls, tutor, matchToTutor, MFCCdist empirical and
% regular dist empirical
% NOTE: this is obsoleted by compareTutorJuvenile.m

%% load DRsylls
birdID = 'Db113'; 
clustDir = ['data' filesep 'cluster-' birdID filesep];
dataDir  = ['data' filesep            birdID filesep];
params = defaultParams;
%% load syllables for a given age
thisAge = 64;
fprintf('Loading syllables for bird %s, age %s', birdID, thisAge)
[juvSylls, juvSessions, juvSpectra, juvFeatures] = loadAgeSylls(birdID, thisAge);
 
%% load acceptedLabels 
clusterIdxs = loadAcceptedLabels(birdID, thisAge);
num2cell(clusterIdxs); [juvSylls.type] = ans{:};

fprintf('Total syllables IDed / total: (%d/%d)\n', sum(~isnan(clusterIdxs)), numel(clusterIdxs));

%% load tutor syllables
tutorFil = [dataDir 'tutor-' birdID '.mat'];
if exist(tutorFil,'file') ~= 2
    error('Tutor file %s for bird %s not found...',tutorFil, birdID);
end
load(tutorFil);

%% load matchToTutor
matchFil = [dataDir filesep 'matchToTutor-age' num2str(thisAge) '.mat'];
if exist(tutorFil,'file') ~= 2
    error('Match file %s for bird %s not found...', matchFil, birdID);
end
load(matchFil);

%% load altClust empirical file (contains global and SAP time warped distances)
% for more information, examine recalcClustersInnerLenWeight
% imports variables 'empMats','empDistrs','distMats', 'clusterIdxs'
[normDistFil, hasNormDist] = getLatestFile([clustDir 'altClustDataAge-' num2str(thisAge) '*.mat']); 
if hasNormDist
    normDists = load(normDistFil); % contain as file
end

%% load altClust empirical file (contains global and SAP time warped distances)
% for more information, examine recalcClustersInnerLenWeight
% imports variables 'empMats','empDistrs','distMats', 'clusterIdxs'
[mfccDistFil, hasMFCCDist] = getLatestFile([clustDir 'mfccClustDataAge-' num2str(thisAge) '*.mat']); 
if hasMFCCDist
    mfccDists = load(mfccDistFil);
end

%% calculate distances between tutorSpecs/tutorFeatures and DRspectra/DRfeatures 
% that have labels only within the correct matchings
% to simplify things: write function that takes in two spectra and 
% features, and the empirical distribution for each.  but what about the
% cosimilarity?  argggggggggggggggggggggggggggggggggggggggggg

% step 0: z-normalize the juvenile syllable statistics
fn = fieldnames(juvFeatures);
juvFeaturesTabled = cellfun(@(x) [juvFeatures.(x)]', fn', 'UniformOutput',false);
juvFeaturesTabled = [juvFeaturesTabled{:}];
[zNormedTable, featMu, featSigma] = zscore(juvFeaturesTabled);

% step 1: get local and global distances to each other syllable
nTutor = numel(tutorSylls);
nJuv = numel(juvSylls);
localDist = struct('raw', zeros(nJuv, nTutor), ...% initialization
                    'emp', zeros(nJuv, nTutor),...
                    'co', zeros(nJuv, nTutor));
globalDist = localDist; fusedDist = rmfield(localDist, 'emp'); 
 localEmp = squareform(normDists.empMats.warpedLocal);
globalEmp = squareform(normDists.empMats.global);
 fusedEmp = squareform(normDists.distMats.cosim);
progressbar(0,0);
for ii = 1:nTutor
    tutorFeats = ((cellfun(@(x) tSFeats(ii).(x), fn') - featMu) ./ featSigma);
    tutorFeats(featSigma == 0) = 0;
    for jj = 1:nJuv
        progressbar([],jj/nJuv);
        localDist.raw(jj, ii) = timeWarpedDistance(tSSpecs(ii), juvSpectra(jj));
    end
    globalDist.raw(:, ii) = pdist2(tutorFeats, zNormedTable, 'euclidean');    
    localDist.emp(:, ii) = interp1(normDists.empDistrs.warpedLocal(2,:), normDists.empDistrs.warpedLocal(1,:), ...
        localDist.raw(:, ii));
    
    % no great way to extrapolate if distance is too far, so just divide by
    % the max for now
    inExcess = find(localDist.raw(:,ii) > normDists.empDistrs.warpedLocal(2,end)); 
    if numel(inExcess) > 0
        localDist.emp(inExcess,ii) = localDist.raw(inExcess,ii) / normDists.empDistrs.warpedLocal(2,end);        
    end
    globalDist.emp(:, ii) = interp1(normDists.empDistrs.global(2,:), normDists.empDistrs.global(1,:), ...
        globalDist.raw(:, ii));
    % no great way to extrapolate if distance is too far, so just divide by
    % the max for now
    inExcess = find(globalDist.raw(:,ii) > normDists.empDistrs.global(2,end)); 
    if numel(inExcess) > 0
        globalDist.emp(inExcess,ii) = globalDist.raw(inExcess,ii) / normDists.empDistrs.global(2,end);        
    end
    if any(isnan(globalDist.emp(:,ii))), keyboard; end;
    localDist.co(:, ii) = pdist2( localDist.emp(:, ii)', localEmp, 'correlation');
    globalDist.co(:, ii) = pdist2(globalDist.emp(:, ii)', globalEmp, 'correlation');
    progressbar(ii/nTutor, 0);
end

fusedDist.raw = sqrt(localDist.co .* globalDist.co);
fusedDist.co  = pdist2(fusedDist.raw', fusedEmp, 'correlation')'; % each tutor syllable against each juvenile syllable

%% for each juvenile syllable, average only the matching tutor syllable distances
distanceScore = NaN(1,nJuv);
for ii = 1:nJuv
   if isnan(juvSylls(ii).type), continue; end
   matchingSylls = strcmp(typeMatch(juvSylls(ii).type), {tutorSylls.type});
   distanceScore(ii) = mean(fusedDist.co(ii, matchingSylls));
end

%% plot/order the syllables by similarity so that they make sense?
% plot in the same window style as checkStereotypy.

