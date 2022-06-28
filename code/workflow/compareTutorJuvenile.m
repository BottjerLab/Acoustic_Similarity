function [distToTut,bestDistScore,bestTutMatch] = compareTutorJuvenile(juvSylls,juvSpectra,juvFeatures,juvSyllClusters,tutSylls,tutSpectra,tutFeatures,params,varargin)
% get tutor song syllables distance to juvenile syllables, taking into
% account both local and global distances.  
% 
% saves distances to tutorSyllableCompare-<session>.mat
if nargin < 8
    params = defaultParams;
end
params = processArgs(params, varargin{:});

%% load altClust empirical file (contains global and SAP time warped distances)
% for more information, examine recalcClustersInnerLenWeight
% imports variables 'empMats','empDistrs','distMats', 'clusterIdxs'
%{
[mfccDistFil, hasMFCCDist] = getLatestFile([clustDir 'mfccClustDataAge-' num2str(age) '*.mat']); 
if hasMFCCDist
    mfccDists = load(mfccDistFil);
end
%}
%% calculate distances between tutorSpecs/tutorFeatures and DRspectra/DRfeatures 
% that have labels only within the correct matchings
% goes through the same steps as in DRcluster

% step 0: z-normalize the juvenile syllable statistics
fn = fieldnames(juvFeatures);
juvFeaturesTabled = cellfun(@(x) [juvFeatures.(x)]', fn', 'UniformOutput',false);
juvFeaturesTabled = [juvFeaturesTabled{:}];
[zNormedTable, featMu, featSigma] = zscore(juvFeaturesTabled);

% step 1: get local and global distances to each other syllable
nTutor = numel(tutSylls);
nJuv = numel(juvSylls);
localDist = struct('raw', zeros(nJuv, nTutor), ...% initialization
                    'emp', zeros(nJuv, nTutor),...
                    'co', zeros(nJuv, nTutor));
globalDist = localDist; fusedDist = rmfield(localDist, 'emp'); 
 localEmp = squareform(juvSyllClusters.empMats.warpedLocal);
globalEmp = squareform(juvSyllClusters.empMats.global);
 fusedEmp = squareform(juvSyllClusters.distMats.cosim);
progressbar('tutor syllables','juvenile syllables');
for ii = 1:nTutor
    tutorFeats = ((cellfun(@(x) tutFeatures(ii).(x), fn') - featMu) ./ featSigma);
    tutorFeats(featSigma == 0) = 0;
    lenTut = tutSylls(ii).stop - tutSylls(ii).start;
    for jj = 1:nJuv
        progressbar([],jj/nJuv);
        lenJuv = juvSylls(jj).stop - juvSylls(jj).start;
        % local distance is normalized by the average length of syllables
        localDist.raw(jj, ii) = timeWarpedDistance(tutSpectra(ii), juvSpectra(jj)) / mean([lenTut, lenJuv]);
    end
    globalDist.raw(:, ii) = pdist2(tutorFeats, zNormedTable, 'euclidean');    
    localDist.emp(:, ii) = interp1(juvSyllClusters.empDistrs.warpedLocal(2,:), juvSyllClusters.empDistrs.warpedLocal(1,:), ...
        localDist.raw(:, ii));
    
    % no great way to extrapolate if distance is too far, so just divide by
    % the max for now
    inExcess = find(localDist.raw(:,ii) > juvSyllClusters.empDistrs.warpedLocal(2,end)); 
    if numel(inExcess) > 0
        localDist.emp(inExcess,ii) = localDist.raw(inExcess,ii) / juvSyllClusters.empDistrs.warpedLocal(2,end);        
    end
    globalDist.emp(:, ii) = interp1(juvSyllClusters.empDistrs.global(2,:), juvSyllClusters.empDistrs.global(1,:), ...
        globalDist.raw(:, ii));
    % no great way to extrapolate if distance is too far, so just divide by
    % the max for now
    inExcess = find(globalDist.raw(:,ii) > juvSyllClusters.empDistrs.global(2,end)); 
    if numel(inExcess) > 0
        globalDist.emp(inExcess,ii) = globalDist.raw(inExcess,ii) / juvSyllClusters.empDistrs.global(2,end);        
    end
    if any(isnan(globalDist.emp(:,ii))), keyboard; end;
    localDist.co(:, ii) = pdist2( localDist.emp(:, ii)', localEmp, 'correlation');
    globalDist.co(:, ii) = pdist2(globalDist.emp(:, ii)', globalEmp, 'correlation');
    progressbar(ii/nTutor, 0);
end

fusedDist.raw = sqrt(localDist.co .* globalDist.co);
fusedDist.co  = pdist2(fusedDist.raw', fusedEmp, 'correlation')'; % each tutor syllable against each juvenile syllable

%% for each juvenile syllable, average the distance for each type of the tutor syllable
% there are multiple instances of the tutor syllable
[tutorTypes,~,matchesTutor] = unique({tutSylls.type});
nTutorTypes = numel(tutorTypes);
distToTut = NaN(nTutorTypes, nJuv);

for ii = 1:nTutorTypes   
   distToTut(ii, :) = mean(fusedDist.co(:, ii == matchesTutor),2)';
end

%% for each juvenile syllable, get different scores to tutor syllable:
% bestDistanceScore <- distance from individual juv syllable to closest tutor syll 
% distanceToConsensus <- distance from individual juv syllable to (tutor syll closest to all
% syllables of that juv syllable's cluster) 
% distanceToCentral <- distance from individual juv syllable to (tutor syll
% closest to the most central syllable of that juv syllable's cluster)
% distanceToHumanMatch <- distance to the matched tutor as described by the
% matchToTutor file

% best distance to any tutor syllable
[bestDistScore, bestInd] = min(distToTut, [],1);
bestTutMatch = tutorTypes(bestInd);

% distanceToConsensus  = NaN(1, nJuv);
% distanceToCentral    = NaN(1, nJuv);
% distanceToHumanMatch = NaN(1, nJuv);

% consensusMatch = NaN(1, nJuv);
% centralMatch   = NaN(1, nJuv);
% humanMatch     = NaN(1, nJuv);

% try to load human matches for a type
% typeMatch = [];
% matchFile = [dataDir 'matchToTutor-age' num2str(age) '.mat'];
% if exist(matchFile, 'file') == 2
%     load(matchFile); % loading typeMatch, cell array of chars;
%     % this is generated by assignClusterToTutor
% end
% if ~isempty(clusterIdxs)
%     % distanceToCluster
%     % first find the closest syllable to each cluster, by summing distances 
%     nClusters = nanmax(clusterIdxs);
%     centralOrder = findMostCentral(normDists.distMats.cosim, clusterIdxs);
%     
%     for ii = 1:nClusters
%         inThisCluster = find(clusterIdxs == ii);
%         thisCentralIdx = centralOrder{ii}(1);
%         
%         distancesToThisCluster = sum(distanceToTutor(:,inThisCluster), 2);
%         [~,closestByConsensus] = min(distancesToThisCluster); 
%         distanceToConsensus(inThisCluster)  = distanceToTutor(closestByConsensus, inThisCluster);
%         consensusMatch(inThisCluster) = closestByConsensus;
%         
%         [~,closestByCentral] = min(distanceToTutor(:,thisCentralIdx));        
%         distanceToCentral(inThisCluster)    = distanceToTutor(closestByCentral, inThisCluster);
%         centralMatch(inThisCluster) = closestByCentral;
%         
%         if ~isempty(typeMatch)
%             closestByHumanMatch = typeMatch{ii} - 'a' + 1;
%             distanceToHumanMatch(inThisCluster) = distanceToTutor(closestByHumanMatch, inThisCluster);
%             humanMatch(inThisCluster) = closestByHumanMatch;
%         end
%     end
%     % get the closest syllable to every cluster       
% end

%% split scores among sessions

%     distToConsensus  = distanceToConsensus(isInSession);
%     distToCentral    = distanceToCentral(isInSession);
%     distToHumanMatch = distanceToHumanMatch(isInSession);
%
%     consMatch = consensusMatch(isInSession);
%     centMatch = centralMatch(isInSession);
%     humMatch  = humanMatch(isInSession);
%     %the juvenile dimension is the last one in theis savefile
% saveFile = [dataDir 'tutorSyllableCompare-' thisSession '.mat'];
% fprintf('Saving %s session data to %s...\n', thisSession, saveFile);
%     save(saveFile, 'distToTutor', 'bestDistScore', 'bestTutMatch',...
%         'distToConsensus','distToCentral','distToHumanMatch',...
%         'consMatch', 'centMatch', 'humMatch');
% save(saveFile, 'distToTutor', 'bestDistScore', 'bestTutMatch');
