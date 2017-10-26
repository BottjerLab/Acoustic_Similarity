%todo: load allNeuronCorrData and call plotTutorJuvenileComparison for
%each neuron/cluster pair that matches the inclusion criteria
 load('data/allNeuronCorrelations.mat');

%% here we load the cluster quality
% get birds and ages first
sessionIDs = {allNeuronCorrData.sessionID};
birdIDs = strtok(sessionIDs, '_');
[uSessions, ~, rIdxSession] = unique(sessionIDs);  % index through ages can go back to sesions
uAges = getAgeOfSession(uSessions);
sessionAges = zeros(size(sessionIDs));
for ii = 1:numel(uAges)
    sessionAges(rIdxSession == ii) = uAges(ii);
end
[sessionQ   , allSubj]  = getClusterQuality(birdIDs, sessionAges, [allNeuronCorrData.syllID]);
[sessionObjQ, allObj ]  = getClusterQuality(birdIDs, sessionAges, [allNeuronCorrData.syllID], true);
foo = num2cell(sessionQ   ); [allNeuronCorrData.clusterQ   ] = foo{:};
foo = num2cell(sessionObjQ); [allNeuronCorrData.clusterObjQ] = foo{:};

% flags
isCore        = [allNeuronCorrData.isCore];
isMUA         = [allNeuronCorrData.isMUA];
isPlastic     = [allNeuronCorrData.isPlastic];
isSignificant = [allNeuronCorrData.sigResponse];
nSylls        = [allNeuronCorrData.nSylls];
%isSignificant = true(1,numel(allNeuronCorrData));
isExcited     = [allNeuronCorrData.isExcited];
isSubjGood    = [allNeuronCorrData.clusterQ]  < 2.5;
isObjGood     = [allNeuronCorrData.clusterObjQ] < 1;

% criteria for neuron inclusion
% must have at least three syllables per quartile
isPresel = isSignificant & ~isMUA & nSylls >= 12 & isSubjGood;

preselData = allNeuronCorrData(isPresel);
for ii = 1:numel(preselData)
    fprintf('Showing correlation for session %s, cluster %d, neuron %d...\n', ...
        preselData(ii).sessionID, ...
        preselData(ii).syllID, preselData(ii).unitNum);
    plotTutorJuvenileComparisonJenny(preselData(ii).sessionID, preselData(ii).unitNum, ...
        preselData(ii).syllID, preselData(ii).isCore, [], 'plot',true,'saveplot', true);
end
%% related: lump syllables into each neuron, on the basis of which
% either which activity is strongest or which correlation is strongest.
[uSessions, ~, sessIdx] = unique({allNeuronCorrData.sessionID});
mostResponsiveData = [];
for ii = 1:numel(uSessions)
    thisSessionData = allNeuronCorrData(sessIdx == ii);
    neurIDs = [thisSessionData.unitNum];
    for jj = 1:max(neurIDs)        
        if any(neurIDs == jj)
            thisUnitData = thisSessionData(neurIDs == jj);
            [~,mostResponsiveIdx] = max([thisUnitData.avgResponse]);
            mostResponsiveData = ...
                [mostResponsiveData thisUnitData(mostResponsiveIdx)];            
        end
    end
end
%%
plotNeuronCorrData(mostResponsiveData,[],'saveplot',true);
% related: also plot rasters in order of similarity

% related: check whether clusters are bimodal
% todo: save plot on subjective / objective clusters

%% todo: rewrite plotting part of allBirdsFeatureAnalysis
allBirdsFeatureAnalysis;
