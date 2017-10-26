% flow of scripts to plot group correlations, starting from syllables
rep = reportOnData;
%% get birds and age of session recordings that contain song
nBirds = numel(rep);
birdIDs = cell(1, nBirds);
ages = cell(1, nBirds);
for ii = 1:nBirds
    birdIDs{ii} = strtok(rep{ii}(1).sessionID, '_');
    agesTmp = getAgeOfSession({rep{ii}.sessionID});
    ages{ii} = unique(agesTmp(~isnan(agesTmp)));
end

%% the commands that produce the clustering
for ii = 1:nBirds
    for jj = 1:numel(ages{ii})
        fullClusterByDay(birdIDs{ii}, ages{ii}(jj));
        % TODO: use an empirical distribution with distances over all sessions        
    end
end
%% optional - edit labels 
% this is an interactive segment but should be done once
% if this is done, clusterIdxs in the corresponding altClustDataAge- file
% will be a structure
for ii = 1:nBirds
    for jj = 1:numel(ages{ii})
        browseAndAccept(birdIDs{ii})
    end
end
%% split the labels from single ages to sessions - this script handles them all
splitAcceptedLabels;
%% get the average firing rates/response strengths for each cluster
writeNeuronStats({'Db113','Dg138','Gy217','Gy242', 'Lb189','Lb277','R204','R247','R288','Y231'});
%% measure amounts of matching and stereotypy - these are independent of the neuronal activity
for ii = 1:nBirds
    for jj = 1:numel(ages{ii})
        fprintf('Analyzing bird %s, age %d...\n', birdIDs{ii}, ages{ii}(jj));
        compareTutorJuvenile(birdIDs{ii}, ages{ii}(jj));
        saveStereotypy(birdIDs{ii},ages{ii}(jj));
    end
end

%% get neural correlation to both matching signals and stereotypy

aNCD = correlateDistanceToFiring(birdIDs);
%%
plotNeuronCorrData(aNCD,[],'saveplot',true);
%%
sessionIDs = {aNCD.sessionID};
manyBirdIDs = strtok(sessionIDs, '_');
[uSessions, ~, rIdxSession] = unique(sessionIDs);  % index through ages can go back to sesions
uAges = getAgeOfSession(uSessions);
sessionAges = zeros(size(sessionIDs));
for ii = 1:numel(uAges)
    sessionAges(rIdxSession == ii) = uAges(ii);
end
[sessionQ   , allSubj]  = getClusterQuality(manyBirdIDs , sessionAges, [aNCD.syllID]);
[sessionObjQ, allObj ]  = getClusterQuality(manyBirdIDs , sessionAges, [aNCD.syllID], true);
foo = num2cell(sessionQ   ); [aNCD.clusterQ   ] = foo{:};
foo = num2cell(sessionObjQ); [aNCD.clusterObjQ] = foo{:};

isCore        = [aNCD.isCore];
isMUA         = [aNCD.isMUA];
%isPlastic     = [allNeuronCorrData.isPlastic];
isSignificant = [aNCD.sigResponse];
nSylls        = [aNCD.nSylls];
%isExcited     = [allNeuronCorrData.isExcited];

isSubjGood    = [aNCD.clusterQ]  < 2.5;

% criteria for unit inclusion
isPresel = isSignificant & ~isMUA & nSylls >= 12 & isSubjGood;

approvedNCD = aNCD(isPresel);
isApprovedCore = isCore(isPresel);
%%
distType = 'central';
nearZ = [approvedNCD.([distType '_nearZ'])];
farZ  = [approvedNCD.([distType '_farZ'])];

cols = [1 0 0; 0.5 0.5 0.5];
plot(nearZ( isApprovedCore), farZ( isApprovedCore), '.', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 16); 
hold on;
plot(nearZ(~isApprovedCore), farZ(~isApprovedCore), '.', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', 16);
plot([-5 5], [-5 5], 'k-');
xlim([-5 5]); ylim([-5 5]);
xlabel('RS z-score for quartile near to tutor');
ylabel('RS z-score for quartile far from tutor');
hold off;
set(gca,'Box', 'off')
set(gcf,'Color', [1 1 1]);

imFile = sprintf('figures/paper/%sNearFarZScoreDistribution.pdf',distType);
fprintf('Writing image to %s...\n',imFile);
export_fig(imFile);

%%
distType = 'inter';
nearZ = [approvedNCD.([distType '_nearZ'])];
farZ  = [approvedNCD.([distType '_farZ'])];

cols = [1 0 0; 0.5 0.5 0.5];
plot(nearZ( isApprovedCore), farZ( isApprovedCore), '.', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 16); 
hold on;
plot(nearZ(~isApprovedCore), farZ(~isApprovedCore), '.', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', 16);
plot([-5 5], [-5 5], 'k-');
xlim([-5 5]); ylim([-5 5]);
xlabel('RS z-score for quartile near to cluster center');
ylabel('RS z-score for quartile far from cluster center');
hold off;
set(gca,'Box', 'off')
set(gcf,'Color', [1 1 1]);

imFile = sprintf('figures/paper/%sNearFarZScoreDistribution.pdf',distType);
fprintf('Writing image to %s...\n',imFile);
export_fig(imFile);
%% optional: plot all members of clusters
for ii = 1:nBirds
    for jj = 1:numel(ages{ii})
        drawClustersUpdated(birdIDs{ii}, ages{ii}(jj));        
    end
end
