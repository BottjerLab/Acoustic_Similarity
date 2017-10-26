%% compileStereotypy - plot figure for summary of response strengths compared to central/peripheral syllables

% include only significant neurons for ?
onlySignificant = false;

reps = reportOnData;
neuronData = [];
progressbar('birds','ages');
nBirds = numel(reps);
for ii = 1:nBirds
    birdID = strtok(reps{ii}(1).sessionID, '_');
    sessionRecords = reps{ii};
    sessionIDs = {sessionRecords.sessionID};
    ages = getAgeOfSession(sessionIDs);
    
    uAges = unique(ages(~isnan(ages)));
    nAges = numel(uAges);
    for jj = 1:numel(uAges)
        try
            tmpData = variableSyllTypes(birdID, uAges(jj));
        catch err
            if strcmp(err.identifier, 'variableSyllTypes:labelsNotFound') || ...
                    strcmp(err.identifier, 'variableSyllTypes:clusterDistssNotFound');
                fprintf('%s\n', err.message);
                continue;
            else
                rethrow(err)
            end
        end
        [tmpData.age] = deal(uAges(jj));
        [tmpData.birdID] = deal(birdID);
        neuronData = [neuronData; tmpData];        
        progressbar([],jj/nAges);
    end
    % if we only want significant neurons, go through
    progressbar(ii/nBirds,0);
end

if onlySignificant
	neuronData([neuronData.allProb]) = [];
end
        
%%
    
%% every single neuron in interaction plot
%{
isPlastic = logical([neuronData.isPlastic]);
isCore    = logical([neuronData.isCore]);
isSUA     =        ~[neuronData.isMUA];
cR = [neuronData.centralRate];
pR = [neuronData.peripheralRate];
for ii = 1:numel(neuronData)
    if isCore(ii),    colType = 'c'; else colType = 'r'; end
    if isPlastic(ii), linType = '-.'; else linType = '-'; end
    plot([1 2], [cR(ii), pR(ii)], [colType linType]);
    hold on;    
end
xlim([0.5 2.5])
set(gca,'XTick', [1 2]);
set(gca,'XTickLabel', {'central', 'peripheral'});
xlabel('syllable type');
ylabel('Firing Rate (Hz)');
%}
%% bar graph: core, SUA
% indices on MUA/SUA/all, shell/core/all, subsong/plastic/all
isPlastic = logical([neuronData.isPlastic]);
isCore    = logical([neuronData.isCore]);
isSUA     =        ~[neuronData.isMUA];
centralMeans = zeros(3,3,3);
periphMeans = zeros(3,3,3);
centralStds = zeros(3,3,3);
periphStds = zeros(3,3,3);
sel = true(size(neuronData))';
for iUnit = 1:3
    if iUnit < 3, sSel = sel & isSUA == (iUnit-1); 
    else sSel = sel; end
    for iPlastic = 1:3
        if iPlastic < 3, s2Sel = sSel & isPlastic == (iPlastic-1); 
        else s2Sel = sSel; end
        for iCore = 1:3
            if iCore < 3, s3Sel = s2Sel & isCore == (iCore - 1); 
            else s3Sel = s2Sel; end
            cData = [neuronData(s3Sel).centralNormRS];
            pData = [neuronData(s3Sel).periphNormRS];
            NData = sum(s3Sel);
            centralMeans(iUnit, iPlastic, iCore) = mean(cData);
             periphMeans(iUnit, iPlastic, iCore) = mean(pData);
            centralStds(iUnit, iPlastic, iCore) = std(cData) / sqrt(NData - 1);
             periphStds(iUnit, iPlastic, iCore) = std(pData) / sqrt(NData - 1);            
        end
    end
end

unitLabels = {'multi-unit', 'single-unit', 'all units'};
bigLabels = {'subsong', 'plastic song', 'all song'};
coreLabels = {'core', 'shell', 'both'};
    
for iUnit = 1:3
    bT = []; sT = [];
    labels = {};
    for iSub = 1:3 % build the data tables
        bT = [bT; [centralMeans(iUnit,:,iSub)', periphMeans(iUnit,:,iSub)']];
        sT = [sT; [centralStds(iUnit,:,iSub)', periphStds(iUnit,:,iSub)']];
        labels = [labels, strcat([bigLabels{iSub} ', '], coreLabels)];
    end
    subplot(3,1,iUnit)
    plotBarError(bT, sT, []);
    set(gca,'XTickLabel', labels);
    xlabel(unitLabels{iUnit});
    ylabel('standard RS');
end
legend({'central', 'peripheral'})
if ~onlySignificant
    export_fig(['figures/centralRSAcrossAllTypes.jpg'])
else
    export_fig(['figures/centralRSAcrossAllTypesSigNeurons.jpg'])
end
