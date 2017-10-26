% script
% test variability measures
birdID = 'Lb277'; % want plastic song with consistent motifs for testing purposes
params = processArgs(defaultParams, birdID);
rep = reportOnData(birdID, '', params, 'verbose','off');

%% todo: align based on auto-labeled syllables

%% loop for cross-correlation
nSessions = numel(rep);

% preallocate the statistics per neuron
nPerFile = zeros(1,nSessions);
for ii = 1:nSessions
    nPerFile(ii) = numel(loadSpikeData(rep(ii).spikeFiles));
end
nPartial = cumsum([0 nPerFile]);
bigCC = NaN(1,sum(nPerFile));
bigShuffled = NaN(1,sum(nPerFile));
mannP = NaN(1,sum(nPerFile));

for ii = 1:nSessions
    % load spike trains and manualMotifs;
    neuronTrains = loadSpikeData(rep(ii).spikeFiles);
    nNeurons = numel(neuronTrains);
    
    sylls = loadFromManifest(rep(ii).manifest, 'approvedSyllables');
    manualMotifs = mergeGaps(sylls, 0.15);
    %manualMtifs = loadFromManifest(rep(ii).manifest, 'manualMotifs');
    latencyAddedMotifs = addPrePost(manualMotifs, params, 'preroll', 40, 'postroll',0);
    nMotifs = numel(manualMotifs);
    
    motifTrains = cell(nNeurons, nMotifs);           
    motifAbsTrains = cell(nNeurons, nMotifs);           
    % for each neuron
    for jj = 1:nNeurons
        bNI = nPartial(ii)+jj; % big neuron index
        [~,motifTrains(jj,:)] = countSpikes(latencyAddedMotifs, neuronTrains{jj},'onset');

        contextMotifs = addPrePost(manualMotifs, params, 'preroll', 200, 'postroll', 200);
        [~,motifAbsTrains(jj,:)] = countSpikes(contextMotifs, neuronTrains{jj},'absolute');
        % scaling would be done here
        
        % shuffled entries correspond to a shuffled control
        subplot(6,1,1:2)
        [bigCC(bNI),ccList, bigShuffled(bNI), shuffledList] = getCC(latencyAddedMotifs, motifTrains(jj,:));        
        if isnan(bigCC(bNI)) || isnan(bigShuffled(bNI)), continue; end;
        
        subplot(6,1,3:4)        
        % raster plot
        for kk = 1:numel(contextMotifs), contextMotifs(kk).children = manualMotifs(kk); end        
        plotRaster(contextMotifs, motifAbsTrains(jj,:));
        
        subplot(6,1,5)
        plotPSTH(contextMotifs, motifAbsTrains(jj,:));
        xl = xlim;
        % calculate MW test for significance of CC
        ccList(isnan(ccList)) = [];
        shuffledList(isnan(shuffledList)) = [];              
        mannP(bNI) = ranksum(ccList, shuffledList);
        fprintf('MW U-test (CC, unshuffled vs. shuffled), neuron% 3d (# =% 6d) :=% 0.5f...\n', bNI, sum(~isnan(ccList)),mannP(bNI))
        
        % calculate Fano factor over time
        [fano, timeBins] = getFano(latencyAddedMotifs, motifTrains(jj,:));
        subplot(6,1,6);
        plot(timeBins, fano, 'b.', xl, [1 1], 'k-'); ylabel('Fano factor');    
        xlim(xl);
    end
end

%% plotting cross-correlation
figure(2)
bigP = ranksum(bigCC, bigShuffled);
mmin = min([bigCC bigShuffled]);
mmax = max([bigCC bigShuffled]);
markerSize = 1-log(mannP); markerSize(isinf(markerSize)) = 100;
scatter(bigCC, bigShuffled, markerSize, 'r')
hold on; plot([mmin mmax], [mmin mmax], 'k-'); hold off;

xlim([mmin mmax]);
ylim([mmin mmax]);
xlabel('real CC');
ylabel('shuffled CC');
title(sprintf('MW U-test, all neurons (CC, unshuffled vs shuffled) := %0.3f',bigP));

%% todo: test getFano