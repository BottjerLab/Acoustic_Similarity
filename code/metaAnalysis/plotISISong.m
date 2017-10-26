function plotISISong(session)

% TODO: REWRITE THIS FUNCTION TO FIT ITS NAME

%session = 'Lb277_3_20_1'; % TODO: batch preparation
birdID = strtok(session, '_');
sessionRecord = reportOnData(birdID, session);
manifest = sessionRecord.manifest;

% don't need values, just the metadata of the song structure
songStruct = loadFromManifest(manifest,'metaStruct');

% get motifs
motifs = loadFromManifest(manifest, 'manualMotifs');

% get spikes
spikes = loadSpikeData(sessionRecord.spikeFiles);

% do analysis
% mark whether or not the spike exists within the song boundaries or not
for ii = 1:numel(spikes) % loop over all the units/neurons
    [motifISIs, nonMotifISIs] = getISIs(spikes{ii}, motifs);
    if isempty(motifISIs) || isempty(nonMotifISIs) %Jenny
        warning ('there are no data so skipping this neuron')
        continue
    end
    
    nBins = 200; %Jenny
    maxLogMotifISI = 0; %Jenny
    minLogMotifISI = []; %Jenny
    maxLogNonMotifISI = 0; %Jenny
    minLogNonMotifISI = []; %Jenny
    
    if maxLogMotifISI == 0 %Jenny
        maxLogMotifISI = max(real(log10(motifISIs)))+1; %Jenny     
    end %Jenny

    if isempty(minLogMotifISI) %Jenny
        minLogMotifISI = floor(min(real(log10(motifISIs)))); %Jenny
    end %Jenny
    
    if maxLogNonMotifISI == 0 %Jenny
        maxLogNonMotifISI = max(real(log10(nonMotifISIs)))+1; %Jenny     
    end %Jenny

    if isempty(minLogNonMotifISI) %Jenny
        minLogNonMotifISI = floor(min(real(log10(nonMotifISIs)))); %Jenny
    end %Jenny
        
    logMotifISIs    = log10(   motifISIs); 
    logNonMotifISIs = log10(nonMotifISIs);
    melogMotif = median(logMotifISIs); %Jenny
    melogNonMotif = median(logNonMotifISIs); %Jenny
    mlogMotif = mean(logMotifISIs); %Jenny
    mlogNonMotif = mean(logNonMotifISIs); %Jenny
    
    binsUsedMotif = logspace(minLogMotifISI,maxLogMotifISI,nBins); %Jenny
    binsUsedNonMotif = logspace(minLogNonMotifISI,maxLogNonMotifISI,nBins); %Jenny
    
    H = ndhist(logMotifISIs', nBins, minLogMotifISI, maxLogMotifISI); %Jenny
    I = ndhist(logNonMotifISIs', nBins, minLogNonMotifISI, maxLogNonMotifISI); %Jenny
    
    %binEdges = [0:0.025:0.5 0.55:0.05:1];Jenny commented out
    
    %minLogISI = min([logMotifISIs; logNonMotifISIs]); Jenny commented out
    %maxLogISI = max([logMotifISIs; logNonMotifISIs]); Jenny commented out
    
    %logStepSize = 0.1; % you can change this if you want; Jenny commented out
    %logBinEdges = minLogISI:logStepSize:maxLogISI; Jenny commented out
    
    figure(ii);
    subplot(211); 
    %bar(logBinEdges, histc(logNonMotifISIs, logBinEdges),'histc'); Jenny
    %commented out
    plot(binsUsedNonMotif, smooth(I,4)); %Jenny
    title(sprintf('ISIs outside song (%d)', numel(logNonMotifISIs)));
    %set(gca, 'XTickLabels', 10.^(logBinEdges)); Jenny commented out
    set(gca, 'XScale', 'log', 'XLim', [10^minLogNonMotifISI 10^maxLogNonMotifISI]); %Jenny
    set(gca, 'YTick', max(I));%Jenny
    key(1) = {['median ISI: ',num2str(melogNonMotif)]}; %Jenny
    key(2) = {['mean ISI: ',num2str(mlogNonMotif)]}; %Jenny
    text(min(binsUsedNonMotif),max(I)/2,key); %Jenny
    subplot(212); 
    %bar(logBinEdges, histc(   logMotifISIs, logBinEdges),'histc'); Jenny
    %commented out
    plot(binsUsedMotif, smooth(H,4)); %Jenny
    title(sprintf('ISIs inside song (%d)' , numel(   logMotifISIs)));
    %set(gca, 'XTickLabels', 10.^(logBinEdges)); Jenny commented out
    set(gca, 'XScale', 'log', 'XLim', [10^minLogNonMotifISI 10^maxLogNonMotifISI]); %Jenny
    set(gca, 'YTick', max(H));%Jenny
    key(1) = {['median ISI: ',num2str(melogMotif)]}; %Jenny
    key(2) = {['mean ISI: ',num2str(mlogMotif)]}; %Jenny
    text(min(binsUsedMotif),max(H)/2,key); %Jenny
end
