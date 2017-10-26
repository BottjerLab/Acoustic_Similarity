%function motifLengthVariation
% script that plots simple length-related metrics of subsong/plastic song;
rep = reportOnData;
birdStats = cell(1,numel(rep));
%%
for ii = 1:numel(rep)
    birdID = strtok(rep{ii}(1).sessionID,'_');
    [vM, mSR, cvM, ages] = motifVariation(birdID);
    birdStats{ii} = struct('ages', ages, 'varMotifLengths', vM,...
        'motifSylRatios', mSR);
end
%% plot based on birdStats
nBirds = numel(birdStats);
cols = jet(nBirds);
birdID = cell(1,nBirds);
for ii = 1:numel(birdStats)
    birdID{ii} = strtok(rep{ii}(1).sessionID, '_');
end
%%
for ii = 1:numel(birdStats)
    subplot(2,1,1)
    plot(birdStats{ii}.ages, birdStats{ii}.varMotifLengths, '-s', 'Color', cols(ii,:)', 'LineWidth', 3);
    
    hold on;
    subplot(2,1,2)
    plot(birdStats{ii}.ages, birdStats{ii}.motifSylRatios, '-s', 'Color', cols(ii,:)', 'LineWidth', 3);
    hold on;
end

for ii = 1:numel(birdStats)
    subplot(2,1,1)
    text(birdStats{ii}.ages(end)-0.1, birdStats{ii}.varMotifLengths(end) - 0.1, birdID{ii});
    subplot(2,1,2)
    text(birdStats{ii}.ages(end)-0.1, birdStats{ii}.motifSylRatios(end) - 0.1, birdID{ii});
end
subplot(2,1,1);
xlabel('Bird Age');
ylabel('Standard deviation in motif lengths (s)');
subplot(2,1,2);
xlabel('Bird Age');
ylabel('Ratio of late/early syllables in motifs');

saveFile = ['figures' filesep 'motifSubsongMeasures.jpg'];
saveCurrFigure(saveFile);
