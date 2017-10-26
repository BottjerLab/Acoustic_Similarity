%% This script examines the syllable durations over sessions

compositeReport = reportOnData;
nBirds = numel(compositeReport);
%%
plasticScores = cell(1, nBirds);
birdID        = cell(1, nBirds);
for ii = 1:nBirds
    report = compositeReport{ii};
    birdID{ii} = strtok(report(1).sessionID,'_');

    % remove sessions that don't have approved syllables
    nSessions = numel(report);
    hasData = false(1,nSessions);
    for jj = 1:nSessions
        hasData(jj) = any(findInManifest(report(jj).manifest, {'approvedSyllables'}));
    end
    
    report = report(hasData);
    nSessions = numel(report);
    
    figure;
    [peakheight, birdAge, peakFitness] = fitBimodalDuration(birdID{ii});    
    saveCurrFigure([pwd filesep 'figures/peakfits' filesep 'dailySyllDurZero-' birdID{ii}]);
    plasticScores{ii} = [birdAge; peakheight; peakFitness]; 
end
%%
figure;
cols = jet(nBirds); 
for ii = 1:nBirds
    plot(plasticScores{ii}(1,:), plasticScores{ii}(2,:), '-s', 'Color', cols(ii,:)', 'LineWidth', 3);     
    xlabel('Age');
    ylabel('Peak height, normalized');
    hold on; 
end
for ii = 1:nBirds
    text(plasticScores{ii}(1,end)-0.1, plasticScores{ii}(2,end) + 0.3, birdID{ii})
end
%legend(birdID)
plot(xlim, 10 * [1 1],'r--', 'LineWidth', 2.5);
hold off
%ylim([0.15 0.4])


saveCurrFigure([pwd filesep 'figures/peakfits' filesep 'allBirdKSZero']);
