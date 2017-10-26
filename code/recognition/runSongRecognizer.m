function scores = runSongRecognizer(songStruct, potentialSyllables, prototypeSong, syllableScorers, pivotSyllable)

% 
% syllableScorers is a set of function handles which are trained to score
% how well a syllable matches a set of features
pivotSyllable = 3; % this should be the loudest syllable, the one most consistently segmented
introSyll = 1:2;


%% computation
potentialOnsets = [potentialSyllables.start];
potentialIdxOnsets = [potentialSyllables.idxStart];

nSongSylls = numel(prototypeSong);
scores = zeros(numel(potentialOnsets),nSongSylls);
progressbar('scoring possible syllables')
nPot = numel(potentialOnsets);
for ii = 1:nPot
    % define the song according to the prototype and see how well it
    % fits...
    trialSong = prototypeSong;
    foo = num2cell([trialSong.start]    + potentialOnsets(ii) - prototypeSong(pivotSyllable).start); [trialSong.start] = foo{:};
    foo = num2cell([trialSong.stop ]    + potentialOnsets(ii) - prototypeSong(pivotSyllable).start); [trialSong.stop ] = foo{:};
    foo = num2cell([trialSong.idxStart] + potentialIdxOnsets(ii) - prototypeSong(pivotSyllable).idxStart); [trialSong.idxStart] = foo{:};
    foo = num2cell([trialSong.idxStop ] + potentialIdxOnsets(ii) - prototypeSong(pivotSyllable).idxStart); [trialSong.idxStop ] = foo{:};

    scores(ii,:) = scoreAsSong(songStruct, trialSong,prototypeSong, syllableScorers);
    progressbar(ii/nPot);
end

