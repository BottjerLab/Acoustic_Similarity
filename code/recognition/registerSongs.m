function registerSongs(songStruct, confirmedSongs, songTemplate, syllFxns, noiseProfile)   
% TODO: figure out if this function is still worthwhile....
if nargin < 5
    noiseProfile = [];
end
params = defaultParams;
    fs = 1/songStruct.interval;
    nSongs = numel(confirmedSongs);
    
    params.preroll = 0; params.postroll = 100; params.fs=fs;
    pivotalSyllable = 3;
    nStandard = numel(syllFxns);
    for ii = 1:nSongs        
        cl = getClipAndProcess(confirmedSongs(ii), songStruct, params, ...
            'doFilterNoise', nargin == 5, noiseFilter', noiseProfile);
        
        params.fine.fs = fs;
        spec = getMTSpectrumStats(cl,params.fine);
        sylls = segmentSyllables(spec);
        
        %% try to match key syllables
        featureReduction = extractFeatures(spec, sylls);
        scores = zeros(numel(sylls), nStandard);
        for kk = 1:nStandard
             [foo, scores(:,kk)] = syllFxns{kk}(featureReduction);
        end
               
        % TODO: still need to map each song syllable to a particular, i.e.
        % give each sylabble its correct label to the normal song
        % training a classifier doesn't seem to be enough
        
        %% plot
        subplot(3,1,1);
        plotWaveform(cl, fs); plotAreaMarks(songTemplate);
        title(sprintf('song #%d, standardized bounds',ii));
        subplot(3,1,2);
        plotWaveform(cl, fs); plotAreaMarks(sylls);
        title(sprintf('song #%d, detected bounds',ii));
        subplot(3,1,3);
        plotDerivGram(spec); plotMarks(sylls,8000);
        title(sprintf('spectrogram'));
        %saveCurrFigure(sprintf('figures/alignment-%02d',ii));
    end
end
