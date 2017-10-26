function scores = scoreIndividualSyllables(songStruct, syllables, syllableScorers)

% score how well a set of intervals (SYLLABLES) matches template syllables
% SYLLABLESCORERS is a set of N function handles that evaluate
% featureReduction
% returns an NxN array of SCORES, where M is the number of syllables

params = defaultParams;
fs = 1/songStruct.interval;
params.fine.fs = fs;

nSylls = numel(syllables);
nSongSylls = numel(syllableScorers);

scores = zeros(nSylls,nSongSylls);
progressbar('scoring possible syllables')
for ii = 1:nSylls
    clip = getClip(syllables(ii), songStruct);
    featureReduction = extractFeatures(getMTSpectrumStats(clip, params.fine));
    for jj = 1:nSongSylls
        [foo, scores(ii,jj)] = syllableScorers{jj}(featureReduction);
        progressbar(((ii-1)*nSongSylls + jj)/(nSylls * nSongSylls));
    end
end

end