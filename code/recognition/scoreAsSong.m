function scores = scoreAsSong(songStruct, candidate, prototypeSong, syllableScorers)

% score how well an interval (CANDIDATE) matches a template song, given the exact
% positioning of sub-syllables within that song (PROTOTYPESONG)
% SYLLABLESCORERS is a set of function handles which is of the same size as
% PROTOTYPESONG
% returns a 1xN array of SCORES, where N is the number of syllables

params = defaultParams;
fs = 1/songStruct.interval;
params.fine.fs = fs;
% prerequisite: the candidate and prototype should have their starts aligned

nSongSylls = numel(syllableScorers);
candidateSylls = adjustTimeStamps(prototypeSong, candidate.start);
scores = zeros(1,nSongSylls);
for ii = 1:nSongSylls
    clip = getClip(candidateSylls(ii), songStruct);
    featureReduction = extractFeatures(getMTSpectrumStats(clip, params.fine));
    [foo, scores(ii)] = syllableScorers{ii}(featureReduction);
end

end