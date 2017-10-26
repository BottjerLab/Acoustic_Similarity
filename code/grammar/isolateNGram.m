function [lenSegment, varSegment] = isolateNGram(events, ngram, gramIdxs, N)
% look at the consistency of the length of each segment
% tests the set of 'events' labeled as n-grams as ngrams, using 
lenSegment = zeros(1,numel(ngram));
varSegment = zeros(1,numel(ngram));
for ii = numel(ngram):-1:1
    candidate = ngram{ii};
    % only test unique ones, non spaced ones??
%    if numel(unique(candidate)) ~= numel(candidate) || any(candidate == ' ')        
%        continue;
%    end
    % find all instances and save them
    candLocs = find(gramIdxs == ii);
    
    % can't check variance if there's only one example
    if numel(candLocs) == 1, continue; end;
    
    % test left and right extensibility
    candSongLE = initEvents;
    candSong = initEvents;
    candSongRE = initEvents;

    for jj = 1:numel(candLocs)
        if candLocs(jj) > 1
            candSongLE(jj) = bookendedClip(events((candLocs(jj)-1):(candLocs(jj)+N-1)));     
        end
        candSong(jj) = bookendedClip(events(candLocs(jj):(candLocs(jj)+N-1)));
        if candLocs(jj) + N <= numel(events)
            candSongRE(jj) = bookendedClip(events(candLocs(jj):(candLocs(jj)+N)));
        end
    end
    lens = [candSong.stop] - [candSong.start];
    lensLE = [candSongLE.stop] - [candSongLE.start];
    lensRE = [candSongRE.stop] - [candSongRE.start];
    lenSegment(ii) = nanmean(lens);
    varSegment(ii) = nanstd(lens);
    
    fprintf('length (in s) of ''%s'' (count = %d): left extend %f +/- %f, of right extend %f +/- %f, of normal %f +/- %f\n',...
        candidate, numel(candLocs),...
        nanmean(lensLE), nanstd(lensLE),...
        nanmean(lensRE), nanstd(lensRE),...
        lenSegment(ii), varSegment(ii));
end
