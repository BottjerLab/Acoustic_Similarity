function [varMotifLengths, motifSylRatios, cvMotifLens, uAges] = motifVariation(birdID)
% function motifVariation
% for each age of a bird, get the variance in motif lengths, motif syllable
% ratios (last two syllables / first two syllables), and coefficient of
% variation of motif lengths

thisBirdRecords = reportOnData(birdID);
ages = getAgeOfSession({thisBirdRecords.sessionID});
[uAges, ~, rIdxAge] = unique(ages);

lastValidAge = find(~isnan(uAges),1,'last');
if lastValidAge < numel(uAges)
    rIdxAge(rIdxAge > lastValidAge) = NaN;
    uAges(lastValidAge+1:end) = [];
end

nAges = numel(uAges);
varMotifLengths = zeros(1,nAges);
motifSylRatios  = zeros(1,nAges);
cvMotifLens     = zeros(1,nAges);

for jj = 1:numel(uAges)
    theseAgeRecords = thisBirdRecords(rIdxAge == jj);    
    % collect two things - length of motifs and ratio of syllables
    % in second half to first half
    motifLens = [];
    sylRatios = [];
    
    for kk = 1:numel(theseAgeRecords)
        if ~findInManifest(theseAgeRecords(kk).manifest, 'derivedMotifs')
            fprintf('derivedMotifs not found, continuing for session %s...\n', theseAgeRecords(kk).sessionID);
            continue;
        end
        mM = loadFromManifest(theseAgeRecords(kk).manifest, 'derivedMotifs');
        
        % take off 50 ms because these syllables all have 50 ms appended to
        % them
        mM = addPrePost(mM, defaultParams, 'preroll', -50, 'postroll', 0);
        if ~findInManifest(theseAgeRecords(kk).manifest, 'approvedSyllables')
            fprintf('approvedSyllables not found, continuing for session %s...\n', theseAgeRecords(kk).sessionID);
            continue;
        end
        
        aS = loadFromManifest(theseAgeRecords(kk).manifest, 'approvedSyllables');
        % lengths of motifs
        mLens = [mM.stop] - [mM.start];
        if size(mLens,2) == 1, mLens = mLens'; end % force row format
        motifLens = [motifLens mLens];
        
        % ratio of lengths of last two syllables / first two syllables motif
        parentIdxs = findParent(mM, aS);
        ratioInMotif = NaN(1,numel(mM));
        for ll = 1:numel(mM)
            % there should be at least 4 syllables in a motif
            childSylls = aS(parentIdxs == ll);             
            if numel(childSylls) < 4, continue; end; 
            firstSylls = childSylls(1:2);
            lastSylls  = childSylls(end-1:end);
            ratioInMotif(ll) = mean([lastSylls.stop] - [lastSylls.start]) / ...
                mean([firstSylls.stop] - [firstSylls.start]);
        end
        sylRatios = [sylRatios ratioInMotif];
    end % end loop over sessions within an age    
    varMotifLengths(jj) = nanstd(motifLens);
    motifSylRatios(jj) = nanmean(sylRatios);
    cvMotifLens(jj) = varMotifLengths(jj) / nanmean(motifLens);
end % end loop over ages
