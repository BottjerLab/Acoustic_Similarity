function newSongs = rescanForCommonString(corrMatrix, stringRep, vocals)

% this takes the recording string, finds the most commong substring of a
% certain stable duration, and looks for other sequences that are similar enough
% to it
% corrMatrix should be full, not triangular
minSylls = [0,0,0,0,1,1,1];
songCandidate = cell(numel(minSylls));
newSongs = cell(0);
newSongsTried = 0;

% assumes corr matrix is full, not triangular
for NSyl = 1:7
    [ngram(NSyl).strs, ngram(NSyl).freqs, ngram(NSyl).idxs] = ...
        mostCommonSubstring(stringRep,NSyl,minSylls(NSyl));
    [lenCand, varCand] = isolateNGram(vocals,ngram(NSyl).strs,ngram(NSyl).idxs,NSyl);
    
    % if the variance of the candidate bout length is smaller than a
    % millisecond and there are at least 3 syllables
    % check the parallel similarity between strings...
    for kk = 1:numel(ngram(NSyl).strs)
        currSongCand = ngram(NSyl).strs{kk};
        if varCand(kk) < 0.1 && NSyl >= 3
            fprintf('Trying song ''%s''', currSongCand);
            songNgramIdx = find(strcmp(ngram(NSyl).strs,currSongCand),1);
            currSelectedSongStarts = find(ngram(NSyl).idxs == songNgramIdx)';
            nSeld = numel(currSelectedSongStarts);
        
        idxs = 0:NSyl-1;
        currRanges = currSelectedSongStarts(:,ones(NSyl,1)) + ...
            idxs(ones(1,nSeld),:);
        
        
        % we use segment similarity weighted by syllable lengths
        sylLengths = zeros(1,NSyl);
        for ii = 1:NSyl
            vocIdxs = (ngram(1).idxs == find(strcmp(ngram(1).strs, currSongCand(ii))));
            sylLengths(ii) = mean([vocals(vocIdxs).stop] - ...
                [vocals(vocIdxs).start]);
        end
        sylWeights = sylLengths / nansum(sylLengths);
        
        totalSim = zeros(1,numel(vocals) - NSyl + 1);
        errorSim = zeros(1,numel(vocals) - NSyl + 1);
        for ii = 1:numel(vocals) - NSyl + 1
            testIdxs = ii:ii+NSyl - 1;
            testTable = testIdxs(ones(1,nSeld),:);
            reducedSimMatrix = corrMatrix(sub2ind(size(corrMatrix),testTable,currRanges));
            
            avgSim = nansum(reducedSimMatrix .* sylWeights(ones(1,nSeld),:),2);
            totalSim(ii) = nanmean(avgSim);
            errorSim(ii) = nanstd(avgSim);
        end
        
        songOnsets = [vocals(1:numel(totalSim)).start];
        
        
        % split into song vs. no song with a mixture of gaussians model
        gmData = [totalSim' errorSim'];
        twoClusterModel = gmdistribution.fit(gmData, 2);
        % pick the cluster with the largest similarity between the two
        [foo, songClusterIdx] = max(squeeze(twoClusterModel.mu(:,1))); 
        [potClusterIdxs, foo, pvals] = cluster(twoClusterModel,gmData);
    
        % plot the song distribution
        newSongCands = (potClusterIdxs == songClusterIdx);
        clf
        plot(songOnsets,totalSim,'r.');
        hold on;
        plot(songOnsets(currSelectedSongStarts),totalSim(currSelectedSongStarts),'g.');
        plot(songOnsets(newSongCands),totalSim(newSongCands),'b.');
        title(sprintf('%d new song candidates', sum(newSongCands)));
        % set correct axes
        endOfSylls = vocals(end).stop;
        xlim([0, endOfSylls + 1]);
        set(gca,'XTick',0:30:endOfSylls);
        hold off
        
        pause;
        
        % accept all new song cands
        newSongsTried = newSongsTried + 1;
        newSongIdxs = find(newSongCands);
        newSongs{newSongsTried} = initEvents;
        for ii = 1:numel(newSongIdxs)
            newSongs{newSongsTried}(ii) = bookendedClip(vocals(newSongIdxs(ii):newSongIdxs(ii)+NSyl-1));
        end      
%         break;
        end
    end
end
