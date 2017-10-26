% make spectrograms

% scan for each day

rep = reportOnData;
nBirds = numel(rep);
%%
params = processArgs(defaultParams, 'fine.freqBands', ...
    linspace(1,10240,params.inter.NfreqBands));
for ii = 3:nBirds
    birdRecords = rep{ii};
    nSessions = numel(birdRecords);
    birdSessions = {birdRecords.sessionID};
    birdID = strtok(birdSessions{1},'_');
    
    ages = getAgeOfSession(birdSessions);
    
    DRmotifs = initEmptyStructArray({'type','start','stop','idxStart','idxStop','file','age'},0);
    hasSylls = true(1,nSessions);
    % find manualMotifs
    
    for jj = 1:nSessions
        thisManifest = birdRecords(jj).manifest;        
        if ~findInManifest(thisManifest, 'approvedSyllables')
            hasSylls(jj) = false;
            continue; 
        end
        isMotifIn = findInManifest(thisManifest, 'manualMotifs');
        
        if isMotifIn
            mM = loadFromManifest(thisManifest, 'manualMotifs');
        else
            aS = loadFromManifest(thisManifest, 'approvedSyllables');
            mM = mergeGaps(aS, 0.5);
        end
        nM = numel(mM);
        if nM < 6,
            continue;
        else
            range = floor(nM/2)-2:floor(nM/2)+2;
        end
        motifs = mM(range);
        
        [~, rootFile] = findInManifest(thisManifest, 'sAudio');         
        [motifs.file] = deal(rootFile);
        [motifs.age ] = deal(ages(jj));
        if isfield(motifs,'length')
            motifs = rmfield(motifs, 'length');
        end
        if size(motifs,2) > 1
            motifs = motifs';
        end
        DRmotifs = [DRmotifs; motifs];
    end
    % max length
    maxLength = max([DRmotifs.stop] - [DRmotifs.start]);
    
    uAges = unique([DRmotifs.age]);
    for jj = 1:numel(uAges)
        thisAge = uAges(jj);
        thisAgeMotifs = DRmotifs([DRmotifs.age] == thisAge);
        for kk = 1:numel(thisAgeMotifs)
            motifLen = thisAgeMotifs(kk).stop - thisAgeMotifs(kk).start;
            [cl, fs] = getClipAndProcess([],thisAgeMotifs(kk));
            
            specFeats = getfield(params, 'fine'); specFeats.fs = fs; 
            spec = getMTSpectrumStats(cl, specFeats);
            him = plotDerivGram(spec,[],'dgram.minContrast', 1e-11);            
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
            set(gca, 'Box', 'off');
            set(gcf, 'Color', [1 1 1]);
            
            pastUnits = get(gcf,'Units');
            set(gcf,'Units','normalized');
            hh = get(gca,'Position'); 
            hh(1) = 0; hh(3) = 1;
            set(gca,'Position', hh);
            set(gca,'YTick', []);
            jj = get(gcf, 'Position');
            jj(1) = 0; jj(3) = motifLen / maxLength;
            set(gcf,'Position', jj);
            set(gcf,'Units', pastUnits);
            
            figPath = [pwd filesep 'figures' filesep 'exampleMotifs' filesep];
            export_fig(sprintf('%s%s-age%d-%d.jpg',figPath, birdID, thisAge, kk));
        end
    end    
end

