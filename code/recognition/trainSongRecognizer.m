function [prototypeSong, syllableFxns] = trainSongRecognizer(songStruct, alignedSongs, songSylls, params, varargin)

if nargin < 4
    params = defaultParams;
end;
params = processArgs(params, varargin{:});

fs = 1/songStruct.interval;
% apply pre/postroll to all songs first
nSongs = numel(alignedSongs);
for ii = 1:nSongs
    alignedSongs(ii) = addPrePost(alignedSongs(ii),params);
end
% given a set of aligned songs and borders,
% define the intervals for a typical song, and a

nTotalSylls = numel(songSylls);
%% first vote on the syllable boundaries
songStartVotes = zeros(1,nTotalSylls);
songStopVotes = zeros(1,nTotalSylls);
for ii = 1:nSongs
    [theseSubSylls, theseIdxs] = getSubEvents(alignedSongs(ii),songSylls);
    zeroAdjustedSylls = adjustTimeStamps(theseSubSylls, -alignedSongs(ii).start);
    songStartVotes(theseIdxs) = [zeroAdjustedSylls.start];
    songStopVotes(theseIdxs)= [zeroAdjustedSylls.stop];
end

% how much variance is there in syllable onsets?
kernelWidth = 0.005; %s

% time resolution parameter
widthRes = 5e-4; %s
songLen = alignedSongs(1).stop - alignedSongs(1).start;
songIdxLen = alignedSongs(1).idxStop - alignedSongs(1).idxStart + 1;
timeEv = 0:widthRes:songLen;

%% estimate the start/stop positions of each event,
% renormalizing so that each syllable has unit 'vote' at its location

kernelParams = {'width',kernelWidth,'kernel','epanechnikov'};

% instead of doing the math, we just find the height empirically
unnormedHeight = max(ksdensity((1:nTotalSylls) * 15 * kernelWidth,timeEv,kernelParams{:}));
startPeakVote = ksdensity(songStartVotes,timeEv,kernelParams{:}) / unnormedHeight;
stopPeakVote = ksdensity(songStopVotes,timeEv,kernelParams{:})  / unnormedHeight;

%% display

% get average clip for display purposes
clipAvg = averageWaveform(songStruct,alignedSongs,songIdxLen);

subplot(2,1,1);
plot((1:songIdxLen)/fs, clipAvg, 'b-');
title('Average Waveform');
xlim([0 songIdxLen/fs]);
subplot(2,1,2);
plot(timeEv, startPeakVote,'c-', timeEv, stopPeakVote,'m-'  )
title('Finding boundaries for songs');
xlim([0 timeEv(end)]);

%% here we use thresholds on the vote to segment syllables
% what percent needs to vote on a syllable?
voteFraction = 0.35;
voteThresh = voteFraction * nSongs;

peakDist = kernelWidth / widthRes;
peakParams = {'minPeakHeight', voteThresh,'minPeakDistance', peakDist};
[startY, stdSyllStarts] = findpeaks(startPeakVote, peakParams{:});
[stopY , stdSyllStops ] = findpeaks(stopPeakVote , peakParams{:});
nStandardSylls = numel(stdSyllStarts);
stdSyllStarts = timeEv(stdSyllStarts);
stdSyllStops  = timeEv(stdSyllStops );
stdSyllLengths = stopY - startY;

%% check (a) are these syllables interwoven at the right intervals?
if all(stdSyllLengths > 0) && all(stopY(1:end-1) - startY(2:end) > 0)
    error('BadIntervals','Intervals are not correct'); % needs to be fixed
end

%% (b) are they symbols representative of real syllables?
% find if any syllables actually match these syllables
sylMatchIndex = NaN(1,nTotalSylls);
for ii = 1:nSongs
    [theseSubSylls, theseIdxs] = getSubEvents(alignedSongs(ii),songSylls);
    zeroAdjustedSylls = adjustTimeStamps(theseSubSylls, -alignedSongs(ii).start);
    theseStarts = [zeroAdjustedSylls.start];
    theseStops  = [zeroAdjustedSylls.stop ];
    
    for jj = 1:numel(theseIdxs)
        matchSyll = abs(theseStarts(jj) - stdSyllStarts) < 2 * kernelWidth & ...
            abs(theseStops(jj) - stdSyllStops) < 2 * kernelWidth;
        if sum(matchSyll) == 1
            sylMatchIndex(theseIdxs(jj)) = find(matchSyll);
        elseif sum(matchSyll) > 1
            warning('matchTooLoose','match criteria are too loose');
        end
    end
end


hold on;
plot(stdSyllStarts, startY, 'c^', 'markerfacecolor', [0 0 1]);
plot(stdSyllStops , stopY , 'm^', 'markerfacecolor', [1 0 0]);
hold off;

%% (c) are the syllables they model similar?
if params.plot
    for ii = 1:nStandardSylls
        theseSylls = songSylls(ii == sylMatchIndex);
        nTheseSylls = numel(theseSylls);
        hl=figure(); clf
        fprintf('Showing spectrograms for syllable #%d/%d (%d)...\n',ii,nStandardSylls,nTheseSylls);
        subplot(ceil(nTheseSylls/4),4,1);
        progressbar(0);
        for jj = 1:nTheseSylls
            figure(hl);
            subplot(ceil(nTheseSylls/4),4,jj);
            waveform = getClip(theseSylls(jj), songStruct);
            
            if params.playsample  && jj == nTheseSylls, playSound(waveform,fs,true); end;
            plotWaveform(waveform,fs);
            set(gca,'XTick',0:0.01:max(get(gca,'XLim')));
            %params.fine.fs = fs;
            %spec = getMTSpectrumStats(waveform,params.fine);
            %plotDerivGram(spec);
            %ylabel(''); set(gca,'YTick',[]); xlabel(''); set(gca,'XTick',[]);
            progressbar(jj/nTheseSylls);
            if jj == 1; title(sprintf('Syllable #%d (%d)',ii,nTheseSylls)); end;
        end
        colormap(gray)
        set(hl,'visible','on');
        saveCurrFigure(sprintf('figures/stdSylls-%02d.tif',ii));
    end
end

%% prepare output
prototypeSong = initEvents(nStandardSylls);
foo = num2cell(1:nStandardSylls); [prototypeSong.type] = foo{:};
foo = num2cell(stdSyllStarts - params.preroll/1000); [prototypeSong.start] = foo{:};
foo = num2cell(stdSyllStops  - params.preroll/1000); [prototypeSong.stop ] = foo{:};
foo = num2cell(floor(fs * (stdSyllStarts - params.preroll/1000))); [prototypeSong.idxStart] = foo{:};
foo = num2cell(floor(fs * (stdSyllStops  - params.preroll/1000))); [prototypeSong.idxStop ] = foo{:};

% now train the classifiers to recognize these syllables
% start by extracting all the features
if nargout == 2
    progressbar('extracting features');
    for ii = nTotalSylls:-1:1 % loop run backward so that preallocation happens
        params.fine.fs = fs;
        clip = getClip(songStruct, songSylls(ii));
        featureReduction(ii) = extractFeatures(getMTSpectrumStats(clip, params.fine));
        progressbar(1-(ii-1)/nTotalSylls);
    end
    
    for ii = 1:nStandardSylls
        isThisSyll = (ii == sylMatchIndex);
        fprintf('Training classifier on syllable #%d (%d)...\n',ii,sum(isThisSyll));
        classer{ii} = trainClassifier(songStruct,songSylls,featureReduction, isThisSyll);
        
        syllableFxns{ii} = @(feature) testClass(classer{ii}, feature);
    end
end
end

function ret = averageWaveform(songStruct, clips, clipLen)
if nargin < 3,
    clipLen = NaN([clips.idxStop] - [clips.idxStart]);
end
allClips = zeros(numel(clips), clipLen);
for ii = 1:numel(clips)
    songLen = clips(ii).idxStop - clips(ii).idxStart + 1;
    allClips(ii,1:songLen) = getClip(songStruct, clips(ii))';
end;
ret = nanmean(allClips,1);
end
