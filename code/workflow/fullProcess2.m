function [bouts, newSongs, boutCorrMatrix] = fullProcess2(matFile, params, varargin)
%% startup - not needed if running in cell mode
if nargin < 2
    params = defaultParams;
end
params = processArgs(params, varargin{:});

%% load recording from spike6 file
% is this a spike 6 matfile?
if numel(who('-file',matFile)) > 1,
    error('fullProcess2:notMATFile','Not a matfile exported from spike 6');
end


disp('Loading...');
songStruct=load(matFile);
fld=fieldnames(songStruct);
songStruct=songStruct.(fld{1});

neededFields = {'values','length','interval'};
%spike6Flds = {'title','comment','interval','scale','offset',...
%    'units','start','length','values'}';

if ~numel(intersect(fieldnames(songStruct), neededFields)) == numel(neededFields)
    %error('notASpike6File','Not a matfile exported from spike 6');
    error('fullProcess2:notASpike6File', 'Not a structure with the right data');
end

clear fld neededFields;
fs = 1 / songStruct.interval;

%% run rough spectrogram - dividing into regions
tic
disp('Separating into regions...');
regions = stepSpectrogram(songStruct,params); % for quiet BOS
toc

%% run fine spectrogram - dividing into syllables using noise reduction
%%%%%%%%%% TODO: derivative marking of boundaries... more effective than solid
%%%%%%%%%% thresholds
fprintf('%d regions found.\n',numel(regions));
tic
disp('Separating into syllables...');
[candidateSyl, features, featureReduction] = parseRegionsIntoSyllables(songStruct, regions, ...
    params); 
toc
%% save some work - partitions
save(prependForSave('divMT-',matFile),'regions','candidateSyl');

%% save some work - features
save(prependForSave('featMT-',matFile),'features', 'featureReduction');

%% Filtering out noise - optional step 
tic
disp('Filtering out noise...');
addpath([pwd '/boosting']);
%load('noiseClassifierMT','forestMT');
load('newNoiseClassifierMT','forestMT');

%isNoise = testClass(forestMT, featureReduction);
isNoise = false(size(candidateSyl));

[candidateSyl(isNoise).type] = deal(3);
vocalizations = candidateSyl(~isNoise);
fprintf('%d syllables found, %d vocalizations.\n',numel(featureReduction), numel(vocalizations));
toc
%% save some optional work - 
save(prependForSave('divMT-',matFile),'-append','vocalizations');

%% arrange into motifs, pretty simple so we can run a one liner for it
tic
disp('Arranging into bouts...');
bouts = arrangeBouts(vocalizations);
toc

%% constructing similarity matrix - 
% this may take a long time, depending on
% how similarity is calculated
tic;
disp('Calculating similarity...');

% not sure which one is better

% slower, syllable matching by scanning - also runs into possible memory
% issues
profile on;
[boutCorrMatrixmean, boutCorrMatrix] = syllableAllCross(songStruct, vocalizations);
profile viewer
% this method is faster but may be sensitive to boundary marking
if ~exist('isNoise')
    isNoise = false(size(candidateSyl));
end

%boutCorrMatrix = 1 - featureDistance(featureReduction(~isNoise));
toc;
% save some more work
% this guy should now return a full matrix
save(prependForSave('corrMT-', matFile), 'boutCorrMatrixmean','boutCorrMatrix');
%% construct clusters of letters
tic;
disp('Hierarchical clustering and alphabet creation...');
[stringRep, clusterIdxs] = createAlphabet(vocalizations, boutCorrMatrix, songStruct, ...
    params,'fs',fs, 'playsample', false, 'plot', true);
strCelled = cellstr(stringRep');
[vocalizations.type] = strCelled{:};

% TODO TODO TODO : peel out the alphabet browser (trying to check
% consistency of each syllable label)
%% arrange into bouts
% mark the beginning of a new bout
followsSilence = [false ([vocalizations(2:end).start] - [vocalizations(1:end-1).stop] > params.silentPeriod)];

stringRep = interlaceVectors(stringRep, followsSilence, ' ');
clusterIdxs = interlaceVectors(clusterIdxs, followsSilence, NaN);
boutCorrMatrix2 = interlaceMatVec(boutCorrMatrix2,followsSilence,NaN);
spacedVocals = interlaceVectors(vocalizations, followsSilence, initEvents(1));
toc;

%% taking most common song and finding new songs
tic; 
disp('Finding most common subword and expanding on it w/ similarity...');

% assumes one triangle of the matrix is zeros or NaNs;
mirroredMatrix = boutCorrMatrix + boutCorrMatrix' - diag(diag(boutCorrMatrix));
% get the longest song

% The song should meet certain criteria: it is the most repeated set of
% syllables and it itself contains distinct syllables.  It becomes the
% template for a song, if it is long enough (>= 3 syllables).
newSongs = rescanForCommonString(mirroredMatrix, stringRep, spacedVocals);
[~,mostSongCounts] = max(cellfun('prodofsize',newSongs));

newSongs = newSongs(mostSongCounts);
save(prependForSave('songMT-',matFile), 'newSongs');
toc;

%% now align
tic
songOffsets = calculateAlignmentOffsets(songStruct,newSongs, params, 'preroll',100,'postroll',100);
[alignedSongs, songSyllables] = alignSegments(songStruct, newSongs, ...
    songOffsets, params, 'postroll',100,'preroll', 100);

%% train a standard song
[prototypeSong, syllableFxns] = trainSongRecognizer(songStruct, alignedSongs, songSyllables, params, 'postroll',100,'preroll', 100);

%% get the best exemplar for that song 
stdizeScore = zeros(numel(alignedSongs), numel(syllableFxns));
progressbar('Finding standardized song scores')
for ii = 1:numel(alignedSongs); 
    progressbar(ii/numel(alignedSongs));
    stdizeScore(ii,:) = scoreAsSong(songStruct,alignedSongs(ii),prototypeSong,syllableFxns); 
end

[foo, exemplarSongIdx] = max(sum(stdizeScore,2));
fprintf('Best song example is #%d\n', exemplarSongIdx);
exemplarSong = alignedSongs(exemplarSongIdx);
plotWaveform(getClip(exemplarSong, songStruct));
toc;

%% save data
save('-v7.3',prependForSave('alignedSongs-', matfile),'exemplarSong','alignedSongs','songSyllables', 'prototypeSong', 'syllableFxns');


%% relabel each syllable according to the standardized score
[foo,bestLabel] = max(stdizeScore,2);
% TODO: investigate how scores could be so low when classification results
% were confident
    
%% now identify as many songs as possible using a rough-spectrogram based match scan
exemplarClip = getClip(exemplarSong, songStruct);
params.rough.fs = fs;
exemplarSpec = getMTSpectrumStats(exemplarClip, params.rough);

fieldsToRemove = intersect(fieldnames(exemplarSpec),{'psd','spectrum','freqs','deriv'});    
exemplarSpec = rmfield(exemplarSpec, fieldsToRemove);
overlapSplit = 1.2 * (exemplarSong.stop - exemplarSong.start);
LEN_OVERLAP = params.fs * overlapSplit / 1000;
LEN_SEGMENT = (songStruct.length - LEN_OVERLAP) / params.Nsplits + ...
    LEN_OVERLAP;

for ii = 15:params.Nsplits
    roughWindow.idxStart = floor(1 + (LEN_SEGMENT - LEN_OVERLAP) * (ii-1));
    roughWindow.idxStop = floor(roughWindow.idxStart + LEN_SEGMENT);
    roughWindow.start = roughWindow.idxStart / fs;
    roughWindow.stop  = roughWindow.idxStop  / fs;
    
    fprintf('Searching on clip %d...\n',ii);
    sample = getClip(roughWindow, songStruct);
    spec = getMTSpectrumStats(sample, params.rough);
    spec = rmfield(spec,fieldsToRemove);
    
    fprintf('Running similarity scan...\n');
    simScan = similarityScan(spec, exemplarSpec);
    plot(linspace(roughWindow.start, roughWindow.stop, numel(simScan)), simScan); xlim([roughWindow.start roughWindow.stop]); pause(0.5);
    
end

%% done!
end
