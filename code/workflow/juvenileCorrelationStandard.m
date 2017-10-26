% JUVENILECORRELATION run this to correlate juvenile-tutor similarity with
% firing rate
% cell by cell if you're patient (Ctrl-Enter on Windows)
% each cell is the light-yellow colored code section

%% This script 
% (1) isolates songs in a juvenile recording, 
% (2) defines syllables in a tutor song,
% (3) isolates syllables in the juvenile recording, 
% (4) labels the syllables according to the tutor song, 
% (5) calculates the syllable based similarity between the tutor syllables
% and the juvenile syllables, and 
% (6) tests the correlation between similarity 
% to the tutor syllables and firing rate.

% Note 1: For this to work, all of Z:\john\parsing\code needs to be in the
% matlab path.  The easiest way to do this is to 
% > cd Z:\john\parsing\
% > startup
% This automatically adds paths and suppresses some annoying warnings.
% Alternately, just highlight this cell and type Ctrl-Enter.

cd Z:\john\parsing
startup
disp('Started up.')
%% Note 2: you can always uiload any of the files that are "uisaved"
uiload
% For silly reasons, the song file (without the prefix) should be loaded
% with the cell below.
%% prep step: load song file
% define files to be operated over
[matFile, matpath] = uigetfile('*.mat','Please choose the song Spike2 file','data');
% load file
songStruct = loadSpikeAudioStruct([matpath matFile]);
songStruct.title = matFile;
fs = 1/songStruct.interval;
fprintf('Loaded spike audio %s.\n',matFile);
%% prep step: set parameters from an alias
params = processArgs(defaultParams,'R204');
params.fs = fs; 
%% (0.5) isolate sounds first
sounds = stepSpectrogram(songStruct, params,'Nsplits',400,'plot',true);
%% arrange into bouts
[bouts, nSyllsInBout] = arrangeBouts(sounds,params,'silentPeriod',0.8);
%% reporting
fprintf('%0.2f s in bouts, %0.1f%% of total recording...\n',...
    sum([bouts.stop] - [bouts.start]), ...
    100*sum([bouts.stop] - [bouts.start]) / (songStruct.interval * songStruct.length));
%% filtering by time (based on the pre-listen) - CAREFUL
% MANUAL STEP: specify
starterRegion = getWholeSongAsRegion(songStruct);
bStart = [870 911 915 928.5 935   943   948   968   974.8 980 985.5 1010   1014 1026   1031 1033.7 1064.5 1117 1147.8 1153   1156 1162.8 1178 1186 1207  1211  1235.5 1290 1366 1490   1492.5 1505];
bStop =  [875 913 926.5 932 941.5 945.5 951.5 973.6 978  983 1006.3 1013.3 1016 1030.5 1033 1064   1112 1144.5 1151.5 1155.5 1162.5 1172 1185 1204 1210.6 1223.7 1271 1336 1483 1492.3 1501   1508] ;
tStart = [200 875 926 936.5 938.2 941 952 973.5 982.5 1006.3 1011 1022 1027 1030 1062   1112   1144 1149 1155   1160 1171.5 1183.5 1203.7 1210 1212.8 1223.7 1234 1241   1244.5 1269   1296  1347  1488 1496   1503 1630];
tStop =  [380 915 933 938   940   943 957 976.8 985.5 1010.5 1016 1025 1029 1034 1064.5 1117.2 1148 1152 1156.5 1163 1174 1186.5 1207 1212.5 1214.5  1233.2 1236  1243.5 1245.5 1273.5 1299.5 1351 1493 1500.5 1507 1686];

bManRange = eventFromTimes(bStart, bStop,fs);
tManRange = eventFromTimes(tStart, tStop,fs);

tbManRange = findIntersection(bManRange, tManRange);
[bManRange, tManRange] = findUniqueAreas(bManRange, tManRange);
%%
starterRegion = returnedMotifs;
hasBParent = findParent(bManRange, starterRegion); 
bosOnlyBouts = returnedMotifs(~isnan(hasBParent));
if ~isempty(bosOnlyBouts), [bosOnlyBouts.type] = deal('BOS'); end;

hasTBParent = findParent(tbManRange, starterRegion); 
tutBosBouts = returnedMotifs(~isnan(hasTBParent));
if ~isempty(tutBosBouts), [tutBosBouts.type] = deal('TUT+BOS'); end;

hasTParent = findParent(tbManRange, starterRegion); 
tutOnlyBouts = returnedMotifs(~isnan(hasTParent));
if ~isempty(tutOnlyBouts), [tutOnlyBouts.type] = deal('TUT'); end;

manualMotifs = sortBy([bosOnlyBouts; tutOnlyBouts; tutBosBouts],'start')

%% refine sound boundaries
ROIs = splitIntoOverlap(bosOnlyBouts,100,0.05);
len = [ROIs.stop] - [ROIs.start];
idx = 1:numel(ROIs);
%[candidateSyl, features, featureReduction] = ...
%[candidateSyl, features] = ...
candidateSyl = ...
    parseRegionsIntoSyllables(songStruct, ROIs, params, ...
    ...%'plot',true, 'quiet', false, 
    'parseSpecType','rough',...
    ...%'syllable.minPower', 1.5e-7, 'syllable.minRiseAmp', 8, 'syllable.minLength', 15, ...
    'plot', false, 'pause', false, 'postroll', 100,'dgram.minContrast',1e-10,...
    'syllable.flatFactor', 8); 

%% re-arrange into bouts
[bouts, nSyllsInBout] = arrangeBouts(candidateSyl,params,'silentPeriod',0.234);
%% re report
fprintf('%0.2f s in bouts, %0.1f%% of total recording...\n',...
    sum([bouts.stop] - [bouts.start]), ...
    100*sum([bouts.stop] - [bouts.start]) / (songStruct.interval * songStruct.length));

%% do a quick filter: bouts shorter than 200 ms probably do not have motifs/subsong
minBoutLength = 0.2; % in seconds;
bouts(([bouts.stop]-[bouts.start]) < minBoutLength) = [];
%% save some work - automatically found partitions
uisave({'candidateSyl', 'bouts'},...
[matpath prependForSave('divMT-',matFile)]);

%% (x.8) find noise 
% define some noise with space between bouts
% NB: since most recordings/recording settings are continuous,
% this can be done once per main session 
% (not necessarily for each 30-minute subsession, but check quality)
notNoise = markedSongs;
approved = false;
while ~approved
    candidateNoise = autodetectNoise(songStruct, notNoise, params);
    fprintf('Is this clip pure noise? Mark if yes...\n'); 
    approved = markRegions(songStruct,candidateNoise);
    % just black this section out and try again
    if ~approved
        notNoise = [notNoise; eventFromTimes(0,candidateNoise.stop,fs)];
    end
end

noiseMask = noiseAnalysis(songStruct, candidateNoise);

% save noise mask
uisave('noiseMask',[matpath prependForSave('noiseMask-', matFile)]);
%%
params.noiseFilter = noiseMask;
%% (1-eps) alternate manual edit: send interval times as text file to spike2
% load with import and procEvts
birdID = strtok(matFile, '_');

sendStarts = returnStarts;
sendStops  = returnStops;
dlmwrite([pwd '\data\' birdID '\boutTimes-' matFile '.txt'],[sendStarts; sendStops]','delimiter',' ','precision',6);
touch([pwd '\data\' birdID '\motifReturn-' matFile '.txt']);
%% (1-eps) load back times and inherit labels
birdID = strtok(matFile, '_');
returnFil = [pwd '\data\' birdID '\motifReturn-' matFile '.txt']
returnFid = fopen(returnFil);

% pick up header and dead line - first column vector is starts
fgetl(returnFid); fgetl(returnFid);
returnStarts = textscan(returnFid,'%f'); returnStarts = returnStarts{1};

% pick up header and dead line - second column vector is stops 
fgetl(returnFid); fgetl(returnFid);
returnStops = textscan(returnFid,'%f'); returnStops = returnStops{1};
    
fclose(returnFid);

returnSIndices = [zeros(numel(returnStarts),1); ones(numel(returnStops), 1)];
[fusedTimes,sIdxs] = sort([returnStarts; returnStops]);
returnSIndices = returnSIndices(sIdxs);

% make sure alternating
if any(diff(returnSIndices) == 0)
    isBadIndices = [diff(returnSIndices) == 0; sIdxs(end)==0];
    fusedTimes(isBadIndices)
end
assert(numel(returnStarts) == numel(returnStops));
assert(all(returnStarts < returnStops) && all(returnStarts(2:end) > returnStops(1:end-1)));

returnedMotifs = eventFromTimes(returnStarts, returnStops, fs);

%% (1) isolate the songs - mark only bout regions that actually contain song
disp('Manual adjustment to find bouts...');

ROIs = bouts;
%parentage = findParent(bouts, candidateSyl);
containsSong = false(1,numel(ROIs));
num2cell([ROIs.stop] - [ROIs.start]); [ROIs.length]=ans{:};
    
sortedROIs = sortBy(ROIs,'length','descend');

% don't look at bouts that are too short - those can't contain any song
%if boutLen < 0.1, continue; end;
   
fprintf('Do these clips contain song? Mark if yes...\n');
containsSong = markRegions(songStruct, sortedROIs);
%[songType, markedROIs] = multiMark(songStruct, sortedROIs);
%%
manualBouts = sortBy(sortedROIs(~containsSong),'start');
%% save first pass on marked songs
uisave({'bouts','manualBouts'},...  
    [matpath prependForSave('manualBouts-',matFile)]);
%% (1.3) try to fix sounds by placing cuts correctly
ROIs = manualBouts;

%note: minPower MUST BE ADJUSTED for different settings
newROIs = fixBoundaries(songStruct, ROIs, params, 'plot', true, 'syllable.minPower', 5e-6, 'boundaryAdvance', 0.2);
%% (1.5 auto) chop bouts into real motifs 

% TODO: include learner training in here
prePostParams = processArgs(params,'preroll',500,'postroll',500);

% make sure there aren't any overlaps in the pre/post rolled bouts
if ~isempty(findOverlaps(addPrePost(manualBouts, prePostParams)))
    warning('Some bouts overlap');
end
motifs = initEvents(0);

sortedBouts = manualBouts;
num2cell([sortedBouts.stop] - [sortedBouts.start]); 
[sortedBouts.length] = ans{:};
% initialize learner
%songLearner = Learner(songStruct, [], params);

ROIs = sortedBouts;
for ii = 1:numel(ROIs)
    thisROI = ROIs(ii);
    fprintf('Song %d/%d (length %0.2f s)... ', ii, numel(ROIs), thisROI.stop - thisROI.start);
    % initialize learner
    if ~exist('songLearner'), songLearner = Learner(songStruct, [], params); end
    if ~exist('callLearner'), callLearner = Learner(songStruct, [], params); end
    
    if songLearner.nExamples < 5 % minimum number of songs for robust classification
        % keep manually classifying until we get 5 songs
        tmpSongs = plotAndAdjust(songStruct, thisROI, ...
        addPrePost(thisROI, prePostParams), params,...
        'editSpecType','inter',...
        'inter.features', {'totalPower','wienerEntropy','deriv','waveform'});
        if ~isempty(tmpSongs)
    %        songLearner.addSpectra(songStruct, tmpSongs);
    %        songLearner.adaptThreshold;
        end
    else
        tic
  %      fprintf('Song %d/%d (length %0.2f s)... ', ii, numel(markedSongs), markedSongs(ii).stop - markedSongs(ii).start);
        tmpSongs = scanForMatches(songStruct, addPrePost(thisROI, prePostParams), songLearner);
        toc
    end
    motifs = [motifs tmpSongs];
end

%% 
uisave('songLearner',[matpath prependForSave('songTemplate-',matFile)]);
uisave({'motifs'}, [matpath prependForSave('motifs-',matFile)]);

%% (1.5 manually) refine/subparse motif boundaries MANUALLY, 
%extending if necessary (define range below)
prePostParams = processArgs(params,'preroll',500,'postroll',500);
newROIs = initEvents(0);

% initialize input
ROIs = mM(isMixed);

%findUniqueAreas(manualBouts, manualMotifs);

% make sure there aren't any overlaps in the pre/post rolled bouts
fprintf('Maximum bout length is %0.1fs...\n', max([ROIs.stop] - [ROIs.start]));
if ~isempty(findOverlaps(addPrePost(ROIs, prePostParams)))
    % find closest approach
    sortROIs = sortBy(ROIs, 'start');
    closestApproach = min([sortROIs(2:end).start] - [sortROIs(1:end-1).stop]);
    
    warning('Some bouts are closer than %dms apart, closest approach is %0.2fms...\n', ...
        prePostParams.preroll + prePostParams.postroll, closestApproach*1000);
    
end

%ROIs = sortBy(manualMotifs, 'length','descend'); %autoMotifs(strcmp({autoMotifs.type},'juvenile'));
nB = numel(ROIs);
for ii = 1:nB %NOTE: define range here;
    thisBout = ROIs(ii);
    fprintf('Bout %d/%d (length %0.2f s)... ', ii, nB, thisBout.stop - thisBout.start);
    if strcmp('L',thisBout.type)
        tmpSongs = thisBout;    
        fprintf('locked, continuing...\n');
    else
        %{
        fprintf('Are we perfect as is?\n');
        if false && markRegions(songStruct, thisBout), 
            tmpSongs = thisBout;
            newROIs = [newROIs; basicEvents(tmpSongs)];
            continue;
        end
        %}
        contextBout = addPrePost(thisBout, prePostParams);
        tmpSongs = plotAndAdjust(songStruct, thisBout, ...
            contextBout, params,...
            'editSpecType','inter', 'adjustLabels', true,...
            'doNoiseFilter',true,'noiseFilter',noiseMask,...
            'inter.features', {'totalPower','wienerEntropy','deriv','waveform'},...
            'dgram.minContrast',1e-10);
        %'noiseFilter', noiseMask);
    end
    newROIs = [newROIs; tmpSongs];
end
%% don't forget!
%manualMotifs = newROIs;
%uisave({'motifs','manualMotifs'}, [matpath prependForSave('motifs-',matFile)]);
uisave({'manualMotifs'}, [matpath prependForSave('motifs-',matFile)]);

%% auto correction: merge all overlapping segments
olaps = findOverlaps(manualMotifs);
adjustedMarkedSongs = manualMotifs; 
markDel = false(numel(manualMotifs),1);
for ii = 1:size(olaps,1)
    adjustedMarkedSongs(olaps(ii,1)).idxStart = manualMotifs(olaps(ii,2)).idxStart; 
    adjustedMarkedSongs(olaps(ii,1)).start = manualMotifs(olaps(ii,2)).start; 
    markDel(olaps(ii,1)) = true;
end
manualMotifs = adjustedMarkedSongs(~markDel); 

%% (1.6) auto correction: where there is an overlap, give overlapped portion to the later segment
olaps = findOverlaps(manualMotifs);
adjustedMarkedSongs = manualMotifs; 
for ii = 1:numel(olaps)/2
    adjustedMarkedSongs(olaps(ii,1)).idxStop = manualMotifs(olaps(ii,2)).idxStart; 
    adjustedMarkedSongs(olaps(ii,1)).stop = manualMotifs(olaps(ii,2)).start; 
end

manualMotifs = adjustedMarkedSongs; 
%% (2) parse juvenile motifs into syllables
ROIs = manualMotifs;%labeledBouts(strcmp('BOS',{labeledBouts.type}));
%ROIs = ROIs(randperm(numel(ROIs)));
params = processArgs(defaultParams,'fs',1/songStruct.interval, 'preroll', 30);
bosSyllables = ...
    parseRegionsIntoSyllables(songStruct, ROIs, params, 'Gy217',...
    'doFilterNoise',true,...
    'noiseFilter', noiseMask,'nps.reduction',-18, ...
    'dgram.minContrast',1e-10,'minCenterFreq', 800,...
    'plot', false, 'pause', false);

%%
uisave({'manualMotifs', 'syllables'}, [matpath prependForSave('syllables-',matFile)]);

%% (2.5) refine (and LABEL) parse with manual work
% do a manual refinement
% tutorParams = params;
% tutorParams.fine.fs = 1/tutorSong.interval; tutorParams.fs = 1/tutorSong.interval;
% tutorSpec = getMTSpectrumStats(tutorSong.values, tutorParams.fine);
% figure(1);
%plotAllFigures(tutorSpec, tutorSylls, tutorParams, 'showLabels',true);
%params = processArgs(defaultParams,'fs',1/songStruct.interval);
params = processArgs(defaultParams,'fs',1/songStruct.interval);

%define the specific regions to operate on (like function arguments)
ROIs = mM;
subROIs = juvSylls;

% resets new results if you want to start over
tmpSyllables = initEvents;

for ii = 1:numel(ROIs) % if you want to continue, change this range
    thisEv = ROIs(ii);
    fprintf('Editing %d/%d...\n', ii, numel(ROIs));
    revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
        'editSpecType', 'inter', 'adjustLabels',true, 'dgram.minContrast',3e-11,...
        'doFilterNoise',true, 'noiseFilter', noiseMask,'nps.reduction',-12);
    tmpSyllables = [tmpSyllables revisedSylls']; 
end

% moves new results to permanent location
manualSyllables = tmpSyllables;
clear tmpSyllables

%% (save) do this after you work on (2.5) %%%%
uisave({'syllables', 'manualSyllables'}, [matpath prependForSave('syllables-',matFile)]);

%% (2.7 optional) unsupervised labeling via 'agglomerative clustering'
% get distances between juvenile syllables 
ROIs = syllables;
profile on;
[boutDistMatrixMean, boutDistMatrix] = syllableAllCross(songStruct, ROIs);
profile viewer
profile off;
uisave({'boutDistMatrixMean','boutDistMatrix'},prependForSave('corrMT-', matFile));

% construct labels for letters according to clusters
nClusters = 5;
disp('Hierarchical clustering and alphabet creation...');
[stringRep, clusterIdxs] = createAlphabet(ROIs, boutDistMatrix, songStruct, ...
    [],'fs',fs, 'playsample', true, 'nClusters', nClusters,'clusterMethod', 'ward', 'plot',true);

clusteredROIs = ROIs;
strCelled = cellstr(stringRep');
[clusteredROIs.type] = strCelled{:};
[ROIs.clusterType] = strCelled{:};
%%
uisave({'syllables', 'manualSyllables', 'clusterSyllables'}, [matpath prependForSave('syllables-',matFile)]);
%% (3-eps) load up tutor songs... (ALTERNATELY, the tutor files usually exist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Here we study the tutor song, similarly parsing out the syllables and
%%%% computing features, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% who's the bird's daddy? what's his song like?
%birdTutor = inputdlg('What is the ID of the tutor of the bird? ', 'Tutor name',1,{'Y130')};
%birdTutor = [birdTutor '_TUT'];
%%% TODO: find a way to get to the song folder instead of assuming the file
%%% is in wavs directory
tutorDir = uigetdir('','Where are the tutor recordings?');

tutorFiles = dir([tutorDir filesep '*.wav']);
tutorFiles = strcat([tutorDir filesep],{tutorFiles.name});
%% (3.pre) set tutor parameters
tutorParams = processArgs(defaultParams,'fs',1/tutorStruct.interval,'preroll',0,'postroll',0);
%% (3) Pick the tutor song and study it
whichIsTutor = 0;
approved = false;
for ii = 1:numel(tutorFiles);
    tutorStruct = loadWavFile(tutorFiles{ii});
    
    tutorParams = processArgs(defaultParams,'fs',1/tutorStruct.interval,'preroll',0,'postroll',0);
    wholeRegion = getWholeSongAsRegion(tutorStruct);
    
    % pre highPass and amplify
    tutorStruct.values = 3*highPassSample(getClip(wholeRegion, tutorStruct), tutorParams);

    approved = markRegions(tutorStruct,wholeRegion,tutorParams,'preroll',0,'postroll',0);
    if(approved)
        whichIsTutor = ii; break; 
    end;
end

% do an automatic parse
[tutorSylls, features] = parseRegionsIntoSyllables(tutorStruct, wholeRegion,tutorParams,...
    'doFilterNoise',false,'syllable.comboLength',5);

%% do a manual refinement on syllable level
wholeRegion = getWholeSongAsRegion(tutorStruct);
tutorParams = processArgs(defaultParams,'fs',1/tutorStruct.interval,'preroll',0,'postroll',0);
tutorSylls = plotAndAdjust(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8,...
    'optGraphs',{'waveform', 'deriv','wienerEntropy','pitchGoodness', 'fundamentalFreq'});

%% do a manual refinement on song level
tutorSongs = plotAndAdjust(tutorStruct,[],wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8,...
    'optGraphs',{'waveform','spectrogram'});

%% (3.5) characterize the tutor syllables / get features
[tSFeats,tSSpecs]=getFeatures(tutorStruct,tutorSylls,tutorParams,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',false); 

%% save the files pertaining to tutor song
uisave({'tutorStruct','tutorSylls','tSFeats','tSSpecs'},...
    [matpath prependForSave('tutor-',matFile)]);

%% (3.7 optional) REVIEW tutor syllables
wholeRegion = getWholeSongAsRegion(tutorStruct);
plotForReview(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'dgram.minContrast',1e-8,...
    'optGraphs',{'waveform', 'fundamentalFreq', 'AM', 'FM', 'wienerEntropy', 'pitchGoodness'});
%% (3.7 optional) REVIEW tutor syllables with spectrogram
currFigure = gcf;
figure
wholeRegion = getWholeSongAsRegion(tutorStruct);
plotForReview(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'dgram.minContrast',1e-9,...
    'optGraphs',{'waveform', 'deriv'});
figure(currFigure)
%% (3.8) load other tutor songs and check similarity metrics
% example
tutorParams = params;
tutorParams.fine.fs = 1/tutorStruct.interval; tutorParams.fs = 1/tutorStruct.interval;
tutorSpec = getMTSpectrumStats(tutorStruct.values, tutorParams.fine);
figure(1);
plotAllFigures(tutorSpec, tutorSylls, tutorParams, 'showLabels',true,'dgram.minContrast',1e-10,...
    'optGraphs',{'waveform', 'fundamentalFreq', 'AM', 'FM', 'wienerEntropy', 'pitchGoodness'});

for ii = 1:numel(tutorFiles)
    otherTutorSong(ii) = loadWavFile(tutorFiles{ii+1});
    wholeTutorRegion(ii) = getWholeSongAsRegion(otherTutorSong(ii));
    % pre highPass and amplify
    otherTutorSong(ii).values = 3*highPassSample(getClip(wholeTutorRegion(ii), otherTutorSong(ii)), params);
    % automatic parse
    otherTutorSylls{ii} = parseRegionsIntoSyllables(otherTutorSong(ii), wholeTutorRegion(ii),...
        params,'doFilterNoise',false,'preroll',0,'postroll',0,...
        'syllable.comboLength',0,'syllable.flatFactor',20);
    figure(2)
    otherTutorSylls{ii} = plotAndAdjust(otherTutorSong(ii), otherTutorSylls{ii}, wholeTutorRegion(ii), params,...
        'editSpecType', 'fine', 'adjustLabels',false,'dgram.minContrast',3e-10,...
    'optGraphs',{'waveform', 'fundamentalFreq', 'AM', 'FM', 'wienerEntropy', 'pitchGoodness'});
    [otSFeats, otSSpecs] = getFeatures(otherTutorSong(ii), ...
        otherTutorSylls{ii}, params, 'doFilterNoise', false);
    
    weightTables = cell(numel(otherTutorSylls{ii}), numel(tutorSylls));
    distVector = cell(numel(otherTutorSylls{ii}), numel(tutorSylls));
    for jj = 1:numel(otherTutorSylls{ii})
        for kk = 1:numel(tutorSylls)
            [distVector{jj,kk} weightTables{jj,kk}]  = standardDistance(otSSpecs(jj),tSSpecs(kk),params);
        end
    end
    
    distScore = cellfun(@(x) mean(x), distVector);
    [bestScore,bestMatch] = min(distScore ,[],2);
    tutorTypes = {tutorSylls.type};
    guessTypes = tutorTypes(bestMatch);
    [otherTutorSylls{ii}.type] = guessTypes{:};
    
    newLabeledSylls{ii} = plotAndAdjust(otherTutorSong(ii), otherTutorSylls{ii}, wholeTutorRegion(ii), params,...
        'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',3e-10,...
    'optGraphs',{'waveform', 'fundamentalFreq', 'AM', 'FM', 'wienerEntropy', 'pitchGoodness'});
    
    % TODO: change weights to features to make similarity score better
    % also investigate more pitch envelope-type techniques
    % manual check
    
    perfOfGT = [];
    perfOfSim = [];
    for jj = 1:numel(otherTutorSylls{ii})
        GTMatch = find(strcmp(tutorTypes,newLabeledSylls{ii}(jj).type),1);    
%        isSimMatch = strcmp(tutorTypes,otherTutorSylls{ii}(jj).type);
        
        if size(weightTables{jj,GTMatch},1) > 1
            matchComparisons = mean(weightTables{jj,GTMatch},1);
        else
            matchComparisons = weightTables{jj,GTMatch};
        end
        %nonmatchComparisons = vertcat(ans{:});
        nonmatchComparisons = vertcat(weightTables{jj,1:numel(tutorSylls) ~= GTMatch});
        
        if size(weightTables{jj,GTMatch},1) > 1
            simMatchComparisons = mean(weightTables{jj,bestMatch(jj)},1);
        else
            simMatchComparisons = weightTables{jj,GTMatch};
        end
        nonsimMatchComparisons = vertcat(weightTables{jj,1:numel(tutorSylls) ~= bestMatch(jj)});
        perfOfGT = vertcat(perfOfGT, mean(matchComparisons(ones(1,size(nonmatchComparisons,1)),:) < nonmatchComparisons));
        perfOfSim = vertcat(perfOfSim, mean(simMatchComparisons(ones(1,size(nonsimMatchComparisons,1)),:) < nonsimMatchComparisons));
    end
end

%% (4.a) calculate features of syllables
% NOTE: there are two types of similarity: sequence-independent and
% sequence-dependent.  
% Sequence-independent similarity summarizes each syllable as a set of
% central tendencies

ROIs = manualSyllables;
nROI = numel(ROIs);

stdDistance = zeros(nROI, numel(tutorSylls));

progressbar('extracting features');
for ii = 1:nROI
    [jSFeats(ii),jSSpec]=getFeatures(songStruct,ROIs(ii),params,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',true, ...
        'noiseFilter', noiseMask,'nps.reduction',-12);   
    progressbar(ii/nROI);
end

% dimensionality fixing
% current policy: z-normalize it all
featNames = fieldnames(jSFeats);
muFeat = initEmptyStructArray(featNames,1);
sigmaFeat = initEmptyStructArray(featNames,1);
progressbar('z-normalizing');
for jj = 1:numel(featNames)
    thisFeat = featNames{jj};
    [znormed, muFeat.(thisFeat), sigmaFeat.(thisFeat)] = zscore([jSFeats.(thisFeat)]');

    % overwrite
    foo = num2cell(znormed); [jSFeats.(thisFeat)] = foo{:};
    progressbar(jj/numel(featNames));
end

% policy: try to remove dominant features
stable = false;
while ~stable
    jsFVs = cellfun(@(x) [jSFeats.(x)], featNames, 'UniformOutput', false);
    jsFVs = vertcat(jsFVs{:});
    tsFV = zeros(numel(featNames), numel(tutorSylls));
    
    progressbar('sequence independent syllable comparison','individual syllables');
    for jj = 1:numel(tutorSylls)
        tsFV(:,jj) = cellfun(@(x) (tSFeats(jj).(x) - muFeat.(x)) / sigmaFeat.(x), featNames);
    end
    
    for ii = 1:nROI
        for jj = 1:numel(tutorSylls);
            % sequence independent distance
            % distance metric: simple L2 for now
            stdDistance(ii, jj) = norm(jsFVs(:,ii) - tsFV(:,jj));
            progressbar([],jj/numel(tutorSylls));
        end
        progressbar(ii/nROI,[])
    end
    
    % outlier test: ratio of variances - the distances are often skewed by
    % one or two most variable features, so we try to minimize this
    % somewhat by removing features (we can also scale them)
    if mean(std(stdDistance,[],2)) / mean(std(stdDistance,[],1)) > 100
        % take out any of the variables which have abnormally high contribution
        % to distance
        featContrib = zeros(numel(featNames), 1);
        for ii = 1:numel(featNames)
            propFeatDist = abs(jsFVs(ii,:)' * ones(1,numel(tutorSylls)) - ones(nROI,1) * tsFV(ii,:)) ...
                ./ stdDistance;
            featContrib(ii) = mean(propFeatDist(:));
        end
        [maxContrib,iMax] = max(featContrib);
        if maxContrib > 0.4
            featNames(iMax) = [];
        else
            stable = true;
        end
    else
        stable = true;
    end
end
%% (4.b) calculate similarity of juvenile sylls w/ true labels to tutor syllables
% using a sequence dependent distance

% index 1: ROI, index 2: tutor syllable type, index 3: type of
% reduction/syllable instance?
syllDistances = zeros(nROI, numel(tutorSylls),numel(reductions));
        
progressbar('getting juvie syllable similarity', 'each tutor syllable');
reductions = {@nanmean, @max, @min, @median, @geomean, @(x) sqrt(nanmean(power(x,2)))};
for ii = 1:nROI
    [jSFeats(ii),jSSpec]=getFeatures(songStruct,ROIs(ii),params,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',true, ...
        'noiseFilter', noiseMask,'nps.reduction',-12);   

    for jj = 1:numel(tutorSylls)
        % this correlation vector is the score for each possible alignment
        % of two syllables of different lengths
        distVector = standardDistance(jSSpec, tSSpecs(jj),params);
        
        % these are various reductions 
        for kk = 1:numel(reductions) %index 
            syllDistances(ii,jj,kk) = reductions{kk}(distVector);
        end
        progressbar([],jj/numel(tutorSylls));   
    end
    progressbar(ii/numel(ROIs),[]);   
end

% take the minimum for the standard accepted score for now
stdDistance = syllDistances(:,:,3);

% save work on tutor similarity
%%
uisave({'stdDistance','tutorStruct','tutorSylls','manualMotifs','manualSyllables','tSFeats','tSSpecs'},...
    prependForSave('tutorSimilarity-',matFile));

%% (4, alternate) calculate dynamic time-warped similarity of juvenile sylls w/ true labels to tutor syllables
progressbar('getting juvie syllable similarity', 'each tutor syllable');

ROIs = manualSyllables;
nF = 5; % number of features
syllFeatDistances = zeros(numel(ROIs), numel(tutorSylls),nF);
syllDistances = zeros(numel(ROIs), numel(tutorSylls));
for ii = 1:numel(ROIs)
    [jSFeats(ii),jSSpec]=getFeatures(songStruct,ROIs(ii),params,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',true, ...
        'noiseFilter', noiseMask,'nps.reduction',-12); 
    
    for jj = 1:numel(tutorSylls)
        % this correlation vector is the score for each possible alignment
        % of two syllables of different lengths
        [syllDistances(ii, jj), syllFeatDistances(ii,jj,:)] = timeWarpedDistance(jSSpec, tSSpecs(jj),params);
        progressbar([],jj/numel(tutorSylls));
    end
    progressbar(ii/numel(ROIs),[]);   
end

% save work on tutor similarity

uisave({'syllFeatDistances','syllDistances','tutorStruct','tutorSylls','manualMotifs','manualSyllables','jSFeats'},...
    prependForSave('TWtutorSimilarity-',matFile));
%% (optional) ASSIGN LABELS AUTOMATICALLY to each juvenile syllable according to its similarity 

% z-standardize each score
% for jj = 1:numel(tutorSylls)
%     stdScore(:,jj) = (stdScore(:,jj) - mean(stdScore(:,jj)))/std(stdScore(:,jj));
% end
% score is the mean across each syllable

% map scores to unique labels
tutorLabels = {tutorSylls.type};
[uLabels, foo, labelsrIdx] = unique(tutorLabels);
nULabels = numel(uLabels);
labelDist = zeros(numel(manualSyllables), nULabels);
for ii = 1:nULabels
%    labelDist(:,ii) = mean(stdDist(:,labelsrIdx == ii),2);
% average over (in time warped case)
    labelDist(:,ii) = mean(syllDistances(:,labelsrIdx == ii),2);
end
[minDistance, bestIdx] = min(labelDist,[],2);
guessedLabels = tutorLabels(bestIdx);

[manualSyllables.type] = guessedLabels{:};

%% (4.5) optional, REVIEW juvenile syllable labels
ROIs = manualMotifs;
subROIs = syllables;
for ii = 1:numel(ROIs)
    plotForReview(songStruct,...
        subROIs,...
        ROIs(ii), ...
        params, ...
        'editSpecType', 'ultrafine', 'dgram.minContrast',3e-10,...
        'doFilterNoise',false, 'noiseFilter', noiseMask,'nps.reduction',-12, ...
        'optGraphs',{'waveform', 'deriv','fundamentalFreq','wienerEntropy'});
    pause;
   % newJuvenileSylls = [newJuvenileSylls revisedSylls]; 
end

%% (optional) REFINE labels with MANUAL work 
% do a manual refinement
%tutorParams = params;
%tutorParams.fine.fs = 1/tutorSong.interval; tutorParams.fs = 1/tutorSong.interval;
%tutorSpec = getMTSpectrumStats(tutorSong.values, tutorParams.fine);
%figure(1);
%plotAllFigures(tutorSpec, tutorSylls, tutorParams, 'showLabels',true);
newJuvenileSylls = initEvents;
for ii = 1:numel(manualMotifs)
    figure(2);
    revisedSylls = plotAndAdjust(songStruct,juvenileSylls,manualMotifs(ii), params, ...
        'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',3e-10,...
        'doFilterNoise',true, 'noiseFilter', noiseMask,'nps.reduction',-12);
    newJuvenileSylls = [revisedSylls newJuvenileSylls]; 
end
% jenny added this here for after manual stuff
uisave({'markedSongs', 'juvenileSylls'}, [matpath prependForSave('markedSongs-',matFile)]);
%juvenileSylls = newJuvenileSylls;
%clear newJuvenileSylls;
%% (optional) CHECK SIMILARITY for a given syllable
selType = input('Which syllable do you want to check? > ','s');
selFilter = strcmp({juvenileSylls.type}, selType);
whichTutorSyll = find(strcmp({tutorSylls.type},selType));
distances = stdDistance(selFilter, whichTutorSyll);
selectedSylls = assignAsField(juvenileSylls(selFilter), 'distance', distances);
selectedSylls = sortBy(selectedSylls, 'distance');
fprintf('Playing original several times...');
clip = getClip(tutorSylls(whichTutorSyll), tutorStruct);
for ii = 1:10
    pause(0.1);
    playSound(clip, 1/tutorStruct.interval,true);
end
fprintf('Done.\n');

goodInput = true;
while goodInput
    sel = input('Which clip do you want to hear? > ', 's');
    sel = str2double(sel);
    goodInput = ~isnan(sel) || sel >= 0;  
    if ~goodInput, continue; end;
    sel = floor(sel);
    if sel == 0
        fprintf('Playing original tutor clip ...');
        clip = getClip(tutorSylls(whichTutorSyll), tutorStruct);
        playSound(clip, 1/tutorStruct.interval,true);
    else
        fprintf('Playing clip %d/%d (distance %f)...', ...
            sel, numel(selectedSylls), selectedSylls(sel).distance);
        clip = getClip(selectedSylls(sel), songStruct);
        playSound(clip, 1/songStruct.interval,true);
    end
    fprintf('Done\n');
end


%% (5) LOAD the SPIKE-SORTED-FILE and calculate firing rate during syllables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% All processing up to this point has been solely on the song file 
%%%% To continue you should have: 
%%%% (a) some measure of distance between syllables and different tutor
%%%% syllables (which we call stdDist),  
%%%% (b) definitions of the tutor syllables (tutorSylls), and
%%%% (c) juvenileSylls with defined boundaries
%%%% AND types that correspond to the tutor syllables' types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is where we LOAD the SPIKE-SORTED-FILE
[matFiles, matSpikePath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file(s)','data', 'MultiSelect', 'on');
spikes = loadSpikeData(strcat(matSpikePath, matFiles));
fprintf('%d neurons loaded...\n', numel(spikes));
% calculate spike data w.r.t segmented regions
%% assess spike counts/rates
ROIs = manualSyllables;
syllLengths = [ROIs.stop] - [ROIs.start];
for ii = 1:numel(spikes)
    [spikeCounts{ii}, spikeInternalTimes{ii}] = countSpikes(ROIs, spikes{ii},'onset');
    spikeRates{ii} = spikeCounts{ii} ./ syllLengths; 
end
totalSpikeRates = sum(cellfun(@sum,spikeRates));
%spikeTimes = sort(vertcat(spikes{:}));

%% (5.x) test if spiking is song related?
neuronData = spikingByCategory(songStruct, manualMotifs, spikes);

%% measure interspike intervals

%% (6, revised) test the syllable-wise correlation between similarity to tutor and firing rate

% get juvenile syllable type
ROIs = manualSyllables;
nR = numel(ROIs);
% syllDistances should be the distance between each ROI and each
% tutorSyllable

% labelDist becomes the distance between each ROI and each
% unique tutorSyllable type
tutorLabels = {tutorSylls.type};
[uLabels, foo, labelsrIdx] = unique(tutorLabels);
nULabels = numel(uLabels);
labelDist = zeros(nR, nULabels);
for ii = 1:nULabels
    % problem: if there's only one column for which labelsrIdx==ii, then it won't work
    isTutorSyll = (labelsrIdx == ii);
    if sum(isTutorSyll) > 1
        labelDist(:,ii) = mean(stdDistance(:,isTutorSyll),2);
    else
        labelDist(:,ii) = stdDistance(:,isTutorSyll);
    end
end
simNames = strcat('sim_',uLabels);
isNames = strcat('isa_', uLabels);
statFields = [simNames isNames];
syllableStats = initEmptyStructArray(statFields,nR);
for ii = 1:nULabels
    foo = num2cell(labelDist(:,ii));
    [syllableStats.(simNames{ii})] = foo{:};
    foo = num2cell(strcmp(uLabels(ii),{ROIs.type}));
    [syllableStats.(isNames{ii})] = foo{:}; 
end

nNeurons = numel(spikes);
corrSimToNeuron = zeros(nULabels * 2, nNeurons);
for ii = 1:nNeurons
    corrSimToNeuron(:,ii) = findSpikeCorrelates(ROIs, syllableStats,spikes{ii});
    pause;
    pause(0.5);
end
figure
hTable = uitable('Data', corrSimToNeuron', 'ColumnName', strcat('pval_',statFields));
set(hTable, 'Position',get(hTable,'Extent'));
%% (6.5) plot PSTH for each syllable type (TODO: plot with baseline)

% get juvenile syllable type, which may be clustered
juvieTypes = {juvenileSylls.type};
[uJTypes, foo, rJIdx] = unique(juvieTypes);
psthPreParams = processArgs(params, 'preroll', 20, 'postroll', 20);

% get tutor syllable type
tutorLabels = {tutorSylls.type};
[labelTypes, foo, rTutorIdx] = unique(tutorLabels);
nLabelTypes = numel(labelTypes);

% set display stuff
sizeFromLength = ceil((syllLengths - min(syllLengths)) / (max(syllLengths) - min(syllLengths)) * 40 + 5);
% set color map
colors = hsv(nLabelTypes);

for ii = 1:numel(uJTypes)
    thisType = (rJIdx == ii);
    extendedTypeSylls = addPrePost(juvenileSylls(thisType), psthPreParams);
    num2cell([extendedTypeSylls.stop] - [extendedTypeSylls.start]);
    [extendedTypeSylls.length] = ans{:};
    
    plotPSTH(sortBy(extendedTypeSylls,'length'),vertcat(spikes{:}),juvenileSylls(thisType));
    title(sprintf('PSTH for syllable %s (%d)', uJTypes{ii}, sum(thisType)));
    
    %figure(2)
   % 
   % theseDistances = vertcat(juvenileSylls(thisType).tutorSim);
   % [foo, bestSyll] = min(theseDistances, [], 2); repSyll = mode(bestSyll);
    %for tutorSyllSel = 1:nLabelTypes;
    %     scatter(theseDistances(:,repSyll), totalSpikeRates(thisType), ...
    %        sizeFromLength(thisType), colors(repSyll,:),'filled');
    %    hold on;
    
    %end
    pause
end
%% (6.8) plot PSTH for all marked songs, no syllables 

ROIs = manualMotifs(strcmp('b',{manualMotifs.type}));
num2cell([ROIs.stop] - [ROIs.start]);
[ROIs.length] = ans{:};
[sortedROIs, indexByLength] = sortBy(ROIs, 'length');
prePostParams = processArgs(params,'preroll',2000,'postroll',2000);
for ii = 1:numel(spikes)
    figure(ii)
     plotPSTH(addPrePost(sortedROIs, prePostParams)...
         , spikes{ii}, sortedROIs ...
         );%... 
         %juvenileSylls,{'a'}); % can change the spikes
%     plotAlignedPSTH(newMarkedSortedSongs,...
%         spikes{ii},...
%         juvenileSylls, {'b','c','d','e','f'}); % can change the spikes
%     title(sprintf('neuron %d',ii));
    pause
end
%% (6.8) plot PSTH for all marked songs, syllables marked

ROIs = manualMotifs(strcmp('juvenile',{manualMotifs.type}));
subROIs = syllables;
neuronList = spikes(highJuv);
num2cell([ROIs.stop] - [ROIs.start]);
[ROIs.length] = ans{:};
[sortedROIs, indexByLength] = sortBy(ROIs, 'length');
prePostParams = processArgs(params,'preroll',500,'postroll',500);
for ii = 1:numel(neuronList)
    figure(ii)
     plotPSTH(ROIs,...
         neuronList{ii},...
         subROIs,  ...
     {}, ... %alignment
     prePostParams); % can change the spikes
     title(sprintf('neuron %d',ii));
    pause
end
%% (6.9) what is similarity for spiking vs non spiking songs?

% flag bouts by spiking / non spiking
spikesInSongs = countSpikes(manualMotifs, vertcat(spikes{:}));

% sum up distance for each song according to its syllable distances, 
% not according to song difference
tutorTypes = unique({tutorSylls.type});
rSyllIdx = zeros(1,numel(juvenileSylls));
for ii = 1:numel(tutorTypes)
    rSyllIdx(strcmp({juvenileSylls.type}, tutorTypes(ii))) = ii;
end

songScore = zeros(1,numel(manualMotifs));
for ii = 1:numel(manualMotifs)
    [subSylls,ind] = getSubEvents(manualMotifs(ii), juvenileSylls);
    % sum all syllable score differences
    songScore(ii) = mean(stdDistance(sub2ind(size(stdDistance),ind, rSyllIdx(ind))));
end

% find alternate distance according to song / song comparison, across the
% same features

% plot nonspiking / spiking on X and distance on Y with error
uNspikes = unique(spikesInSongs);
for ii = 1:numel(uNspikes)
    yy(ii) = nanmean(songScore(spikesInSongs == uNspikes(ii)));
    ey(ii) = nanstd(songScore(spikesInSongs == uNspikes(ii)));
end

subplot(211)
errorbar(uNspikes, yy, ey, 'b.')
xlabel('Number of spikes')
ylabel('Mean Similarity to best syllable in each syllable');
subplot(212)
hist(spikesInSongs,min(uNspikes):max(uNspikes));
xlabel('Number of spikes')
ylabel('Count');

%% bar chart of spikes for each syllable
tutorTypes = {tutorSylls.type};
rSyllIdx = zeros(1,numel(juvenileSylls));
averageSpikesPerType = zeros(1,numel(tutorTypes));
semSpikesPerType = zeros(1,numel(tutorTypes));

for ii = 1:numel(tutorTypes)
    subSel = strcmp({juvenileSylls.type}, tutorTypes(ii));
    rSyllIdx(subSel) = ii;
    for jj = 1:numel(spikeCounts)
        averageSpikesPerType(ii,jj) = mean(spikeRates{jj}(subSel));
        semSpikesPerType(ii,jj) = std(spikeRates{jj}(subSel))/sqrt(sum(subSel));
    end
end
subplot(111)
barweb(averageSpikesPerType',semSpikesPerType');
xlabel('neurons')
ylabel('spike rates (Hz)')
legend(strcat('Syllable:', tutorTypes))
% split with low distance / hi distance 

% plot low distance / hi distance on X and spiking on Y, baseline as
% sidebar

% what are the structures we need for correlating firing rates within syllable?
% syllable boundaries, feature reductions, and spikes
% what do we need for correlating fireing rates within a song? 
% song boundaries, (subsyllables, their sub features), and spikes

%% evaluate correlates to song in syllable-independent fashion
tutorFrame = addPrePost(bookendedClip(tutorSylls), params);
tutorState = stateVector(tutorFrame,tutorSylls);
silDistance = zeros(1,numel(manualMotifs));

% similarity in terms of spacing (hamming distance)
parentage = findParent(addPrePost(manualMotifs, params), manualSyllables);
for ii = 1:numel(manualMotifs)
    silDistance(ii) = minEditDistance(...
        stateVector(addPrePost(manualMotifs(ii), params), ...
        manualSyllables(parentage==ii)), ...
        tutorState);
end

% similarity to length included but never shows correlation
songScores = struct('length',num2cell([manualMotifs.stop] - [manualMotifs.start]),...
    'silDistance', num2cell(silDistance));
for ii = 1:numel(spikes)
    [basalRate(ii), spikeRate(ii), rateSEM(ii), rate_tp(ii)] = findSpikeRates(manualMotifs, spikes{ii});
    silenceEditTp = findSpikeCorrelates(manualMotifs, songScores,spikes{ii});
end
%uisave({'basalRate', 'spikeRate', 'rateSEM', 'rate_tp'},prependForSave('spikeRates-', matFile));
