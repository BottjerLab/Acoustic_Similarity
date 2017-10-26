% JUVENILECORRELATION run this to correlate juvenile-tutor similarity with
% firing rate
% cell by cell if you're patient (Ctrl-Enter on Windows)

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

cd Z:\john\parsing
startup

%% Note 2: you can always uiload any of the files that are "uisaved"
uiload
% For silly reasons, the song file (without the prefix) should be loaded
% with the cell below.
%% prep step: load song file
% define files to be operated over
[matFile, matpath] = uigetfile('*.mat','Please choose the song Spike2 file','data');
% load file
songStruct = load([matpath matFile]);
% trick to get the main struct into a standard name, if there's only one
% variable in the file
fld=fieldnames(songStruct);
songStruct=songStruct.(fld{1});
fs = 1/songStruct.interval;

%% prep step: set parameters from an alias
params = processArgs(defaultParams,'R204');

%% (0.5) isolate sounds first
sounds = stepSpectrogram(songStruct, params,'Nsplits',400,'plot',true);
%% arrange into bouts
[bouts, nSyllsInBout] = arrangeBouts(sounds,params,'silentPeriod',0.5);
%% refine sound boundaries
%[candidateSyl, features, featureReduction] = ...
candidateSyl = ...
    parseRegionsIntoSyllables(songStruct, bouts(15:end), params, 'plot', false, 'quiet',true,'pause', false,...
'syllable.minPower',5e-7,'syllable.borderRise',8e-6, 'noiseFilter', noiseMask, 'dgram.minContrast', 1e-12); 

% TODO: clean out obvious spikes/wingflaps?

%% arrange into bouts
[bouts, nSyllsInBout] = arrangeBouts(candidateSyl,params,'silentPeriod',0.5);

% do a quick filter: bouts shorter than 200 ms probably do not have bouts
minBoutLength = 0.2; % in seconds;
bouts(([bouts.stop]-[bouts.start]) < minBoutLength) = [];
%% save some work - automatically found partitions
uisave({'sounds', 'candidateSyl', 'bouts'},...
    [matpath prependForSave('divMT-',matFile)]);

%% (x.8) find noise 
% define some noise with space between bouts
notNoise = bouts;
approved = false;
while ~approved
    candidateNoise = autodetectNoise(songStruct, notNoise);
    fprintf('Is this clip pure noise? Mark if yes...\n'); 
    approved = markRegions(songStruct,candidateNoise);
    % just black this section out and try again
    if ~approved
        notNoise = [notNoise; candidateNoise];
    end
end

noiseMask = noiseAnalysis(songStruct, candidateNoise);

% save noise mask
uisave('noiseMask',[matpath prependForSave('noiseMask-', matFile)]);
%%
params.noiseFilter = noiseMask;
%% (1) isolate the songs
disp('Manual adjustment to find bouts...');

%parentage = findParent(bouts, candidateSyl);
containsSong = false(1,numel(bouts));
num2cell([bouts.stop] - [bouts.start]); [bouts.length]=ans{:};
    
sortedBouts = sortBy(bouts,'length','descend');
% don't look at bouts that are too short - those can't contain any song
%if boutLen < 0.1, continue; end;
   
fprintf('Do these clips contain song? Mark if yes...\n');
containsSong = markRegions(songStruct, sortedBouts);

manualBouts = sortedBouts(containsSong);   
%% save first pass on marked songs
uisave('manualBouts',[matpath prependForSave('manualBouts-',matFile)]);

%% (1.5 auto) refine/subparse song boundaries and chop 'songs' which are here bouts 
% into real songs 

% TODO: include learner training in here
prePostParams = processArgs(params,'preroll',500,'postroll',500);

newMarkedSongs = initEvents(0);

sortedMarkedSongs = markedBouts;
num2cell([sortedMarkedSongs.stop] - [sortedMarkedSongs.start]); 
[sortedMarkedSongs.length] = ans{:};
sortedMarkedSongs = sortBy(sortedMarkedSongs,'length','descend');
% initialize learner
%songLearner = Learner(songStruct, [], params);

ROIs = sortedMarkedSongs;
for ii = 2:numel(ROIs)
    thisROI = ROIs(ii);
    fprintf('Song %d/%d (length %0.2f s)... ', ii, numel(ROIs), thisROI.stop - thisROI.start);
    if songLearner.nExamples < 5 % minimum number of songs for robust classification
        % keep manually classifying until we get 5 songs
        tmpSongs = plotAndAdjust(songStruct, thisROI, ...
        addPrePost(thisROI, prePostParams), params,...
        'editSpecType','inter',...
        'inter.features', {'totalPower','wienerEntropy','deriv','waveform'});
        if ~isempty(tmpSongs)
            songLearner.addSpectra(songStruct, tmpSongs);
            songLearner.adaptThreshold;
        end
    else
        tic
  %      fprintf('Song %d/%d (length %0.2f s)... ', ii, numel(markedSongs), markedSongs(ii).stop - markedSongs(ii).start);
        tmpSongs = scanForMatches(songStruct, addPrePost(thisROI, prePostParams), songLearner);
        toc
    end
    newMarkedSongs = [newMarkedSongs tmpSongs];
end

%% 
uisave('songLearner',[matpath prependForSave('songTemplate-',matFile)]);

%% (1.5 manually) refine/subparse song boundaries MANUALLY, 
%extending if necessary (define range below)
prePostParams = processArgs(params,'preroll',500,'postroll',500);
newROIs = initEvents(0);

% initialize learner
ROIs = sortedMarkedSongs;
songLearner = Learner(songStruct, [], params);
nB = numel(ROIs);
for ii = 1:nB %NOTE: define range here;
    thisBout = ROIs(ii);
    fprintf('Bout %d/%d (length %0.2f s)... ', ii, nB, thisBout.stop - thisBout.start);
    if strcmp('L',thisBout.type)
        tmpSongs = thisBout; 
        fprintf('locked, continuing...\n');
    else
        contextBout = addPrePost(thisBout, prePostParams);
        tmpSongs = plotAndAdjust(songStruct, thisBout, ...
            contextBout, params,...
            'editSpecType','inter', 'adjustLabels', true,...
            'inter.features', {'totalPower','wienerEntropy','deriv','waveform'});
        %'noiseFilter', noiseMask);
    end
    newROIs = [newROIs; tmpSongs];
end
%%
manualMotifs = newROIs;

%%
uisave({'markedMotifs','newMarkedSongs'}, [matpath prependForSave('manualMarkedSongs-',matFile)]);

%%
uisave({'manualMotifs'}, [matpath prependForSave('motifs-',matFile)]);

%% temporary
for ii = 1:size(olaps,1); fprintf('%Olap #d: \n',ii); 
    cl1 = getClipAndProcess(songStruct, newMarkedSongs(olaps(ii,1)),defaultParams); 
    cl2 = getClipAndProcess(songStruct, newMarkedSongs(olaps(ii,2)),defaultParams);
    clo = getClipAndProcess(songStruct, olapEvents(ii),defaultParams); 
    playSound(cl1, fs); playSound(0.3,fs); playSound(cl2, fs); playSound(0.3,fs); playSound(clo, fs);
end;

%% auto correction: merge all overlapping segments
olaps = findOverlaps(newMarkedSongs);
adjustedMarkedSongs = newMarkedSongs; 
markDel = false(numel(newMarkedSongs),1);
for ii = 1:size(olaps,1)
    adjustedMarkedSongs(olaps(ii,1)).idxStart = newMarkedSongs(olaps(ii,2)).idxStart; 
    adjustedMarkedSongs(olaps(ii,1)).start = newMarkedSongs(olaps(ii,2)).start; 
    markDel(olaps(ii,1)) = true;
end
newMarkedSongs = adjustedMarkedSongs(~markDel); 

%% (1.6) auto correction: where there is an overlap, give overlapped portion to the later segment
olaps = findOverlaps(newMarkedSongs);
adjustedMarkedSongs = newMarkedSongs; 
for ii = 1:numel(olaps)/2
    adjustedMarkedSongs(olaps(ii,1)).idxStop = newMarkedSongs(olaps(ii,2)).idxStart; 
    adjustedMarkedSongs(olaps(ii,1)).stop = newMarkedSongs(olaps(ii,2)).start; 
end
% verification
%{
for ii = 1:size(olaps,1)
    fprintf('>>>>>>>>>>>>>>Olap #%d: \n',ii);
    cl1 = getClipAndProcess(songStruct, newMarkedSongs(olaps(ii,1)),defaultParams); 
    cl2 = getClipAndProcess(songStruct, adjustedMarkedSongs(olaps(ii,1)),defaultParams);
    newMarkedSongs(olaps(ii,1))
    playSound(cl1, fs, true); pause(0.5); 
    adjustedMarkedSongs(olaps(ii,1))
    playSound(cl2, fs, true);
end
%}
newMarkedSongs = adjustedMarkedSongs; 
%% (2) parse juvenile song into syllables
ROIs = manualBouts;
params = processArgs(defaultParams,'fs',1/songStruct.interval, 'preroll', 50);
syllables = parseRegionsIntoSyllables(songStruct, ROIs, params, ...
    'noiseFilter', noiseMask,'nps.reduction',-12,'plot',false, 'dgram.minContrast',1e-11, ...
    'syllable.comboLength',3,... % gap size in ms
    'syllable.borderRise',1.2e-4,'syllable.minPower',7e-6,'pause',false);
%%
uisave({'manualMotifs', 'syllables'}, [matpath prependForSave('syllables-',matFile)]);
%% (2.5) refine parse with manual work
% do a manual refinement
% tutorParams = params;
% tutorParams.fine.fs = 1/tutorSong.interval; tutorParams.fs = 1/tutorSong.interval;
% tutorSpec = getMTSpectrumStats(tutorSong.values, tutorParams.fine);
% figure(1);
%plotAllFigures(tutorSpec, tutorSylls, tutorParams, 'showLabels',true);
params = processArgs(defaultParams,'fs',1/songStruct.interval);

ROIs = manualMotifs;
subROIs = syllables;
% resets new results
tmpSyllables = initEvents;
for ii = 1:numel(ROIs)
    thisEv = ROIs(ii);
    fprintf('Editing %d/%d...\n', ii, numel(ROIs));
    revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
        'editSpecType', 'inter', 'adjustLabels',true, 'dgram.minContrast',3e-12,...
        'doFilterNoise',false, 'noiseFilter', [],'nps.reduction',-12);
    tmpSyllables = [tmpSyllables revisedSylls']; 
end

% moves new results to permanent location
syllables = tmpSyllables;
clear tmpSyllables

%% (3-eps) load up tutor songs... 
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
    'doFilterNoise',false,'syllable.comboLength',18);

%% do a manual refinement on syllable level
tutorSylls = plotAndAdjust(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8,...
    'optGraphs',{'waveform', 'deriv','FM','AM', 'mTD'});

%% do a manual refinement on song level
tutorSongs = plotAndAdjust(tutorStruct,[],wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8,...
    'optGraphs',{'waveform', 'deriv','FM','AM', 'mTD'});

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
wholeRegion = getWholeSongAsRegion(tutorStruct);
plotForReview(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'dgram.minContrast',1e-9,...
    'optGraphs',{'waveform', 'deriv'});

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



%% (4) calculate similarity of juvenile sylls w/ true labels to tutor syllables
progressbar('getting juvie syllable similarity', 'each tutor syllable');

reductions = {@nanmean, @max, @min, @median, @geomean, @(x) sqrt(nanmean(power(x,2)))};
syllDistances = zeros(numel(juvenileSylls), numel(tutorSylls),numel(reductions));
for ii = 1:numel(juvenileSylls)
    [jSFeats(ii),jSSpec]=getFeatures(songStruct,juvenileSylls(ii),params,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',true, ...
        'noiseFilter', noiseMask,'nps.reduction',-12); 
    
    for jj = 1:numel(tutorSylls)
        % this correlation vector is the score for each possible alignment
        % of two syllables of different lengths
        distVector = standardDistance(jSSpec, tSSpecs(jj),params);
        
        % these are various reductions 
        
        for kk = 1:numel(reductions)
            syllDistances(ii,jj,kk) = reductions{kk}(distVector);
        end
        progressbar([],jj/numel(tutorSylls));   
    end
    progressbar(ii/numel(juvenileSylls),[]);   
end

% take the minimum for the standard accepted score for now
stdDist = syllDistances(:,:,3);

% save work on tutor similarity
%%
uisave({'syllDistances', 'stdDist','tutorStruct','tutorSylls','manualMotifs','juvenileSylls','tSFeats','tSSpecs'},...
    prependForSave('tutorSimilarity-',matFile));

%% (4, alternate) calculate dynamic time-warped similarity of juvenile sylls w/ true labels to tutor syllables
progressbar('getting juvie syllable similarity', 'each tutor syllable');

nF = 5; % number of features
syllFeatDistances = zeros(numel(juvenileSylls), numel(tutorSylls),nF);
syllDistances = zeros(numel(juvenileSylls), numel(tutorSylls));
for ii = 1:numel(juvenileSylls)
    [jSFeats(ii),jSSpec]=getFeatures(songStruct,juvenileSylls(ii),params,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',true, ...
        'noiseFilter', noiseMask,'nps.reduction',-12); 
    
    for jj = 1:numel(tutorSylls)
        % this correlation vector is the score for each possible alignment
        % of two syllables of different lengths
        [syllDistances(ii, jj), syllFeatDistances(ii,jj,:)] = timeWarpedDistance(jSSpec, tSSpecs(jj),params);
        progressbar([],jj/numel(tutorSylls));
    end
    progressbar(ii/numel(juvenileSylls),[]);   
end

% take the minimum for the standard accepted score for now
% save work on tutor similarity
uisave({'syllFeatDistances','syllDistances','tutorSong','tutorSylls','manualMotifs','juvenileSylls','tSFeats','tSSpecs'},...
    prependForSave('TWtutorSimilarity-',matFile));
%% (optional) ASSIGN LABELS AUTOMATICALLY to each juvenile syllable according to its similarity 

% % z-standardize each score
% for jj = 1:numel(tutorSylls)
%     stdScore(:,jj) = (stdScore(:,jj) - mean(stdScore(:,jj)))/std(stdScore(:,jj));
% end
% map scores to unique labels
tutorLabels = {tutorSylls.type};
[uLabels, foo, labelsrIdx] = unique(tutorLabels);
nULabels = numel(uLabels);
labelDist = zeros(numel(juvenileSylls), nULabels);
for ii = 1:nULabels
    labelDist(:,ii) = mean(stdDist(:,labelsrIdx == ii),2);
end
[minDistance, bestIdx] = min(labelDist,[],2);
guessedLabels = tutorLabels(bestIdx);

[juvenileSylls.type] = guessedLabels{:};

%% (4.5) optional, REVIEW juvenile syllable labels
for ii = 1:numel(newMarkedSongs)
    plotForReview(songStruct,juvenileSylls,newMarkedSongs(ii), params, ...
        'editSpecType', 'fine', 'dgram.minContrast',3e-10,...
        'doFilterNoise',true, 'noiseFilter', noiseMask,'nps.reduction',-12, ...
        'optGraphs',{'waveform', 'deriv','fundamentalFreq'});
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
for ii = 1:numel(newMarkedSongs)
    figure(2);
    revisedSylls = plotAndAdjust(songStruct,juvenileSylls,newMarkedSongs(ii), params, ...
        'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',3e-10,...
        'doFilterNoise',true, 'noiseFilter', noiseMask,'nps.reduction',-12);
    newJuvenileSylls = [revisedSylls newJuvenileSylls]; 
end
%juvenileSylls = newJuvenileSylls;
%clear newJuvenileSylls;
%% (optional) CHECK SIMILARITY for a given syllable
selType = input('Which syllable do you want to check? > ','s');
selFilter = strcmp({juvenileSylls.type}, selType);
whichTutorSyll = find(strcmp({tutorSylls.type},selType));
distances = stdDist(selFilter, whichTutorSyll);
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

%% jenny added this here for after manual stuff
uisave({'markedSongs', 'juvenileSylls'}, [matpath prependForSave('markedSongs-',matFile)]);
%% (optional) get distances between juvenile syllables for hierarchical unsupervised clustering
profile on;
[boutDistMatrixMean, boutDistMatrix] = syllableAllCross(songStruct, juvenileSylls);
profile viewer
profile off;
uisave({'boutDistMatrixMean','boutDistMatrix'},prependForSave('corrMT-', matFile));

%% (optional) construct CLUSTERS of letters
disp('Hierarchical clustering and alphabet creation...');
[stringRep, clusterIdxs] = createAlphabet(juvenileSylls, boutDistMatrix, songStruct, ...
    [],'fs',fs, 'playsample', true, 'nClusters',15,'clusterMethod', 'centroid', 'plot',true);
strCelled = cellstr(stringRep');
[juvenileSylls.type] = strCelled{:};

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
[matFile, matSpikePath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data');
spikes = loadSpikeData([matSpikePath matFile]);
% calculate spike data w.r.t segmented regions
%% assess spike counts/rates
ROIs = juvenileSylls;
syllLengths = [ROIs.stop] - [ROIs.start];
for ii = 1:numel(spikes)
    [spikeCounts{ii}, spikeInternalTimes{ii}] = countSpikes(ROIs, spikes{ii},'onset');
    spikeRates{ii} = spikeCounts{ii} ./ syllLengths; 
end
totalSpikeRates = sum(cellfun(@sum,spikeRates));
%spikeTimes = sort(vertcat(spikes{:}));

%% (5.x) option: load ANOTHER spike file (make sure (5) is run first)
%%syllLengths = [juvenileSylls.stop] - [juvenileSylls.start];
while(strcmp(questdlg('Load another spike file?', 'Load more spikes', 'OK','No', 'No'), 'OK'))
    [matFile, matpath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data');
    spikeData = load([matpath matFile]); 
    clusterFields = fieldnames(spikeData);
    
    % append to spikes cell array and recount
    newSpikes = loadSpikeData([matpath matFile]);
    spikes = [spikes newSpikes];
% for ii = 1:numel(newSpikes)
%         [spikeCounts{end+1}, spikeInternalTimes{end+1}] = countSpikes(juvenileSylls, newSpikes{ii},'onset');
%         spikeRates{end+1} = spikeCounts{end} ./ syllLengths;
%     end
%     clear spikeData clusterFields
end
   
% if exist('spikeRates'), totalSpikeRates = sum(vertcat(spikeRates{:})); end;
%[spikeCounts, spikeInternalTimes] = countSpikes(juvenileSylls, spikeTimes,'onset');
%% (5.x) test if spiking is song related?
% get non song firing rates

ROIs = manualMotifs;
syllLengths = [ROIs.stop] - [ROIs.start];
for ii = 1:numel(spikes)
    [spikeCounts{ii}, spikeInternalTimes{ii}] = countSpikes(ROIs, spikes{ii},'onset');
    spikeRates{ii} = spikeCounts{ii} ./ syllLengths; 
end
totalSpikeRates = sum(cellfun(@sum,spikeRates));

totalTime = songStruct.interval * songStruct.length;
isAdult = false(numel(ROIs),1);
%isAdult = strcmp({ROIs.type},'adult');

% get local baseline rate ([-3, -2], [2, 3] seconds before, that do not
% overlap with the song
localBaseEvents = eventFromTimes(...
    [max(0, [ROIs.start] - 3) min(songStruct.length/fs, [ROIs.stop] + 2)], ...
    [max(0, [ROIs.start] - 2), min(songStruct.length/fs, [ROIs.stop] + 3)],  fs);

[localBaseNotSong, ~, localBaseInSong] = findIntersection(localBaseEvents, ROIs);

fullDur = sum([localBaseEvents.stop] - [localBaseEvents.start]);
notSongDur = sum([localBaseNotSong.stop] - [localBaseNotSong.start]);
isSongDur = sum([localBaseInSong.stop] - [localBaseInSong.start]);

fprintf('Baseline samples (possibly redundantly) %fs, %f%% valid...\n', ...
    fullDur, notSongDur/fullDur * 100);
fprintf('Neuron | Global | Local | (p)   | Adult | (p)   | Juvenile | (p)  \n');
for ii = 1:numel(spikes)
    % get song firing rates   
    adultSongs = ROIs(isAdult);
    [spikesInAdultSong, ~, ratesInAdultSong] = countSpikes(adultSongs,spikes{ii});
    adultSongLengths = [adultSongs.stop] - [adultSongs.start];
    
    juvieSongs = ROIs(~isAdult);
    [spikesInJuvieSong, ~, ratesInJuvieSong] = countSpikes(juvieSongs,spikes{ii});
    juvieSongLengths = [juvieSongs.stop] - [juvieSongs.start];
    
    % get global (non-song) baseline rate
    baseRate = (numel(spikes{ii}) - sum(spikesInAdultSong) - sum(spikesInJuvieSong)) / ...
        (totalTime - sum(adultSongLengths) - sum(juvieSongLengths));

    % get local baseline rate (may count some times twice)
    [spikesInValidBase,~, ratesInValidBase] = countSpikes(localBaseNotSong,spikes{ii});
    localValidBaseRate = sum(spikesInValidBase) / notSongDur; 
    % overlapping with sound

    [h,p,ci,stats]=ttest(ratesInValidBase - baseRate);
    [h2,p2,ci2,stats2]=ttest(ratesInAdultSong - localValidBaseRate);
    [h3,p3,ci3,stats3]=ttest(ratesInJuvieSong - localValidBaseRate);
    
    
    fprintf('%6d | %6.3f | %4.3f | %4.3f | %4.3f | %4.3f | %8.3f | %4.3f \n', ...
        ii, baseRate, ...
        localValidBaseRate, p,  ...
        mean(ratesInAdultSong),p2, ...
        mean(ratesInJuvieSong),p3);
end
%% (6) test the correlation between similarity 
% to the tutor syllables and firing rate.

% get juvenile syllable type
syllLabels = {juvenileSylls.type};

% get tutor syllable type
tutorLabels = {tutorSylls.type};
[labelTypes, foo, rTutorIdx] = unique(tutorLabels);
nLabelTypes = numel(labelTypes);

% get minimum distance for each syllable type (b/c tutor song has multiple
% examples of each syllable)
for ii = 1:nLabelTypes
    stdDistMin(:,ii) = min(stdDist(:,rTutorIdx == ii),[],2);
    rLabelIdx(strcmp(syllLabels,labelTypes{ii})) = ii; % what is the most descriptive label?
end

% assign similarity to the syllables
mat2cell(stdDistMin, ones(numel(juvenileSylls),1));
[juvenileSylls.tutorSim] = ans{:};

% set color map
colors = winter(nLabelTypes);

[minDistance, bestIdx] = min(stdDistMin,[],2);
minDistance = minDistance';
%meanDistance = mean(stdDist, 2);
%totalSpikeCounts = {sum(vertcat(spikeCounts{:}))};
totalSpikeCounts = spikeCounts;

nUnits = numel(totalSpikeCounts);

sizeFromLength = ceil((syllLengths - min(syllLengths)) / (max(syllLengths) - min(syllLengths)) * 40 + 5);
sizeFromSim    = ceil((minDistance -  min(minDistance)) / (max(minDistance) - min(minDistance)) * 40 + 5);

for ii = 1:nUnits
    figure(ii+1)
    if nUnits > 1, subplot(ceil(nUnits/2), 2, ii); end;
    for jj = 1:nLabelTypes % for each tutor syllable type
    
        isThisSyll = (rLabelIdx == jj);
        if ~any(isThisSyll), continue; end;
        
        subplot(ceil(nLabelTypes/2), 2, jj)
        %xxdat = minDistance(isThisSyll);
        %xxdat = stdDist(isThisSyll,jj);
        xxdat = syllDistances(isThisSyll);
        yydat = totalSpikeCounts{ii}(isThisSyll); %./ syllLengths(isThisSyll);
        sizeDat = sizeFromLength(isThisSyll);
        hh = scatter(xxdat, yydat, sizeDat, ... 
            colors(jj,:),'filled');
        xlabel('Distance (arbitrary)');
        ylabel('Spike rates (Hz)');
        title(sprintf('Syllable %s, Unit %d, %d syllables',labelTypes{jj}, ii, sum(isThisSyll)));
        % fit a trend line
        hold on;
        [linfit, ~,~,~, fitStats] = regress(yydat', [ones(numel(xxdat),1) xxdat']);
        fprintf('\ttutor syllable [%s] (# = %d), fit for trend, r^2 = %0.3g, F = %0.3g, p = %0.3g\n', ...
            labelTypes{jj}, sum(rLabelIdx == jj), fitStats(1), fitStats(2),fitStats(3));
        legend(sprintf('fit for trend, r^2 = %0.3g, F = %0.3g, p = %0.3g\n',...
            fitStats(1), fitStats(2),fitStats(3)));
        plot(xxdat, linfit(1) + xxdat * linfit(2), '--','Color',colors(jj,:),...
            'HandleVisibility', 'off');
        hold off;
        pause;
        % get and plot averages for each juvenile labeled group
%        meanDistance = mean(minDistance(rLabelIdx == jj));
%        meanSpikeRate = mean(spikeRates{ii}(rLabelIdx == jj));
%        meanSize = mean(sizeFromLength(rLabelIdx==jj));
%        plot(meanDistance, meanSpikeRate, 's', 'Color', colors(jj,:), 'MarkerSize',meanSize);
        % compare to each possible syllable 
        %xx = stdDistMin(:,jj); %xx = xx(:);
        %yy = spikeRates{ii};
        %sz = sizeFromLength;
        %hh = scatter(xx,yy, ...
        %    sz, colors(jj,:),'filled');
        %hold on;
    end
    hold off;
    xlabel('distance to nearest syllable');
    ylabel('spike rates (#)');
    title(sprintf('Spiking rates against similarity, unit %d',ii));
    %if ii == 1
    %    legend(labelTypes)
    %end
end

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
%% (6.8) plot PSTH for all marked songs (TODO: plot with baseline)

ROIs = juvenileSongs;
num2cell([ROIs.stop] - [ROIs.start]);
[ROIs.length] = ans{:};
[sortedROIs, indexByLength] = sortBy(ROIs, 'length');
prePostParams = processArgs(params,'preroll',500,'postroll',500);
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
%% (6.9) what is similarity for spiking vs non spiking songs?

% flag bouts by spiking / non spiking
spikesInSongs = countSpikes(newMarkedSongs, vertcat(spikes{:}));

% sum up distance for each song according to its syllable distances, 
% not according to song difference
tutorTypes = unique({tutorSylls.type});
rSyllIdx = zeros(1,numel(juvenileSylls));
for ii = 1:numel(tutorTypes)
    rSyllIdx(strcmp({juvenileSylls.type}, tutorTypes(ii))) = ii;
end

songScore = zeros(1,numel(newMarkedSongs));
for ii = 1:numel(newMarkedSongs)
    [subSylls,ind] = getSubEvents(newMarkedSongs(ii), juvenileSylls);
    % sum all syllable score differences
    songScore(ii) = mean(stdDist(sub2ind(size(stdDist),ind, rSyllIdx(ind))));
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

%% evaluate correlates in syllable-independent fashion
tutorFrame = addPrePost(bookendedClip(tutorSylls), params);
tutorState = stateVector(tutorFrame,tutorSylls);
silDistance = zeros(1,numel(newMarkedSongs));

% similarity in terms of spacing (hamming distance)
parentage = findParent(addPrePost(newMarkedSongs, params), juvenileSylls);
for ii = 1:numel(newMarkedSongs)
    silDistance(ii) = minEditDistance(...
        stateVector(addPrePost(newMarkedSongs(ii), params), ...
        juvenileSylls(parentage==ii)), ...
        tutorState);
end

% similarity to length included but never shows correlation
songScores = struct('length',num2cell([newMarkedSongs.stop] - [newMarkedSongs.start]),...
    'silDistance', num2cell(silDistance));
%{
% similarity in terms of resemblance to tutor song
[tutorSummary, tutorSpec] = getFeatures(tutorSong,tutorFrame,params,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',false);

features = {'fundamentalFreq','FM','AM','wienerEntropy','aperiodicity','length'};

progressbar('Scoring songs on multiple similarity measures');
for ii = 1:numel(markedSongs)
    [songSummary, songSpec] = getFeatures(songStruct, ...
        markedSongs(ii), params, ...
        'plot',false,'verbose',true,'playsample',false, ...
        'doFilterNoise',true, 'noiseFilter', noiseMask,'nps.reduction',-12);
    [timeSongDist, indivFeatureDists] = standardDistance(songSpec, tutorSpec);
    [songScores(ii).songDist, minIndexDist] = min(timeSongDist);
    for jj = 1:numel(features)
        songScores(ii).([features{jj} '_dist']) = indivFeatureDists(minIndexDist, jj);
    end
    progressbar(ii/numel(markedSongs));
end
    %}
for ii = 1:numel(spikes)
    [basalRate(ii), spikeRate(ii), rateSEM(ii), rate_tp(ii)] = findSpikeCorrelates(newMarkedSongs, songScores,spikes{ii});
end
uisave({'basalRate', 'spikeRate', 'rateSEM', 'rate_tp'},prependForSave('spikeRates-', matFile));
