%%%% Here we study the tutor song, by similarly parsing out the syllables and
%%%% computing features, etc.  Also we plot spectrograms with the parsed
%%%% syllables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note: look at storeTutorInfo if you want to concatenate samples

%birdTutor = inputdlg('What is the ID of the tutor of the bird? ', 'Tutor name',1,{'Y130')};
%birdTutor = [birdTutor '_TUT'];
%%% TODO: find a way to get to the song folder instead of assuming the file
%%% is in wavs directory
birdID = inputdlg('Which bird are you looking for?', 'Which bird?', 1, {''});
birdID = birdID{1};
[tutorID tutorDir] = getBirdTutor(birdID);

tutorFiles = dir([tutorDir filesep '*.wav']);
tutorFiles = strcat([tutorDir filesep],{tutorFiles.name});
nFiles = numel(tutorFiles);
fprintf('Found %d samples of tutor %s for bird %s.\n', nFiles, tutorID, birdID);
%% Pick a good representation of the tutor song

fprintf('Mark if you want this file to represent tutor song...\n');
isGood = false;
ctr = 0;
while ~isGood
    ctr = ctr + 1;
    tutorStruct = loadWavFile(tutorFiles{ctr});
    tutorParams = processArgs(defaultParams,'fs',1/tutorStruct.interval,'noroll');
    tutorParams.fine.features = [tutorParams.fine.features 'harmonicPitch'];
    wholeRegion = getWholeSongAsRegion(tutorStruct);
    
    % playback the song and see if it's typical    
    isGood = markRegions(tutorStruct, wholeRegion, tutorParams);
end

%% do an automatic parse
tutorSylls = parseRegionsIntoSyllables(tutorStruct, wholeRegion, ...
    tutorParams, 'doFilterNoise', false, 'syllable.comboLength', 5, 'syllable.minPower', 3e-8,'syllable.borderRise', 1e-4);

%% refine syllables and adjust labels at syllable level
tutorSylls = plotAndAdjust(tutorStruct, tutorSylls, wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',3e-11,...
    'optGraphs',{'waveform', 'deriv','wienerEntropy','pitchGoodness', 'fundamentalFreq'});
[tutorSylls.file] = deal(tutorFiles{ctr});

%% characterize the tutor syllables / get features
tutorParams = defaultParams;
tutorParams.fine.features = [tutorParams.fine.features 'harmonicPitch'];
tutorParams.reduceFeatures = {'totalPower' 'wienerEntropy' 'mTD', 'mFD', 'FM',...
    'AM', 'rawAM', 'centerFreq', 'pitchGoodness', 'harmonicPitch', 'aperiodicity', 'fundamentalFreq'};

specFieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'};

[tSFeats,tSSpecs]=getFeatures(tutorStruct,tutorSylls,tutorParams,...
    'plot',false,'verbose',true,'playsample',false,'doFilterNoise',false);

tSSpecs = [tSSpecs{:}];
elimFields = setdiff(fieldnames(tSSpecs), specFieldsToKeep);
tSSpecs = rmfield(tSSpecs, elimFields);
%% save the files pertaining to tutor song
matpath = ['data' filesep birdID filesep];
uisave({'tutorStruct','tutorSylls','tSFeats','tSSpecs'},...
    [matpath 'tutor-' birdID '.mat']);

%% (3.7 optional) REVIEW tutor syllables
%{
wholeRegion = getWholeSongAsRegion(tutorStruct);
plotForReview(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'dgram.minContrast',1e-8,...
    'optGraphs',{'waveform', 'fundamentalFreq', 'AM', 'FM', 'wienerEntropy', 'pitchGoodness'});
%}
%% (3.7 optional) REVIEW tutor syllables with spectrogram
figure
wholeRegion = getWholeSongAsRegion(tutorStruct);
plotForReview(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'dgram.minContrast',1e-9,...
    'showDgramRegionStyle', 'fancy',...
    'optGraphs',{'deriv'});

%%
set(gca,'Position', [0 0 1 1]);
set(gca,'Box', 'off');
set(gcf,'Units', 'normalized');
set(gcf, 'Position', [0 1/2 1 1/4]);

%%
figDir = ['figures' filesep];
export_fig([figDir 'tutor-' birdID '.jpg']);