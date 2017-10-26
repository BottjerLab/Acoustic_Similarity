diary on
clock

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
params = processArgs(defaultParams,'Y231early');

%% (0.5) isolate sounds first
sounds = stepSpectrogram(songStruct, params,'Nsplits',400,'plot',true);

%% refine sound boundaries
[candidateSyl, features, featureReduction] = ...
    parseRegionsIntoSyllables(songStruct, sounds, params, 'plot',false,...
    'syllable.minPower', 2.5e-6, 'syllable.minLength', 15, 'pause', false); 

% TODO: clean out obvious spikes/wingflaps?
%% arrange into bouts
[bouts, nSyllsInBout] = arrangeBouts(candidateSyl);
%% save some work - automatically found partitions
save([matpath prependForSave('divMT-',matFile)], 'sounds', 'candidateSyl', 'bouts', 'features', 'featureReduction');

%% (x.8) find noise 
% define some noise with space between bouts
notNoise = bouts;
approved = false;
while ~approved
    candidateNoise = autodetectNoise(songStruct, notNoise); % todo: improve/verify noise finding
    fprintf('Is this clip pure noise? Mark if yes...\n'); 
    approved = true;
    disp('auto yes');
    %    approved = markRegions(songStruct,candidateNoise);
    % just black this section out and try again
%    if ~approved
%        notNoise = [notNoise candidateNoise];
%    end
end
noiseMask = noiseAnalysis(songStruct, candidateNoise);

% save noise mask
save([matpath prependForSave('noiseMask-', matFile)], 'noiseMask');

params.noiseFilter = noiseMask;
%markedSongs(insertCtr:end) = [];    

clock

%% substitute them for now
markedBouts = bouts;
%% (1.5 auto) refine/subparse song boundaries and chop 'songs' which are here bouts 
% into real songs 

% TODO: include learner training in here
prePostParams = processArgs(params,'preroll',500,'postroll',500);

newMarkedSongs = initEvents(0);

% initialize learner
%songLearner = Learner(songStruct, [], params);

sortedMarkedBouts = markedBouts;
num2cell([markedBouts.stop] - [markedBouts.start]); [sortedMarkedBouts.length] = ans{:};

sortedMarkedBouts = sortBy(sortedMarkedBouts,'length');
for ii = 174:174%numel(sortedMarkedBouts) % or 218
    profile on
    thisBout = sortedMarkedBouts(ii);
    fprintf('Bout %d/%d (length %0.2f s)... ', ii, numel(sortedMarkedBouts), thisBout.length);
        tic
    %if songLearner.nExamples < 5 % minimum number of songs for robust classification
        % keep manually classifying until we get 5 songs
        tmpSongs = plotAndAdjust(songStruct, thisBout, ...
            addPrePost(thisBout, prePostParams), params,...
            'editSpecType','inter',...
            'inter.features', {'totalPower','wienerEntropy','deriv','waveform'});
     %   if ~isempty(tmpSongs)
    %        songLearner.addSpectra(songStruct, tmpSongs);
    %        songLearner.adaptThreshold;
    %    end
    %else
    %    tmpSongs = scanForMatches(songStruct, addPrePost(thisBout, prePostParams), songLearner);
    %    fprintf('%d song%s found ... ', numel(tmpSongs), char('s' * (1~=numel(tmpSongs))));
    %end
    toc % generates new line
    %newMarkedSongs = [newMarkedSongs tmpSongs];
    profile viewer
end

%% auto correction: merge all overlapping segments
olaps = findOverlaps(newMarkedSongs);
adjustedMarkedSongs = newMarkedSongs; 
markDel = false(numel(newMarkedSongs),1);
for ii = 1:size(olaps,1)
    adjustedMarkedSongs(olaps(ii,1)).idxStart = newMarkedSongs(olaps(ii,2)).idxStart; 
    adjustedMarkedSongs(olaps(ii,1)).start = newMarkedSongs(olaps(ii,2)).start; 
    markDel(olaps(ii,1)) = true;
end
adjustedMarkedSongs = adjustedMarkedSongs(~markDel); 

%clear adjustedMarkedSongs
%%
uisave({'markedBouts','newMarkedSongs'},[matpath prependForSave('autoMarkedSongs-',matFile)]);
diary off