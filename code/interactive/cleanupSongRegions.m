function [newSongSyllables,modelSylls] = cleanupSongRegions(songStruct, songs, syllables, exampleSong, params, varargin)
%CLEANUPSONGREGIONS relabel songs syllables based on an example in the recording
%
% relabel individual syllables within songs
% it is recommended to split this job into parts in case something
% goes wrong

% TODO: allow labels for each interval
if nargin < 5 || isempty(params)
    params = defaultParams;
end

params = processArgs(params,varargin{:});
params.inter.fs = params.fs;
MAXSONGS = 300;
%% adjust the new boundaries for each song

%% adjust and label just the model syllables
params.adjustLabels = true; % now we are changing the labels
params.optGraphs = {'waveform','totalPower','wienerEntropy','centerFreq','deriv'};
disp('Setting the model syllables...');
extendedSong = addPrePost(exampleSong, processArgs(params,'preroll',50,'postroll',50));
modelSylls = plotAndAdjust(songStruct,...
    getSubEvents(exampleSong, syllables), extendedSong, params);

% show model
modelCl = getClip(extendedSong, songStruct);
modelSylls = adjustTimeStamps(modelSylls, -extendedSong.start);


modelSpec = getMTSpectrumStats(modelCl, params.inter);

figure(1);
plotAllFigures(modelSpec,modelSylls,params,'showLabels',true);
title('Model syllables');

% preallocate events
insertCtr = 1;
newSongSyllables = initEvents(MAXSONGS);

figure(2)
for ii=1:numel(songs)
    %% adjust and label each song's syllables
    
    fprintf('Adjusting syllables for song %d...\n',ii);
    % adjust the syllables in a song, using the model song as a visual template
    tempSylls = plotAndAdjust(songStruct,...,
        getSubEvents(songs(ii),syllables),songs(ii), params);
    
    % if there are new syllables, add them
    nNew = numel(tempSylls);
    if nNew > 0
        newSongSyllables(insertCtr:insertCtr+nNew-1) = tempSylls;
        insertCtr = insertCtr + nNew;
        if insertCtr > numel(newSongSyllables)
            newSongSyllables = [newSongSyllables initEvents(MAXSONGS)];
        end
    end
    
    % clear figure and callbacks
    clf('reset');
end
newSongSyllables(insertCtr:end) = [];