%% prep step: load song file
% define files to be operated over
[matFile, matpath] = uigetfile('*.mat','Please choose the song Spike2 file','data');
% load file
sessionID = strrep(matFile,'.mat','')
birdID = strtok(sessionID, '_')
info = reportOnData(birdID, sessionID, [],'verbose','false');

%% (x.8) find noise 
% define some noise with space between bouts
% NB: since most recordings/recording settings are continuous,
% this can be done once per main session 
% (not necessarily for each 30-minute subsession, but check quality)
songStruct = loadFromManifest(info.manifest, 'songStruct');
%%
notNoise = loadFromManifest(info.manifest, 'manualMotifs');
notNoise = [notNoise; eventFromTimes(0,43, fs)];
fs = 1/songStruct.interval;
approved = false;
while ~approved
    candidateNoise = autodetectNoise(songStruct, notNoise, params, 'noiseDynamicRange' ,1e-4, 'noiseLength', 2);
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