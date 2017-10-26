function newSongs = cleanupSongs(songStruct, songs, bouts, exampleSong, params, varargin)
% manually adjust boundaries for songs
% it is recommended to split this job into parts in case something
% goes wrong

% TODO: allow labels for each interval 
if nargin < 5 || isempty(params)
    params = defaultParams;
end

params = processArgs(params,varargin{:});
params.inter.fs = params.fs;

% add a small buffer for each bout
params.adjustLabels = false;
boutParams = processArgs(params,'preroll',params.silentPeriod * 1000,...
    'postroll',params.silentPeriod * 1000);
bouts = addPrePost(bouts,boutParams);

% calculate stats for the model song 
contextSong = addPrePost(exampleSong, params);
songInContext = adjustTimeStamps(exampleSong, -contextSong.start);
contextSpec = getMTSpectrumStats(getClip(contextSong, songStruct), params.inter);
% show the model song
figure(1); 
plotAllFigures(contextSpec,songInContext,params,'optGraphs',{'waveform','totalPower'},'showLabels',true);
title('Model song');

MAXSONGS = 300;
newSongs = initEvents(MAXSONGS);
ctr=1;
% adjust the new boundaries for each song
for ii = 1:numel(bouts)
    fprintf('Editing in bout %d/%d...\n',ii,numel(bouts));
    figure(2)
    % plot the editor for the edited song: just the waveform should be enough to detect song
     tempSongs = plotAndAdjust(songStruct,...
        getSubEvents(bouts(ii),songs),bouts(ii), params,'optGraphs',{'waveform'});

    nNew = numel(tempSongs);
    if nNew > 0
        newSongs(ctr:ctr+nNew-1) = tempSongs;
        ctr = ctr + nNew; % new insert point
        if ctr > numel(newSongs)
            newSongs = [newSongs initEvents(MAXSONGS)];
        end
    end
    clf('reset');
end

% clip the empty songs off the array
newSongs(ctr:end) = [];

end