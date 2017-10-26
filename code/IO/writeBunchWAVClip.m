function writeBunchWAVClip(songStruct, events, fileStem, resample_fs)
% WRITEBUNCHWAVCLIP writes clip of a song structure to file
% note: use eventFromTimes (with seconds input) to get a manually denoted
% event

nEv = numel(events);
[cl{1},fs] = getClip(events(1), songStruct);
for ii = 2:numel(events)
    cl{ii} = getClip(events(ii), songStruct);
end

minClipLen = 1.0;  %seconds

clLens = [events.stop]-[events.start];
cumLens = [0 cumsum(clLens)];


maxClips = ceil(cumLens(end)/minClipLen);
refRanges = cell(1,maxClips);

ii = 1; startEntry = 1;
cumLens(end) = Inf; % to provide termination
while startEntry < nEv
    refRanges{ii} = startEntry:find(cumLens > minClipLen + cumLens(startEntry), 1) - 1;
    startEntry = refRanges{ii}(end);
    ii = ii+1;
end
if ii > 1
    refRanges(ii:end) = [];
end

concatCl = cell(1,numel(refRanges));
for ii = 1:numel(refRanges)
    % concatenate clip
    % insert silence in between
    idxRange = refRanges{ii};    
    clipsToStick(1:2:numel(idxRange) * 2 - 1) = cl(idxRange);
    
    silClip = zeros(0.2*fs,1);
    for jj = 2:2:numel(idxRange)*2 - 2; clipsToStick{jj} = silClip; end
    
    concatCl{ii} = vertcat(clipsToStick{:});
end

% file dialog to get a file name
if nargin < 3
    [filename pathname] = uiputfile;
    filename = [pathname filesep filename];
end
% enforce that file ends in '.wav'

for ii = 1:numel(concatCl)
    thisFile = sprintf('%s%03d-%03d.wav',fileStem, min(refRanges{ii}), max(refRanges{ii}));
    fprintf('Writing to [%s]...\n', thisFile);
    if nargin < 4        
        wavwrite(concatCl{ii}, fs, thisFile);
    else
        wavwrite(resample(concatCl{ii}, resample_fs, floor(fs/10)*10),resample_fs, thisFile)
    end
end