function songStruct = loadManyWavFile(filenames, silLength, fs)
% puts many wav files together
% optional arguments silLength introduces silence (in seconds)
% between each clip
% fs forces the fs rate, otherwise it will be the first frame rate we see
% if there is no specified fs, the fs should be equal across all videos
if nargin < 2, silLength = 0; end;
if nargin < 3, fs = NaN; end

songStruct = loadWavFile(filenames{1});
for ii = 2:numel(filenames)
    fil = filenames{ii};
    filInfo = dir(fil);
    
    [y, fsTemp, nbits] = wavread(fil,'native');
    fprintf('Reading file [%s], size %.1f MB, sample rate %0.1f kHz, encoding depth %d...\n',...
        fil, filInfo.bytes/(2^20),fsTemp/1000, nbits);
    if isnan(fs), fs = fsTemp; end;
    
    if silLength > 0 
        y = [zeros(silLength * fs, 1); y];
    end
    
    % populate some fields
    songStruct.values = [songStruct.values; y];
    songData = whos('songStruct');
    fprintf('Size is now %0.3f MB\n', songData.bytes/(2^20));
end
songStruct.interval = 1/fs;
songStruct.length = numel(songStruct.values);
end