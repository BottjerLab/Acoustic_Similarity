function writeWAVClip(songStruct, event, filename, resample_fs)
% WRITEWAVCLIP writes clip of a song structure to file
% note: use eventFromTimes (with seconds input) to get a manually denoted
% event
if nargin < 2
    event = getWholeSongAsRegion(songStruct);
end
[cl,fs] = getClip(event, songStruct);

% file dialog to get a file name
if nargin < 3
    [filename pathname] = uiputfile;
    filename = [pathname filesep filename];
end
% enforce that file ends in '.wav'
if strcmp(filename(end-4:end), '.wav') ~= 0
    filename = [filename '.wav'];
end

if nargin < 4
    wavwrite(cl, fs, filename);
else
    wavwrite(resample(cl, resample_fs, floor(fs/10)*10),resample_fs, filename)
end