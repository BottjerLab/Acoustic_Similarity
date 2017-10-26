function events = readAudacityLabelFile(fs, filename)

if nargin < 1
    error('Please give fs as first argument');
end

if nargin < 2
	[filename filepath] = uigetfile(pwd);
	filename = [filepath filesep filename];
end


labeltimes = [];
fid = fopen(filename);
while ~feof(fid)
    str = fgetl(fid);
    nums = sscanf(str,'%f'); 
    labeltimes(end+1) = nums(1);
end
fclose(fid);
%labeltimes = labeltimes{1};
if mod(numel(labeltimes),2) == 1
	warning('Labeltimes are uneven');
end

if numel(labeltimes) == 0
	disp('abort soon, checking with keyboard');
	keyboard
end

songStarts = labeltimes(1:2:end-1);
songStops  = labeltimes(2:2:end);

events = eventFromTimes(songStarts, songStops, fs);

