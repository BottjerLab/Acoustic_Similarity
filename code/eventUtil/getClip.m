function [ret, fs] = getClip(region, songStruct)
%GETCLIP extracts usable clip from sound structure
%  
%   getClip(SONGSTRUCT, REGION) extracts the waveform from the values field
%   of a region, given the index bounds of the region, converting to double.
%   This is useful for (a) plotting the waveform and (b) doing spectral 
%   and subsequent analysis on the sample.  For operations on intervals, 
%   however, we prefer to pass the region structure and song structure 
%   whenever possible.
%   Note for programmers: the songStructure is basically treated as a
%   memory-hogging, copy-on-write variable in MATLAB, while other vectors 
%   (often with redundant data) may often be copied.

%   revision, 7/30/13: Matlab v7.3 and later supports partial disk reading to
%   memory, so we use getClipDR with a named file reference.
%   
%   In cases where the precision of the original audio file is less than
%   64-bit, the promotion is automatically adjusted.

if nargin == 1 || isempty(songStruct)
    % do we have the correct file, and does it exist?
    if ~isfield(region,'file')
        error('getClip:badFileReference', 'region does not have file reference');
    end
    
    [fp fstem fext] = fileparts(region.file);
    metaFile = [fp filesep 'meta-' fstem fext];

    if ~exist(region.file,'file')
        error('getClip:badFileReference', 'region file reference does not exist');
    end
    if ~exist(metaFile, 'file')
        error('getClip:badFileReference', 'region meta file reference does not exist');
    end
    fileObj = matfile(region.file);
    shell = load(metaFile); audioLen = shell.metaStruct.length;
    
    if nargout == 2
        fs = 1/shell.metaStruct.interval;
    end
else
    audioLen = songStruct.length;
    if nargout == 2
        fs = 1/songStruct.interval;
    end
end

% bounds checking
if floor(region.idxStart) < 1 
    warning('getClip:outOfBounds', ...
            'region begins before songStruct start boundary, clipping');
    region.idxStart = 1;
end
if floor(region.idxStop) > audioLen;
    warning('getClip:outOfBounds', ...
            'region begins after songStruct end boundary. clipping');
    region.idxStop = audioLen;
end

if nargin == 1 || isempty(songStruct)
    rawVals = fileObj.sAudio(floor(region.idxStart):floor(region.idxStop),1);
else
    rawVals = songStruct.values(floor(region.idxStart):floor(region.idxStop));
end

% conversion to double, based on information in 'doc wavread'
switch class(rawVals)
    case 'char'
        ret = double(rawVals)/(2^7) - 0.5;
    case 'int16'
        ret = double(rawVals)/(2^15);
    case 'int32'
        ret = double(rawVals)/(2^23);
    case 'single'
        ret = double(rawVals);
    case 'double'
        ret = rawVals;
end