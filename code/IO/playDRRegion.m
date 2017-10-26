function playDRRegion(events, doBlock)
% thin wrapper around playSound for disk-stored events
if nargin < 2, doBlock = true; end;

for ii = 1:numel(events)
    clip = getClip(events(ii));
    
    % get the fs of that recording...
    [fp,fs,fext] = fileparts(events(ii).file);
    metaFile = [fp filesep 'meta-' fs fext];
    metaStruct=[]; load(metaFile);
    fs = 1/metaStruct.interval;
    
    playSound(clip,fs,doBlock);
end