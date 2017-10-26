function [metaStruct, sAudio] = reformatAsDR(filename)

if nargin < 1
    [fname fp] = uigetfile('Multiselect', 'on');
    filename = strcat(fp, fname);
end
if isstr(filename), filename = {filename}; end;

fileExt = '.mat';    
for ii = 1:numel(filename)
    %% error check: make sure it's a session file
    [fp sessionName] = fileparts(filename{ii});
    seps = strfind(fp,filesep);
    birdID = fp(seps(end)+1:end);
    assert(strcmp(birdID, strtok(sessionName, '_')));

    metaFile = [fp filesep 'meta-' sessionName fileExt];
    if exist(metaFile, 'file') == 2,
        fprintf('Already have 7.3 version for session %s, continuing...\n',...
            sessionName);
        continue;
    end
    sessData = reportOnData(birdID, sessionName);
    [sS, sFile] = loadFromManifest(sessData.manifest,'songStruct');
    if isempty(sS),
        error('missing audio from session %s...\n', sessionName);
    end
    sAudio = sS.values;
    DRfile = [fp filesep sessionName fileExt];
    fprintf('Writing 7.3 compatible file for session %s...\n', sessionName);
    save(DRfile, 'sAudio', '-v7.3');
    metaStruct = rmfield(sS,'values');
    fprintf('Writing meta for session %s to file %s...\n', sessionName, metaFile);
    save(metaFile, 'metaStruct');
end
