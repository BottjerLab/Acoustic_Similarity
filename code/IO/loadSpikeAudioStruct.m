function songStruct = loadSpikeAudioStruct(file, opt)
[fp, fs, fext] = fileparts(file);
metaFile = [fp filesep 'meta-' fs fext];

if exist(metaFile, 'file')
    metaStruct = []; load(metaFile);
    
    songStruct = metaStruct;
    if nargin == 2 && strcmp(opt, 'novalues')
        return;
    else
        % trick to get the main struct into a standard name, if there's only one
        % variable in the file
        foo = load(file);
        fld=fieldnames(foo);
        songStruct.values = foo.(fld{1});
    end
else
    foo = load(file);
    fld=fieldnames(foo);
    songStruct = foo.(fld{1});    
end

end