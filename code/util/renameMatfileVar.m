function renameMatfileVar(matfile, varsOld, varsNew)

assert(numel(varsOld) == numel(varsNew))
WS = load(matfile);
for ii = 1:numel(varsOld)
    varOld = varsOld{ii};
    varNew = varsNew{ii};
    if strcmp(varOld, varNew), continue; end;
    if ~exist(matfile, 'file')
        error('file [%s] does not exist', matfile);
    end
    if ~isfield(WS, varOld)
        warning('variable [%s] not found in file [%s]', varOld, matfile);
        return;
    end
    if isfield(WS, varNew)
        error('variable [%s] already in [%s]', varNew, matfile);
    end
    
    WS.(varNew) = WS.(varOld);
    WS = rmfield(WS, varOld);
end
save(matfile, '-struct','WS');
