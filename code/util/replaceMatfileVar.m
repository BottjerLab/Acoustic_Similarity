function replaceMatfileVar(matfile, varsOld, newValues)
% rewrites value varsOld (a cell array of strings) 
% with newValues (a cell array of variables)
assert(numel(varsOld) == numel(newValues))
WS = load(matfile);
for ii = 1:numel(varsOld)
    varOld = varsOld{ii};
    varNew = newValues{ii};
    if strcmp(varOld, varNew), continue; end;
    if ~exist(matfile, 'file')
        error('file [%s] does not exist', matfile);
    end
    if ~isfield(WS, varOld)
        warning('variable [%s] not found in file [%s]', varOld, matfile);
        return;
    end
    WS.(varOld) = varNew;
end
save(matfile, '-struct','WS');
