function [isInManifest, keyFile] = findInManifest(manifest, varName)

% varName should be a character string or cell array of character strings
uncell = false;
if ~iscell(varName), uncell = true; varName = {varName}; end;
nEntries = numel(varName);

isInManifest = false(size(varName)); keyFile = cell(size(varName));
% manifest itself is empty
if isempty(manifest), return; end;


% the manifest shouldn't be created if the original session audio file
% weren't there
for ii = 1:numel(varName)
    if strcmp(varName{ii},'songStruct')
        [~,keyEntry]=min(cellfun('length',{manifest.originalFile}));
        keyFile{ii} = manifest(keyEntry).originalFile;
        isInManifest(ii) = true;
        return;
    end
    % todo: pick the last updated rather than last listed, more robust
    maniIndex = find(strcmp(varName{ii},{manifest.name}), 1, 'last');
    if ~isempty(maniIndex)
        isInManifest(ii) = true;
        manifestEntry = manifest(maniIndex);
        keyFile{ii} = manifestEntry.originalFile;
    end
end
if uncell, keyFile = keyFile{1}; end;