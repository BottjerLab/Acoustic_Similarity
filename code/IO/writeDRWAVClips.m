function fNames = writeDRWAVClips(DRevents, filestem, params, varargin)
if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});
% options for future: preroll context at low volume

warnstate = warning('MATLAB:audiovideo:wavwrite:dataClipped','off');
[fdir, fstem] = fileparts(filestem);

[s, mess] = mkdir(fdir);
if isempty(mess)
    fprintf('Created directory %s ...\n',fdir);
end

nEv = numel(DRevents);
fNames = cell(1,nEv);
fieldWidth = 4;
progressbar(sprintf('Writing wave files (%d)', nEv));
for ii = 1:nEv
    [cl, fs] = getClip(DRevents(ii));
    fNames{ii} = [fdir filesep fstem '-' sprintf(sprintf('%%0%dd',fieldWidth),ii)];
    wavwrite(cl, fs, fNames{ii});
    progressbar(ii/nEv);
end

warnstate = warning('MATLAB:audiovideo:wavwrite:dataClipped',warnstate);