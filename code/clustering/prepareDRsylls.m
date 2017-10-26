function prepareDRsylls(birdID, syllName, params, varargin)
if nargin < 2
    syllName = 'approvedSyllables';
end

if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});

birdPath = [pwd filesep 'data' filesep birdID filesep];
if ~exist(birdPath, 'dir'),
    error('Data for %s does not exist, aborting...', birdID);
end

% compile syllables
rep = reportOnData(birdID, [], params, 'rejectNoNeuronSessions', 'true');
DRsylls = struct([]);
sessionNum = [];
for ii = 1:numel(rep)
    iSession = rep(ii).sessionID;
    iMani = rep(ii).manifest;
    
    if ~findInManifest(iMani, syllName),
        warning('Variable %s not found assigned for session %s, skipping session...', ...
            syllName, iSession);
        continue;
    end
    sylls = loadFromManifest(iMani, syllName);
    [sylls.file] = deal([birdPath filesep iSession '.mat' ]);
    [sylls.age]  = deal(getAgeOfSession(iSession));
    DRsylls = [DRsylls sylls];
    sessionNum = [sessionNum ii * ones(1,numel(sylls))];
end

% get sampling rate
ii = 1;
while ~findInManifest(rep(ii).manifest, 'metaStruct'), ii = ii+1; end;
fs = readSamplingRate(rep(ii).manifest);
params.fine.fs = fs;

% remove any syllables that are too short
isTooShort = (params.fine.windowSize / 1000 > [DRsylls.stop] - [DRsylls.start]);
DRsylls(isTooShort) = [];
sessionNum(isTooShort) = [];
fprintf('Removing %d syllables that are too short...\n', sum(isTooShort));
%%
% get spectra and feature tables...
% TODO: if necessary, split them into HD space and then compile... how?
N = numel(DRsylls);
fieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'};
spectra = initEmptyStructArray(fieldsToKeep, N);
featureTable = cell(1,N);
progressbar(sprintf('Calculating spectra & features for regions (# = %d)',N));
for ii = 1:N
    % get noisemask
    if ii==1 || sessionNum(ii-1) ~= sessionNum(ii)
        [nMExist, nMFile] = findInManifest(rep(sessionNum(ii)).manifest, 'noiseMask');
        if ~nMExist
            warning('Noise mask not available for session %s, please load or exit... ', ...
                rep(sessionNum(ii)).sessionID);
            keyboard
        else
            fprintf('Loading noise mask from %s...\n',nMFile);
            noiseMask = loadFromManifest(rep(sessionNum(ii)).manifest, 'noiseMask');
        end
        if ~exist('noiseMask', 'var'), error('Noise mask not loaded'); end
    end
    
    cl = getClipAndProcess([],DRsylls(ii), params, 'noroll', ...
        'doFilterNoise', true, 'noiseFilter', noiseMask);
    pF = params.fine;
    pF.features = [pF.features 'harmonicPitch'];
    tmpSpec = getMTSpectrumStats(cl, pF);
    featureTable{ii} = extractFeatures(tmpSpec);
    progressbar(ii/N);
    for jj = 1:numel(fieldsToKeep)
        spectra(ii).(fieldsToKeep{jj}) = tmpSpec.(fieldsToKeep{jj});
    end
    % TODO: intermediary saves
end
featureTable = [featureTable{:}];

saveFile = [birdPath 'allSpecs-' birdID '.mat'];
fprintf('Saving %s features and spectra to %s...\n',birdID, saveFile);
save(saveFile,'DRsylls','spectra','featureTable');