function [spectra,featureTable] = prepareDRsylls(songStruct,motifs,noiseMask)
% Function to calculate the spectra and features of a given number of
% syllables.
%
% Requires the following inputs:
%   songStruct: the song channel data/struct exported from Spike2
%   motifs: the manualMotifs struct generated from segmentAndCluster.m
%   noiseMask: the noiseMask struct generated from segmentAndCluster.m
%
% Edited by EL 2021

%% CONSTANTS: Change if desired
FOLDER_NAME = 'syllable files'; % name of folder to save all files in

%%
% if nargin < 2
%     syllName = 'approvedSyllables';
% end
% 
% if nargin < 3
    params = defaultParams;
% end
% params = processArgs(params, varargin{:});
% 
% birdPath = [pwd filesep 'data' filesep birdID filesep];
% if ~exist(birdPath, 'dir'),
%     error('Data for %s does not exist, aborting...', birdID);
% end
% 
% % compile syllables
% rep = reportOnData(birdID, [], params, 'rejectNoNeuronSessions', 'true');
% DRsylls = struct([]);
% sessionNum = [];
% for ii = 1:numel(rep)
%     iSession = rep(ii).sessionID;
%     iMani = rep(ii).manifest;
%     
%     if ~findInManifest(iMani, syllName),
%         warning('Variable %s not found assigned for session %s, skipping session...', ...
%             syllName, iSession);
%         continue;
%     end
%     sylls = loadFromManifest(iMani, syllName);
%     [sylls.file] = deal([birdPath filesep iSession '.mat' ]);
%     [sylls.age]  = deal(getAgeOfSession(iSession));
%     DRsylls = [DRsylls sylls];
%     sessionNum = [sessionNum ii * ones(1,numel(sylls))];
% end
% 
% % get sampling rate
% ii = 1;
% while ~findInManifest(rep(ii).manifest, 'metaStruct'), ii = ii+1; end;
% fs = readSamplingRate(rep(ii).manifest);
% params.fine.fs = fs;

% remove any syllables that are too short
isTooShort = (params.fine.windowSize / 1000 > [motifs.stop] - [motifs.start]);
motifs(isTooShort) = [];
fprintf('Removing %d syllables that are too short...\n', sum(isTooShort));
%%
% get spectra and feature tables...
% TODO: if necessary, split them into HD space and then compile... how?
N = numel(motifs);
fieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'};
spectra = initEmptyStructArray(fieldsToKeep, N);
featureTable = cell(1,N);
progressbar(sprintf('Calculating spectra & features for regions (# = %d)',N));
for ii = 1:N
    % get noisemask
%     if ii==1 || sessionNum(ii-1) ~= sessionNum(ii)
%         [nMExist, nMFile] = findInManifest(rep(sessionNum(ii)).manifest, 'noiseMask');
%         if ~nMExist
%             warning('Noise mask not available for session %s, please load or exit... ', ...
%                 rep(sessionNum(ii)).sessionID);
%             keyboard
%         else
%             fprintf('Loading noise mask from %s...\n',nMFile);
%             noiseMask = loadFromManifest(rep(sessionNum(ii)).manifest, 'noiseMask');
%         end
%         if ~exist('noiseMask', 'var'), error('Noise mask not loaded'); end
%     end
    
    cl = getClipAndProcess(songStruct,motifs(ii), params, 'noroll', ...
        'doFilterNoise', true, 'noiseFilter', noiseMask);
    pF = params.fine;
    pF.features = [pF.features 'harmonicPitch'];
    pF.fs = 1/songStruct.interval;
    tmpSpec = getMTSpectrumStats(cl, pF);
    featureTable{ii} = extractFeatures(tmpSpec);
    progressbar(ii/N);
    for jj = 1:numel(fieldsToKeep)
        spectra(ii).(fieldsToKeep{jj}) = tmpSpec.(fieldsToKeep{jj});
    end
    % TODO: intermediary saves
end
featureTable = [featureTable{:}];

% fprintf('Saving %s features and spectra to %s...\n',birdID, saveFile);
% if not(isfolder(FOLDER_NAME))
%     mkdir(FOLDER_NAME)
% end
% save([FOLDER_NAME,'\allspecs-'],'spectra','featureTable');