function clusterIdxs = loadAcceptedLabels(birdID, age, session)
% this loads all accepted labels for a given bird / age
% NOTE: session is an optional argument
% note 2: eventual goal is to make more syllables labeled
dataDir  = ['data' filesep birdID filesep];
clusterIdxs = [];

if nargin < 3
    acceptedLabelFile = [dataDir 'acceptedLabels-' birdID '-age' num2str(age) '.mat'];
    if ~exist(acceptedLabelFile, 'file')
        error('loadAcceptedLabels:fileNotFound','No acceptedLabels file found for %s age %d', birdID, age);
    end
    load(acceptedLabelFile);
    clTypes = fieldnames(clusterIdxs);
    if numel(clTypes) == 3,
        % choose the regrouped among approved, regrouped, augmented (see browseAndAccept)
        prefType = clTypes{2};
    else
        prefType = clTypes{1}; % choose the approved
    end
    clusterIdxs = clusterIdxs.(prefType);
elseif nargin == 3
    [ageSylls, syllSessions] = loadAgeSylls(birdID, age);
    isInSession = strcmp(session, syllSessions);
    clusterIdxs = clusterIdxs(isInSession);
end