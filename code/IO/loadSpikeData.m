function [spikes, nNeuronsPerFile, isMUA] = loadSpikeData(matFiles, MUAdata)
%LOADSPIKEDATA give spike times in a cell array of vectors, one per neuron
% 
% spikes = loadSpikeData(matFiles) can take any of a number of files
% with no argument, prompts user for files with input dialog

% returns: spikes, a cell array of the number of neurons
% nNeuronsPerFile, the count of neurons in a vector of the number of files
% isMUA, flags for multiunitness, is a vector of the number of neurons
if nargin < 1 
    [matSpikeFiles, matSpikePath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data','MultiSelect','on');
    matFiles = strcat(matSpikePath, matSpikeFiles);
end
spikes = []; nNeuronsPerFile = []; isMUA = [];
if isempty(matFiles)
    return;
end

if ~iscell(matFiles)
    matFiles = {matFiles};
end
nFiles = numel(matFiles);
annoyingWarning = 'MATLAB:unknownElementsNowStruc';
warnState = warning('query', annoyingWarning);
warning('off', annoyingWarning);

spikes = [];
nNeuronsPerFile = zeros(1, nFiles);
for ii = 1:nFiles
    if exist(matFiles{ii}, 'file') ~= 2
        spikes = [spikes {}];
        continue;
    end
    spikeData = load(matFiles{ii}); 
    if ~isfield(spikeData, 'MClust_FeatureTimestamps') % just naked times
        clusterFields = fieldnames(spikeData);
        nNeuronsPerFile(ii) = numel(clusterFields);
        subSpikes = cell(1, nNeuronsPerFile(ii));
        % put the structure data into a cell array, one cell for each unit found
        for jj = 1:nNeuronsPerFile(ii)
            subSpikes{jj} = spikeData.(clusterFields{jj});
        end
    else % times in the featureTimestamps, indices in the clusters -> myPoints
        nNeuronsPerFile(ii) = numel(spikeData.MClust_Clusters);
        subSpikes = cell(1, nNeuronsPerFile(ii));
        for jj = 1:nNeuronsPerFile(ii)
            subSpikes{jj} = spikeData.MClust_FeatureTimestamps(spikeData.MClust_Clusters{jj}.myPoints);
            % spikeData.MClust_Clusters{jj}.myPoints 
        end
    end
    spikes = [spikes subSpikes];
end
warning(warnState.state, annoyingWarning);

if nargout >= 3
    if nargin < 2
        MUAdata = getSpikeMUAData;
    end
    isMUA = getMUAStatus(matFiles, nNeuronsPerFile, MUAdata);
end
end
