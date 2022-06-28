function [noiseMask,manualMotifs,syllables] = segmentAndCluster(songStruct,motifData)
% Function to identifty individual syllables from a singing channel 
% exported from Spike2 in the form of a struct.
% This struct (denoted with the precursor approvedSyllables) includes info
% on the following:
% - The start/stop time of a syllable in seconds
% - The start/stop index of a syllable in the raw voltage data
%
% Also saves/returns a struct representing the noise mask (denoted with the 
% precursor noiseMask).
%
% Assumes the Spike2-exported file contains your motif onsets/offsets (as 
% a level channel; with 'motif' somewhere in the title) and the song 
% channel.
%
% Edited by EL 2021

%% CONSTANTS: Change if desired
% FOLDER_NAME = 'syllable files'; % name of folder to save all files in

%% prep step: load song file (the new way, disk readable files)
fs = 1/songStruct.interval;
params = processArgs(defaultParams,'Gy242');
params.fs = fs; 

%% load motif onset/offset times
returnStarts = motifData.times(logical(motifData.level));
returnStops = motifData.times(~logical(motifData.level));
    
returnSIndices = [zeros(numel(returnStarts),1); ones(numel(returnStops), 1)];
[fusedTimes,sIdxs] = sort([returnStarts; returnStops]);
returnSIndices = returnSIndices(sIdxs);

% make sure alternating
if any(diff(returnSIndices) == 0)
    isBadIndices = [diff(returnSIndices) == 0; sIdxs(end)==0];
    fusedTimes(isBadIndices)
end
assert(numel(returnStarts) == numel(returnStops));
assert(all(returnStarts < returnStops) && all(returnStarts(2:end) > returnStops(1:end-1)));

manualMotifs = eventFromTimes(returnStarts, returnStops, fs);
%% (1.2) find noise 
% define some noise with space between bouts
% NB: since most recordings/recording settings are continuous,
% this can be done once per main session 
% (not necessarily for each 30-minute subsession, but check quality)
notNoise = manualMotifs;
approved = false;
while ~approved
    candidateNoise = autodetectNoise(songStruct, notNoise, params);
    fprintf('Is this clip pure noise? Mark if yes...\n'); 
    approved = markRegions(songStruct,candidateNoise);
    % just black this section out and try again
    if ~approved
        notNoise = [notNoise; eventFromTimes(0,candidateNoise.stop,fs)];
    end
end

noiseMask = noiseAnalysis(songStruct, candidateNoise);

% save noise mask
% if not(isfolder(FOLDER_NAME))
%     mkdir(FOLDER_NAME)
% end
% save([FOLDER_NAME,'\noiseMask-',matFile],'noiseMask');

%% (2) parse juvenile motifs into syllables
ROIs = manualMotifs;
params = processArgs(defaultParams,'fs',1/songStruct.interval, 'preroll', 30, ...
    'inter.freqBands', linspace(1,10240,params.inter.NfreqBands));
syllables = ...
    parseRegionsIntoSyllables(songStruct, ROIs, params, 'Gy242',...
    'doFilterNoise',true,...
    'noiseFilter', noiseMask,'nps.reduction',-18, ...
    'dgram.minContrast',3e-10,'minCenterFreq', 800,...
    'syllable.minLength', 10,...
    'plot', false, 'pause', false);

% save([FOLDER_NAME,'\syllables-',matFile],'manualMotifs','syllables');

%% (2.5) refine (and LABEL) parse with manual work
% do a manual refinement
params = processArgs(defaultParams,'fs',1/songStruct.interval,...
    'editSpecType', 'inter',...
    'inter.freqBands', linspace(1,10240,params.inter.NfreqBands));

%define the specific regions to operate on (like function arguments)
ROIs = manualMotifs;
[~,subROIs] = checkRegions(syllables,fs);

% resets new results if you want to start over
tmpSyllables = initEvents;

for ii = 1:numel(ROIs) % if you want to continue, change this range
    thisEv = ROIs(ii);
    fprintf('Editing %d/%d...\n', ii, numel(ROIs));
    revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
        'editSpecType', 'inter', 'adjustLabels',true, 'dgram.minContrast',1e-12,...
        'doFilterNoise',true, 'noiseFilter', noiseMask,'nps.reduction',-12);
    tmpSyllables = [tmpSyllables revisedSylls']; 
end

% moves new results to permanent location
syllables = tmpSyllables;
clear tmpSyllables

%% (save) do this after you work on (2.5) %%%%
% approvedSyllables = manualSyllables;
% uisave({'approvedSyllables'}, [matpath filesep prependForSave('approvedSyllables-',matFile)]);
% 
% % %% (2.7 optional) unsupervised labeling via 'agglomerative clustering'
% % % get distances between juvenile syllables 
% % ROIs = syllables;
% % profile on;
% % [boutDistMatrixMean, boutDistMatrix] = syllableAllCross(songStruct, ROIs);
% % profile viewer
% % profile off;
% % uisave({'boutDistMatrixMean','boutDistMatrix'},prependForSave('corrMT-', matFile));
% % 
% % % construct labels for letters according to clusters
% % nClusters = 5;
% % disp('Hierarchical clustering and alphabet creation...');
% % [stringRep, clusterIdxs] = createAlphabet(ROIs, boutDistMatrix, songStruct, ...
% %     [],'fs',fs, 'playsample', true, 'nClusters', nClusters,'clusterMethod', 'ward', 'plot',true);
% % 
% % clusteredROIs = ROIs;
% % strCelled = cellstr(stringRep');
% % [clusteredROIs.type] = strCelled{:};
% % [ROIs.clusterType] = strCelled{:};