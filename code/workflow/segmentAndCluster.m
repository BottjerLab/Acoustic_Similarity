
%% prep step: load song file (the new way, disk readable files)
% define files to be operated over
[matFile, matpath] = uigetfile('*.mat','Please choose the song Spike2 file','data');
% load file
songStruct = loadSpikeAudioStruct([matpath matFile]);
%%songStruct.title = matFile;
fs = 1/songStruct.interval;
fprintf('Loaded spike audio %s.\n',matFile);
%% prep step: set parameters from an alias
params = processArgs(defaultParams,'Gy242');
params.fs = fs; 

%% (1-eps) load motif onset/offset times from Spike
birdID = strtok(matFile, '_');
returnFil = [matpath '\motifReturn-' matFile '.txt'];
returnFil = strrep(returnFil, '.mat','')
returnFid = fopen(returnFil);

% pick up header and dead line - first column vector is starts
fgetl(returnFid); fgetl(returnFid);
returnStarts = textscan(returnFid,'%f'); returnStarts = returnStarts{1};

% pick up header and dead line - second column vector is stops 
fgetl(returnFid); fgetl(returnFid);
returnStops = textscan(returnFid,'%f'); returnStops = returnStops{1};
    
fclose(returnFid);

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
uisave('noiseMask',[matpath prependForSave('noiseMask-', matFile)]);

%% (2) parse juvenile motifs into syllables
ROIs = manualMotifs;%labeledBouts(strcmp('BOS',{labeledBouts.type}));
%ROIs = ROIs(randperm(numel(ROIs)));
params = processArgs(defaultParams,'fs',1/songStruct.interval, 'preroll', 30, ...
    'inter.freqBands', linspace(1,10240,params.inter.NfreqBands));
syllables = ...
    parseRegionsIntoSyllables(songStruct, ROIs, params, 'Gy242',...
    'doFilterNoise',true,...
    'noiseFilter', noiseMask,'nps.reduction',-18, ...
    'dgram.minContrast',3e-10,'minCenterFreq', 800,...
    'syllable.minLength', 10,...
    'plot', false, 'pause', false);

uisave({'manualMotifs', 'syllables'}, [matpath prependForSave('syllables-',matFile)]);

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
manualSyllables = tmpSyllables;
clear tmpSyllables

%% (save) do this after you work on (2.5) %%%%
approvedSyllables = manualSyllables;
uisave({'approvedSyllables'}, [matpath filesep prependForSave('approvedSyllables-',matFile)]);

%% (2.7 optional) unsupervised labeling via 'agglomerative clustering'
% get distances between juvenile syllables 
ROIs = syllables;
profile on;
[boutDistMatrixMean, boutDistMatrix] = syllableAllCross(songStruct, ROIs);
profile viewer
profile off;
uisave({'boutDistMatrixMean','boutDistMatrix'},prependForSave('corrMT-', matFile));

% construct labels for letters according to clusters
nClusters = 5;
disp('Hierarchical clustering and alphabet creation...');
[stringRep, clusterIdxs] = createAlphabet(ROIs, boutDistMatrix, songStruct, ...
    [],'fs',fs, 'playsample', true, 'nClusters', nClusters,'clusterMethod', 'ward', 'plot',true);

clusteredROIs = ROIs;
strCelled = cellstr(stringRep');
[clusteredROIs.type] = strCelled{:};
[ROIs.clusterType] = strCelled{:};