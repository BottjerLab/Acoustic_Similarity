%startup a current session
function loadSessionDataFiles(session)
dataDir = 'data\';

%session = 'Lb277_3_27';
birdID = strtok(session, '_');
dataSubdir = [dataDir dataSubdir filesep];
% look at what files are there
dir([dataSubdir '*' session '*.mat']);

fil=[dataSubdir session '_voice.mat'];
if exist(fil,'file'), load(fil), else fprintf('%s not found...\n', fil); end

fil = [dataSubdir 'markedSongs-' session '_voice.mat']  ;
if exist(fil,'file'), load(fil), else fprintf('%s not found...\n', fil); end

fil = [dataSubdir 'spikeRates-' session '_ch13-16_times.mat'];
if exist(fil,'file'), load(fil), else fprintf('%s not found...\n', fil); end

fil = [dataSubdir 'noiseMask-' session '_voice.mat'];
if exist(fil,'file'), load(fil), else fprintf('%s not found...\n', fil); end
    
fil = [dataSubdir 'tutorSimilarity-' session '_voice.mat'];
if exist(fil,'file'), load(fil), else fprintf('%s not found...\n', fil); end

% rename the songStruct (it's the biggest one) with a little dirty eval
vars = whos; 
[~,bigVarIdx] = max([vars.bytes]);
bigVarName = vars(bigVarIdx).name;
eval(sprintf('songStruct = %s; clear %s', bigVarName, bigVarName));

%% (5) LOAD the SPIKE-SORTED-FILE and calculate firing rate during syllables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% All processing up to this point has been solely on the song file 
%%%% To continue you should have: 
%%%% (a) some measure of distance between syllables and different tutor
%%%% syllables (which we call stdDist),  
%%%% (b) definitions of the tutor syllables (tutorSylls), and
%%%% (c) juvenileSylls with defined boundaries
%%%% AND types that correspond to the tutor syllables' types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is where we LOAD the SPIKE-SORTED-FILE
[matFile, matSpikePath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data');
spikes = loadSpikeData([matSpikePath matFile]);
% calculate spike data w.r.t segmented regions
syllLengths = [juvenileSylls.stop] - [juvenileSylls.start];
for ii = 1:numel(spikes)
    [spikeCounts{ii}, spikeInternalTimes{ii}] = countSpikes(juvenileSylls, spikes{ii},'onset');
    spikeRates{ii} = spikeCounts{ii} ./ syllLengths; 
end
totalSpikeRates = sum(vertcat(spikeRates{:}));
%spikeTimes = sort(vertcat(spikes{:}));

%% (5.x) option: load ANOTHER spike file (make sure (5) is run first)
syllLengths = [juvenileSylls.stop] - [juvenileSylls.start];
while(strcmp(questdlg('Load another spike file?', 'Load more spikes', 'OK','No', 'No'), 'OK'))
    [matFile, matpath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data');
    spikeData = load([matpath matFile]); 
    clusterFields = fieldnames(spikeData);
    
    % append to spikes cell array and recount
    newSpikes = loadSpikeData([matpath matFile]);
    spikes = [spikes newSpikes];
    for ii = 1:numel(newSpikes)
        [spikeCounts{end+1}, spikeInternalTimes{end+1}] = countSpikes(juvenileSylls, newSpikes{ii},'onset');
        spikeRates{end+1} = spikeCounts{end} ./ syllLengths;
    end
    clear spikeData clusterFields
end
   
if exist('spikeRates'), totalSpikeRates = sum(vertcat(spikeRates{:})); end;
%[spikeCounts, spikeInternalTimes] = countSpikes(juvenileSylls, spikeTimes,'onset');