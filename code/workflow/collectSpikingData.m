[matFile matpath] = uigetfile('data\','Loading bout boundaries...');
manifest=whos('-file',[matpath matFile]);
if ~any(strcmp({manifest.name},'newMarkedSongs'))
    whos('-file',[matpath matFile]);
    warning('File does not contain newMarkedSongs, please correct...'); 
    keyboard
    localManifest = who;
    if ~any(strcmp(localManifest,'newMarkedSongs'))
        error('Still no newMarkedSongs...');
    end
end
load([matpath matFile])
%%
% this is where we LOAD the SPIKE-SORTED-FILE
[matSpikeFiles, matSpikePath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data','MultiSelect','on');

spikes = loadSpikeData(strcat(matSpikePath, matSpikeFiles));
areTheyShell = strcmp(questdlg('Is this a file of shell neurons?', 'Shell?', 'Yes','No','Yes'),'Yes');
isShell = areTheyShell * ones(numel(spikes),1);
% calculate spike data w.r.t segmented regions

% load additional spike files 
while(strcmp(questdlg('Load another spike file?', 'Load more spikes', 'OK','No', 'No'), 'OK'))
    [matSpikeFiles, matpath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data','MultiSelect','on');
    
    % append to spikes cell array and recount
    newSpikes = loadSpikeData(strcat(matSpikePath, matSpikeFiles));
    if ~isempty(newSpikes)
        spikes = [spikes newSpikes];
        areTheyShell = strcmp(questdlg('Is this a shell file?', 'Shell?', 'Yes','No', 'Yes'),'Yes');
        isShell = [isShell; areTheyShell * ones(numel(newSpikes),1)];
    end
end

%% load original file just to pull out length;
[matFile, matpath] = uigetfile('*.mat','Please choose the song Spike2 file','data');
% load file
songStruct = load([matpath matFile]); 
fld=fieldnames(songStruct);
songStruct=songStruct.(fld{1});
songStruct = rmfield(songStruct,'values');
totalTime = songStruct.interval * songStruct.length; % time in seconds
%% test if spiking is song related?
% spike analysis

ROIs = newMarkedSongs;


isAdult = strcmp({ROIs.type},'adult');

% get local baseline rate ([-3, -2], [2, 3] seconds before, that do not
% overlap with the song
localBaseEvents = eventFromTimes(...
    [max(0, [ROIs.start] - 3) min(totalTime, [ROIs.stop] + 2)], ...
    [max(0, [ROIs.start] - 2), min(totalTime, [ROIs.stop] + 3)], 1/songStruct.interval);

% clear the empties out
localBaseEvents([localBaseEvents.start]==[localBaseEvents.stop])=[];

localBaseNotSong = findUniqueAreas(localBaseEvents, ROIs);
localBaseInSong = findIntersection(localBaseEvents, ROIs);

% clear the empties out
localBaseNotSong([localBaseNotSong.start]==[localBaseNotSong.stop])=[];
localBaseInSong([localBaseInSong.start]==[localBaseInSong.stop])=[];

fullDur = sum([localBaseEvents.stop] - [localBaseEvents.start]);
notSongDur = sum([localBaseNotSong.stop] - [localBaseNotSong.start]);
isSongDur = sum([localBaseInSong.stop] - [localBaseInSong.start]);

fprintf('Baseline samples (possibly redundantly) %0.2fs, %f%% valid...\n', ...
    fullDur, notSongDur/fullDur * 100);
fprintf('Neuron | Global | Local | (p)   | Adult | (p)   | Juvenile | (p)  \n');

for ii = 1:numel(spikes)
    % get song firing rates   
    % may have proportionality problem in calc'ing average of spiking rates, 
    % not average spiking rate
    adultSongs = ROIs(isAdult);
    [spikesInAdultSong, ~, stats(ii).ratesInAdultSong] = countSpikes(adultSongs,spikes{ii});
    stats(ii).adultSongLengths = [adultSongs.stop] - [adultSongs.start];
    stats(ii).adultRate = mean(stats(ii).ratesInAdultSong);
    if ~isempty(adultSongs)
        stats(ii).adultRateSEM = std(stats(ii).ratesInAdultSong) / sqrt(numel(adultSongs)-1);
    else
        stats(ii).adultRateSEM = NaN;
    end
    
    juvieSongs = ROIs(~isAdult);
    [spikesInJuvieSong, ~, stats(ii).ratesInJuvieSong] = countSpikes(juvieSongs,spikes{ii});
    stats(ii).juvieSongLengths = [juvieSongs.stop] - [juvieSongs.start];
    stats(ii).juvieRate = mean(stats(ii).ratesInJuvieSong);
    if ~isempty(juvieSongs)
        stats(ii).juvieRateSEM = std(stats(ii).ratesInJuvieSong) / sqrt(numel(juvieSongs)-1);
    else
        stats(ii).juvieRateSEM = NaN;
    end
    
    % get global (non-song) baseline rate
    stats(ii).baseRate = (numel(spikes{ii}) - sum(spikesInAdultSong) - sum(spikesInJuvieSong)) / ...
        (totalTime - sum(stats(ii).adultSongLengths) - sum(stats(ii).juvieSongLengths));

    % get local baseline rate (may count some times twice) 
    [spikesInValidBase,~, stats(ii).ratesInValidBase] = countSpikes(localBaseNotSong,spikes{ii});
    stats(ii).localValidBaseRate = mean(stats(ii).ratesInValidBase); 
    stats(ii).localValidBaseRateSEM = std(spikesInValidBase) / notSongDur ...
        / sqrt(numel(spikesInValidBase)-1);
    % overlapping with sound

    [h1,p1,ci1,stats1]=ttest(stats(ii).ratesInValidBase - stats(ii).baseRate);
    [h2,p2,ci2,stats2]=ttest(stats(ii).ratesInAdultSong - stats(ii).localValidBaseRate);
    [h3,p3,ci3,stats3]=ttest(stats(ii).ratesInJuvieSong - stats(ii).localValidBaseRate);
    
    fprintf('%6d | %6.3f | %4.3f | %4.3f | %4.3f | %4.3f | %8.3f | %4.3f \n', ...
        ii,stats(ii). baseRate, ...
        stats(ii).localValidBaseRate, p1,  ...
        stats(ii).adultRate,p2, ...
        stats(ii).juvieRate,p3);
end

num2cell(isShell); [stats.isShell]=ans{:};
clear isShell
%% save results
uisave('stats',prependForSave('spikeSongStats-', matFile));

%% plot results
ploterr([stats.localValidBaseRate], [stats.juvieRate],[stats.localValidBaseRateSEM],[stats.juvieRateSEM], 'b.')
hold on; plot([0 25],[0 25], 'k-'); hold off;
clear stats