%% Pick the tutor song and study it
birdID = 'Lb277';
tutorDir = uigetdir('',sprintf('Where are the tutor recordings for %s?',birdID));

tutorFiles = dir([tutorDir filesep '*.wav']);
tutorFiles = strcat([tutorDir filesep],{tutorFiles.name});

whichIsTutor = 0;
approved = false;
for ii = 1:numel(tutorFiles);
    tutorStruct = loadWavFile(tutorFiles{ii});
    
    tutorParams = processArgs(defaultParams,'fs',1/tutorStruct.interval,'preroll',0,'postroll',0);
    tutorParams.fine.features = [tutorParams.fine.features, 'harmonicPitch', 'mfcc'];
    tutorParams.reduceFeatures = [tutorParams.reduceFeatures 'mfcc'];
    wholeRegion = getWholeSongAsRegion(tutorStruct);
    
    % pre highPass and amplify
    tutorStruct.values = 3*highPassSample(getClip(wholeRegion, tutorStruct), tutorParams);

    approved = markRegions(tutorStruct,wholeRegion,tutorParams,'preroll',0,'postroll',0);
    if(approved)
        whichIsTutor = ii; break; 
    end
end

% do an automatic parse
[tutorSylls, features] = parseRegionsIntoSyllables(tutorStruct, wholeRegion,tutorParams,...
    'doFilterNoise',false,'syllable.comboLength',18);

% do a manual refinement on syllable level
tutorSylls = plotAndAdjust(tutorStruct,tutorSylls,wholeRegion, tutorParams, ...
    'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8,...
    'optGraphs',{'waveform', 'deriv','FM','AM', 'mTD'});

[tSFeats,tSSpecs]=getFeatures(tutorStruct,tutorSylls,tutorParams,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',false); 
%% load DRsylls
load([pwd filesep 'data' filesep birdID filesep 'allSpecs-' birdID '.mat']); % load drsylls

sForSylls = cell(1,numel(DRsylls));
for ii = 1:numel(DRsylls)
    [~,sForSylls{ii}] = fileparts(DRsylls(ii).file);
end
[uSessions, ~, sessIdx] = unique(sForSylls);

% get syllables for a certain session
seld = listdlg('ListString',uSessions,'SelectionMode','single','Name','Which session?');
seldSession = uSessions{seld};
juvenileSylls = DRsylls(sessIdx == seld);
%% (5) LOAD the SPIKE-SORTED-FILE and calculate firing rate during syllables
[matFile, matSpikePath] = uigetfile('*.mat','Please choose the SPIKING Spike2 file','data');
spikes = loadSpikeData([matSpikePath matFile]);
% calculate spike data w.r.t segmented regions
for ii = 1:numel(spikes)
    [spikeCounts{ii}, ~, spikeRates{ii}] = countSpikes(juvenileSylls, spikes{ii},'onset');    
end
%%
progressbar('getting juvie syllable similarity', 'each tutor syllable');

fs = readSamplingRate(getfield(reportOnData(birdID, seldSession),'manifest'));
reductions = {@nanmean, @max, @min, @median, @geomean, @(x) sqrt(nanmean(power(x,2)))};
nRed = numel(reductions);
syllDistances = zeros(numel(juvenileSylls), numel(tutorSylls), nRed + 3);

juvParams = processArgs(defaultParams,...
        'fs', fs, 'plot',false,'verbose',true,'playsample',false,'doFilterNoise',true);
juvParams.fine.features = [juvParams.fine.features, 'harmonicPitch', 'mfcc'];
juvParams.reduceFeatures = [juvParams.reduceFeatures 'mfcc'];

for ii = 1:numel(juvenileSylls)
    [jSFeats,jSSpec]=getFeatures([],juvenileSylls(ii),juvParams); 
    jSSpec = jSSpec{1};
    
    for jj = 1:numel(tutorSylls)
        % this correlation vector is the score for each possible alignment
        % of two syllables of different lengths
        distVector = standardDistance(jSSpec, tSSpecs{jj},defaultParams);
        
        % these are various reductions 
        
        for kk = 1:numel(reductions)
            syllDistances(ii,jj,kk) = reductions{kk}(distVector);
        end
        
        syllDistances(ii, jj, nRed + 1) = timeWarpedDistance(jSSpec, tSSpecs{jj});
        syllDistances(ii, jj, nRed + 2) = timeWarpedDistanceMFCC(jSSpec, tSSpecs{jj});
        
        commFeatures = intersect(fieldnames(tSFeats), fieldnames(jSFeats));
        fV1 = cellfun(@(x) tSFeats(jj).(x), commFeatures);
        fV2 = cellfun(@(x)     jSFeats.(x), commFeatures);
        syllDistances(ii, jj, nRed + 3) = sqrt(dot(fV1-fV2, fV1-fV2));        
    end
    progressbar(ii/numel(juvenileSylls),[]);   
end

%% (6) test the correlation between similarity
% to the tutor syllables and firing rate.

% get juvenile syllable type
while ~exist('syllLabels', 'var')
    fprintf('Import me the syllable labels for bird %s, session %s, and place in syllLabels as a numeric vector, any way you can... ', ...
        birdID, seldSession);
    keyboard
end

% get tutor syllable type
tutorLabels = str2double({tutorSylls.type});
[labelTypes, ~, rTutorIdx] = unique(tutorLabels);
nLabelTypes = numel(labelTypes);

% get minimum distance for each syllable type (b/c tutor song has multiple
% examples of each syllable)
stdDist = syllDistances(:,:,2);
stdDistMin = zeros(size(stdDist,1),nLabelTypes);
for ii = 1:nLabelTypes
    stdDistMin(:,ii) = min(stdDist(:,rTutorIdx == ii),[],2);
end

% assign similarity to the syllables
mat2cell(stdDistMin, ones(numel(juvenileSylls),1));
[juvenileSylls.tutorSim] = ans{:};
%%
% set color map
colors = winter(nLabelTypes);

syllLengths = [juvenileSylls.stop] - [juvenileSylls.start];
[minDistance, bestIdx] = min(stdDistMin,[],2);
minDistance = minDistance';

nUnits = numel(totalSpikeCounts);

sizeFromLength = ceil((syllLengths - min(syllLengths)) / (max(syllLengths) - min(syllLengths)) * 40 + 5);
sizeFromSim    = ceil((minDistance -  min(minDistance)) / (max(minDistance) - min(minDistance)) * 40 + 5);

for ii = 3:3%nUnits
    if nUnits > 1, subplot(ceil(nUnits/2), 2, ii); end;
    for jj = 4:4 % for each tutor syllable type
        figure
        for kk = 1:max(syllLabels) % for each empirical label
            isThisSyll = (syllLabels == kk);
            if ~any(isThisSyll), continue; end;

            subplot(1, max(syllLabels), kk)
            xxdat = stdDistMin(isThisSyll, jj);
            yydat = spikeRates{ii}(isThisSyll); %./ syllLengths(isThisSyll);
            sizeDat = sizeFromLength(isThisSyll);
            hh = scatter(xxdat, yydat, sizeDat, ... 
                colors(jj,:),'filled');
            xlabel('Distance (arbitrary)');
            ylabel('Spike rates (Hz)');
            title(sprintf('Syllable %s, Unit %d, %d syllables',labelTypes(jj)+('a'-1), ii, sum(isThisSyll)));
            % fit a trend line
            hold on;
            [linfit, ~,~,~, fitStats] = regress(yydat', [ones(numel(xxdat),1) xxdat]);
            fprintf('\ttutor syllable [%s] (# = %d), fit for trend, r^2 = %0.3g, F = %0.3g, p = %0.3g\n', ...
                labelTypes(jj)+('a'-1), sum(syllLabels == jj), fitStats(1), fitStats(2),fitStats(3));
            legend(sprintf('fit for trend, r^2 = %0.3g, F = %0.3g, p = %0.3g\n',...
                fitStats(1), fitStats(2),fitStats(3)));
            plot(xxdat, linfit(1) + xxdat * linfit(2), '--','Color',colors(jj,:),...
                'HandleVisibility', 'off');
            hold off;
        end
    end
    hold off;
    ylabel('spike rates (#)');
%    title(sprintf('Spiking rates against similarity, unit %d',ii));

end