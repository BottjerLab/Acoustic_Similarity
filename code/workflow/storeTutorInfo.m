% make tutor spectrograms/diagrams to study
% note: for studying one bird, one wav only, use tutorDisplay
dataDir = [pwd filesep 'data' filesep];
%{
allReports = reportOnData;
allReports(cellfun('isempty',allReports)) = [];
birdID = cell(1,numel(allReports));
for ii = 1:numel(allReports)
    birdID{ii} = strtok(allReports{ii}(1).sessionID,'_');
end
%}
birdID = {'Y231'};
tutorID = getBirdTutor(birdID);

% 
for ii = 1:numel(birdID)
    birdDir = [dataDir birdID{ii}];
    mkdir(birdDir, 'tutor');
    
    % find the tutor song directory
    tutorSongDir = '..\..\Tutor Song';
    tutorExtraDir = [tutorSongDir filesep 'Tutors no longer in aviary']; 
    % check for existence of tutor song directory first
    if ~exist(tutorSongDir,'dir')
        error('%s does not exist as directory', tutorSongDir);
    end
    indivTutorDirs = dir(tutorSongDir); % these are only 
    findIdx = ~cellfun('isempty', strfind({indivTutorDirs.name}, tutorID{ii}));
    inExtra = false;
    thisTutorDir = [];
    if any(findIdx)
        thisTutorDir = [tutorSongDir filesep indivTutorDirs(findIdx).name];
    else
        indivTutorDirs = dir(tutorSongDir); % these are only
        findIdx = ~cellfun('isempty', strfind({indivTutorDirs.name}, tutorID{ii}));
        if any(findIdx)
            inExtra = true;
            thisTutorDir = [tutorExtraDir filesep indivTutorDirs(findIdx).name];            
        end
    end
    
    if isempty(thisTutorDir)
        warning('storeTutorInfo:nosong','Can''t find tutor %s''s song for bird %s', tutorID{ii}, birdID{ii});
        continue;
    end
    
    fprintf('Found directory thisTutorDir for tutor info for bird %s...\n', thisTutorDir, birdID{ii});
    
    % store the following:
    % one representative, concatenated in .mat format equivalent to spike2
    % imports
    tutorWavs = dir([thisTutorDir filesep '*.wav']);
    tutorWavs = {tutorWavs.name};
    % pick representative clips by manual audio marking
    nRegions = numel(tutorWavs);
    isMarked = false(1,nRegions);
    ii = 1; %#ok<FXSET>
    helpAsked = false;
    % method adapted from markRegions, with the exception that audio source
    % is individual clips
    fprintf('Mark clips with ''y'' that you would like to include as part of the tutor song analysis, press ! to quit at any time.\n');
    while ii <= nRegions
        % get clip
        songStruct = loadWavFile([thisTutorDir filesep tutorWavs{ii}]);
        clip = getClip(getWholeSongAsRegion(songStruct), songStruct);
        
        % get player object
        fs = 1/songStruct.interval;
        player = playSound(clip, fs, false);
        
        % request and parse input
        timeInfo = sprintf('%0.2fs', numel(clip)/fs);
        
        helpStr = '? for help';
        fullHelpStr = 'y to mark, u to undo, r to replay, , g to jump, ? for help, ! to abort';
        if helpAsked
            helpStr = fullHelpStr;
        end
        prompt = sprintf('Mark #%d/%d [%s] (%s, default is no mark) ? ', ii, nRegions, timeInfo, helpStr);
        char = input(prompt,'s');
        if isempty(char)
            char = 'n';
        end
        if length(char) > 1
            char = char(1);
        end
        helpAsked = false;
        % kill the soundplayer once input is received to allow shortcutting behavior
        stop(player);
        
        switch lower(char)
            case 'u'
                if ii > 1, ii = ii - 1; end; thisParams = params; continue;
            case 'r'
                continue;
            case '?'
                helpAsked = true;
                continue;
            case '!'
                fprintf('\nExiting...\n');
                break;
            case 'y'
                isMarked(ii) = true;
            case 'g'
                val = floor(input('Which clip to jump to? '));
                assert(val > 1 && val <= nRegions);
                ii = val;
                continue;
        end
        ii = ii + 1;
    end    
    if ~any(isMarked), 
        fprintf('No clips selected as examples of the tutor');
        continue;
    end
    
    % todo: truncate long silence?
    markedWavs = strcat([thisTutorDir filesep], tutorWavs(isMarked));
    nMarked = numel(markedWavs);
    tutorSylls = cell(1,nMarked);
    tSSpecs = cell(1,nMarked); tSFeats = cell(1,nMarked);
    durs = zeros(1,nMarked);
    for jj = 1:nMarked
        tmpStruct = loadWavFile(markedWavs{jj});
        wholeRegion = getWholeSongAsRegion(tmpStruct);
        fs = 1/tmpStruct.interval;
        % construct syllable boundaries, semi-automatically parsed and
        % labeled
        tutorParams = processArgs(defaultParams,'fs',fs,'preroll',0,'postroll',0);
        
        tmpSylls = parseRegionsIntoSyllables(tmpStruct, wholeRegion,tutorParams,...
            'doFilterNoise',false,'syllable.comboLength',18);
        
        tutorSylls{jj} = plotAndAdjust(tmpStruct,tmpSylls,wholeRegion, tutorParams, ...
            'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8,...
            'optGraphs',{'waveform', 'deriv','FM','AM', 'mTD'});
        
        [tSFeats{jj},tSSpecs{jj}]=getFeatures(tmpStruct,tutorSylls,tutorParams,...
            'plot',false,'verbose',true,'playsample',false,'doFilterNoise',false);
    
        durs(jj) = tmpStruct.length * tmpStruct.interval;
        keyboard
    end
    
    % concatenate the larger structures
    tutorStruct = loadManyWavFile(markedWavs, 0, fs);
    cumDurs = cumsum(durs);
    
    for jj = 2:nMarked
        tutorSylls{jj} = adjustTimeStamps(tutorSylls{jj}, cumDurs(jj), fs);
    end
    tutorSylls = vertcat(tutorSylls{:});
    tSFeats = vertcat(tSFeats{:});
    tSSpecs = vertcat(tSSpecs{:});
    keyboard
    %uisave({'tutorStruct','tutorSylls','tSFeats','tSSpecs'},...
    %    [matpath prependForSave('tutor-',matFile)]);
    
end