% splitAcceptedLabels
% take acceptedLabels files and split them by session
compositeReport = reportOnData('','',defaultParams,'quiet',true);
nBirds = numel(compositeReport);

%%
ageRecords = [];
for ii = 1:nBirds
    % get the whole report for one bird
    allBirdSessions = sortBy(compositeReport{ii}, 'sessionID');
    thisBirdID = strtok(allBirdSessions(1).sessionID, '_');
    % we want all syllables for subsong/plastic song statistics and label
    % clustering
   
    % simplify this script - no need to break down by sessions
    % also add subsong/plastic song distinctions and columns in excel for
    % comparison

    % get all the syllables for the bird
    DRsylls = [];
    allSyllFileName = ['data' filesep thisBirdID filesep 'allSpecs-' thisBirdID '.mat'];
    if exist(allSyllFileName, 'file') == 2
        load(allSyllFileName, 'DRsylls');
    else
        fprintf('Can''t find syllables for bird %s\n', thisBirdID);
        continue;
    end
    
    % what is the age of the sessions in the main file?
    [syllSessions, ~, syllSessIdxs] = unique({DRsylls.file});    
    for jj = 1:numel(syllSessions)
        [~,syllSessions{jj}] = fileparts(syllSessions{jj});
    end
    dataFilAges = getAgeOfSession(syllSessions);        
    
    uAges = unique(dataFilAges);
    for jj = 1:numel(uAges)    
        thisAge = uAges(jj);        
        % get syllables that correspond to this age
        sessIdxsOfThisAge = find(dataFilAges == thisAge);
        isThisAge = ismember(syllSessIdxs, sessIdxsOfThisAge);
        theseSylls = DRsylls(isThisAge);
        theseSyllIdxs = syllSessIdxs(isThisAge);
        % get the labels corresponding to this age
        try
            theseSyllLabels = loadAcceptedLabels(thisBirdID, thisAge);
        catch ME
            if strcmp(ME.identifier, 'loadAcceptedLabels:fileNotFound')
                fprintf('%s.\n',ME.message);
                continue;
            else
                rethrow(ME)
            end
        end
        if numel(theseSylls) ~= numel(theseSyllLabels)
            warning('Mismatch between label number and syllable number...\n');
            continue;
        end
        % split the labels into sessions here
        for kk = 1:numel(sessIdxsOfThisAge)
            thisSessionIdx = sessIdxsOfThisAge(kk);
            thisSession = syllSessions{thisSessionIdx};
            acceptedLabels = theseSyllLabels(theseSyllIdxs == thisSessionIdx);
            fileName = ['data' filesep thisBirdID filesep 'acceptedLabels-' thisSession '.mat'];
            fprintf('Saving accepted labels for %s to %s...\n', thisSession, fileName);
            save(fileName, 'acceptedLabels');
        end
    end
end