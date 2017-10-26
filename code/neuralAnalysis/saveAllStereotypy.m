% plots stereotypical and unusual syllables for all sessions
reps = reportOnData;
neuronData = [];
nBirds = numel(reps);
%%
for ii = 1:nBirds
    birdID = strtok(reps{ii}(1).sessionID, '_');
    sessions = reps{ii};
    sessionIDs = {sessions.sessionID};
    ages = getAgeOfSession(sessionIDs);
    
    uAges = unique(ages(~isnan(ages)));
    nAges = numel(uAges);
    for jj = 1:numel(uAges)
        try
            saveStereotypy(birdID, uAges(jj));
            
        catch err
            if strncmp(err.identifier, 'saveStereotypy:', 15)
                fprintf('%s\n', err.message);
                continue;
            else
                rethrow(err)
            end
        end        
    end
end
