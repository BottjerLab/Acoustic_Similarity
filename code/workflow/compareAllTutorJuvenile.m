birdRep = reportOnData;
for ii = 1:numel(birdRep)
    birdID = strtok(birdRep{ii}(1).sessionID,'_');
    sessions = {birdRep{ii}.sessionID};
    ages = getAgeOfSession(sessions);
    ages = unique(ages(~isnan(ages)));
    for jj = 1:numel(ages)
        try
            fprintf('Writing for %s, age %d...\n', birdID, ages(jj));
            compareTutorJuvenile(birdID, ages(jj), defaultParams, 'plot',false);
        catch err
            fprintf('ERROR: %s\n',err.message);
        end        
    end
end