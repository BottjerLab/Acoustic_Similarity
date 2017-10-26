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
            checkStereotypy(birdID, uAges(jj));
            set(gcf, 'Units', 'normalized','Position', [0 0 1 1])
            saveFile = sprintf('figures/clusterExamples/%s-%d.jpg',birdID, uAges(jj));
            fprintf('Saving to %s...', saveFile);
            export_fig(saveFile);
            fprintf('done.\n');
            close;
        catch err
            if strncmp(err.identifier, 'checkStereotypy:', 16)
                fprintf('%s\n', err.message);
                continue;
            else
                rethrow(err)
            end
        end        
    end
end
