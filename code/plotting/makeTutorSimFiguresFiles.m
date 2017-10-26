rep = reportOnData;
%%
for ii = 10:numel(rep) %start from Gy242
    sessions = {rep{ii}.sessionID};
    birdID = strtok(sessions{1}, '_');
    ages = getAgeOfSession(sessions);
    uAges = unique(ages(~isnan(ages)));
    startJ = 1;
    if ii == 4, startJ = 3; end;
    for jj = startJ:numel(uAges)
        compareTutorJuvenile(birdID, uAges(jj));
        close;
    end
end
