
%% temporary script to rewrite clusters with consistent distances in cluster folders, not subfolders
rep = reportOnData; 
for ii = 10:numel(rep)
    birdID = strtok(rep{ii}(1).sessionID,'_')
    ages = getAgeOfSession({rep{ii}.sessionID});
    ages = ages(~isnan(ages));
    ages = unique(ages);
    for jj = 5:numel(ages)
        if strcmp('Gy242', birdID) && jj == 2, continue; end %too many syllables in this case
        fullClusterByDay(birdID, ages(jj)); % problem with this is sometimes folder gets locked out for writing
    end
end