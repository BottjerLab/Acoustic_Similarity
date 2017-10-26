function [clusterQuality, allQualities] = getClusterQuality(queryBirdID, queryAge, queryClusterNum, isObjective)
if nargin < 4
    isObjective = false;
end
sheet = 'Cluster Quality';
if isObjective,
    sheet = 'Objective Cluster Quality';
end
    xlFile = 'birdSummaries.xlsx';
    [~,~,rawData] = xlsread(xlFile, sheet, 'A3:S38');
    
    keyBirdID = rawData(:,1);
    keyBirdAges = [rawData{:,2}]';
    allQualities = cell2mat(rawData(:,4:end)); % is the answer in the column Y?
    
    % cell-ify single bird ID strings
    if isstr(queryBirdID), queryBirdID = {queryBirdID}; end;
    
    clusterQuality = NaN(size(queryBirdID));
    for ii = 1:numel(queryBirdID)
        isRightSession = strcmp(queryBirdID{ii}, keyBirdID);
        isRightAge = (queryAge(ii) == keyBirdAges);
        matchedRecord = find(isRightSession & isRightAge);
        if numel(matchedRecord) == 1 % should be unique            
            clusterQuality(ii) = allQualities(matchedRecord, queryClusterNum(ii));
        elseif isempty(matchedRecord)
            warning('getClusterQuality:noMatch','Unmatched record for %s/age %d', ...
                queryBirdID{ii}, queryAge(ii));
        else
            error('getClusterQuality:multiMatch','Duplicated record');
        end
    end
end