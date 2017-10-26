function [thisAge, isPlastic, plasticScore] = getAgeOfSession(thisSession)
% maybe consider using persistent variable pattern - but sync will be issue
    xlFile = 'birdSummaries.xlsx';
        [~,~,rawData] = xlsread(xlFile,'Session Records','B2:E100');
    
    if nargout == 3
        [~,~,allPlasticTable] = xlsread(xlFile,'Day Records','A2:I100');
        allPlasticTable(:,3:end-1) = [];
        birdKeys = allPlasticTable(:,1);
        ageKeys = [allPlasticTable{:,2}]';
        allScores = [allPlasticTable{:,3}]';
    end
    sessionIDs = rawData(:,1);
    ages = rawData(:,3);
    plasticness = logical(strcmp(rawData(:,4),'Y')); % is the answer in the column Y?
    % cell-ify single strings
    if isstr(thisSession), thisSession = {thisSession}; end;
    
    thisAge   = NaN(size(thisSession));
    isPlastic = NaN(size(thisSession));
    plasticScore = NaN(size(thisSession));
    for ii = 1:numel(thisSession)
        isRightSession = strcmp(thisSession{ii}, sessionIDs);
        if sum(isRightSession) == 1 % should be unique
            thisAge(ii) = ages{isRightSession};
            isPlastic(ii) = plasticness(isRightSession);
            if nargout == 3
                thisBird = strtok(thisSession, '_');
                isRightBirdAge = (strcmp(birdKeys, thisBird) & ageKeys == thisAge(ii));
                if sum(isRightBirdAge) == 1
                    plasticScore(ii) = allScores(isRightBirdAge);
                end
            end
        end
    end
end