function [tutors, tutorFolders] = getBirdTutor(lookupBird)
% returns tutor IDs and tutor folder names for a collection of birdIDs 
% maybe consider using persistent variable pattern - but sync will be issue
    xlFile = 'birdSummaries.xlsx';
    [~,~,rawData] = xlsread(xlFile,'Bird Records');
    
    % remove headers: 
    rawData(1,:) = [];
    
    birdIDs = rawData(:,1); % always first column
    tutorIDs = rawData(:,5); % always fifth column
    tutorDests = rawData(:,7); % always seventh column 
    % cell-ify single strings
    isSingle = false;
    if ischar(lookupBird), isSingle = true; lookupBird = {lookupBird}; end
    
    tutors = cell(size(lookupBird));
    for ii = 1:numel(lookupBird)
        isRightBird = strcmp(lookupBird{ii}, birdIDs);
        if sum(isRightBird) == 1 % should be unique
            tutors{ii} = tutorIDs{isRightBird};
            tutorFolders{ii} = strcat('..\..\Tutor Song\', tutorDests{isRightBird});
        end        
    end
    if isSingle
        tutors = tutors{1};
        tutorFolders = tutorFolders{1};
    end
end