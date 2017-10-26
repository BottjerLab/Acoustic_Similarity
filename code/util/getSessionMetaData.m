function metaData = getSessionMetaData
    xlFile = 'birdSummaries.xlsx';
    [~,~,rawData]=xlsread(xlFile,'Session Records');

    % remove first row, which is header
    rawData(1,:) = [];
    
    % input columns: birdID, sessionID, recording date, age @ recording,
    % extra unnecessary columns...
    
    % output columns: birdID, sessionID, age @ recording
    metaData = struct('birdID',rawData(:,1), 'sessionID', rawData(:,2), 'age', rawData(:,4));
end