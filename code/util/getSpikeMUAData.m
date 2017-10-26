function metaData = getSpikeMUAData
    xlFile = 'birdSummaries.xlsx';
    %[~,~,rawData]=xlsread(xlFile,'Multi-unit Clusters','','basic'); % <--
    %necessary for excel 2013
    [~,~,rawData]=xlsread(xlFile,'Multi-unit Clusters',''); % <-- only works with excel 2010
    % workaround would be to get excel version by checking the 
    % Excel = actxserver('excel.application'); version property.
    % remove first row, which is header
    rawData(1,:) = [];
    
    % get rid of NaN rows
    % output columns: birdID, sessionID+channels, MUA clusters 
    
    nMUA = size(rawData,1);
    isRep = false(1,nMUA);
    pastPoint = 1;
    % note: this assumes some sorted/connected order in the excel sheet
    % (i.e. same entries are adjacent)    
    for ii = 1:nMUA
        if numel(rawData{ii,1})==1 && isnan(rawData{ii,1}) 
            % this is where the data ends
            rawData(ii:end,:) = [];
            break;
        end
        rawData{ii,2} = [rawData{ii,1} '_' rawData{ii,2}];
        sessionID{ii} = rawData{ii,2}(1:find(rawData{ii,2}=='_',1,'last')-1);
        if ii > 1 && strcmp(rawData{ii-1,2}, rawData{ii,2})
            isRep(ii) = true;
        elseif ii > 1
            rawData{pastPoint,3} = [rawData{pastPoint:(ii-1),3}];
            pastPoint = ii;
        end
    end
    % remove repetitions
    
    rawData(isRep,:) = [];
    sessionID(isRep) = []; 
 
    % import into struct array format for named fields
    metaData = struct('birdID',rawData(:,1), 'sessionID', sessionID', ...
        'trodeID', rawData(:,2), 'MUAunits', rawData(:,3));
    
end