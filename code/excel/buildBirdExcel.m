% buids summary of all data for excel visualization

% note: should be a script b/c xlswrite1 uses Excel from the base but we
% don't want to keep it open 
% get the manifests
compositeReport = reportOnData('','',defaultParams,'quiet',true);
nBirds = numel(compositeReport);

% get MUA info
if ~exist('MUAinfo','var')
    MUAinfo = getSpikeMUAData;
end

%% open Excel server
% todo: make sure excel server is open 
Excel = actxserver('Excel.Application');

% open the specific file (note: the file must exist)
xlFilename = [pwd filesep 'birdSummaries.xlsx'];
if ~exist(xlFilename,'file')
    error('File %s does not exist...', xlFilename);
end
%Workbook = invoke(Excel.Workbooks,'Open', xlName);
Workbook = Excel.Workbooks.Open(xlFilename);

% example - get the sheet names
nSheets = Workbook.Sheets.Count;
sheetNames = cell(1,nSheets);
% seems like array indexing has to be done these ways instead of more
% naturally with operator()
for ii = 1:nSheets
    sheetNames{ii} = get(get(Workbook.Sheets, 'Item', ii),'Name');
end
%%
sessSheet = Item(Workbook.Sheets, find(strcmp('Session Records', sheetNames)));
sessSheet.Activate;

%%
% put the useful info in 'Session Records' worksheet 
% Columns = 
% Bird ID, Session ID, Date of Recording, Age, #syllables, song time, 
% #neurons and core/shell SUA/MUA

% pointer counter for the rows
% each Row = one session
rowCount = 2;

for ii = 1:nBirds
    % get the whole report for one bird
    allBirdSessions = compositeReport{ii};
    
    % sort the manifests by date first...    
    allBirdSessions = sortBy(allBirdSessions, 'sessionID');
    
    % always the same bird ID in this loop
    thisBirdID = strtok(allBirdSessions(1).sessionID, '_');
    for jj = 1:numel(allBirdSessions)
        thisSession = allBirdSessions(jj).sessionID;
        thisManifest = allBirdSessions(jj).manifest;
        theseSpikeFiles = allBirdSessions(jj).spikeFiles;
        nUnits=0; isMUAUnit = false(0); isCoreUnit = false(0);
        
        % only list with neurons
%        if isempty(theseSpikeFiles), continue; end;
 
        % only list with manual approved syllables
        mS = loadFromManifest(thisManifest,'approvedSyllables');
               
        if isempty(mS), continue; end;
        
        % how many neurons do we have? what kind? (core/shell, SUA/MUA)
        if ~isempty(theseSpikeFiles)
            isCoreFile = ~cellfun('isempty',strfind(theseSpikeFiles, 'core'));
            isShellFile = ~cellfun('isempty',strfind(theseSpikeFiles, 'shell'));

            fprintf('session = %s... \n',thisSession);
            spikes = {}; nNPerFile = []; isMUAUnit = [];
            [spikes, nNPerFile, isMUAUnit] = loadSpikeData(theseSpikeFiles, MUAinfo);
            nUnits = numel(spikes);
 
            isCoreUnit = false(nUnits, 1);
            
            cumIndex = [0 cumsum(nNPerFile)];
            for kk = 1:numel(theseSpikeFiles)
                unitsFromFile = (cumIndex(kk)+1):cumIndex(kk+1);
                
                % core/shell ness depends on the file folder
                isCoreUnit(unitsFromFile) = ~isempty(strfind(theseSpikeFiles{kk}, 'core'));
            end
        end
        
        
        singLength = sum([mS.stop] - [mS.start]);
        
        % now write to the excel file: columns A and B    
        xlswrite1(xlFilename, {thisBirdID, thisSession}, 'Session Records', ...
            sprintf('A%d:B%d',rowCount,rowCount));

        % now write to the excel file: columns F-P
        statLine = [numel(mS), singLength, ...
            nUnits, sum(isCoreUnit & ~isMUAUnit), sum(~isCoreUnit & ~isMUAUnit),...
                    sum(isCoreUnit &  isMUAUnit), sum(~isCoreUnit & isMUAUnit)];
        xlswrite1(xlFilename, statLine, ...
            'Session Records', sprintf('F%d:%c%d',rowCount,'F' + numel(statLine) - 1, rowCount));

        % increment after each session ...
        rowCount = rowCount + 1;            
    end    
end
%% formatting:
% autofit columns
for ii = 1:15 % todo: auto find # columns
    invoke(Item(sessSheet.Columns,ii),'Autofit');
end

%% wrap up
fprintf('Please review and save AFTER unpausing to accept changes to [%s] (paused)...\n', xlFilename);
Excel.visible = true;
pause;
%% close up shop

Workbook.Close;
Excel.Quit;

%% clean up after ourselves in base workspace
clear Workbook rangeStr rowCount mS mB mM birdID sessionID ...
    sortIdx allBirdSessions sessionIDs compositeReport nSheets sheetNames nBirds ...
    sessionSheet
