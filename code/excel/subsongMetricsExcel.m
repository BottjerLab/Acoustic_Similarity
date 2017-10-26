% buids summary of the birds by day for excel visualization
% syllable counts and metrics about plastic/subsong 

% note: this file is a script b/c xlswrite1 uses Excel from the base but we
% don't want to keep it open
% get the manifests
compositeReport = reportOnData('','',defaultParams,'quiet',true);
nBirds = numel(compositeReport);

% get MUA info
if ~exist('MUAinfo','var')
    MUAinfo = getSpikeMUAData;
end

%%
ageRecords = [];
for ii = 1:nBirds
    % get the whole report for one bird
    allBirdSessions = sortBy(compositeReport{ii}, 'sessionID');
    thisBirdID = strtok(allBirdSessions(1).sessionID, '_');
    % we want all syllables for subsong/plastic song statistics and label
    % clustering
   
    % simplify this script - no need to break down by sessions
    % also add subsong/plastic song distinctions and columns in excel for
    % comparison

    DRsylls = [];
    allSyllFileName = ['data' filesep thisBirdID filesep 'allSpecs-' thisBirdID '.mat'];
    if exist(allSyllFileName, 'file') == 2
        % get the syllables for the whole bird
        load(allSyllFileName, 'DRsylls');
    else
        fprintf('Can''t find syllables for bird %s\n', thisBirdID);
        continue;
    end
    
    % what is the age of the sessions in the main file?
    [syllSessions, ~, syllSessIdxs] = unique({DRsylls.file});    
    for kk = 1:numel(syllSessions)
        [~,syllSessions{kk}] = fileparts(syllSessions{kk});
    end
    dataFilAges = getAgeOfSession(syllSessions);        
    
    uAges = unique(dataFilAges);

    thisAgeRecord = initEmptyStructArray(...
        {'birdID', 'age', 'nSylls', 'tSylls', 'nLabeled', 'tLabeled',...
        'ksfit', 'peakHeight', 'motifVar', 'syllMotifRatio'}, numel(uAges));
    for jj = 1:numel(uAges)    
        thisAge = uAges(jj);        
        % get syllables that correspond to this age
        sessIdxsOfThisAge = find(dataFilAges == thisAge);
        theseSylls = DRsylls(ismember(syllSessIdxs, sessIdxsOfThisAge));
        thisAgeRecord(jj).birdID = thisBirdID;
        thisAgeRecord(jj).age = thisAge;
        thisAgeRecord(jj).nSylls = numel(theseSylls);
        thisAgeRecord(jj).tSylls = sum([theseSylls.stop] - [theseSylls.start]);
        % get the labels corresponding to this age
        try
            theseSyllLabels = loadAcceptedLabels(thisBirdID, thisAge);
        catch ME
            if strcmp(ME.identifier, 'loadAcceptedLabels:fileNotFound')
                fprintf('%s.\n',ME.message);
                thisAgeRecord(jj).nLabeled = 0;
                thisAgeRecord(jj).tLabeled = 0;
                continue;
            else
                rethrow(ME)
            end
        end
        if numel(theseSylls) ~= numel(theseSyllLabels)
            warning('Mismatch between label number and syllable number...\n');
            continue;
        end
        thisAgeRecord(jj).nLabeled = sum(~isnan(theseSyllLabels));
        thisAgeRecord(jj).tLabeled = sum([theseSylls(~isnan(theseSyllLabels)).stop] - ...
            [theseSylls(~isnan(theseSyllLabels)).start]);
    end    
    % now metrics:

    % exponential fit metric
    kstat = fitExponentialDuration(thisBirdID);
    num2cell(kstat); [thisAgeRecord.ksfit] = ans{:};
    
    % peak above exponential metric (NB: should implement multi-peak)
    peakHeight = fitBimodalDuration(thisBirdID); % <- rename this too
    num2cell(peakHeight); [thisAgeRecord.peakHeight] = ans{:};
    
    % motif length variance
    % also average ratio of syllables between last and first two in motif
    % and coefficient of variation of motif length
	[vM, mSR, cvM] = motifVariation(thisBirdID);
    
    num2cell(vM);  [thisAgeRecord.motifVar]       = ans{:};
    num2cell(mSR); [thisAgeRecord.syllMotifRatio] = ans{:};
    num2cell(cvM); [thisAgeRecord.motifCV] = ans{:};  

        % TODO: these two metrics look into approvedSyllables, which in
    % the case of Gy242 may not be all the syllables in DRsylls.
    % this is because the clustering for Gy242 takes a long time -
    % I'll wait until a weekend to redo the clusters / correct the
    % discrepancy but it most cases this won't change

    ageRecords = [ageRecords; thisAgeRecord];
end

save('data/birdAgeRecords.mat','ageRecords');
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

xlswrite1(xlFilename, {'Bird ID', 'Age', '# sylls', 'Elapsed sylls', '# labeled',...
    'Elapsed labeled', '% # labeled', '% time labeled', ...
    'Exponential distribution fit', 'Gaussian peak height', ...
    'Var motif length (s)', 'CV of motif length','Ratio of late/early sylls'}, 'Day Records', ...
    'A1:M1');

for ii = 1:numel(ageRecords)
    thisRecord = ageRecords(ii);
    pctLabeled = thisRecord.nLabeled/thisRecord.nSylls;
    if isnan(pctLabeled), pctLabeled = 0; end;
    pctTLabeled = thisRecord.tLabeled/thisRecord.tSylls;
    if isnan(pctTLabeled), pctTLabeled = 0; end;
    
    % write data to the 'Day Records' sheet
    xlswrite1(xlFilename, {thisRecord.birdID, thisRecord.age, thisRecord.nSylls, ...
        thisRecord.tSylls, thisRecord.nLabeled, thisRecord.tLabeled, ...
        pctLabeled, pctTLabeled, thisRecord.ksfit, thisRecord.peakHeight,...
        thisRecord.motifVar, thisRecord.motifCV, thisRecord.syllMotifRatio}, 'Day Records', ...
        sprintf('A%d:M%d',ii+1,ii+1));    
end
%{
for ii = 1:numel(compositeReport)
    if ~isfield(compositeReport{ii},'age'), continue; end;
    thisBirdID = strtok(compositeReport{ii}(1).sessionID, '_');
    % sum up by ages
    uAges = unique([compositeReport{ii}.age]);
    uAges(isnan(uAges)) = [];
    for jj = 1:numel(uAges)
        thisAge = uAges(jj);
        
        rightAges = compositeReport{ii}([compositeReport{ii}.age] == thisAge);
        flds = {'nSylls','lenSylls','nLabeled','lenLabeled'};
        summary = initEmptyStructArray(flds,1);
        for kk = 1:numel(flds);
            summary.(flds{kk}) = sum([rightAges.(flds{kk})]);
        end
        
        % write values in a row to the excel file
        
        % increment after each session ...
        rowCount = rowCount + 1;
    end
end
%}
%% wrap up
fprintf('Please review and save AFTER unpausing to accept changes to [%s] (paused)...\n', xlFilename);
Excel.visible = true;
pause;
%% close up shop

Workbook.Close;
Excel.Quit;

%% clean up after ourselves in base workspace
clear Workbook Excel
