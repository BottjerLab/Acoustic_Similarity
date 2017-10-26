%% This script examines the syllable durations over sessions
function [ksamples, sampNumbers] = oversampleTest(birdID, selAge)
%birdID = 'Gy242';
report = reportOnData(birdID,'',[],'verbose',false);

% only look for sessions with approved syllables
nSessions = numel(report);
hasData = false(1,nSessions);
for ii = 1:nSessions
    hasData(ii) = any(findInManifest(report(ii).manifest, {'approvedSyllables'})); %,'bosSyllables','manualSyllables','syllables'}));
end
report = report(hasData);
nSessions = numel(report);
%% bookkeeping
dateStr = cell(size(report));
for ii = 1:nSessions
    % get the date for each session
    sess = report(ii).sessionID;
    delimIdx = strfind(sess, '_');
    dateStr{ii} = sess((delimIdx(1)+1):(delimIdx(3) - 1));
end

% organize data - not every session has a date (in  the excel file), 
% but every session has a recording date embedded in the name
possAges = getAgeOfSession({report.sessionID}); 
[uDate, ~, dateIdx] = unique(dateStr);

nAges = numel(uDate);
% syllable Records
sR = initEmptyStructArray({'sessions', 'age', 'kstat', 'nIncluded', 'tPrime'}, nAges);
for ii = 1:nAges
    foo = possAges(dateIdx == ii); foo(isnan(foo)) = [];
    if isempty(foo), error('sylDurationByBird:noAgeForSession', 'Age not found...'); end;
    sR(ii).age = foo(1);
end

%% gather syllable length empirical distributions 
fprintf('Session %s, %d sessions with data...\n', birdID, nSessions)
syllLensPerDate = cell(1,nAges);
for ii = 1:nSessions
    sylls = loadFromManifest(report(ii).manifest, 'approvedSyllables');
    syllLengths = [sylls.stop]-[sylls.start];
    syllLensPerDate{dateIdx(ii)} = [syllLensPerDate{dateIdx(ii)} syllLengths]; 
    sR(dateIdx(ii)).sessions = [sR(dateIdx(ii)).sessions report(ii).sessionID];
end
%% select only the age we want
syllLensPerDate = syllLensPerDate{[sR.age] == selAge};
sR([sR.age] ~= selAge) = [];
nAges = 1;
%% oversampling mini-experiment
expRegion = [0.025 0.4]; % seconds
sampNumbers = 100 * 2.^(0:15);
trials = 50;
ksamples = zeros(numel(sampNumbers), trials);
nPop = numel(syllLensPerDate);
progressbar(0)
for jj = 1:trials
    for ii = 1:numel(sampNumbers)
        oversamp = randi(nPop,1,sampNumbers(ii));
        boots = syllLensPerDate(oversamp);
        ksamples(ii,jj) = ksfitstat(boots, expRegion);
    end
    progressbar(jj/trials)
end
