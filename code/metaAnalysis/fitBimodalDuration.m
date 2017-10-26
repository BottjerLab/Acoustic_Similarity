%% This script examines the syllable durations over sessions
function [peakSize, ages, fitQ] = fitBimodalDuration(birdID)
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
sR = initEmptyStructArray({'sessions', 'age', 'peakheight', 'peakfit'}, nAges);
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

tau = zeros(1,nAges);
for ii = 1:nAges
    syllLens = syllLensPerDate{ii};
    binWidth = 0.003; % in seconds
    syllBins = 0:binWidth:0.65;
    
    % get the fit from 200-400 milliseconds
    lensDistr = hist(syllLens, syllBins) / numel(syllLens) / binWidth;
    [tau(ii), nFit] = fitExpOnInterval(syllLens, [0.2 0.4]);
    expFit = exp(-syllBins/tau(ii)) / tau(ii);

    % fit a gaussian to the residual    
    residDistr = lensDistr - expFit;
    
    zeroedDistr = residDistr; zeroedDistr(zeroedDistr<0) = 0;
    [cfun, opts] = fit(syllBins', zeroedDistr', 'gauss1');
    sR(ii).peakheight = cfun.a1;
    sR(ii).peakfit = opts.rsquare;
    
    subplot(nAges, 3, 1 + 3*(ii-1));
    plot(syllBins, lensDistr, 'k-', syllBins, expFit, 'r-');    
    ylabel('Prob. density (s^{-1})');
    title(sprintf('%s, age %d (fit on %d)', birdID, sR(ii).age, nFit));
    if ii == nAges, xlabel('Syllable duration (ms)'); end
    set(gca,'Box', 'off');
    ylim([0 15])
    xlim([0 0.5]);

    subplot(nAges, 3, 2 + 3*(ii-1));
    semilogy(syllBins, lensDistr, 'k-', syllBins, expFit, 'r-');
    set(gca,'Box', 'off');
    xlim([0 0.5]);
    set(gca,'YTick',[0.01 0.1 1 10]);
    ylim([0.01 15]);
    
    subplot(nAges, 3, 3*ii)
    plot(syllBins, lensDistr, 'k-', syllBins, expFit'+cfun(syllBins), 'b-');
    title(sprintf('Peak = %0.2f, fit = %0.2f', sR(ii).peakheight, sR(ii).peakfit)); 
    if ii == nAges, xlabel('Syllable duration (ms)'); end
    set(gca,'Box', 'off');
    ylim([0 15])
    xlim([0 0.5]);
    
    % try to fit
end
set(gcf,'Color', [1 1 1]);

peakSize = [sR.peakheight];
ages = [sR.age];
fitQ = [sR.peakfit];
