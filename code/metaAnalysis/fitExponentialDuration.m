%% This script examines the syllable durations over sessions
function [kstatValues, ages] = fitExponentialDuration(birdID)
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
sR = initEmptyStructArray({'sessions', 'age', 'kstat', 'nIncluded', 'tPrime', 'adjKS', 'adjIncl'}, nAges);
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
%% plotting, not always necessary
%{
nR = floor(sqrt(nSessions));
nC = ceil(nSessions / nR);
figure(1)
for ii = 1:nSessions
    subplot(nR, nC, ii);
    syllLengths = plotDurationDistr(sylls); 
    title(report(ii).sessionID,'interpreter', 'none');
    xlim([0 0.5])
end
%}
%% do subsong/plastic song analysis with correction for sample size (exponential, maybe MOG/log-normals)
% distribution is roughly exponential from 25-400 ms
expRegion = [0.025 0.4]; % seconds
%progressbar(sprintf('%s: Ages', birdID), 'Trials');
for ii = 1:nAges
    syllLens = syllLensPerDate{ii};
    nPop = numel(syllLens);
    
    [sR(ii).kstat, sR(ii).nIncluded, sR(ii).tPrime] ...
        = ksfitstat(syllLens, expRegion);
    
    %{
    sampleNumbers = 100 * 2.^(0:10);
    nTrials = 20;
    ksTrial = zeros(nTrials, numel(sampleNumbers));
    tTrial  = zeros(nTrials, numel(sampleNumbers));
    for jj = 1:nTrials
        for kk = 1:numel(sampleNumbers)
            oversamp = randi(nPop,1,sampleNumbers(kk));
            bootstrapSample = syllLens(oversamp);
            [ksTrial(jj,kk), ~, tTrial(jj,kk)] = ...
                ksfitstat(bootstrapSample, expRegion);
        end
        progressbar([],jj/nTrials);
    end
    
    fitParams = polyfit(log(sampleNumbers), log(mean(ksTrial)),1);
    surrogate = 0;
    sR(ii).adjKS = exp(fitParams(2)); %exp(polyval(fitParams, log(surrogate)));
    sR(ii).adjIncl = surrogate;
    %keyboard
    progressbar(ii/nAges,[]);
        %}
end
%%
 % take out day 48
 %{
 if strcmp(birdID,'Y231')
     ageToEliminate = find(birdAge == 48);
     uDate(ageToEliminate) = [];
     syllLensPerDate(ageToEliminate) = [];
     nInInterval(ageToEliminate) = [];
     tPrime(ageToEliminate) = [];
     kstat(ageToEliminate) = [];
     birdAge(ageToEliminate) = [];
 end
 %}
%% plotting
%{

% histogram binning for display purposes
binWidth = 0.003;
syllBins = 0:binWidth:0.65; 

% prep axes
hax = zeros(1,nAges);
nR = nAges; nC = 1;
for ii = 1:nAges   
    hax(ii) = subplot(nR,nC,ii);
    hold off;
    % make axis wider
    pos = get(gca,'Position');
    set(gca, 'Position', [0.1 pos(2) 0.8 pos(4)]);
    
    % x-coordinates
    binCenters = syllBins + binWidth/2;
    
    % histogram
    syllDatePDF = hist(syllLensPerDate{ii},syllBins);
    bar(binCenters, syllDatePDF, 1, 'FaceColor', [0 0 0])
    hold on;
    % trend line
    distrFunc = @(x,t) 1/t * (exp(-expRegion(1)/t) - exp(-expRegion(2)/t)).^(-1) * exp(-x/t); 
    plot(binCenters, sR(ii).nIncluded * binWidth * distrFunc(binCenters, sR(ii).tPrime), 'r-',...
        'LineWidth',2);
    
    % axis labels
    fontN = 'Arial';
    if ii == numel(uDate)
        xlabel('syllable length (ms)','FontName', fontN, 'FontSize',14);
    end
    ylabel({'syllable count',sprintf('age %d dph', sR(ii).age)},'FontName', fontN, 'FontSize',14);
    xlim([0 0.4]);    
    
    % text labels

    text(0.32, ylim * [0.4 0.6]', ...
        sprintf('Score = %0.2f, (# = %d)\nAdj score = %0.2f, (# = %d)', ...
            sR(ii).kstat, sR(ii).nIncluded, sR(ii).adjKS, sR(ii).adjIncl),...
        'HorizontalAlignment', 'right', 'VerticalAlignment','middle',...
        'FontName', fontN, 'FontSize',14);   
    % appearance tweaks
    set(gca,'TickLength',[0 0], 'FontName', fontN, 'FontSize',14, ...
        'Box', 'off');    
    
end
%subplot(hax(2)); ylim([0 150]); %?
set(gcf,'Color', [1 1 1]);
subplot(hax(1));
    %}
%%
kstatValues = [sR.kstat];
ages = [sR.age];
%%+
%saveCurrFigure([pwd filesep 'figures' filesep 'dailySyllDur-' birdID '.pdf']);