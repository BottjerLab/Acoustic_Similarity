function [stats, statFamily] = findSpikeCorrelates(events, statistics, evR)
%function findSpikeCorrelates(events, statistics, spikeTimes)
% events is a struct array of events
% statistics is a struct array of scalars, one for each event
% evSpikeRates gives he 
% old: spikeTimes is a vector
% test if spiking is significantly higher/lower during the event than
% during other times
% evRS is the normalized response after the baseline is removed

% these commands are kept because they correspond to the use case by
% juvenileCorrelation*, but the analysis part of those files is obsoleted.
%maxTime = max(spikeTimes); % an approximation
%evSpikeCounts = countSpikes(events, spikeTimes);
%evLengths = [events.stop] - [events.start];


% do a quick test on spiking within/without 

% todo: replace evSpikeRates with an input argument, response strengths,
% based on proper baselines.
%evSpikeRates = evSpikeCounts ./ evLengths;
iterWarn = warning('off', 'stats:statrobustfit:IterationLimit');
%conditionedWarn = warning('off', 'stats:glmfit:IllConditioned');
%regressWarn  = warning('off','stats:regress:RankDefDesignMat');
evSpikeRates = evR;
if size(evSpikeRates, 1) == 1
    evSpikeRates = evSpikeRates';
end
%basalRate = (numel(spikeTimes) - sum(evSpikeCounts)) / (maxTime - sum(evLengths));
%rateSEM = std(evSpikeRates)/sqrt(numel(events)-1);
%[~,tp]=ttest(evSpikeRates);

%fprintf('Basal rate %0.2f Hz, event rate %0.2f +/- %0.2f Hz, (p = %0.3f)\n', ...
%    basalRate, mean(evSpikeRates), rateSEM, tp);

%%%% test significance for each statistic
statNames = fieldnames(statistics); 
nFeats = numel(statNames); % number of features does not include constant term
dataTable = cellfun(@(x) [statistics.(x)], statNames, 'UniformOutput', false);
dataTable = vertcat(dataTable{:})';

nSylls = numel(evR);
dataTableWO = [dataTable ones(nSylls, 1)]; % with constant ones

% remove linearly dependent columns for linear regression 
[~,jb] = rref(dataTableWO, 1e-5);
indpCols = false(nFeats+1,1);
indpCols(jb) = true;
fprintf('Redundant features removed: %d of %d features retained...\n', ...
    length(jb), nFeats);
statNames = [statNames; 'intercept']; 

% for figuremaking
%nr = round(sqrt(nStats));
%nc = ceil(nStats/nr);
%dataTableWOnes = [ones(nSylls,1) dataTable];

% these statistics aren't corrected for collinearity, so these are kept for
% historical purposes (used by juvenileCorrelation*)
%[loads, ~, fitStats] = glmfit(dataTable(:, indpCols), evSpikeRates, 'normal'); 
%stats.pCorr = (fitStats.p(2:end) .* sign(loads(2:end)))';
%stats.tCorr = (fitStats.t(2:end) .* sign(loads(2:end)))';

init = false;
outliers = false(nSylls,1);
nPrevOutliers = 0;

while ~init || sum(outliers) > nPrevOutliers % any new outliers?

    init = true;
    fullMdl = LinearModel.fit(dataTableWO(:, indpCols), evSpikeRates, eye(length(jb), length(jb)+1), ... % include constant value on own terms
        'RobustOpts', 'on', 'VarNames', [statNames(indpCols); 'RS'], 'Exclude', outliers); % response strength is dependent value
    
    % find outliers if the cooks distance is way out from others
    nPrevOutliers = sum(outliers);
    outliers = outliers | ...
        (fullMdl.Diagnostics.CooksDistance > 8 * nanstd(fullMdl.Diagnostics.CooksDistance) | ...
         isinf(fullMdl.Diagnostics.CooksDistance));
     
    % remove new dependent features
    [~,jb] = rref(dataTableWO(~outliers,:), 1e-5);
    indpCols(:) = false; indpCols(jb) = true;
end
%fprintf('Removing %d outliers...\n',sum(outliers));

stats.pCorr(jb) = fullMdl.Coefficients.pValue .* sign(fullMdl.Coefficients.Estimate); %JMA added in (jb) to pCorr and tCorr and changed .pValue etc from (2:end)
stats.tCorr(jb) = fullMdl.Coefficients.tStat .* sign(fullMdl.Coefficients.Estimate);
stats.RSquare = fullMdl.Rsquared.Ordinary; %JMA added

anovaModel = anova(fullMdl,'summary');
stats.fAll = anovaModel.F(2);
stats.pAll = anovaModel.pValue(2);
stats.dfAll = anovaModel.DF(2:3)';

% these lines manually calculate the total F statistic for regression
% significance
%{
ssResidual = sum(fitStats.resid.^2); glmSSE = ssResidual;
ssTotal = var(evSpikeRates) * (numel(evSpikeRates)-1);
dfResid = fitStats.dfe;
dfRegress = nSylls - 1 - dfResid;
stats.fAll = ((ssTotal - ssResidual) / dfRegress) / (ssResidual / dfResid);
stats.pAll = 1 - fcdf(stats.fAll,dfRegress, dfResid); % these calculations are the same as 
stats.dfAll = [dfRegress, dfResid];
%}
fprintf('Full correlation probability: F(%d,%d) = %0.3f, p = %0.3g...\n', ...
    stats.dfAll(1), stats.dfAll(2), stats.fAll, stats.pAll);

% get some partial F scores for each individual feature by using nested
% models

stats.fPartial = zeros(1, nFeats);
stats.pPartial = zeros(1, nFeats);
for ii = 1:nFeats
    if ~indpCols(ii), continue; end
    
    %retainedCols = indpCols; retainedCols(ii) = false;
    %partialMdl = LinearModel.fit(dataTable(:, retainedCols), evSpikeRates, 'RobustOpts', 'on');
    partialMdl = removeTerms(fullMdl, statNames{ii});
    stats.fPartial(ii) = (partialMdl.SSE - fullMdl.SSE) / (fullMdl.SSE / fullMdl.DFE); % still might have outliers...
    stats.pPartial(ii) = 1 - fcdf(stats.fPartial(ii), 1, partialMdl.DFE);
end

if nargout >= 2
    % get the partial F scores for a set of features - more useful
    % family of features is all features that have a certain prefix
    [families, ~, familyIdx] = unique(strtok(statNames(1:end-1),'_'));
    familyIdx(end+1) = NaN;
    nFamilies = max(familyIdx);
    statFamily = initEmptyStructArray({'feature', 'Fpartial', 'Ppartial', 'df', 'R2','F','P'}, nFamilies);
    for ii = 1:nFamilies
        if all(~indpCols(familyIdx==ii)), continue; end
        
        %noFamilyCols = indpCols; noFamilyCols(familyIdx==ii) = false;
        
        %familyMdl = LinearModel.fit(dataTable(:, noFamilyCols), evSpikeRates, 'RobustOpts', 'on');
        familyStr = strcat(statNames(indpCols & familyIdx==ii), '+'); 
        familyStr = [familyStr{:}];
        familyStr(end) = [];
        familyOutMdl = removeTerms(fullMdl, familyStr);
        statFamily(ii).feature = families{ii};
        
        %[~,~,tmpStats] = glmfit(dataTable(:,noFamilyCols), evSpikeRates, 'normal');
%        ssPartialModel = sum(tmpStats.resid.^2);
        
        df = [sum(indpCols(familyIdx==ii)), familyOutMdl.DFE];    
        statFamily(ii).df = df;
        statFamily(ii).Fpartial = ((familyOutMdl.SSE - fullMdl.SSE) / df(1)) / (fullMdl.SSE / df(2));
        statFamily(ii).Ppartial = 1 - fcdf(statFamily(ii).Fpartial, df(1), df(2));
%        statFamily(ii).Fpartial = ((ssPartialModel - ssResidual) / df(1)) / (ssResidual / df(2));
%        statFamily(ii).Ppartial = 1 - fcdf(statFamily(ii).Fpartial, sum(familyIdx==ii), dfResid);
        
%        familyData = [ones(nSylls,1) dataTable(:, familyIdx == ii)];
        
        % 
        familyOnlyMdl = LinearModel.fit(dataTable(:, indpCols(familyIdx==ii)), evSpikeRates);
        %[~,~,~,~,familyOnlyStat] = regress(evSpikeRates, familyData);
        familyOnlyAnova = anova(familyOnlyMdl,'summary');
        
        statFamily(ii).R2 = familyOnlyMdl.Rsquared;
        statFamily(ii).F  = familyOnlyAnova.F(2); %jenny changed to (2) from (end)because was always NaN
        statFamily(ii).P  = familyOnlyAnova.pValue(2); %jenny changed to (2) from (end)because was always NaN
    end
end
warning(iterWarn);
%warning(regressWarn);
%warning(conditionedWarn);

%         figure;
%{
for ii = 1:nStats
    %fprintf('\tTesting significance on %s... ', statNames{ii});
    iName = statNames{ii};
    
    % compiled stat is the generic, feature or syllable-based value
    compiledStat = [statistics.(iName)];
    
    % we regress on the spike counts or spike rates 
    regressedStat = evSpikeRates;
    if ~all(size(compiledStat) == size(regressedStat)),
        fprintf('\nWarning: %s is not acceptable - size mismatch.\n', iName);
    end
    
    % logistic regression on binary
    if all(compiledStat == 0 | compiledStat == 1)
        if all(compiledStat == 0) || all(compiledStat == 1)
            pCorr(ii) = NaN;
            continue;
        end
        
        % note: the test statistics that come from this logistic fit are the same as those
        % that come out of a simple t-test on (is vs. is not decisions)
        
        % this is reversed because compiledStat is a binary variable
        [logitfit, ~, fitstats] = glmfit(regressedStat', compiledStat', 'binomial','link','logit');
        
        % signed p-value
        pCorr(ii) = fitstats.p(2) * sign(logitfit(2)); 
        

%         subplot(nr, nc, ii);
        
        % plot spike count vs. statistic
%         xModel = linspace(min(regressedStat), max(regressedStat), 500);
%         yModel = glmval(logitfit, xModel, 'logit');
%         plot(regressedStat, compiledStat, 'r.', ...,
%             xModel, yModel, 'b--','HandleVisibility', 'off');
%         hold on;
%         ylabel(sprintf('%s', iName))
%         xlabel('spike rate (Hz)')
%         plot(basalRate * [1 1], ylim, 'g-');
%         hold off;
    else % linear regression on multi-valued discrete/continuous variable
        [linfit, ~, fitstats] = glmfit(compiledStat', regressedStat'); 
        
        % signed p-value for significance
        pCorr(ii) = fitstats.p(2) * sign(linfit(2));
        
%         subplot(nr, nc, ii);
        
        % plot spike count vs statistic 
%         plot(compiledStat, regressedStat, 'r.', ...,
%             compiledStat, linfit(1) + compiledStat * linfit(2), ...
%             'b--','HandleVisibility', 'off');
%         hold on;
%         xlabel(sprintf('%s', iName))
%         ylabel('spike rate (Hz)')
%         plot(xlim, basalRate * [1 1], 'g-');
%         hold off;
    end
%     fprintf('\n\t\tStatistic [%s], fit for trend, t = %0.3f, p = %0.3f\n', ...
%         iName, fitstats.t(2), fitstats.p(2));
    
%    drawnow
end
%}

