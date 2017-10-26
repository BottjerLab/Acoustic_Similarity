function [pCorr, tCorr] = findCorrelates(statistics, DV)
%function findCorrelates(statistics, DV)
% massive fishing expeditions...

% take off the warning
pastWarn = warning('off','stats:glmfit:IllConditioned');
DV_SEM = std(DV)/sqrt(numel(DV)-1);

fprintf('DV (# = %d): %0.2f +/- %0.2f\n', ...
    numel(DV), mean(DV), DV_SEM);

% test significance for each statistic
statNames = fieldnames(statistics);
nStats = numel(statNames);
pCorr = zeros(nStats,1);

for ii = 1:nStats
    %fprintf('\tTesting significance on %s... ', statNames{ii});
    iName = statNames{ii};
    
    % compiled stat is the generic, feature or syllable-based value
    compiledStat = [statistics.(iName)];
    
    % we regress on the spike counts or spike rates
    regressedStat = DV;
    if ~all(size(compiledStat) == size(regressedStat)),
        fprintf('\nWarning: %s is not acceptable - size mismatch.\n', iName);
    end
    
    [linfit, ~, fitstats] = glmfit(compiledStat', regressedStat');

    % signed p-value for significance
    pCorr(ii) = fitstats.p(2) * sign(linfit(2));
    tCorr(ii) = fitstats.t(2);

    fprintf('\t\tStatistic [%s], fit for trend, t = %0.3f, p = %0.3f\n', ...
        iName, fitstats.t(2), fitstats.p(2));
    
    %drawnow
end

% undo the warning state
warning(pastWarn.state, pastWarn.identifier);