% fit syllable distributions
% based on veit/aronov/fee papers on breaths/etc

function [kstat, nIncluded, tPrime] = ksfitstat(X, fitRange)
% default? fitRange = [0.02 0.5]; % seconds
%warning('off', 'stats:lillietest:OutOfRangePLow');

MLEfunc = @(t, sampMean) t + fitRange(1) - sampMean - ...
    (diff(fitRange) * exp(-(diff(fitRange)/t))) / (1 - exp(-(diff(fitRange)/t)));

isWithin = X > fitRange(1) & X < fitRange(2);
sample = X(isWithin);
nIncluded = sum(isWithin);
tBounds = [0 10*diff(fitRange)]; % seconds
tPrime = fzero(@(x) MLEfunc(x, mean(sample)), tBounds); % the sufficient statistic
if versionNumber > 8
    [~,p,kstat] = lillietest(sample, 0.05, 'exp', 'MCTol', 1e-6);
else
    [~,p,kstat] = lillietest(sample, 0.05, 'exp', 1e-3);
end
kstat = kstat * sqrt(numel(sample)); 
end