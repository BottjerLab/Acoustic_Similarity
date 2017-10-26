function [Tau, nIncl] = fitExpOnInterval(values, interval)
    a = interval(1); b = interval(2);
    incl = values > a & values < b;
    nIncl = sum(incl);
    meanV = mean(values(incl));
    mleF = @(tau) tau + a - meanV - ((b - a) * exp(-(b - a)/tau)) / (1 - exp(-(b - a)/tau));
    Tau = fzero(mleF, mean(interval));
end