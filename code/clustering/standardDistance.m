function [totalDists, indivDists] = standardDistance(spectrum1, spectrum2, params, varargin)
%STANDARDDISTANCE Running, sample by sample scan of two spectral feature sets
%
%  COVAR = standardDistance(SPECTRUM1, SPECTRUM2) returns the similarity of two
%  sounds, sample by sample, according to the standardized difference found between
%  their features.  The correlations are taken along ever calculated
%  feature of the sample (power, entropy, average pitch, etc.) and is a
%  VECTOR of values, one for each possible offset of the spectrum relative
%  to the other one.
%
%  See also TIMEWARPEDDISTANCE
%% argument handling
if nargin < 3
    params = defaultParams;
end
params = processArgs(params,varargin{:});

%% order by length
% convert lengths
len1 = length(spectrum1.times);
len2 = length(spectrum2.times);
if(len1 > len2)
    temp = spectrum1; 
    spectrum1 = spectrum2; 
    spectrum2 = temp; 
    len1 = length(spectrum1.times);
    len2 = length(spectrum2.times);
end

%% name all features that should be compared
featCatalog = params.featureCatalog;

featureSel = false(1,numel(featCatalog));
for ii = 1:numel(featCatalog)
    featureSel(ii) = any(strcmp(fieldnames(spectrum1), featCatalog(ii).name)) && ...
        any(strcmp(fieldnames(spectrum2), featCatalog(ii).name));
end
featCatalog = featCatalog(featureSel);
nF = numel(featCatalog);

%% convert structure to array to speed indexing
% specArray1 = zeros(len1,nF);
% specArray2 = zeros(len2,nF);
% for ii = 1:nF
%     if featCatalog(ii).doLog
%         natVel1 = log(spectrum1.(featCatalog(ii).name));
%         natVel2 = log(spectrum2.(featCatalog(ii).name));
%     else
%         natVel1 = spectrum1.(featCatalog(ii).name);
%         natVel2 = spectrum2.(featCatalog(ii).name);
%     end
%     specArray1(:,ii) = (natVel1 - featCatalog(ii).median)/featCatalog(ii).MAD;
%     specArray2(:,ii) = (natVel2 - featCatalog(ii).median)/featCatalog(ii).MAD;
% end
% clear natVel1 natVel2
%% score differences
% determine the time intervals that we want to run similarity
nSkipMs = 0.4; % in ms
nSkipSamples = ceil(nSkipMs/1000 / (spectrum1.times(2) - spectrum1.times(1)));
nComparisons = numel(1:nSkipSamples:len2 - len1 + 1);

indivDists = zeros(nComparisons, nF);
for jj = 1:nF % loop over the features
    % get standardized vectors
    fv1 = getStdVector(spectrum1, featCatalog(jj));
    fullfv2 = getStdVector(spectrum2,featCatalog(jj));

    % loop over every nSkipSamples of the longer sound
    % (restrict to 1 for completeness)    
    ctr = 1;
    for ii = 1:nSkipSamples:len2 - len1 + 1
        indivDists(ctr,jj) = norm(fv1-fullfv2(ii:(ii + len1 - 1)))/len1;
        ctr = ctr+1;
    end
end

%% duration distance - append this to the distances
durDev = 40; % this is an empirical number based on Lb189_4_25_5 
durationDist = abs(spectrum1.times(end) - spectrum2.times(end)) / (durDev / 1000);
indivDists(:,nF+1) = durationDist;

%% combine scores along features - where weighting can be directly controlled
weightVec = [params.featureCatalog.weights];  
weightVec(end+1) = 0.5; %length is discouraged a little
weightVec = weightVec/sum(weightVec);
totalDists = weightVec * indivDists';

end

function fv = getStdVector(spectrum, featEntry, sampStart, sampEnd)
if nargin == 2
    sampStart = 1;
    sampEnd = numel(spectrum.(featEntry.name));
end
fv = spectrum.(featEntry.name)(sampStart:sampEnd);
if featEntry.doLog
    fv = log(fv);
end
fv = (fv - featEntry.median) ./ featEntry.MAD;
end
