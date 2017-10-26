function [totalDist, featureDists, matchInterval] = timeWarpedDistance(spectrum1, spectrum2, params, varargin)
%TIMEWARPEDDISTANCE Dynamic time-warped distance of two spectral feature sets
%
%  [totalDist, featureDists] = timeWarpedDistance(SPECTRUM1, SPECTRUM2) 
%  operates as standardDistance does, but takes into account some local time
%  warping measure as well, which allows two sequences do have different time courses.
%  
%  A warping cost can be implemented, as well as a cost for the difference in duration.
%
% totalDist is the distance (L2) over all the features, and featureDists are the
% component (L1 difference) distance.

% matchInterval returns the times (in sec) of the shorter interval (with relation to
% that interval) that best matches the other spectrum


%% argument handling
if nargin < 3
    params = defaultParams;
end
params = processArgs(params,varargin{:});

margin = params.maxWarpAllowed;
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
len1 = length(spectrum1.times);
len2 = length(spectrum2.times);

specArray1 = zeros(nF,len1);
specArray2 = zeros(nF,len2); % rows are different features, columns are time points
for ii = 1:nF
    if featCatalog(ii).doLog
        specArray1(ii,:) = log(spectrum1.(featCatalog(ii).name)) / featCatalog(ii).MAD;
        specArray2(ii,:) = log(spectrum2.(featCatalog(ii).name)) / featCatalog(ii).MAD;
    else
        specArray1(ii,:) = spectrum1.(featCatalog(ii).name) / featCatalog(ii).MAD;
        specArray2(ii,:) = spectrum2.(featCatalog(ii).name) / featCatalog(ii).MAD;
    end
end
[totalDist, featureDists, matchInterval] = tWDistance(specArray1, specArray2);
