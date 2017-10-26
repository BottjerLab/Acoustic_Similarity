function [covar, indivDists] = similarityScan(spectrum1, spectrum2, params, varargin)
%SIMILARITYSCAN Running, sample by sample scan of two spectral feature sets

%  COVAR = similarityScan(SPECTRUM1, SPECTRUM2) returns the similarity of two
%  sounds, sample by sample, according to the correlation found between
%  their features.  The correlations are taken along ever calculated
%  feature of the sample (power, entropy, average pitch, etc.) and is a
%  VECTOR of values, one for each possible offset of the spectrum relative
%  to the other one.
%  
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
	[covar, indivDists] = similarityScan(spectrum2, spectrum1, params);
    return;
end

%% find all features that can be correlated
features = whichFeatures(spectrum1);
nF = numel(features);

%% look at 	x) / var(y)) over the different fields, for each matching window
% TODO: make a better way to weight features
%weightVec =  [24.7557  629.1870   28.8207   24.3587  133.1410  679.5239   27.8380   63.5714   41.6093   55.1736];
% TODO: define similarity more robustly, with real euclidean distances and
% more weight on frequency measures
weightVec = ones(1,nF) / nF; % uniform weight, can be changed
%weightVec = weightVec/sum(weightVec);
% determine the time intervals that we want to run similarity
nSkipMs = 0.2; % in ms
nSkipSamples = floor(nSkipMs/1000 * params.fs); 
nComparisons = numel(1:nSkipSamples:len2 - len1 + 1);
covar =  zeros(1,nComparisons);
ctr = 1;
indivDists = zeros(nComparisons, nF);
for ii = 1:nSkipSamples:len2 - len1 + 1
  segStart = ii;
  segEnd = ii + len1 - 1;
  for jj = 1:nF % loop over the features
    fv1 = spectrum1.(features{jj});
    fv2 = spectrum2.(features{jj})(segStart:segEnd);
    indivDists(ctr,jj) = 1 - pdist([fv1;fv2],'cosine');
    %corr=cc(2,1);
    % this will happen when one of the vectors is constant
%    if isnan(corr), corr = 0; end
    covar(ctr) = covar(ctr) + indivDists(ctr,jj) * weightVec(jj);
  end
  ctr = ctr+1;
end

function ret = featsToArray(featSummary,features)
    nF = numel(features);
    nS = numel(featsToArray.times);
    ret = zeros(nF,nS);
    for ii = 1:nF
        ret(ii,:) = featSummary.(features{ii});
    end