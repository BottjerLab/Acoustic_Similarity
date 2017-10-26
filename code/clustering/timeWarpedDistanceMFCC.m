function [totalDist, featureDists, matchInterval] = timeWarpedDistanceMFCC(spectrum1, spectrum2, params, varargin)
%% argument handling
if nargin < 3
    params = defaultParams;
end
params = processArgs(params,varargin{:});

% check if the 
assert(isfield(spectrum1,'mfcc') && isfield(spectrum2, 'mfcc'))
margin = params.maxWarpAllowed;
mf1 = spectrum1.mfcc;
mf2 = spectrum2.mfcc;
[totalDist, featureDists, matchInterval] = tWDistance(mf1, mf2);
