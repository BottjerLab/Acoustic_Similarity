function sim = similarityRegions(songStruct, keyreg, otherRegs, noiseGate, params, varargin)
%% parameter handling
if nargin < 5; params = defaultParams; end
params = processArgs(params, varargin{:});

fs = 1/songStruct.interval;

%% preprocessing
% filter out noise first
[keyclip,keyspec] = cleanClip(keyreg);


for ii = 1:numel(otherRegs)
  [iclip,ispec] = cleanClip(otherRegs(ii));

  sim{ii} = similarityScan(keyspec, ispec);

end
function [clip,spec] = cleanClip(region)

region = addPrePost(region,params);
if nargin >= 3
    clip = noiseGate(songStruct, region, noiseProfile);
else
    clip = getClip(region, songStruct);
end

if ~params.quiet, playSound(clip,fs); end

% high pass the result
params.fs=fs;
clip = highPassSample(clip,params);

%% spectral analysis
params.fine.fs = fs;
spec = getSpectrumStats(clip, params.fine);

end
end