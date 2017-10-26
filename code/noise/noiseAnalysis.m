function noiseFFT = noiseAnalysis(songStruct, region, params,varargin)
%NOISEANALYSIS find maximal noise FFT for noise reduction
%  noiseFFT = noiseAnalysis(songStruct, region) takes a particular region
%  of the recording and finds the maximum level of noise within each
%  frequency band.  It outputs various volume thresholds of noise in a N x 1
%  as a function of frequency within a given region.  
%  This does not work well if the region has any sounds, clicks, etc.
if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});

fs = 1/songStruct.interval;
clip = getClip(region, songStruct);

params.noiseReduce.fs = fs;
spec = getMTSpectrumStats(clip, params.noiseReduce);

%plotAllFigures(clip, spec, [], params);
%title(sprintf('PSD from %s--%s with centroid freq.', ...
%    sampleToTimeString(sampleStart), ...
%    sampleToTimeString(sampleEnd)));
%drawnow;

if params.playsample
    fprintf('Playing what should be a noise sample...');
    playSound(clip,fs,false);
    fprintf('Done.\n');
end
ws = params.nps.sampleWindow;

noiseFFT = zeros(1,numel(spec.freqs));
% get the absolute value for each
absPSD = abs(spec.psd);
% get the maximum noise level for each

%TODO: could still be optimized...
for ifreq = 1:numel(spec.freqs)
    smoothedFreqSpec = zeros(1,numel(spec.times) - ws + 1);
    ptr = 1;
    while ptr <= numel(spec.times) - ws + 1
        [minval, idx] = min(absPSD(ifreq,ptr:ptr+ws-1));
        smoothedFreqSpec(ptr:ptr+idx-1) = minval;
        ptr = ptr + idx;
    end
    noiseFFT(ifreq) = max(smoothedFreqSpec);
end
end