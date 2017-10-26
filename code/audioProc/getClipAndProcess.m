function [clip, fs] = getClipAndProcess(songStruct, region, params, varargin)
% GETCLIPANDPROCESS add pre/post roll, bandpass, and filter noise
%
%   function clip = getClipAndProcess(songStruct, region, params) returns the
%   clean clip belonging to a specific region, as preprocessed by the
%   default params.
% 
%   revision note: an empty songStruct first argument can be given if the
%   region structure has a 'file' field which will access partial read data
if nargin < 3
    params = defaultParams;
end

params = processArgs(params, varargin{:});

region = addPrePost(region, params);
% either filter noise or get raw values
[clip, fs] = getClip(region, songStruct);

params.fs = fs;

minFiltLength = (params.noiseReduce.windowSize * 2 - ...
    params.noiseReduce.nOverlap) / 1000 + params.nps.attack; % minimum length in seconds
if params.doFilterNoise && length(clip) / fs > minFiltLength
    if isempty(songStruct) && isempty(params.noiseFilter) && isfield(region,'file');
        [filedir filestem] = fileparts(region.file);
        % TODO: check here
        manifest = struct('name','noiseMask','originalFile',[filedir filesep 'noiseMask-' filestem '.mat']);
        params.noiseFilter = loadFromManifest(manifest,'noiseMask');
    end
    if isempty(params.noiseFilter)
        warning('getClipAndProcess:noNoiseFilter','Cannot find noise filter, skipping filtering for clip...');
    else
        clip = noiseGate(clip, params.noiseFilter);   
    end
end

% band pass the result
clip = highPassSample(clip, params);
clip =  lowPassSample(clip, params);

function reconstructClip = noiseGate(clip, noiseProfile)
%NOISEGATE  Extracts a noise filtered version of a clip
%   clip = noiseGate(songStruct, region, noiseProfile) returns a clip which
%   has been cleaned of its noise, following certain parameters
%   The spectral noise filter is equivalent to a bank of narrow band-passed
%   limiters for each frequency band.  This has the general advantage of 
%   precisely defining syllable boundaries in noisier recording conditions.
%   
%   noiseProfile is the 1D reduced FFT of noise given by noiseAnalysis.m

params.noiseReduce.fs = fs;  
spec = getMTSpectrumStats(clip, params.noiseReduce);

nps = params.nps;

% define regions where gating is required -todo, apply hysteresis
gain = ones(size(spec.spectrum)) * nps.reduction;
for ifreq = 1:numel(spec.freqs)
    gain(ifreq,spec.psd(ifreq,:) >= noiseProfile(ifreq)) = 0.0;
end

% smooth in frequency space
if nps.freqSmooth > 0.0
    dF = spec.freqs(2) - spec.freqs(1);
    fWindowSize = nps.freqSmooth * 8/dF; % 4 half widths in either direction
    fWindow = (-fWindowSize / 2 : fWindowSize / 2) * dF;
    fWindow = exp(- (fWindow / nps.freqSmooth).^2);
    fWindow = fWindow ./ sum(fWindow) ; % normalize
    gain = conv2(fWindow,[1],gain,'same');
end
% apply attack, hold, release
dT = spec.times(2) - spec.times(1);
attackRate  = -nps.reduction * dT / nps.attack;
releaseRate = -nps.reduction * dT / nps.release;
holdSamples = floor(nps.hold / dT);

for ifreq = 1:numel(spec.freqs)
    holdClock = 0;
    for iTp = 1:numel(spec.times)-1
        dG = gain(ifreq, iTp+1) - gain(ifreq, iTp);     
        if dG > attackRate % attack is on
            gain(ifreq, iTp+1) = gain(ifreq, iTp) + attackRate;
        elseif gain(ifreq,iTp + 1)== 0.0 % gate is open, hold is ready
            holdClock = holdSamples;
        elseif holdClock > 0 % hold is being used up
            gain(ifreq,iTp+1) = 0.0;
            holdClock = holdClock - 1;
        elseif dG < -releaseRate
            gain(ifreq, iTp+1) = gain(ifreq, iTp) - releaseRate;
        end       
    end
end

% reduce volume in those bands
adjustedSpec = spec.spectrum .* ...
    (10.^(gain / 10));

% the inverse spectrogram does not always reproduce the exact volume, but
% this is pretty close to about +/- %1
scaleFac = 2.195; % empirically found by running spectrogram and inverse
winSS = floor(params.noiseReduce.windowSize * params.noiseReduce.fs/1000);
overlapSS = floor(params.noiseReduce.nOverlap * params.noiseReduce.fs/1000);
reconstructClip = invspecgram(adjustedSpec, params.noiseReduce.NfreqBands, ...
    params.noiseReduce.fs, winSS, overlapSS);
reconstructClip = reconstructClip / scaleFac;

% debug function to plot either the spectrogram or the shape of the
% spectral noise filter
    function hndl = plotGram(mat)
        mat = mat + eps;
        hndl = surf(spec.times, spec.freqs, 10*log10(abs(mat)),'EdgeColor','none');
        view(0,90);
        xlim([min(spec.times) max(spec.times)]); xlabel('Time (s)');
        ylim([min(spec.freqs) max(spec.freqs)]); ylabel('Frequency (Hz)');
        colorbar;
        title(sprintf('max=%0.3f,min=%0.3f',10*log10(abs(max(mat(:)))),...
            10*log10(abs(min(mat(:))))));
    end
end
end