function candidateNoise = autodetectNoise(songStruct, soundRegions, params, varargin)

if nargin < 3, params = defaultParams; end;
params = processArgs(params, varargin{:});
fs = 1/songStruct.interval;
params.fs = fs;
%% find a noise region

candidateNoise = initEvents(1);
NOISE_LEN = params.noiseLength; % s
TIME_STEP = 0.5;
isQuiet = false;

regionStart = [soundRegions.start];
regionStop =  [soundRegions.stop];
tStamp = 1/fs;
if params.verbose, fprintf('Autodetecting noise ... '); end;
while ~isQuiet && NOISE_LEN > 0.5
    % are we colliding with an actual event?
    while any(regionStart < (NOISE_LEN + tStamp) & regionStop > tStamp)
        tStamp = tStamp + TIME_STEP;
        if tStamp + NOISE_LEN > songStruct.length * songStruct.interval %did we overrun the clip?
            NOISE_LEN = NOISE_LEN * 0.8; % shorten gradually
            TIME_STEP = TIME_STEP * 0.8;
            tStamp = 1/fs;
        end
    end
    
    candidateNoise.start = tStamp;
    candidateNoise.stop = candidateNoise.start + NOISE_LEN;
    candidateNoise.idxStart = floor(candidateNoise.start * params.fs); % sample indices
    candidateNoise.idxStop  = floor(candidateNoise.stop  * params.fs); % sample indices
    
    sample = getClip(candidateNoise, songStruct);

    % test noise
    params.rough.fs=params.fs;
    spectrum = getMTSpectrumStats(sample, params.rough);
    
    if params.plot, plotAllFigures(spectrum,[],params); 
        if params.pause, pause; end; 
    end;
    
    range = max(spectrum.totalPower) - min(spectrum.totalPower);
    variation = std(spectrum.totalPower);
    fprintf('Candidate [%0.1f-%0.1f]Current range is %0.2g, variability is %0.2g\n', ...
        candidateNoise.start, candidateNoise.stop, range, variation);
    
    % does it satisfy our criteria for noise, i.e.
    % small dynamic range and small variation in power?
    isQuiet = range < params.noiseDynamicRange && ...
        variation < params.noiseVariation;    
    tStamp = tStamp + TIME_STEP;
end

% play the sample
if ~params.quiet
    fprintf('Playing what should be all noise...\n');
    playSound(getClip(candidateNoise, songStruct), fs, true);
end

if params.plot,
    title('Current noise profile')
    if params.pause
        pause
    end
end
