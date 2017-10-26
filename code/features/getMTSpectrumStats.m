 function spectrumData = getMTSpectrumStats(sample, specParams)
%GETMTSPECTRUMSTATS calculates spectrogram and associated audio metrics
%
%     function spectrumData = getMTSpectrumStats(sample, specParams) takes
%     in a sample (as a linear array) and specParams, and returns a struct
%     array of data about the sample clip, based on the specParams.
%     The specParams structure should include the following fields:
%             fs: the sample rate of the clip
%     windowSize: the length of the spectrogram window (in ms)
%       nOverlap: the length of the spectrogram window overlap (in ms)
%      freqBands: (optional) the locations of the frequency bands (useful
%                 for limiting the frequency range)
%     NfreqBands: if freqBands is not defined, the number of frequency bands
%     multitaper: [t/f] whether or not to use multitaper methods
%pitchValueIfNaN: the pitch value if YIN gives no value / a NaN
%       features: (optional) determines which features should be calculated
%                 and returned.  psd, freqs and times are mandatory fields
%
%     In defaultParams, specParams are substructures of the parameter
%     structure, called rough, inter, noiseReduce, fine, etc.
%
%     The structure returned from spectrumData contains the following
%     fields:
%                 psd: the short-time spectrogram data in a 2D array,
%                      whose size is numel(freqs) x numel(times)
%               freqs: the set of frequencies over which the spectrogram is
%                      calculated
%               times: the set of times over which the spectrogram windows
%                      begin (or are centered)
%     totalPower (opt.): the volume (in normalized units) of the signal
%  wienerEntropy (opt.): the Wiener entropy of the frequency content
%          deriv (opt.): the short-time derivative spectrogram data in a 2D
%                      array, whose size is numel(freqs) x numel(times)
%   The rest of the features require that deriv be included.
%     centerFreq (opt.): the amplitude-weighted instantaneous mean frequency
%fundamentalFreq (opt.): an estimate of the fundamental frequency (from YIN)
%   aperiodicity (opt.): an estimate of the amount of non-periodic energy in
%                      the signal (from YIN)
%  harmonicPitch (opt.): an estimate of the fundamental frequency (from
%                      cepstral methods)
%  pitchGoodness (opt.): an estimate of the amount of periodic energy in the
%                      signal (from cepstral methods)
%           mfcc (opt.): mel-frequency cepstral coefficients
%       waveform (opt.): the actual waveform

% check if chronux is in path
ourPaths = path;
if specParams.multitaper && all(isempty(strfind(ourPaths,'chronux')))
    error('ChronuxNotIncluded',['Chronux toolbox is not included in path,' ...
        'visit http://chronux.org/'])
end;

% default for features parameter
if ~isfield(specParams,'features')
    specParams.features = {'totalPower','deriv','wienerEntropy','mfcc','centerFreq','fundamentalFreq','harmonicPitch'};
    %specParams.features = {'totalPower','deriv'};
end


% calculate window sizes
win_s = specParams.windowSize / 1000;
overlap_s = specParams.nOverlap / 1000;

if win_s < numel(sample) / specParams.fs && ...
           numel(sample) / specParams.fs < win_s + (win_s - overlap_s)
    overlap_s_new = 1.001 * (2 * win_s - numel(sample) / specParams.fs); % small adjustment in overlaps
    warning('Changing overlap from %0.2f ms to %0.2f ms', overlap_s * 1000, overlap_s_new * 1000);
    overlap_s = overlap_s_new;
end
if specParams.multitaper
    nTapers = 1; timeBandwidth = 1;
    mtparams = struct('tapers',[timeBandwidth, nTapers],'Fs',specParams.fs,'err',0,'pad',4);
    if isfield(specParams,'freqBands'),
        mtparams.fpass = [min(specParams.freqBands) max(specParams.freqBands)];
    else
        %TODO: implement Nfreqbands
    end
    
    
    [spectrumData.psd, spectrumData.times, spectrumData.freqs] = ...
        mtspecgramc(sample, [win_s win_s - overlap_s],mtparams);
    
    if isempty(spectrumData.times), 
        warning('getMTSpectrumStats::tooshort', 'sample is too short');
        spectrumData = [];
        return;
    end;
    spectrumData.psd = spectrumData.psd';
    spectrumData.freqs = spectrumData.freqs';
else % no multi-taper
    winSS = floor(win_s * specParams.fs);
    overlapSS = floor(overlap_s * specParams.fs);
    if isfield(specParams,'freqBands')
        [spectrumData.spectrum, spectrumData.freqs, spectrumData.times, spectrumData.psd] = ...
            spectrogram(sample, winSS, overlapSS, specParams.freqBands, specParams.fs);
    else
        [spectrumData.spectrum, spectrumData.freqs, spectrumData.times, spectrumData.psd] = ...
            spectrogram(sample, winSS, overlapSS, specParams.NfreqBands, specParams.fs);
    end
end

%% get total power
if any(strcmp(specParams.features,'totalPower'))
    spectrumData.totalPower = getTotalPower(spectrumData);
end

%% get wiener entropy: high = noise, low = pure tone, mid = stack/chirp ramp
if any(strcmp(specParams.features,'wienerEntropy'))
    spectrumData.wienerEntropy = getEntropy(spectrumData);
end

%% get mfcc - cepstral coefficients (note: these are vectored for each time point)
if any(strcmp(specParams.features,'mfcc')) 
    % kicking out to external function
    spectrumData.mfcc = getMFCC(spectrumData, specParams);
end

%% get derivative spectrogram
if any(strcmp(specParams.features,'deriv'))
    [spectrumData.deriv spectrumData.mTD, spectrumData.mFD, spectrumData.FM, spectrumData.AM spectrumData.rawAM] ...
        = getSpectralDerivative(sample);
end

%% get center frequency - simply weighted average by spectrogram
if any(strcmp(specParams.features,'centerFreq'))
    spectrumData.centerFreq = getCentroidFreq(spectrumData);
end

%% get pitch goodness
if any(strcmp(specParams.features,'harmonicPitch'))
    [spectrumData.pitchGoodness, spectrumData.harmonicPitch] = ...
        getHarmonicPitch(spectrumData, specParams);
end

%% get fundamental frequency
if any(strcmp(specParams.features,'fundamentalFreq'))
    [spectrumData.aperiodicity, spectrumData.fundamentalFreq] = ...
        getYINPitch(sample, spectrumData, specParams);
end

%% pass along waveform
if any(strcmp(specParams.features,'waveform'))
    spectrumData.waveform = sample;
end


%% function definitions
    function power = getTotalPower(spectrum)
        dFreq = spectrum.freqs(2) - spectrum.freqs(1);
        power = sqrt(dot(spectrum.psd,spectrum.psd)) * dFreq; % amp factor
    end

    function centerfreq = getCentroidFreq(spectrum)
        % weighted average of frequency
        centerfreq = spectrum.freqs' * abs(spectrum.deriv) ./ sum(abs(spectrum.deriv));
    end

    function entropy = getEntropy(spectrum)
        % Wiener entropy, defined as the log ratio of GM to AM of the power
        % pure white noise means that Wiener entropy should be 0, and pure tone
        % should have -inf Wiener entropy
        
        AMpsd = mean(spectrum.psd,1);
        logGMpsd = mean(log(spectrum.psd),1);
        entropy = logGMpsd - log(AMpsd);
    end

% get mag gradient of spectral derivative, as well as frequency Modulation
    function [deriv, mTD, mFD, freqMod, ampMod, rawAmpMod] = getSpectralDerivative(sample)
        % use multitaper estimates, if available
        if specParams.multitaper
            nTapers = 3; timeBandwidth = 2;
            mtparams.tapers = [timeBandwidth, nTapers];
            [specDerivs, foo, bar] = ...
                mtdspecgramc(sample, [win_s win_s - overlap_s],[0 pi/2],mtparams); %#ok<NASGU>
            
            timeDeriv = squeeze(specDerivs(1,:,:))';
            freqDeriv = squeeze(specDerivs(2,:,:))';
        else
            % FIXME: what happens if spectrumData.times is only one
            % element?
            % just take the discrete derivative of the spectrum
            dt = spectrumData.times(2) - spectrumData.times(1);
            df = spectrumData.freqs(2) - spectrumData.freqs(1);
            freqDeriv = conv2([1 0 -1], 1, spectrumData.psd, 'same') / df;
            timeDeriv = conv2(1, [1 0 -1], spectrumData.psd, 'same') / dt;
        end
        % get modulation at each time point
        mTD = max(abs(timeDeriv.^2),[],1);
        mFD = max(abs(freqDeriv.^2),[],1);
        %mTD = max(abs(timeDeriv),[],1);
        %mFD = max(abs(freqDeriv),[],1);
        
        % get maximum amplitude of derivative at each time point and get
        % freqModulation
        derivSqrd = timeDeriv.^2 + freqDeriv.^2;
        [foo, imax] = max(abs(derivSqrd), [], 1);
        %freqMod = atan2(abs(timeDeriv(imax)), abs(freqDeriv(imax)));
        freqMod = atan2(mTD, mFD);
        % get fractional parts 
        tPart = sin(freqMod);
        fPart = cos(freqMod);
        freqMod = freqMod * 180 / pi;
        
        % calculate amplitude modulation
        
         dt = spectrumData.times(2) - spectrumData.times(1);
        ampMod = -sum(timeDeriv,1) ./ getTotalPower(spectrumData) / (dt*1000);
        rawAmpMod = -sum(timeDeriv,1) / (dt*1000);
        
        % compute derivative as mixture of time and frequency components
        deriv = timeDeriv .* (ones(numel(spectrumData.freqs),1) * tPart) + ...
            freqDeriv .* (ones(numel(spectrumData.freqs),1) * fPart);
    end

    function [goodness, pitch] = getHarmonicPitch(spectrum, params)
        
        % harmonic pitch is estimated from deriv-cepstrum.
        % doesn't seem to work too well
        dcepstrum = fft(spectrum.deriv,[],1);
        % take out the non-symmetric part 
        dcepstrum = dcepstrum(1:ceil(end/2),:);
        
        % pitch goodness is just the amount of power in the frequency band,
        % unscaled
        [goodness, pitch] = max(abs(dcepstrum),[],1);
        pitch = params.fs./pitch;
    end

% fundamental frequency, may have NaNs
    function [aperiodicity, pitch] = getYINPitch(sample, spectrum, params)
        % get pitch estimate
        pitchOut = yin(sample,params.fs);
        
        % FIXME: the addend is an empirical fudge factor to get pitch estimates to
        % align timewise, because of the different windowing that YIN does
        addedOffset = pitchOut.wsize/3 / params.fs;
        estTimes = (1:numel(pitchOut.f0))/(params.fs/pitchOut.hop) + addedOffset;
      
        % interpolation
        pitch = interp1(estTimes,440 * power(2,pitchOut.f0),spectrum.times);
        aperiodicity = interp1(estTimes,pitchOut.ap0,spectrum.times);
        
        % define NaN fillin - useful for later statistics sometimes
        aperiodicity(isnan(aperiodicity)) = 1;
        pitch(isnan(pitch)) = params.pitchValueIfNaN;
    end
end